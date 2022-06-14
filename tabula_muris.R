args=commandArgs(TRUE)

# LOADING LIBRARIES
library(stats)
library(dplyr)
library(Seurat)
library(ontologyIndex)
library(SnowballC)
library(lsa)
library(corrplot)
library(scran)
library(SCnorm)
library(scuttle)
library(DESeq2)
library(rliger)
library(harmony)
library(ggplot2)
library(pdfCluster)
library(R.filesets)
library(doParallel)
library(iterators)
library(foreach)


if (length(args) != 8) {
  stop("Error: 8 arguments required")
}

# ARGUMENTS
# path to 10X data directory
input <- args[1]
# path to Cell Ontology obo file
cell_ont <- args[2]
# path to annotation file
annotations <- args[3]
# CLUSTERING PARAMETERS
# integer with end of dimension
dimension_end <- as.numeric(args[4])
# get dimension with dimension end -> 1:dimension_end
dimension <- 1:as.numeric(dimension_end)
# integer with resoltion value
resolution <- as.double(args[5])
# string of which assay is used (e.g 'RNA', 'originalexp', ..)
which_assay <- args[7]
# string of which data is used (e.g 'counts', 'data', ..)
counts_data <- args[8]


# -----------------------------------------functions----------------------------------------------

# load and prepare data
load_pbmc <- function(input){
  pbmc.data <- Read10X(data.dir = input)
  project_name <- strsplit(basename(input), "-", 2)[[1]][2]
  # init seurat obect with non-normalized data
  pbmc <- CreateSeuratObject(counts=pbmc.data, project = project_name, min.cells = 3, min.features = 350)
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  # scaling the data
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  return(pbmc)
}

# SCRAN 
scran_normalize <- function(gene_expression){
  # transform to single cell experiment object
  scran_gene_expr <- SingleCellExperiment::SingleCellExperiment(assays = list('counts' = gene_expression))
  # cursory quality control to remove cells with low total counts
  qcstats <- perCellQCMetrics(scran_gene_expr)
  qcfilter <- quickPerCellQC(qcstats)
  #scran_gene_expr <- scran_gene_expr[,!qcfilter$discard]
  # clustering before normalization 
  clusters <- quickCluster(scran_gene_expr)
  scran_gene_expr <- computeSumFactors(scran_gene_expr, cluster = clusters)
  # normalization
  scran_gene_expr_normalized <- logNormCounts(scran_gene_expr)
  # get counts
  counts <- logcounts(scran_gene_expr_normalized)
  return(as.data.frame(counts))
}

# SCNORM NOCH NICHT FUNKTIONARE
scnorm_normalize <- function(gene_expression){
  # transform data frame to matrix
  matrix_gene_expr <- as.matrix(gene_expression)
  # transform to single cell experiment object
  scnorm_gene_expr <- SingleCellExperiment::SingleCellExperiment(assays = list('data' = matrix_gene_expr))
  conditions = rep(c(1), each= ncol(matrix_gene_expr))
  # check count-depth relationship
}

# DESeq2
deseq_normalize <- function(gene_expression){
  # create meta data for samples
  samples <- colnames(gene_expression)
  meta <- data.frame(samples)
  rownames(meta) <- meta$samples
  # add pseudo count to raw expression
  dds <- DESeqDataSetFromMatrix(countData = round(gene_expression), 
                                colData = meta,
                                design = ~1)
  # normalization
  dds <- estimateSizeFactors(dds)
  # get counts
  normalized <- counts(dds, normalized = T)
  return(as.data.frame(normalized))
}

# SEURAT
seurat_normalize <- function(seurat_obj, gene_expr){
  pbmc <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
  counts_norm <- as.data.frame(pbmc[[which_assay]]@data)
  colnames(counts_norm) <- colnames(gene_expr)
  return(counts_norm)
}

# prepare all batches and put in one list
prepare_batches <- function(batch_list){
  # create new list for all batches
  ob_list <- c()
  for (batch in batch_list) {
    # load and prepare pbmc data
    batch_pbmc <- load_pbmc(input = batch)
    batch_name <- basename(batch)
    batch_pbmc[["batch"]] <- batch_name
    # sub sample data to speed up calculations??
    ob_list <- c(ob_list, batch_pbmc)
  }
  return(ob_list)
}

# LIGER
liger_batch_correction <- function(seurat_obj, dimension){
  # isolate RNA count matrices from each batch and format in a list
  seurat_list <- sapply(unique(seurat_obj@meta.data[["Library.Batch"]]), 
                        function(batch) GetAssayData(subset(seurat_obj, 
                                                            subset = Library.Batch == batch), 
                                                     slot=counts_data, assay=which_assay))
  # create liger object
  liger_obj <- createLiger(seurat_list)
  # normalize
  liger_obj <- rliger::normalize(liger_obj)
  # select variable genes
  liger_obj <- selectGenes(liger_obj, do.plot = T)
  # scale
  liger_obj <- scaleNotCenter(liger_obj)
  # find optimal k ~20 minutes
  #k.suggest <- suggestK(liger_obj, num.cores=4, gen.new=T, plot.log2=F)
  kval = 30
  # find optimal lamda ~ 
  #lambda.suggest <- suggestLambda(liger_obj, k=kval)
  lamda_val = 5
  # quantile align the factor loadings across datasets
  liger_obj <- optimizeALS(liger_obj, k = kval, lambda = lamda_val)
  liger_obj <- quantileAlignSNF(liger_obj)
  # add iNMF to seurat obj
  seurat_obj@reductions$iNMF <- CreateDimReducObject(
    embeddings = liger_obj@H.norm[colnames(seurat_obj),],
    loadings = t(liger_obj@W ),
    key= "iNMF",
    assay=which_assay
  )
  # UMAP of liger data
  liger_umap <- RunUMAP(seurat_obj, reduction='iNMF', dims = dimension)
  liger_umap <- FindNeighbors(liger_umap, reduction='iNMF')
  liger_umap <- FindClusters(liger_umap, resolution = 0.5)
  # plot cluster and batch with UMAP
  DimPlot(liger_umap, group.by='seurat_clusters', reduction='umap')
  DimPlot(liger_umap, group.by='Library.Batch', reduction='umap')
  return(liger_umap)
}

# Harmony
harmony_batch_correctioin <- function(seurat_obj, dimension){
  # perform PCA
  seurat_obj <- run_pca(seurat_obj)
  # run harmony on RNA assay 
  # plot_convergence -> make sure harmony function gets better each round
  seurat_obj <- seurat_obj %>% RunHarmony("Library.Batch", assay.use=which_assay, plot_convergence = T)
  # UMAP of harmony data
  harmony_umap <- RunUMAP(seurat_obj, reduction='harmony', dims = dimension)
  harmony_umap <- FindNeighbors(harmony_umap, reduction='harmony')
  harmony_umap <- FindClusters(harmony_umap, resolution = 0.5)
  # plot cluster and batch with UMAP
  DimPlot(harmony_umap, group.by='seurat_clusters', reduction='umap')
  DimPlot(harmony_umap, group.by='Library.Batch', reduction='umap')
  return(harmony_umap)
}

# Seurat Batch Effect Correction 
seurat_batch_correction <- function(batch_list, dimension){
  # prepare all batches
  ob_list <- prepare_batches(batch_list)
  # identify anchors, commonly shared variable genes and integrate samples
  seurat_obj <- FindIntegrationAnchors(object.list = ob_list)
  seurat_obj <- IntegrateData(anchorset = seurat_obj)
  DefaultAssay(seurat_obj) <- "integrated"
  seurat_obj <- ScaleData(seurat_obj, do.center = T, do.scale = F)
  seurat_obj <- RunPCA(seurat_obj)
  DimPlot(seurat_obj, dims = c(1, 2), reduction = "pca", split.by = "batch")
  # UMAP of seurat data
  seurat_cluster <- RunUMAP(seurat_obj, reduction='pca', dims = dimension)
  seurat_cluster <- FindNeighbors(seurat_cluster, reduction='pca')
  seurat_cluster <- FindClusters(seurat_cluster, resolution = 0.5)
  # plot cluster and batch with UMAP
  DimPlot(seurat_cluster, reduction = "umap", group.by = "seurat_clusters")
  DimPlot(seurat_cluster, reduction = "umap", group.by = "batch")
  return(seurat_cluster)
}

# cluster data with dimension and resolution as args
cluster_pbmc <- function(pbmc, dim, res){
  pbmc <- run_pca(pbmc)
  # finding clusters
  pbmc <- FindNeighbors(pbmc, dims = dim)
  pbmc <- FindClusters(pbmc, resolution = res)
  return(pbmc)
}

# run PCA
run_pca <- function(pbmc){
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  ElbowPlot(pbmc)
  return(pbmc)
}

# create gene expression matrix with one col and delete after adding all cols
build_cluster_matrix <- function(mat, cluster, idents){
  cluster_matrix <- data.frame(matrix(ncol = 1, nrow = nrow(mat)))
  # rownames are genes and colnames clustered barcodes
  rownames(cluster_matrix) <- rownames(mat)
  cols <- rownames(idents)[idents[1] == cluster]
  # add col from expression matrix
  for (col in cols) {
    cluster_matrix[col] <- select(mat, col)[[1]]
  }
  # remove first empty col
  cluster_matrix[1] <- NULL
  # get row wise sum
  return(rowSums(cluster_matrix))
}

# create dataframe for sums
sum_clustered <- function(mat, idents){
  sum_cluster_matrix <- data.frame(matrix(nrow = nrow(mat), ncol = length(unique(idents[[1]]))))
  # rownames same as in gene expression matrix
  rownames(sum_cluster_matrix) <- rownames(mat)
  for (col in 1:ncol(sum_cluster_matrix)) {
    # from every cluster matrix filter sum column and put in matrix
    extract_sum <- build_cluster_matrix(mat, col-1, idents)
    sum_cluster_matrix[col] <- extract_sum
  }
  return(sum_cluster_matrix)
}

create_big_matrix <- function(expression_matrix){
  big_pears <- data.frame(matrix(nrow = length(colnames(expression_matrix)), ncol = length(colnames(expression_matrix))))
  # names need to be same as in idents
  names <- paste0("cluster.", 1:ncol(expression_matrix)-1)
  colnames(big_pears) <- names
  rownames(big_pears) <- names
  return(big_pears)
}

# fill observed matrix with spearman
fill_spear <- function(expression_matrix, big_pears){
  # loop through every row and col and fill matrix
  for (i in 1:nrow(big_pears)){
    # only half matrix
    for (j in i:nrow(big_pears)) {
      # calc spearman for two columns
      spear <- cor.test(expression_matrix[[i]], expression_matrix[[j]], method = "spearman", exact = FALSE)
      # put spearman value in big general spearman matrix
      big_pears[i,j] <- round(spear$estimate[["rho"]], digits = 2)
    }
  }
  # fill symmetric matrix
  big_pears[lower.tri(big_pears)] <- t(big_pears)[lower.tri(big_pears)]
  return(big_pears)
}

# fill observed matrix with pearson
fill_pears <- function(expression_matrix, big_pears){
  # loop through every row and col and fill matrix
  for (i in 1:nrow(big_pears)){
    # only half matrix
    for (j in i:nrow(big_pears)) {
      # calc pearson for two columns
      spear <- cor.test(expression_matrix[[i]], expression_matrix[[j]], method = "pearson", exact = FALSE)
      # put pearson value in big general pearson matrix
      big_pears[i,j] <- round(spear$estimate[["cor"]], digits = 2)
    }
  }
  # fill symmetric matrix
  big_pears[lower.tri(big_pears)] <- t(big_pears)[lower.tri(big_pears)]
  return(big_pears)
}

# fill table with cell ids of data frame
fill_cl_ids <- function(table, id_df){
  table$cell_ontology_id <- id_df[table$cell_type,]
  return(table)
}

# get most occurring cl id for one cluster
get_most <- function(cluster, idents, anno_table){
  # filter for one cluster -> -1 clusters start at 0
  idents <- idents[idents[1] == cluster-1,]
  barcodes <- paste(idents[,"orig.name"], sapply(strsplit(rownames(idents), "-"), "[", 1), sep = "_")
  # barcodes to cl ids
  cl_ids <- as.vector(anno_table[barcodes, "cell_ontology_id"])
  cl_ids <- cl_ids[cl_ids!="" & !is.na(cl_ids)]
  return(names(sort(table(cl_ids), decreasing = T)[1]))
}


# get expected similarity with most occurring cl id from jac or cos matrix
get_expected <- function(jac_or_cos, observed, idents, anno_table){
  # number of clusters
  n_clusters <- length(unique(idents[[1]]))
  exp_cluster <- data.frame(matrix(nrow = n_clusters, ncol = n_clusters, dimnames = list(rownames(observed), colnames(observed))))
  for (cluster1 in 1:n_clusters){
    # get cl id which occurs most in cluster
    most_cl1 <- get_most(cluster1, idents, anno_table)
    for (cluster2 in cluster1:n_clusters) {
      most_cl2 <- get_most(cluster2, idents, anno_table)
      # get value from jaccard or cosine table
      # if most_cl is empty (some barcodes don't have a cl id)
      if (length(most_cl1)==0 | length(most_cl2)==0){
        exp_cluster[cluster1, cluster2] <- NA
      } else {
        exp_cluster[cluster1, cluster2] <- jac_or_cos[most_cl1, most_cl2]
      }
    }
  }
  # fill matrix
  exp_cluster[lower.tri(exp_cluster)] <- t(exp_cluster)[lower.tri(exp_cluster)]
  return(exp_cluster)
}

# get expected matrix for all clusters
expected_sim <- function(jac_or_cos, observed, idents, anno_table) {
  # number of clusters
  n_clusters <- nrow(observed)
  exp_cluster <- data.frame(matrix(nrow = n_clusters, ncol = n_clusters, dimnames = list(rownames(observed), colnames(observed))))
  for (cluster1 in 1:n_clusters) {
    # filter for one cluster -> -1 clusters start at 0
    barcodes1 <- rownames(idents)[idents[1] == cluster1-1]
    # barcodes to cl ids
    cl_ids1 <- as.vector(sapply(X = barcodes1, function(x) cell_to_cl(x, anno_table)))
    cl_ids1 <- cl_ids1[cl_ids1!="" & !is.na(cl_ids1)]
    # number of cl ids
    n <- as.numeric(length(cl_ids1))
    # get table with occurrences for each unique cl id
    occurences1 <- as.data.frame(table(cl_ids1))
    # number of unique cl ids
    uniques1 <- unique(cl_ids1)
    # get values from Jaccard or Cosine similarity matrix
    tmp1 <- jac_or_cos[uniques1, uniques1, drop = F]
    score1 <- quantify(tmp1, occx = occurences1, n)
    # half matrix
    for (cluster2 in cluster1:n_clusters) {
      # filter for second cluster
      barcodes2 <- rownames(idents)[idents[1] == cluster2-1]
      cl_ids2 <- as.vector(sapply(X = barcodes2, function(x) cell_to_cl(x, anno_table)))
      cl_ids2 <- cl_ids2[cl_ids2!="" & !is.na(cl_ids2)]
      m <- as.numeric(length(cl_ids2))
      # table with occurrences for each unique cl id
      occurences2 <- as.data.frame(table(cl_ids2))
      uniques2 <- unique(cl_ids2)
      cat(cluster1, cluster2, "\n")
      # select subset
      tmp2 <- jac_or_cos[uniques2, uniques2, drop = F]
      score2 <- quantify(tmp2, occx = occurences2, m)
      # fill expected matrix with distance
      exp_cluster[cluster1, cluster2] <- 1-(abs(score1-score2))
    }
  }
  # fill symmetric matrix
  exp_cluster[lower.tri(exp_cluster)] <- t(exp_cluster)[lower.tri(exp_cluster)]
  return(exp_cluster)
}

# get quantified value for occurrence table based on frequences 
quantify <- function(tablou, occx, n){
  for (x in 1:nrow(occx)) {
    x_val <- grep(occx[x,1], rownames(tablou))
    for (y in 1:nrow(occx)) {
      y_val <- grep(occx[y,1], colnames(tablou))
      tablou[x_val, y_val] <- tablou[x_val, y_val]*occx[x,"Freq"]*occx[y,"Freq"]
    }
  }
  return(round(sum(tablou)/n^2, 2))
}

# converts barcode to cl id from annotation table
cell_to_cl <- function(cell, anno_table){
  return(anno_table[anno_table$cell == cell, "cell_ontology_id"])
}

# create empty matrix for jaccard and cosine similarity matrix
create_all_matrix <- function(){
  # get all cl ids
  cl_ids <- ontology$id
  # create new empty matrix and change col and row names
  all_matrix <- data.frame(matrix(nrow = length(cl_ids), ncol = length(cl_ids)))
  colnames(all_matrix) <- cl_ids
  rownames(all_matrix) <- cl_ids
  return(all_matrix)
}

# fill matrix for jaccard
fill_matrix_jac <- function(){
  jac_matrix <- create_all_matrix()
  # loop through every row and col and fill matrix
  for (i in 1:nrow(jac_matrix)){
    # only half matrix
    for (j in i:nrow(jac_matrix)) {
      jac_matrix[i,j] <- compute_2_cl(i,j)
    }
  }
  # fill matrix
  jac_matrix[lower.tri(jac_matrix)] <- t(jac_matrix)[lower.tri(jac_matrix)]
  return(jac_matrix)
}

# fill matrix for cosine
fill_matrix_cos <- function(){
  cos_matrix <- create_all_matrix()
  all_cl <- colnames(cos_matrix)
  # loop through every row and col and fill matrix
  for (i in 1:nrow(cos_matrix)){
    for (j in i:nrow(cos_matrix)) {
      # calc cosine similarity
      cos_matrix[i,j] <- round(cosine(ancestor_representation(ontology$ancestors[[i]], all_cl), ancestor_representation(ontology$ancestors[[j]], all_cl)), digits = 2)
    }
  }
  cos_matrix[lower.tri(cos_matrix)] <- t(cos_matrix)[lower.tri(cos_matrix)]
  return(cos_matrix)
}


# JACCARD SIMILARITY
# get the intersection between two cell ontologies
compute_2_cl <- function(cl1, cl2){
  # get all ancestors for one cl id
  list1 <- ontology$ancestors[[cl1]]
  list2 <- ontology$ancestors[[cl2]]
  return(calculate_jaccard(list1, list2))
}

# compute jaccard intersection
calculate_jaccard <- function(list1, list2){
  inter <- intersect(list1, list2)
  percentage <- length(inter) / (length(list1) + length(list2) - length(inter))
  return(round(percentage, digits = 2))
}

# get mean coordinates for one cluster
mean_coordinates <- function(cluster, coordinates, idents){
  # extract barcodes and coordinates for one cluster
  barcodes <- rownames(idents)[idents[1] == cluster-1]
  cluster_coordinates <- subset(coordinates, rownames(coordinates) %in% barcodes)
  # get mean value for x and y axe
  x <- mean(as.matrix(cluster_coordinates[1]))
  y <- mean(as.matrix(cluster_coordinates[2]))
  return(c(x,y))
}

# euclidean distance
euclidean <- function(a, b) sqrt(sum((a - b)^2))

# transform barcode to cell ontology id
transform_to_id <- function(i, healthy_anno, cell_list, barcodes){
  barcode <- as.numeric(which(rownames(healthy_anno) == barcodes[i]))
  cell_type <- healthy_anno[barcode,]$cell_type
  cell_id <- as.vector(cell_list[cell_type])
  return(cell_id)
}

# calculate distance matrix for all clusters
distance_matrix_cluster <- function(coordinates, observed, idents){
  n_cluster <- as.numeric(nrow(observed))
  distance_matrix <- data.frame(matrix(nrow = n_cluster, ncol = n_cluster, dimnames = list(rownames(observed), rownames(observed))))
  for (cluster in 1:n_cluster) {
    # get mean position of cluster
    coordinates_a <- mean_coordinates(cluster, coordinates, idents)
    print(coordinates_a)
    for (cluster2 in cluster:n_cluster) {
      coordinates_b <- mean_coordinates(cluster2, coordinates, idents)
      # compute euclidean distance between two clusters
      distance_matrix[cluster, cluster2] <- euclidean(coordinates_a, coordinates_b)
    }
  }
  distance_matrix[lower.tri(distance_matrix)] <- t(distance_matrix)[lower.tri(distance_matrix)]
  return(distance_matrix)
}



# COSINE SIMILARITY
# create binary list with 0 and ancestors 1
ancestor_representation <- function(ancestors, all){
  # list with length of CL ids
  list <- numeric(length(all))
  # for every ancestor set to 1
  list[match(ancestors, all)] <- 1
  return(list)
}

# grid run to find optimal values for clustering 
grid_run <- function(samples, dimensions, resolutions, jac_or_cos, gene_expression){
  runs <- data.frame(matrix(ncol = 3))
  colnames(runs) <- c("dimension", "resolution", "score")
  for (dim in dimensions) {
    for (res in resolutions) {
      score <- run_single(all_samples, dim, res, jac_or_cos, gene_expression)
      d <- paste(c(dim[1], dim[length(dim)]), collapse = ":")
      runs <- rbind(runs, c(d, res, score))
      cat(dim, res, score, "\n")
    }
  }
  write.table(runs, file = "results.tsv", sep = "\t")
  return(runs)
}

# prepare sample list 
prepare_samples <- function(batch_samples){
  sample_list <- c()
  for (sample in 1:nrow(batch_samples)) {
    sample_list <- c(sample_list, batch_samples[[2]][sample])
  }
  # rename 
  names(sample_list) <- batch_samples[,1]
  return(sample_list)
}

# run all batch effection
run_all_batch <- function(seurat_obj, batch_list, dimension, resolution, jac_or_cos, gene_expression, anno_table, name){
  # run Seurat
  print(Sys.time())
  cat("starting Seurat Batch Correction of ", name)
  seurat_umap <- seurat_batch_correction(batch_list, dimension)
  seurat_correlation <- run_single(seurat_umap, dimension, resolution, jac_or_cos, gene_expression, anno_table, name)
  seurat_correlation <- round(seurat_correlation, digits = 5)
  rm(seurat_umap)
  print(Sys.time())
  cat("finished Seurat Batch Correction of", name)
  # run LIGER
  print(Sys.time())
  print("starting LIGER Batch Correction...")
  liger_umap <- liger_batch_correction(seurat_obj, dimension)
  liger_correlation <- run_single(liger_umap, dimension, resolution, jac_or_cos, gene_expression, anno_table, name)
  liger_correlation <- round(liger_correlation, digits = 5)
  rm(liger_umap)
  print(Sys.time())
  cat("finished LIGER Batch Correction of ", name)
  # run Harmony
  print(Sys.time())
  print("starting Harmony Batch Correction...")
  harmony_umap <- harmony_batch_correctioin(seurat_obj, dimension)
  harmony_correlation <- run_single(harmony_umap, dimension, resolution, jac_or_cos, gene_expression, anno_table, name)
  harmony_correlation <- round(harmony_correlation, digits = 5)
  rm(harmony_umap)
  print(Sys.time())
  cat("finished Harmony Batch Correction of " , name)
  # run with no batch correlation
  # cluster all samples
  seurat_obj <- cluster_pbmc(seurat_obj, dimension, resolution)
  original_correlation <- run_single(seurat_obj, dimension, resolution, jac_or_cos, gene_expression, anno_table, name)
  original_correlation <- round(original_correlation, digits = 5)
  # print all correlation values
  cat(name, " Original correlation value: ", original_correlation)
  cat(name, " Batch Effect Correction Methods", "\n")
  cat(name, " LIGER correlation value: ", liger_correlation)
  cat(name, " Harmony correlation value: ", harmony_correlation)
  cat(name, " Seurat correlation value: ", seurat_correlation)
  list <- c(original_correlation, liger_correlation, harmony_correlation, seurat_correlation)
  return(list)
}

# run all normalization methods
run_all_normalization <- function(gene_expression, seurat_obj, dimension, resolution, jac_or_cos, anno_table, name){
  #seurat_obj <- cluster_pbmc(seurat_obj, dimension, resolution)
  # run Seurat v3
  print(Sys.time())
  cat("starting Seurat Normalization of ", name)
  seurat_normalized <- seurat_normalize(seurat_obj, gene_expression)
  seurat_correlation <- run_single(seurat_obj, dimension, resolution, jac_or_cos, seurat_normalized, anno_table, name)
  seurat_correlation <- round(seurat_correlation, digits = 5)
  rm(seurat_normalized)
  print(Sys.time())
  cat("finished Seurat Normalization of ", name)
  # run SCRAN
  print(Sys.time())
  cat("starting SCRAN Normalization of", name)
  scran_normalized <- scran_normalize(gene_expression)
  scran_correlation <- run_single(seurat_obj, dimension, resolution, jac_or_cos, scran_normalized, anno_table, name)
  scran_correlation <- round(scran_correlation, digits = 5)
  rm(scran_normalized)
  print(Sys.time())
  cat("finished SCRAN Normalization of", name)
  # run DESeq2
  print(Sys.time())
  cat("starting DESeq2 Normalization of ", name)
  deseq_normalized <- deseq_normalize(gene_expression)
  deseq_correlation <- run_single(seurat_obj, dimension, resolution, jac_or_cos, deseq_normalized, anno_table, name)
  deseq_correlation <- round(deseq_correlation, digits = 5)
  rm(deseq_normalized)
  print(Sys.time())
  cat("finished DESeq2 Normalization of", name)
  # print all correlation values
  cat(name, " Normalization Methods", "\n")
  cat(name, " SCRAN correlation value: ", scran_correlation, "\n")
  cat(name, " DESeq2 correlation value: ", deseq_correlation, "\n")
  cat(name, " Seurat correlation value: ", seurat_correlation, "\n")
  list <- c(scran_correlation, deseq_correlation, seurat_correlation)
  return(list)
}

run_single <- function(seurat_obj, dimension, resolution, jac_or_cos, gene_expression, anno_table, name){
  # create output directory
  #out = paste("cluster", max(dimension), resolution, sep = "_")
  #out = ("cluster_3_0.5_cosine")
  #dir.create(out, showWarnings = F)
  # UMAP
  #umap <- RunUMAP(seurat_obj, dimension)
  #png(file.path(out, "UMAP.png"))
  #plot(DimPlot(umap, reduction = "umap", label = T)
  # idents reform
  idents <- as.data.frame(Idents(seurat_obj))
  idents$orig.name <- as.data.frame(seurat_obj$orig.ident)[[1]]
  idents$barcodes <- sapply(strsplit(rownames(idents), "-"), "[", 1)
  colnames(gene_expression) <- rownames(idents)
  cat("print idents",nrow(idents))
  cat("print gene expression", nrow(gene_expression))
  #rownames(idents) <- paste(idents$orig.name, idents$barcodes, sep = "_")
  # OBSERVED -> gene expression + spearman correlation
  cluster_expression <- sum_clustered(gene_expression, idents)
  cat("calculating observed similarities of ", name)
  spear <- fill_pears(cluster_expression, create_big_matrix(cluster_expression))
  cat("spear finished:", spear[[2]][1], " ", spear[[3]][1])
  #write.table(spear, file.path(out, "observed.tsv"), sep = "\t")
  # EXPECTED -> CL ids + jaccard/cos similarity
  cat("\n jacorcos ", jac_or_cos[[2]][1], " spear ", spear[[2]][1], " idents ", idents[[1]][1], " anno table ", anno_table[[1]][1])
  cat("calculation expected similarities of ", name)
  expected <- get_expected(jac_or_cos, spear, idents, anno_table)
  cat("expected finished:", expected[[2]][1], " ", expected[[3]][1])
  #seurat_obj <- RunTSNE(seurat_obj)
  #seurat_obj <- RunUMAP(seurat_obj, dims = dimension)
  #coordinates <- as.data.frame(seurat_obj@reductions[["tsne"]]@cell.embeddings)
  #spear_new <- distance_matrix_cluster(coordinates, expected, idents)
  #write.table(expected, file.path(out, "expected.tsv"), sep = "\t")
  expected <- expected[rowSums(expected, na.rm = T) != 0,colSums(expected, na.rm = T) != 0]
  print("fehler?")
  spear <- spear[rownames(expected),colnames(expected)]
  # CORRELATION
  cat("spear", spear[[2]][1], "\n", "expected", expected[[2]][1])
  correlation <- cor(spear, expected)
  # CORRELATION HEATMAP
  #cor_heatmap <- corrplot(correlation, method = "color")
  #png(file.path(out, "correlation_plot.png"))
  #plot(cor_heatmap)
  # CORRELATION VALUE
  # mean of absolute values from diagonal from correlation matrix
  cat(name, as.numeric(diag(correlation)), nrow(correlation))
  cor_score <- sum(abs(as.numeric(diag(correlation)))) / nrow(correlation)
  return(cor_score)
}

# for parallel, number of cores:25
my_cluster <- makeForkCluster(25, outfile = ".tabula.log")
registerDoParallel(my_cluster)

# for for all batches, batch effect correction and normalization methods
run_all_for_all <- function(sample_list, dimension, resolution, jac_or_cos, anno_table){
  # loop through every sample
  df <- foreach(i = seq_along(sample_list), .combine = rbind) %dopar% {
    name <- names(sample_list[i])
    print(name)
    sample <- sample_list[[name]]
    print(sample)
    # if sample has batches
    if (length(sample) > 1) {
      # PARSE ALL SAMPLES IN GIVEN DIRECTORY AND MERGE TO ONE BIG SEURAT OBJECT
      print("parse all samples...")
      all_samples <- NULL
      for (one in sample) {
        cat("processing: ", one, "\n")
        pbmc <- load_pbmc(one)
        pbmc[["Library.Batch"]] = basename(one)
        # first one
        if (is.null(all_samples)) {
          all_samples <- pbmc
        } else {
          # merge with last sample
          all_samples <- merge(all_samples, pbmc, project = "pairwise comparison")
        }
      }
      all_samples <- cluster_pbmc(all_samples, dimension, resolution)
      # get idents
      idents <- as.data.frame(Idents(all_samples))
      idents$orig.name <- as.data.frame(all_samples$orig.ident)[[1]]
      idents$barcodes <- sapply(strsplit(rownames(idents), "-"), "[", 1)
      rownames(idents) <- paste(idents$orig.name, idents$barcodes, sep = "_")
      # add pseudo count to gene expression
      # get gene_expr
      print("combining gene expression for clusters...")
      gene_expr <- as.data.frame(all_samples@assays[["RNA"]]@counts+1)
      colnames(gene_expr) <- rownames(idents)
      print("HAAAALT STOP")
      cat("first ident:", idents[[1]][1], " first gene:", gene_expr[[1]][1])
      cat("erstes idents", nrow(idents), "erstes gene expr", nrow(gene_expr))
      cat("\t first ident:", idents[[1]][1], " first gene:", gene_expr[[1]][1])
      # get all batch corrections and normalizations in lists
      batch_list <- run_all_batch(all_samples, sample, dimension, resolution, jac_or_cos, gene_expr, anno_table, name)
      normalization_list <- run_all_normalization(gene_expr, all_samples, dimension, resolution, jac_or_cos, anno_table, name)
      # fill row in data frame
      return(c(as.character(name), batch_list[1], batch_list[2], batch_list[3], batch_list[4], 
               normalization_list[1], normalization_list[2], normalization_list[3]))
    }
  }
  print("bis hier gut")
  stopCluster(my_cluster)
  return(df)
}

# -------------------------------------------------------------------------------------------------
# load CL reference data
ontology <- get_ontology(cell_ont)
anno_table <- read.csv(annotations)
rownames(anno_table) <- anno_table$cell


# CREATE JACCARD AND COSINE SIMILARITY IF THEY ARE NOT INPUT
if(file.exists(args[6])){
  print("load jac and cos matrices...")
  jac_or_cos <- read.table(args[6], sep = "\t", header = T)
  # CL.0000006 must be CL:0000006
  colnames(jac_or_cos) <- sapply(gsub("\\.", ":", colnames(jac_or_cos)), "[", 1)
}else {
  # create and fill matrices with Cosine and Jaccard similarity
  if (args[6] == "jac") {
    print("create and fill jaccard similarity matrix...")
    jac_or_cos <- fill_matrix_jac()
  } else {
    print("create and fill cosine similarity matrix...")
    jac_or_cos <- fill_matrix_cos()
  }
}

batch_samples <- as.data.frame(list.dirs(input, full.names = T, recursive = F))
colnames(batch_samples) <- "path"
batch_samples$sample <- basename(batch_samples[,1])
batch_samples$organ <- sapply(strsplit(batch_samples$sample, "-"), "[", 1)


batch_samples <- aggregate(batch_samples$path, by = list(organ=batch_samples$organ), FUN=c)
batch_samples <- prepare_samples(batch_samples)

df <- run_all_for_all(batch_samples, dimension, resolution, jac_or_cos, anno_table)

save.Object(df, "result.RData")

# dataset name to rownames
df <- data.frame(df, row.names = 1)
# change col names of df
colnames(df) <- c("original", "LIGER", "Harmony", "Seurat", "SCRAN", "DESeq2", "Seurat")

name <- paste("tabula_muris_neu", dimension_end, resolution, ".tsv", sep = "-")

write.table(df, file=name, sep = "\t", quote = F)

# # cluster all samples
# all_samples <- cluster_pbmc(all_samples, dimension, resolution)
# 
# # get gene_expr
# print("combining gene expression for clusters...")
# gene_expr <- as.data.frame(all_samples@assays[[which_assay]]@counts)
# colnames(gene_expr) <- rownames(idents)





# RUN ANALYSIS
#score <- run_single(all_samples, dimension, resolution, jac_or_cos, gene_expression)
#cat("Run with: ", "Dim: ", dimension, " ;Res: ", resolution, "\n",
    #"had correlation score of: ", score)

# RUN GRID OPTIONS
#results <- grid_run(samples = all_samples, 
#                    dimensions = list(1:5, 1:6, 1:7, 1:8, 1:9),
#                   resolutions = c(0.75, 1, 1.25, 1.5),
#                   jac_matrix, gene_expression)

# save.image(file.path(out, "all.Rdata"))
print("FINISH")