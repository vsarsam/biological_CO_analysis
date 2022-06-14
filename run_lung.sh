#!/bin/bash
#
#SBATCH --time="UNLIMITED"
#SBATCH -J run_lung
#SBATCH -c 40

srun Rscript lung_atlas.R lung_cell_atlas/seurat_lung.RDS cll.obo annotations_droplet.csv $1 data/jac_matrix.tsv "originalexp" "data" data/cell_type_list.RDS
