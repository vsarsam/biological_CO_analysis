#!/bin/bash
#
#SBATCH --time="UNLIMITED"
#SBATCH -J run_Tabula_Muris
#SBATCH -c 25
#SBATCH --mem 100008M

srun Rscript old.R droplet/ cll.obo annotations_droplet.csv $1 data/jac_matrix.tsv "RNA" "counts"