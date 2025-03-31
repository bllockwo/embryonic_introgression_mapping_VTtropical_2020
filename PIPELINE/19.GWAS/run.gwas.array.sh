#!/usr/bin/env bash  
#  
#SBATCH -J gwas  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=2-117

module load Rtidyverse 

Rscript --vanilla fullgrm.GWAS.R ${SLURM_ARRAY_TASK_ID}


