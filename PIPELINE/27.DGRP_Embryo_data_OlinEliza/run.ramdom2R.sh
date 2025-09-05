#!/usr/bin/env bash  
#  
#SBATCH -J ran2r  
#SBATCH -c 6  
#SBATCH -N 1 # on one node  
#SBATCH -t 8:00:00   
#SBATCH --mem 20G   
#SBATCH -o ./slurmOutput/%x.%A_%a.out  
#SBATCH -p general  
#SBATCH --array=1-53

module load Rtidyverse 

Rscript --vanilla Randomization_test_fix2R.R ${SLURM_ARRAY_TASK_ID}


