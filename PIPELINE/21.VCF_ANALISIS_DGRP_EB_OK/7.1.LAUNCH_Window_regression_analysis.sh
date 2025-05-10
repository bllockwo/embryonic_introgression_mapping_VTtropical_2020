#!/bin/sh
#
#SBATCH -J pcacor
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=15G
#SBATCH --time=3:00:00
#SBATCH -o ./slurmOut/pcacor.%A_%a.out # Standard output
#SBATCH -p general
#SBATCH --array=1-560

module load Rgeospatial

Rscript --vanilla 7.Window_regression_analysis.R \
${SLURM_ARRAY_TASK_ID} 
