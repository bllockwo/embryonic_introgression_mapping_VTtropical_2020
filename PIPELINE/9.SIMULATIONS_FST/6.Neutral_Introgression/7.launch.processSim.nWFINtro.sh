#!/usr/bin/env bash
#
#SBATCH -J nWF.ag # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fst.%A_%a.out # Standard output
#SBATCH --array=2-608%50

module load Rtidyverse

Rscript \
--vanilla \
6.proces.Simulations.nWFIntro.R \
${SLURM_ARRAY_TASK_ID}


