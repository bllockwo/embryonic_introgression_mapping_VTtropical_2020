#!/usr/bin/env bash
#
#SBATCH -J collecD # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fst.%A_%a.out # Standard output

module load Rtidyverse

Rscript \
--vanilla \
2.collect_drift.R

date
echo "done"