#!/usr/bin/env bash
#
#SBATCH -J collecO # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 90G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fst.%A_%a.out # Standard output

module load Rtidyverse

Rscript \
--vanilla \
3.create_joint_object.R

date
echo "done"