#!/usr/bin/env bash
#
#SBATCH -J collect # A single job name for the array
#SBATCH -c 1
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 80G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/collect.%A_%a.out # Standard output

module load Rtidyverse

Rscript \
--vanilla \
collect_Nunez24_mod.R \
$1 $2 $3

##Example

##LAUNCH_collect_Nunez24_mod.sh 8 "temp.max" "5.Cville"

#mod.i= as.numeric(args[1]) --> 8
#variable.i=args[2] --> temp.max
#cluster.i=args[3] --> 5.Cville


