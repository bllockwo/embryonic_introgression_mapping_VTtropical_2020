#!/usr/bin/env bash
#
#SBATCH -J rnpv_mods # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon
#SBATCH --array=1-7

######
module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

Rscript \
--vanilla \
1.rnvp.test.models.R \
${SLURM_ARRAY_TASK_ID}

date
echo "done"
