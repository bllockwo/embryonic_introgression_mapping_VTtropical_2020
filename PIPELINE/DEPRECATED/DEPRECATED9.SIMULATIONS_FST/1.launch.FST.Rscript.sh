#!/usr/bin/env bash
#
#SBATCH -J fst.calc # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fst.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/fst.%A_%a.err # Standard error
#SBATCH --array=1-231

module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

Rscript \
--vanilla \
1.run.FST.window.R \
${SLURM_ARRAY_TASK_ID}

date
echo "done"