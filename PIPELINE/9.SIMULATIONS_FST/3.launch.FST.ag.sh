#!/usr/bin/env bash
#
#SBATCH -J fst.ag # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon
#SBATCH -o ./slurmOutput/fst.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/fst.%A_%a.err # Standard error
#SBATCH --array=1-60

module load spack/spack-0.18.1
spack load r@4.2.1 r-sf

met=master.file.fst.txt

s1=$( cat $met  | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $1 }' )
s2=$( cat $met| sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $2 }' )
sn=$( cat $met | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $3 }' )
ch=$( cat $met | sed '1d' | sed "${SLURM_ARRAY_TASK_ID}q;d" | awk '{ print $4 }' )

echo $s1 $s2 $sn $ch

Rscript \
--vanilla \
2.collect.Fsts.R $s1 $s2 $sn $ch

date
echo "done"