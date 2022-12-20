#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 50
#SBATCH -N 1 # on one node
#SBATCH -t 24:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon

pwd
echo $SLURM_CPUS_PER_TASK 

### modules
  spack load singularity

###################################
# Part  1. Get Sample information #
###################################
  #SLURM_ARRAY_TASK_ID=1

  pop=$( cat $4  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  numFlies=30

  echo $pop

 cp $2/${pop}_1.fq.gz $1
 cp $2/${pop}_2.fq.gz $1

###################################
# Part  2. Run Docker             #
###################################

  singularity run \
  --home $1 \
  $1/dest_freeze1_latest.sif \
  $1/${pop}_1.fq.gz \
  $1/${pop}_2.fq.gz \
  ${pop} \
  $3 \
  --cores $SLURM_CPUS_PER_TASK \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --num-flies ${numFlies} \
  --do_poolsnp \
  --do_snape

#clean up
rm $1/${pop}_1.fq.gz $1/${pop}_2.fq.gz 

echo "done"
date

