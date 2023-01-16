#!/usr/bin/env bash
#
#SBATCH -J panamPool # A single job name for the array
#SBATCH -c 55
#SBATCH -N 1 # on one node
#SBATCH -t 30:00:00 #<= this may depend on your resources
#SBATCH --mem 55G #<= this may depend on your resources
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

  numFlies=50
  
  #sampleid=$( cat $4  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
  #read_file=$( cat $4  | sed '1d' | awk '{ print $9 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

  #echo $sampleid
  #echo $read_file
	
  #SLURM_ARRAY_TASK_ID=6
  #pop=$( cat /netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/sampleId_meta_data.txt | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

	sampleid=PanamaPool

  cp /netfiles02/lockwood_lab/IntrogressionProject/Wild_samps_data/Panama/SRR3053560_1.fastq $1
  cp /netfiles02/lockwood_lab/IntrogressionProject/Wild_samps_data/Panama/SRR3053560_2.fastq $1


###################################
# Part  2. Run Docker             #
###################################

  singularity run \
  --home $1 \
  $1/dest_freeze1_latest.sif \
  $1/SRR3053560_1.fastq \
  $1/SRR3053560_2.fastq  \
  ${sampleid} \
  $3 \
  --cores $SLURM_CPUS_PER_TASK \
  --max-cov 0.95 \
  --min-cov 4 \
  --base-quality-threshold 25 \
  --num-flies ${numFlies} \
  --do_poolsnp \
  --do_snape

#clean up
rm $1/SRR3053560_1.fastq $1/SRR3053560_2.fastq 

echo "done"
date

