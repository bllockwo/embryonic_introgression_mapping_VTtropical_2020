#!/usr/bin/env bash
#
#SBATCH -J estPI # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon
#SBATCH --array=1-20


sample_metadat=/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data_Jan5.2023.txt

### Load mdoules
module load spack/spack-0.18.1
spack load samtools

#declare npstats
npstat=/users/j/n/jnunez2/software/npstat/npstat

#extract samples
sample=$( cat $sample_metadat  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )
echo $sample

#decalre the bam file
bam=/netfiles02/lockwood_lab/IntrogressionProject/Mapped_bams/$sample/$sample.mel.bam

### index the bam file
#samtools \
#index \
#$bam

### Prepare loop

for chr in 2L 2R 3L 3R X
do
echo $chr
##exctract slice
samtools mpileup -r $chr $bam > $sample.$chr.pileup
##run npstat
$npstat -n 30 \
-l 1000000 \
-nolowfreq 3 \
-mincov 10 \
$sample.$chr.pileup
##remove slice
rm $sample.$chr.pileup
done

echo "done"
date
