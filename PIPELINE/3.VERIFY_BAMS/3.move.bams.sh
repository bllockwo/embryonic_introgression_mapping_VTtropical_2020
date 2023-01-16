#!/usr/bin/env bash
#
#SBATCH -J move_files # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon
#SBATCH --array=1-20


sample_metadat=/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data_Jan5.2023.txt

sample=$( cat $sample_metadat  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

echo $sample

file_in=~/scratch/test_DEST/pipeline_output/$sample

cp -r $file_in /netfiles02/lockwood_lab/IntrogressionProject/Mapped_bams

date
echo "done"
