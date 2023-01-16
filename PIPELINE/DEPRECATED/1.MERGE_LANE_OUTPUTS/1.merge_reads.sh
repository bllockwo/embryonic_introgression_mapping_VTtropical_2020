#!/usr/bin/env bash
#
#SBATCH -J mergeReads # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 12:00:00 #<= this may depend on your resources
#SBATCH --mem 40G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/merge.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/merge.%A_%a.err # Standard error
#SBATCH -p bluemoon
#SBATCH --array=1


## This script will merge reads from different runs of the same individual
## Dec 19, 2022
## Jcb Nunez, PhD - jnunez2@uvm.edu

seq_metadat=/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data.txt
sample_metadat=/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/sampleId_meta_data.txt
reads_loc=/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData

sample=$( cat $sample_metadat  | sed '1d' | awk '{ print $1 }' | sed "${SLURM_ARRAY_TASK_ID}q;d" )

mkdir MERGED_READS_fixed

mkdir joint_$sample
mkdir joint_$sample/processed

### Forward
count=0
for i in $(grep "${sample}" $seq_metadat | awk -v add=$reads_loc '{ print add"/"$3"/"$4"_1.fq.gz" }' | tr '\n' ' ')
do
echo $i
(( count++ ))
echo $count
cp $i joint_$sample
gzip -d joint_$sample/*.fq.gz
ls joint_$sample/*.fq
mv joint_$sample/*.fq joint_$sample/$sample.$count.1.fq
mv joint_$sample/$sample.$count.1.fq joint_$sample/processed
done

cat joint_$sample/processed/*.1.fq > MERGED_READS_fixed/${sample}_1.fq
gzip MERGED_READS_fixed/${sample}_1.fq

### Reverse
count=0
for i in $(grep "${sample}" $seq_metadat | awk -v add=$reads_loc '{ print add"/"$3"/"$4"_2.fq.gz" }' | tr '\n' ' ')
do
echo $i
(( count++ ))
echo $count
cp $i joint_$sample
gzip -d joint_$sample/*.fq.gz
ls joint_$sample/*.fq
mv joint_$sample/*.fq joint_$sample/$sample.$count.2.fq
mv joint_$sample/$sample.$count.2.fq joint_$sample/processed
done

cat joint_$sample/processed/*.2.fq > MERGED_READS_fixed/${sample}_2.fq
gzip MERGED_READS_fixed/${sample}_2.fq



echo "done"
date
