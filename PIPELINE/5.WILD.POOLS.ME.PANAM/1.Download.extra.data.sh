#!/usr/bin/env bash
#
#SBATCH -J download # A single job name for the array
#SBATCH -c 40
#SBATCH -N 1 # on one node
#SBATCH -t 6:00:00 #<= this may depend on your resources
#SBATCH --mem 50G #<= this may depend on your resources
#SBATCH -p bluemoon

### wd=/scratch/aob2x/dest
### grep -E "ES_ba_12|AT_gr_12" /scratch/aob2x/dest/DEST/populationInfo/samps.csv | cut -f1,13 -d',' > /scratch/aob2x/fastq/todl.csv
### run as: sbatch --array=1-$( wc -l < /scratch/aob2x/fastq/todl.csv ) /scratch/aob2x/dest/DEST/mappingPipeline/scripts/download_SRA.sh
### sacct -j 18750292

### MANIE

wget http://berglandlab.uvadcos.io/dest_mapped/DESTv1/ME_bo_09_fall.r2/ME_bo_09_fall.r2.masked.sync.gz

### PANAMA
module load spack/spack-0.18.1
spack load sratoolkit@3.0.0

#SLURM_ARRAY_TASK_ID=1

#sranum=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/fastq/todl.csv | cut -f2 -d',' )
#sampName=$( sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/aob2x/fastq/todl.csv | cut -f1 -d',' )

sranum=SRR3053560
out=/netfiles02/lockwood_lab/IntrogressionProject/Wild_samps_data/Panama

#echo $sampName " / " $sranum

prefetch \
-o $out/${sranum}.sra \
${sranum}

fasterq-dump \
--split-files \
--split-3 \
--outfile $out/${sranum} \
-e 10 \
$out/${sranum}.sra

gzip $out/${sranum}_1.fastq
gzip $out/${sranum}_2.fastq

#rm /scratch/aob2x/fastq/${sranum}.sra
