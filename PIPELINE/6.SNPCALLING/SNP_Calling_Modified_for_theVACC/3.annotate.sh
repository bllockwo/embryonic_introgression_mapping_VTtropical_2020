#!/usr/bin/env bash
#
#SBATCH -J dockerMap # A single job name for the array
#SBATCH -c 45
#SBATCH -N 1 # on one node
#SBATCH -t 5:00:00 #<= this may depend on your resources
#SBATCH --mem 45G #<= this may depend on your resources
#SBATCH -o ./slurmOutput/RunDest.%A_%a.out # Standard output
#SBATCH -e ./slurmOutput/RunDest.%A_%a.err # Standard error
#SBATCH -p bluemoon
#!/usr/bin/env bash


#module load htslib bcftools intel/18.0 intelmpi/18.0 parallel R/3.6.3
module load spack/spack-0.18.1
spack load bcftools@1.14
spack load openjdk@11.0.15_10
spack load r@4.2.1

### Load picard
SNPEFF=/users/j/n/jnunez2/software/snpEff/exec/snpeff
PICARD=/users/j/n/jnunez2/software/picard/build/libs/picard-2.27.5-5-gbdfea14-SNAPSHOT-all.jar 
JAVAMEM=49G
PROJNAME=BLockIntro

###
##### params
  popSet="PoolSeq"
  method="PoolSNP"
  maf="001"
  mac=5
  version=test
  script_dir="/users/j/n/jnunez2/scratch/test_DEST/poolSNP/snpCalling/innerscripts"
  wd="/users/j/n/jnunez2/scratch/test_DEST/poolSNP/output"
  syncPath1orig="/users/j/n/jnunez2/scratch/test_DEST/pipeline_output/*/*masked.sync.gz"
  syncPath2orig="/users/j/n/jnunez2/scratch/test_DEST/pipeline_output/*/*masked.sync.gz"
#####

java -Xmx$JAVAMEM \
 -jar $PICARD MergeVcfs \
  I=${wd}/sub_bcf/dest.2L.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.2R.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.3L.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.3R.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.4.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.X.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  I=${wd}/sub_bcf/dest.Y.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  O=${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  D=../holo_dmel_6.12.dict \

#### ANNOTATE  
  $SNPEFF \
  BDGP6.32.105 \
  ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.vcf \
  > ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf

#### FIX HEADER
echo "fix header" #this is now fixed in PoolSNP.py
  sed -i '0,/CHROM/{s/AF,Number=1/AF,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
  sed -i '0,/CHROM/{s/AC,Number=1/AC,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
  sed -i '0,/CHROM/{s/AD,Number=1/AD,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf
  sed -i '0,/CHROM/{s/FREQ,Number=1/FREQ,Number=A/}' ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf

  #bcftools view -h ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/tmp.header

  #bcftools reheader --threads 10 -h ${wd}/tmp.header \
  #-o ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.header.bcf \
  #${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.bcf

####
echo "make GDS"
  Rscript --vanilla \
  /users/j/n/jnunez2/scratch/test_DEST/poolSNP/snpCalling/innerscripts/vcf2gds.R \
  ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf

echo "bgzip & tabix"
  bgzip -c  ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf > ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz
  tabix -p vcf ${wd}/$PROJNAME.${popSet}.${method}.${maf}.${mac}.${version}.ann.vcf.gz
