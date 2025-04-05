##### Create portable VCFs for analysis.

### analysis at --> /users/j/c/jcnunez/scratch/fst_brent/slice_VCFs

module load vcftools/0.1.16
module load gcc/13.3.0-xp3epyt \
bcftools/1.19-iq5mwek \
htslib/1.19.1-6ivqauw


######

gvcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz 
name=top2R_1MB_VA
chr=2R
begin=19551633
end=21551633

vcftools \
--gzvcf $gvcf \
--recode \
--recode-INFO-all \
--out $name \
--chr $chr \
--from-bp $begin \
--to-bp $end

bgzip $name.recode.vcf
tabix $name.recode.vcf.gz

####

gvcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Single_Individuals/Dmel_inds_Taylor.wSNPids.vcf.gz 
name=top2R_1MB_World
chr=2R
begin=19551633
end=21551633

vcftools \
--gzvcf $gvcf \
--recode \
--recode-INFO-all \
--out $name \
--chr $chr \
--from-bp $begin \
--to-bp $end

bgzip $name.recode.vcf
tabix $name.recode.vcf.gz
bcftools query -l $name.recode.vcf.gz


####

gvcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz
name=top2R_1MB_DGRP
chr=2R
begin=15003258
end=16503258


vcftools \
--gzvcf $gvcf \
--recode \
--recode-INFO-all \
--out $name \
--chr $chr \
--from-bp $begin \
--to-bp $end

bgzip $name.recode.vcf
tabix $name.recode.vcf.gz
bcftools query -l $name.recode.vcf.gz
