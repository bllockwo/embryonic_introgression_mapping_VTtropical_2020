module load vcftools
plink=/netfiles/nunezlab/Shared_Resources/Software/plink_linux_x86_64_20240818/plink

DGRP=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz
#2R:16439138
#2R_16439138_SNP


vcftools \
--gzvcf $DGRP \
--chr 2R \
--remove-indels \
--recode --recode-INFO-all \
--out Embryo.snp.dgrp

#--from-bp 15439138 \
#--to-bp 17439138 \


$plink --vcf Embryo.snp.dgrp.recode.vcf \
--aec \
--r2 \
--ld-snp 2R_16439138_SNP \
--ld-window-kb 1000 \
--ld-window 99999  \
--ld-window-r2 0 \


vcftools --vcf Embryo.snp.dgrp.recode.vcf \
--thin 10 \
--geno-r2 \
--ld-window 10000 \
--out ${eF}_file

vcftools --vcf Embryo.snp.dgrp.recode.vcf \
--window-pi 20000 \
--window-pi-step 10000 \
--out Embryo.snp.dgrp

vcftools --vcf Embryo.snp.dgrp.recode.vcf \
--TajimaD  20000 \
--out Embryo.snp.dgrp


######
Cville=/netfiles/nunezlab/Drosophila_resources/Datasets/2023.Nunez_et_al_Supergene_paper/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

vcftools \
--gzvcf $Cville \
--chr 2L \
--from-bp 9582788 \
--to-bp 9582788 \
--remove-indels \
--maf 0.01 \
--recode --recode-INFO-all \
--out ir94e
