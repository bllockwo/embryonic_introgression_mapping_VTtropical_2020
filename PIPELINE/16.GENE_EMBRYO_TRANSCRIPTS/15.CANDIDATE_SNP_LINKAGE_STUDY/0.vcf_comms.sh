module load vcftools
plink=/netfiles/nunezlab/Shared_Resources/Software/plink_linux_x86_64_20240818/plink

DGRP=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz
#2R:16439138
#2R_16439138_SNP


vcftools \
--gzvcf $DGRP \
--chr 2R \
--remove-indels \
--max-missing 0.1 \
--recode --recode-INFO-all \
--out Embryo.snp.dgrp.SFS

vcftools \
--gzvcf $DGRP \
--chr 2R \
--remove-indels \
--max-missing 0.8 \
--maf 0.01 \
--recode --recode-INFO-all \
--out Embryo.snp.dgrp.maf

#--from-bp 15439138 \
#--to-bp 17439138 \


$plink --vcf Embryo.snp.dgrp.maf.recode.vcf \
--aec \
--r2 \
--ld-snp 2R_16439138_SNP \
--ld-window-r2 0 \
--ld-window-kb 999999999 \
--ld-window 999999999  \

####
vcftools --vcf Embryo.snp.dgrp.SFS.recode.vcf \
--fst-window-size 20000 \
--fst-window-step 10000 \
--weir-fst-pop snp0.txt \
--weir-fst-pop snp2.txt \
--out Embryo.FST.snp.dgrp


###

vcftools --vcf Embryo.snp.dgrp.SFS.recode.vcf \
--window-pi 20000 \
--window-pi-step 10000 \
--keep snp0.txt \
--out Embryo0.snp.dgrp

vcftools --vcf Embryo.snp.dgrp.SFS.recode.vcf \
--window-pi 20000 \
--window-pi-step 10000 \
--keep snp2.txt \
--out Embryo2.snp.dgrp

vcftools --vcf Embryo.snp.dgrp.SFS.recode.vcf \
--TajimaD  20000 \
--keep snp0.txt \
--out Embryo0.snp.dgrp

vcftools --vcf Embryo.snp.dgrp.SFS.recode.vcf \
--TajimaD  20000 \
--keep snp2.txt \
--out Embryo2.snp.dgrp


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
