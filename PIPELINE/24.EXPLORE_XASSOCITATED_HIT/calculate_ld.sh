module load vcftools

input_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz

vcftools \
--gzvcf $input_vcf \
--geno-r2 \
--ld-window 50000 \
--chr X \
--from-bp 15587404 \
--to-bp 15687404 \
--out Pi_X_LD

module load vcftools

input_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz

vcftools \
--gzvcf $input_vcf \
--geno-r2 \
--ld-window 50000 \
--chr X \
--from-bp 15587404 \
--to-bp 15687404 \
--out Pi_X_LD_VA

#########
va_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz
#zcat $va_vcf | cut -f-8 | grep "2R_20551633_SNP"
#2R: 20,551,633 and X: 15,743,371

module load plink
plink --ld "2R_20551633_SNP" "X_15743371_SNP" --vcf $va_vcf --allow-extra-chr --double-id

plink --ld "2R_16439138_SNP" "X_15637404_SNP" --vcf $dgrp_vcf --allow-extra-chr --double-id


