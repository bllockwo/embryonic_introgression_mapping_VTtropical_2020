module load vcftools
module load plink


#########
va_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/2024.Nunez_et_al_Genetics/Single_Individuals/CM_pops.AllChrs.Whatshap.shapeit.annot.wSNPids.vcf.gz
dgrp_vcf=/netfiles/nunezlab/D_melanogaster_resources/Datasets/DGRP2/dgrp2.vcf.gz
#zcat $va_vcf | cut -f-8 | grep "2R_20551633_SNP"
#2R: 20,551,633 and X: 15,743,371

plink --ld "2R_20551633_SNP" "X_15602941_SNP" --vcf $va_vcf --allow-extra-chr --double-id

plink --ld "2R_16439138_SNP" "X_15496974_SNP" --vcf $dgrp_vcf --allow-extra-chr --double-id

### Other interesting mutations to explore ... 
#### 2R:20,551,633 and X:15,602,941
plink --ld "X_15602941_SNP" "X_15648874_SNP" --vcf $va_vcf --allow-extra-chr --double-id
plink --ld "X_15496974_SNP" "X_15542907_SNP" --vcf $dgrp_vcf --allow-extra-chr --double-id
