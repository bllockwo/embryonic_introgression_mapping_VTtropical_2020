#### GOwinda analysis

gowinda=/netfiles/nunezlab/Shared_Resources/Software/GOwinda/Gowinda-1.12.jar

all_SNPS=/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_setsall_SNP.txt
all_outs=/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_GENOMEWIDE_outliers_SNP.txt

assocGO=/netfiles/nunezlab/D_melanogaster_resources/Annotations/funcassociate_go_associations_FBGNstyle.txt
DMELgtf=/netfiles/nunezlab/D_melanogaster_resources/Annotations/dmel-all-r6.12.gtf

####
java -jar $gowinda \
--snp-file $all_SNPS \
--candidate-snp-file $all_outs \
--gene-set-file $assocGO \
--annotation-file $DMELgtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition gene \
--threads 4 \
--output-file ALL_outliers_gowinda_results.txt \
--mode gene \
--min-genes 5

####
out2R=/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_2R_window_SNP.txt


java -jar $gowinda \
--snp-file $all_SNPS \
--candidate-snp-file $out2R \
--gene-set-file $assocGO \
--annotation-file $DMELgtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition gene \
--threads 4 \
--output-file chr2R_win_outliers_gowinda_results.txt \
--mode gene \
--min-genes 5

####
outX=/netfiles02/lockwood_lab/IntrogressionProject/GOWINDA_sets/GOWINDA_X_window_SNP.txt


java -jar $gowinda \
--snp-file $all_SNPS \
--candidate-snp-file $outX \
--gene-set-file $assocGO \
--annotation-file $DMELgtf \
--simulations  100000 \
--min-significance 1 \
--gene-definition gene \
--threads 4 \
--output-file chrX_win_outliers_gowinda_results.txt \
--mode gene \
--min-genes 5

