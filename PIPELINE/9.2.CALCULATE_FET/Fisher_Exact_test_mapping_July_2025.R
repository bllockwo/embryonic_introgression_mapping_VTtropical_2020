library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(QTLseqr)
library(ggplot2)
library(cowplot)

# Load snp tables
dat_t <- read.delim(file = "allele_frequencies.txt", sep = "\t", header = TRUE, row.names = 1)

ad_t <- read.delim(file = "allele_depths.txt", sep = "\t", header = TRUE, row.names = 1)

dp_t <- read.delim(file = "total_depths.txt", sep = "\t", header = TRUE, row.names = 1)

#Make data frame

sk_snps <- data.frame(dat_t$SK_0_1, dat_t$VT8_0_1, dat_t$SKF_1_1, dat_t$SKF2_2_1, dat_t$SKF3_3_1, dat_t$VT8F_1_1, dat_t$VT8F2_2_1, dat_t$VT8F3_3_1, ad_t$SK_0_1, ad_t$VT8_0_1, ad_t$SKF_1_1, ad_t$SKF2_2_1, ad_t$SKF3_3_1, ad_t$VT8F_1_1, ad_t$VT8F2_2_1, ad_t$VT8F3_3_1, dp_t$SK_0_1, dp_t$VT8_0_1, dp_t$SKF_1_1, dp_t$SKF2_2_1, dp_t$SKF3_3_1, dp_t$VT8F_1_1, dp_t$VT8F2_2_1, dp_t$VT8F3_3_1, row.names = row.names(dat_t))
names(sk_snps) <- c("sk_af", "vt8_af", "skf1_af", "skf2_af", "skf3_af", "vt8f1_af", "vt8f2_af", "vt8f3_af", "sk_ad", "vt8_ad", "skf1_ad", "skf2_ad", "skf3_ad", "vt8f1_ad", "vt8f2_ad", "vt8f3_ad", "sk_dp", "vt8_dp", "skf1_dp", "skf2_dp", "skf3_dp", "vt8f1_dp", "vt8f2_dp", "vt8f3_dp")
row_names <- row.names(sk_snps)
row_names <- str_split_fixed(row_names, "_", 3)
CHROM <- row_names[,1]
POS <- as.numeric(row_names[,2])
sk_snps <- data.frame(CHROM, POS, sk_snps)

# SKF1

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$skf1_ad[i], sk_snps$skf1_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

Fisher_test_results <- data.frame(sk_snps$CHROM, sk_snps$POS)
colnames(Fisher_test_results)[1:2] <- c("CHROM","POS")

skf1_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, skf1_f.test_pval)

# SKF2

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$skf2_ad[i], sk_snps$skf2_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

skf2_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, skf2_f.test_pval)

# SKF3

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$skf3_ad[i], sk_snps$skf3_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

skf3_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, skf3_f.test_pval)

# VT8F1

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$vt8f1_ad[i], sk_snps$vt8f1_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

vt8f1_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, vt8f1_f.test_pval)

# VT8F2

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$vt8f2_ad[i], sk_snps$vt8f2_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

vt8f2_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, vt8f2_f.test_pval)

# VT8F3

# Run a for loop to calculate Fisher's exact tests across the genome. 
# Indicate the data frame and the crosses you want to analyze.

f.test_p_val <- vector("numeric", nrow(sk_snps))

for (i in 1:nrow(sk_snps)){
  allele_table <- c(sk_snps$vt8f3_ad[i], sk_snps$vt8f3_dp[i], sk_snps$vt8_ad[i], sk_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

vt8f3_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, vt8f3_f.test_pval)

write.table(Fisher_test_results, file="Fisher_test_results_SK_crosses.txt", sep = "\t", row.names = FALSE)

# Combine replicates to look for large-effect loci.
# Do this by adding allele counts and total depths across replicates.
# Do this separately for SK and CH crosses.

sk_all_ad <- sk_snps$skf1_ad + sk_snps$skf2_ad + sk_snps$skf3_ad + sk_snps$vt8f1_ad + sk_snps$vt8f2_ad + sk_snps$vt8f3_ad
  
sk_all_dp <- sk_snps$skf1_dp + sk_snps$skf2_dp + sk_snps$skf3_dp + sk_snps$vt8f1_dp + sk_snps$vt8f2_dp + sk_snps$vt8f3_dp
  
added_snps <- data.frame(sk_snps$CHROM, sk_snps$POS, sk_all_ad, sk_all_dp, sk_snps$vt8_ad, sk_snps$vt8_dp)
names(added_snps) <- c("CHROM", "POS", "sk_all_ad", "sk_all_dp", "vt8_ad", "vt8_dp")

# Fisher's exact test
# SK altogether

f.test_p_val <- vector("numeric", nrow(added_snps))

for (i in 1:nrow(added_snps)){
  allele_table <- c(added_snps$sk_all_ad[i], added_snps$sk_all_dp[i], added_snps$vt8_ad[i], added_snps$vt8_dp[i])
  dim(allele_table) <- c(2,2)
  f.test <- fisher.test(allele_table)
  f.test_p_val[i] <- f.test$p.value
  if (i%%1000 == 0){print(i)}
}

# Add results to data frame of Fisher's Exact test results 

sk_all_f.test_pval <- f.test_p_val
Fisher_test_results <- cbind(Fisher_test_results, sk_all_f.test_pval)

# Now filter out low coverage SNPs
# SK crosses
sk_all_fisher <- data.frame(Fisher_test_results$CHROM, Fisher_test_results$POS, Fisher_test_results$sk_all_f.test_pval, added_snps$sk_all_dp, added_snps$vt8_dp)
names(sk_all_fisher) <- c("CHROM", "POS", "f.test_pval", "sk_all_dp", "vt8_dp")

sk_all_fisher <- sk_all_fisher %>% 
  filter((sk_all_dp + vt8_dp) >= (mean(c(sk_all_dp + vt8_dp))*(210-1))/(mean(c(sk_all_dp + vt8_dp)) + 210))

# The number of low-coverage SNPs that were filtered out:
nrow(added_snps)-nrow(sk_all_fisher)
# [1] 1199

# Negative log p-value calculation

sk_all_fisher <- sk_all_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))

# Number of SNPs that are significant
# Bonferroni
sk_all_fisher %>% filter(negLog10Pval_ftest > -log10(0.05/nrow(sk_all_fisher))) %>% nrow()
# [1] 677

# B-H FDR
sk_all_fisher %>% filter(negLog10Pval_ftest > -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01))) %>% nrow()
# [1] 23044

# Make new columns based on significance threshold (TRUE/FALSE)
# Bonferroni
sk_all_fisher <- sk_all_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(sk_all_fisher)))

# B-H FDR
sk_all_fisher <- sk_all_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)))

# Manhattan Plots
# SK
p_sk_all <- ggplot(sk_all_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_sk_all, 
       file = "SK_all_FET_Bonferroni_030823.png", 
       w = 10, h = 6)

#### Make separate Manhattan plots of Fisher's Exact Test (FET) results for each replicate cross.

# First, find max values of -log10 pvalue in order to set y-axis limits.
max(skf1_fisher$negLog10Pval_ftest, skf2_fisher$negLog10Pval_ftest, skf3_fisher$negLog10Pval_ftest,
    vt8f1_fisher$negLog10Pval_ftest, vt8f2_fisher$negLog10Pval_ftest, vt8f3_fisher$negLog10Pval_ftest,
    chf1_fisher$negLog10Pval_ftest, chf2_fisher$negLog10Pval_ftest, chf3_fisher$negLog10Pval_ftest,
    vt10f1_fisher$negLog10Pval_ftest, vt10f2_fisher$negLog10Pval_ftest, vt10f3_fisher$negLog10Pval_ftest)
# [1] 14.51966

# SKF1
skf1_fisher_filtered <- skf1_fisher_filtered %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
skf1_fisher_filtered <- skf1_fisher_filtered %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(skf1_fisher_filtered)))

p_skf1 <- ggplot(skf1_fisher_filtered, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(skf1_fisher_filtered)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = skf1_fisher_filtered$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_skf1, 
       file = "SKF1_FET_Bonferroni_040623.png", 
       w = 10, h = 6)

# SKF2
skf2_fisher <- skf2_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
skf2_fisher <- skf2_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(skf2_fisher)))

p_skf2 <- ggplot(skf2_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_skf2, 
       file = "SKF2_FET_Bonferroni_030823.png", 
       w = 10, h = 6)

# SKF3
skf3_fisher <- skf3_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
skf3_fisher <- skf3_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(skf3_fisher)))

p_skf3 <- ggplot(skf3_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_skf3, 
       file = "SKF3_FET_Bonferroni_030823.png", 
       w = 10, h = 6)

# VT8F1
vt8f1_fisher <- vt8f1_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
vt8f1_fisher <- vt8f1_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(vt8f1_fisher)))

p_vt8f1 <- ggplot(vt8f1_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_vt8f1, 
       file = "VT8F1_FET_Bonferroni_030823.png", 
       w = 10, h = 6)

# VT8F2
vt8f2_fisher <- vt8f2_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
vt8f2_fisher <- vt8f2_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(vt8f2_fisher)))

p_vt8f2 <- ggplot(vt8f2_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_vt8f2, 
       file = "VT8F2_FET_Bonferroni_030823.png", 
       w = 10, h = 6)

# VT8F3
vt8f3_fisher <- vt8f3_fisher %>% mutate(negLog10Pval_ftest = -log10(f.test_pval))
vt8f3_fisher <- vt8f3_fisher %>% mutate(threshold = negLog10Pval_ftest > -log10(0.05/nrow(vt8f3_fisher)))

p_vt8f3 <- ggplot(vt8f3_fisher, aes(x=as.numeric(POS)/1e6, y=negLog10Pval_ftest)) +
  facet_grid(.~CHROM, scales = "free_x", space = "free_x") +
  geom_point(aes(color=threshold), size = 2, alpha=0.4) +
  scale_color_manual(values = c("#0072B2", "#D55E00")) +
  #geom_line(aes(y=skf1_trop_ftest_smooth), color="black") +
  geom_hline(yintercept = -log10(0.05/nrow(sk_all_fisher)), color = "#D55E00") +
  #geom_hline(yintercept = -log10(QTLseqr::getFDRThreshold(pvalues = sk_all_fisher$f.test_pval, alpha = 0.01)), color = "#D55E00") +
  ylim(0,15) +
  scale_x_continuous(breaks = c(0,5,10,15,20,25,30)) +
  labs(x="Genomic position (Mb)", y=expression(-log[10](italic(p)))) +
  theme_half_open(font_size = 20) +
  panel_border() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(size = 14))

ggsave(p_vt8f3, 
       file = "VT8F3_FET_Bonferroni_030823.png", 
       w = 10, h = 6)
