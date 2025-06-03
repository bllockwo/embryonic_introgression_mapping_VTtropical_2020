library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

FBgn0034518_summary <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/16.GENE_TRANSCRIPTS/Expression_Kylie/SP70_expression_summary.csv")


FBgn0034518_summary$genome <- factor(FBgn0034518_summary$genome, levels = c("SK", "VT8", "introgression"))
ggplot(FBgn0034518_summary, aes(x=treatment, y=count, color=genome, fill = genome, shape = genome)) +
  geom_point(size=4, alpha=0.75) +
  geom_errorbar(aes(ymin=count-se, ymax=count+se), width=.15, alpha=0.75, size = 1) +
  geom_line(aes(color=genome, group=genome, linetype = genome), alpha=0.75, size = 1) +
  scale_color_manual(values = c("SK" = "firebrick3", "VT8" = "steelblue", "introgression" = "blueviolet"),
                     labels = c("SK", "VT", "introgression")) +
  scale_fill_manual(values = c("SK" = "firebrick3", "VT8" = "steelblue", "introgression" = "blueviolet"),
                    labels = c("SK", "VT", "introgression")) +
  scale_shape_manual(values = c("SK" = 17, "VT8" = 15, "introgression" = 16),
                     labels = c("SK", "VT", "introgression")) +
  scale_linetype_manual(values = c("SK" = "solid", "VT8" = "dashed", "introgression" = "twodash"),
                        labels = c("SK", "VT", "introgression")) +
  scale_x_discrete(labels = c("25°C", "34°C")) +
  ylim(400,1200) +
  ylab("Normalized counts") +
  labs(title = "SP70") +
  theme_classic()
