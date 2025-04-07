#### Gene structure and expression
library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)

####
dev <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/16.GENE_TRANSCRIPTS/developmental_expression.csv")
locs <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/16.GENE_TRANSCRIPTS/gene_locations.txt")
names(locs) = c("chr","source","feature", "start", "end", "qual", "strand",
                "shift","annotation", "gene_name"
                ) #10 cols

###
locs %>%
  filter(feature == "exon") %>%
  ggplot(aes(
    x=gene_name,
    y=start,
    yend=end
  )) + 
  geom_segment(linewidth = 2, color = "grey50") +
  geom_hline(yintercept = 20551633, linetype = "dashed") +
  coord_flip() + theme_classic()

locs %>%
  filter(feature %in% c("CDS","exon") ) %>%
  ggplot(aes(
    x=gene_name,
    y=start,
    yend=end,
    color=feature
  )) + 
  geom_segment(linewidth = 2, alpha = 0.5) +
  geom_hline(yintercept = 20551633, linetype = "dashed") +
  coord_flip() + theme_classic()


####
dev %>%
  melt(id = c("Lifestage", "Life_State",   "Short_Name")) %>%
  ggplot(aes(
    x=(Life_State),
    y=value
  )) + geom_bar(stat = "identity")+ 
  facet_grid(variable~.) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
