### Mine the duplication rates
### 

library(tidyverse)
library(vroom)
library(foreach)
library(magrittr)


qcs <-vroom("/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data_Jan5.2023.txt", col_names = T)

dup.rat.head <- read.table( pipe( paste("sed -n 7p ", "~/scratch/test_DEST/pipeline_output/" , 
                                   paste(qcs$Sample_id[1], "/", sep = ""),
                                   paste(qcs$Sample_id[1], ".mark_duplicates_report.txt", sep = ""), 
                                   sep = "")))


dup.rate <- foreach(i=1:20, .combine = "rbind")%do%{
  
  message(i)
  
dup.rat <- read.table( pipe( 
  paste("sed -n 8p ", "~/scratch/test_DEST/pipeline_output/" , 
                    paste(qcs$Sample_id[i], "/", sep = ""),
              paste(qcs$Sample_id[i], ".mark_duplicates_report.txt", sep = ""), 
              sep = "")
  )
)


dup.rat %<>% as.data.frame() %>% mutate(samp = qcs$Sample_id[i])
return(dup.rat)
}

names(dup.rate)[1:10] = c("LIBRARY",
                          "UNPAIRED_READS_EXAMINED",
                          "READ_PAIRS_EXAMINED",
                          "SECONDARY_OR_SUPPLEMENTARY_RDS",
                          "UNMAPPED_READS",
                          "UNPAIRED_READ_DUPLICATES",
                          "READ_PAIR_DUPLICATES",
                          "READ_PAIR_OPTICAL_DUPLICATES",
                          "PERCENT_DUPLICATION",
                          "ESTIMATED_LIBRARY_SIZE"
                          )

#names(qcs) = c("samp","path", "batch")

pass.samps =
  c("CH_0_1",
    #"CH_0_2",
    "CHF_1_1",
    "CHF2_2_1",
    #"CHF3_3_1",
    "CHF3_3_2",
    "SK_0_1",
    "SKF_1_1",
    "SKF2_2_1",
    "SKF3_3_1",
    "VT10_0_1",
    #"VT10_0_2",
    "VT10F_1_1",
    #"VT10F2_2_1",
    "VT10F2_2_2",
    "VT10F3_3_1",
    "VT8_0_1",
    "VT8F_1_1",
    "VT8F2_2_1",
    "VT8F3_3_1")

dup.rate %>%
  filter(samp %in% pass.samps) -> table.seq.inf

table.seq.inf$PERCENT_DUPLICATION %>% mean



#####

dup.rate %>%
  left_join(qcs) -> data.pcr.dups

data.pcr.dups %>%
  ggplot(aes(x=batch, y=PERCENT_DUPLICATION, fill = batch)) + 
  geom_boxplot() ->
  PERCENT_DUPLICATION
ggsave(PERCENT_DUPLICATION, file = "PERCENT_DUPLICATION.png", w = 5, h =4)

t.test(PERCENT_DUPLICATION~batch, data=data.pcr.dups)

