library(data.table)
library(tidyverse)

ld <- fread("Pi_X_LD.geno.ld")

ld %>%
  filter(POS1 == 15637404 | POS2 == 15637404) %>%
  mutate(pos1_corr = ifelse(POS1 == 15637404, POS1, POS2)) %>%
  mutate(pos2_corr = ifelse(POS1 == 15637404, POS2, POS1)) %>%
  ggplot(aes(
    x=pos2_corr,
    y=`R^2`
  )) + geom_point() ->
  ld_x
ggsave(ld_x, file = "ld_x.pdf")