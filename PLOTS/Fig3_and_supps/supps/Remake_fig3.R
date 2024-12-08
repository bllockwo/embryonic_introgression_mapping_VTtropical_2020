### Plot #3
### Replicable plot
library(tidyverse)
library(reshape2)
library(magrittr)

##### plot
real_enr = get(load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/9.SIMULATIONS_FST/plots/FST_enrichment/fst.binomial.test.VT8.Rdata"))
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/9.SIMULATIONS_FST/plots/FST_enrichment/enrich.sim.data.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/9.SIMULATIONS_FST/plots/FST_window/fst.wins.mean.simsIntro.Rdata")
load("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/embryonic_introgression_mapping_VTtropical_2020/PIPELINE/9.SIMULATIONS_FST/plots/FST_window/real_fst.mean.Rdata")


fst.wins %>%
  group_by(chr,winSTA,winEND) %>%
  summarize(
    sdfst = sd(fstm, na.rm = T),
    fst95 = quantile(fstm, 0.95),
    fst05 = quantile(fstm, 0.05),
    fstm = mean(fstm, na.rm = T),
  ) ->
  sim.summaries

real_fst$fstm[real_fst$fstm < 0] = 0

ggplot()+
  geom_line(data=real_fst,
            aes(
              x=( winSTA   + winEND)/2e6,
              y=fstm,
              color = samp2
            ) 
  ) +
  geom_point(data=real_fst,
             aes(
               x=( winSTA   + winEND)/2e6,
               y=fstm,
               color = samp2
             ) 
  ) +
  geom_ribbon(data=sim.summaries,
              alpha = 0.8,linewidth = 1.2,
              aes(
                x=( winSTA   + winEND)/2e6,
                y=fstm,
                ymin=0,
                ymax=fstm+2.576*(sdfst/sqrt(100)),
              )
  )  +
  facet_grid(samp2~chr, scales = "free_x") +
  scale_color_brewer(palette = "Set3") +
  theme_minimal() + theme(legend.position = "bottom")



#####
sim_data %>%
  group_by(winSTA,  winEND, chr) %>%
  summarise(sim95 = quantile(p.binom.adj, 0.95)) ->
  sim.95.envelope

ggplot() +
  geom_line(
    data=real_enr,
    aes(
      x=( winSTA   + winEND)/2e6,
      y=-log10(p.binom.adj),
      color=samp2
    )) +
  geom_line(
    data=sim.95.envelope,
    aes(
      x=( winSTA   + winEND)/2e6,
      y=-log10(sim95)
    ), color = "grey", linewidth = 1.9) + 
  scale_color_brewer(palette = "Set3") +
  theme_minimal() +  
  facet_grid(samp2~chr, scales = "free_x")

#### zoomed in plots
real_enr %>%
  filter(p.binom.adj < 0.001) %>%
  dcast(winSTA+winEND+chr~samp2, value.var = "p.binom.adj",
        fun.aggregate = mean) %>%
  .[complete.cases(.),] %>%
  mutate(o.mean = (SKF_1_1+SKF2_2_1+SKF3_3_1+VT8F_1_1+VT8F2_2_1+VT8F3_3_1)/6 ) %>%
  arrange(o.mean)

#
#1  19115753 20115753  2R  
#2  19615753 20615753  2R  
#3  15123410 16123410   X
