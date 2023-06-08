### Explore Brent's data
### An R script
### Jcbn June 8, 2023

library(tidyverse)
library(magrittr)
library(data.table)
library(foreach)

### where are the models?
root <- "/gpfs2/scratch/jcnunez/Brent_Introgression/rnpv_analysis"
mod_files <- system(paste("ls", root, sep = " "), intern = T)

### load models

mod_o <- foreach(i=mod_files, .combine = "rbind")%do%{

tmp <- get(load(paste(root,i, sep = "/")))

tmp %>%
group_by(perm_type, chr, pos_mean) %>% 
  summarize(uci = quantile(rnp.binom.p, 0.1))  %>%
  mutate(metric = "rnvp") -> rnvp.tmp

tmp %>%
group_by(perm_type, chr, pos_mean) %>%  
  summarize(uci = quantile(rnp.wZa.p, 0.1))  %>%
  mutate(metric = "wza") -> wza.tmp

rbind(rnvp.tmp, wza.tmp ) %>%
mutate(model = i) -> tmp.ag

return(tmp.ag)

} 

### visualize

mod_o %>%
filter(chr == "2R") %>%
mutate(uci_mod = case_when(metric == "rnvp" ~ log10(uci)*-1,
						   metric == "wza" ~ log10(uci))) ->
mod_o_plot

mod_o_plot %>%
 ggplot(
    aes(
      x=pos_mean/1e6,
      y=uci_mod,
      color=perm_type
    )
  ) +
   geom_vline(xintercept = 15391154/1e6) +
   geom_vline(xintercept = 20276334/1e6) +
  geom_line(alpha = 0.7) +
  ylab("P-values") +
  xlab("Genome Position (Mb)") +
  theme_bw() +
  scale_color_manual(values = c("grey","red","grey","purple") ) +
  facet_grid(model~., scales = "free") ->
  mod.2R.plot
	
ggsave(mod.2R.plot, file = "mod.2R.plot.pdf")

#### visualize FET
intro.fet <- fread("/netfiles02/lockwood_lab/IntrogressionProject/QTL_mapping_Brent/Fisher_test_results_All_crosses.txt")
names(intro.fet)[1:2] = c("chr","pos")

intro.fet %>%
filter(chr == "2R") %>%
ggplot(
    aes(
      x=pos/1e6,
      y=-log10(sk_all_f.test_pval),
    )
) +   geom_line(alpha = 0.7) +
  theme_bw() +
   geom_vline(xintercept = 15391154/1e6) +
   geom_vline(xintercept = 20276334/1e6) ->
  mod.2R.plot.fet
ggsave(mod.2R.plot.fet, file = "mod.2R.plot.fet.pdf")


##3 Joint PLOT FET RNPV
intro.fet %>%
dplyr::select(chr, pos_mean=pos, sk_all_f.test_pval ) %>%
mutate(uci_mod = log10(sk_all_f.test_pval)) -> fet.sk.plot

setDT(mod_o_plot)
setDT(fet.sk.plot)

### ==> 
 ggplot() +
    geom_vline(xintercept = 15391154/1e6) +
    geom_vline(xintercept = 20276334/1e6) +
    geom_line(
    data = mod_o_plot[metric == "rnvp"],
    aes(
      x=pos_mean/1e6,
      y=uci_mod,
      color=perm_type
    )) +
   geom_line(
    data = fet.sk.plot[chr == "2R"],
    aes(
      x=pos_mean/1e6,
      y=uci_mod,
    )) +
  ylab("P-values") +
  xlab("Genome Position (Mb)") +
  theme_bw() +
  xlim(14.5,25) +
  facet_grid(model~., scales = "free") ->
  joint.fet.glm.plot
	
ggsave(joint.fet.glm.plot, file = "joint.fet.glm.plot.png")


