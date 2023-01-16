### Pretty report Qulimap
### 

library(tidyverse)
library(magrittr)
library(vroom)
library(foreach)
library(patchwork)
library(reshape2)
## load dat
## 

#### ---> this is a file that contains the information of the files to fastQC, including the address.
qcsamps <-vroom("/netfiles02/lockwood_lab/IntrogressionProject/IntrogressionRawData/mapping_meta_data_Jan5.2023.txt", col_names = T)

#system(paste("ls", paste(qcsamps$X2[1], "raw_data_qualimapReport", sep = "/")),
#       intern = T) -> items.to.load

items.to.load <- c("coverage_across_reference.txt", "insert_size_across_reference.txt", "mapping_quality_across_reference.txt", "mapped_reads_gc-content_distribution.txt")

o.all <- foreach(k=1:length(items.to.load), .combine = "rbind")%do%{
tmp.o <- foreach(j=1:dim(qcsamps)[1], .combine = "rbind")%do%{
  
tmp <- vroom( paste( paste("QUAL_", qcsamps$Sample_id[j], sep = "") , "raw_data_qualimapReport", 
                               items.to.load[k], sep = "/")) %>%
  mutate(samp = qcsamps$Sample_id[j],
         #group = qcsamps$X3[j],
         metric = items.to.load[k])

names(tmp)[1:2] = c("pos","var")

return(tmp)

}

maxt = ifelse(items.to.load[k] %in% c("insert_size_histogram.txt",
                                      "duplication_rate_histogram.txt",
                                      "coverage_across_reference.txt"
                                      )
              , 1e9, 200)

#tmp.o %>%
#  filter(var != 0) %>%
#  filter(var < maxt) %>%
#  ggplot(aes(x=group, y=(var), fill=group)) + 
#  ylab(paste(items.to.load[k])) + 
#  geom_boxplot() +
#  ggtitle("T-test: BAD vs REG", subtitle = t.test(var~ group,data = tmp.o)$p.value)-> tmp.box
  
#tmp.o %>%
#  filter(var > 0) %>%
#  filter(var < maxt) %>%
#  ggplot(aes(x=pos, y=(var), color=group)) + 
#  ylab(paste(items.to.load[k])) + 
# geom_point() -> tmp.line

#ggsave(tmp.box+tmp.line, file = paste(items.to.load[k],"png", sep ="."), h= 4, w=8)

return(tmp.o)
}

23513712+25286936+28110227+32079331+1348131+23542271 -> dmel.genome

l2L = 0
r2L = 23513712

l2R = 23513712
r2R = 23513712+25286936

l3L = 23513712+25286936
r3L = 23513712+25286936+28110227

l3R = 23513712+25286936+28110227
r3R = 23513712+25286936+28110227+32079331

l4 = 23513712+25286936+28110227+32079331
r4 = 23513712+25286936+28110227+32079331+1348131

lX = 23513712+25286936+28110227+32079331+1348131
rX = 23513712+25286936+28110227+32079331+1348131+23542271

o.all$metric = gsub(".txt", "", o.all$metric)

o.all %>% 
  filter(metric %in% c("coverage_across_reference", 
                       "insert_size_across_reference", 
                       "mapping_quality_across_reference")) %>%
  mutate(bam_slice = case_when(pos < dmel.genome ~ "Dmel",
                               pos >= l2L & pos < r2L ~ "2L",
                               pos >= l2R & pos < r2R ~ "2R",
                               pos >= l3L & pos < r3R ~ "3L",
                               pos >= l3R & pos < r3R ~ "3R",
                               pos >= l4 & pos < r4 ~ "4",
                               pos >= lX & pos < rX ~ "X",
                           pos >= dmel.genome ~ "Holo")) %>% 
  group_by(samp, bam_slice, metric)  %>%
  summarise(m = mean(var),
            sd = sd(var)) ->
  o.ag.summ

save(o.ag.summ, file = "qc.stats.pools.Rdata")


#
# CAN START HERE from saved object! 
#

load("qc.stats.pools.Rdata")

##select samples. 
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

o.ag.summ %>%
  filter(samp %in% pass.samps) %>%
  filter(bam_slice == "Dmel") %>%
  group_by(metric) %>%
  summarise(m = mean(m))
###
o.ag.summ %>%
  filter(samp %in% pass.samps) %>%
  filter(bam_slice == "Dmel") %>%
  dcast(samp~metric) -> dcast.tmp

cor.test(dcast.tmp$coverage_across_reference, dcast.tmp$insert_size_across_reference)
cor.test(dcast.tmp$mapping_quality_across_reference, dcast.tmp$insert_size_across_reference)

###
foreach(i = c("coverage_across_reference", 
              "insert_size_across_reference", 
              "mapping_quality_across_reference"))%do%{
                
o.ag.summ %>%
  filter(samp %in% pass.samps) %>%
  separate(samp, into = c("pop", "gen", "Rep"), sep = "_") %>%
  filter(metric == i) %>%
  filter(bam_slice == "Dmel") -> tmp.dat
  
  if(i != "insert_size_across_reference"){
    lim.l = mean(tmp.dat$m)-(3.0*mean(tmp.dat$sd))
    lim.r = mean(tmp.dat$m)+(3.0*mean(tmp.dat$sd))
  } else if(i == "insert_size_across_reference"){
    lim.l = mean(tmp.dat$m)-(25*mean(tmp.dat$sd))
    lim.r = mean(tmp.dat$m)+(25*mean(tmp.dat$sd))
  }

  
  tmp.dat %>%             
  ggplot(aes(
    x=pop,
    y=m,
    ymin = m - sd,
    ymax = m + sd,
    #color=Rep,
  )) + 
    coord_flip() +
  geom_hline(yintercept = mean(tmp.dat$m), linetype = "dashed") +
  geom_errorbar( width = 0.3, position=position_dodge(width=0.5)) +
  geom_point(shape = 21, fill = "white" , size = 1.8, position=position_dodge(width=0.5)) +
  ggtitle(gsub("_"," ",i)) +
  theme_bw() + 
  xlab("Pool") +
  ylab("QC metric (Mean +/- SD)") +
  ylim(  lim.l, lim.r) +
    theme(plot.title = element_text(size = 12)) +
  facet_grid(gen~., scales = "free_y") ->
  metrics.mapping
  
assign(i, metrics.mapping)
  
ggsave(metrics.mapping, file = paste(i, ".pdf", sep = ""), w= 4, h = 4)

}

ggsave(coverage_across_reference+mapping_quality_across_reference+insert_size_across_reference, file = "Figure.QC.pdf",
       h =3.5, w = 9)


