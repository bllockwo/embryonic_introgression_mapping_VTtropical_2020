#### summary statistics pi and D.
#### 

library(tidyverse)
library(magrittr)
library(vroom)
library(foreach)
library(patchwork)
library(reshape2)



###
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

###
stats.fi =  list.files(pattern = "pileup.stats")
data.frame(file = stats.fi,
           file.name = stats.fi) %>%
  separate(file, into =  c("samp", "chr", "samtools", "file"), 
           sep = "\\.") -> 
  guide.files

pi.d.stats =
foreach(i=1:dim(guide.files)[1], .combine = "rbind")%do%{
  
  tmp <- vroom(guide.files$file.name[i]) 
  tmp %>%
    dplyr::select(window, length, read_depth, 
                 S,  Watterson, Pi, Tajima_D) %>%
    mutate(chr = guide.files$chr[i],
           samp =  guide.files$samp[i] ) %>%
    group_by(chr) %>%
    mutate(pos.end = cumsum(length)) %>%
    mutate(pos.begin =  pos.end-length) ->
    o
  
  return(o)
}


#######
#save(pi.d.stats, file = "pi.d.stats.Rdata")
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
load("pi.d.stats.Rdata")
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE


pi.d.stats %>%
filter(samp %in% pass.samps) ->
  pi.d.stats
  
pi.d.stats$samp = gsub("VT8", "V8" , pi.d.stats$samp) 
pi.d.stats %>%
  separate(samp, remove = F, into = c("pop", "bio.rep", "tech.rep"), sep = "_" ) %>%
  mutate(pop.group = substr(pop, start = 1, stop = 2)) ->
  pi.d.stats
pi.d.stats$pop.group = gsub("VT", "VT10" , pi.d.stats$pop.group) 
pi.d.stats$pop.group = gsub("V8", "VT8" , pi.d.stats$pop.group) 

#plot coverage
pi.d.stats %>%
  ggplot(aes(
    x=window-1,
    y=read_depth,
    color = bio.rep,
  )) +
  ggtitle("Mapping Coverage") +
  ylab("Coverage X") +
  xlab("Genome window (Mbp)") +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="Replicate")) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(pop.group~chr, scales = "free_x") ->
  depth.stats

ggsave(depth.stats, file = "depth.stats.pdf", w = 9, h = 5)

######

pi.d.stats %>%
  ggplot(aes(
    x=window-1,
    y=Pi,
    color = bio.rep,
  )) +
  ggtitle("Nucleotide diversity") +
  ylab(expression(pi)) +
  xlab("Genome window (Mbp)") +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(color=guide_legend(title="Replicate")) +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(pop.group~chr, scales = "free_x") ->
  pi.stats

ggsave(pi.stats, file = "pi.stats.pdf", w = 9, h = 5)

####
pi.d.stats %>%
  ggplot(aes(
    x=chr,
    y=log10(Pi),
    fill = bio.rep,
  )) +
  #ggtitle("Nucleotide diversity") +
  ylab(expression(log[10](pi))) +
  xlab("Genome window (Mbp)") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="Replicate")) +
  scale_fill_brewer(palette = "Dark2") +
  facet_grid(.~pop.group, scales = "free_x") ->
  pi.stats.box

ggsave(pi.stats.box, file = "pi.stats.box.pdf", w = 8, h = 3.0)

###
pi.d.stats %>%
  group_by(pop.group, bio.rep) %>%
  summarise(mean.pi = mean(Pi)*100)

pi.d.stats %>%
  group_by( bio.rep) %>%
  summarise(mean.pi = mean(Pi)*100)

pi.d.stats %<>% mutate(parental_t = bio.rep==0)

lm(Pi~parental_t, data =pi.d.stats) %>% summary
lm(Pi~parental_t, data =pi.d.stats) %>% aov %>% summary  

lm(Pi~bio.rep*chr, data =pi.d.stats) %>% summary
lm(Pi~bio.rep*chr, data =pi.d.stats) %>% aov %>% TukeyHSD  
