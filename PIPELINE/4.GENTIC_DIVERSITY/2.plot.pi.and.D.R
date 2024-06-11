#### summary statistics pi and D.
#### 

library(tidyverse)
library(magrittr)
library(vroom)
library(foreach)
library(patchwork)
library(reshape2)
library(data.table)

###
samps <- fread("/netfiles02/lockwood_lab/IntrogressionProject/Population_files/samp.guide.files.introg.txt")
samps %>% 
  filter(Pop %in% c("SK","VT8",
                    "SKF","VT8F",
                    "SKF2","VT8F2",
                    "SKF3","VT8F3",
                    "ME_temp","PAN_trop"
  )) -> samps.flt

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
load("/netfiles02/lockwood_lab/IntrogressionProject/NucleotideDiversity/pi.d.stats.Rdata")
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE
########### <------ LOAD HERE


pi.d.stats %>%
filter(samp %in% samps.flt$Sample_id) ->
  pi.d.stats

pi.d.stats %<>%
  mutate(type = case_when(
    samp %in% c("SK_0_1", "V8_0_1") ~ "Parentals",
    samp %in% c("SKF_1_1", "SKF2_2_1", "SKF3_3_1") ~ "F16-SK-mother",
    samp %in% c("V8F_1_1", "V8F2_2_1", "V8F3_3_1" ) ~ "F16-VT8-mother",
  )) 
pi.d.stats$type = factor(pi.d.stats$type, levels = c("Parentals","F16-SK-mother", "F16-VT8-mother"))

###
pi.d.stats$samp = gsub("VT8", "V8" , pi.d.stats$samp) 
pi.d.stats %>%
  separate(samp, remove = F, into = c("pop", "bio.rep", "tech.rep"), sep = "_" ) %>%
  mutate(pop.group = substr(pop, start = 1, stop = 2)) ->
  pi.d.stats
pi.d.stats$pop.group = gsub("VT", "VT10" , pi.d.stats$pop.group) 
pi.d.stats$pop.group = gsub("V8", "VT8" , pi.d.stats$pop.group) 

##### Plot nuc. Diversity

pi.d.stats %>%
  ggplot(aes(
    x=samp,
    y=(Pi),
    fill = chr,
  )) +
  ylab(expression(pi)) +
  xlab("Chromosome") +
  geom_boxplot(outlier.shape = NA) +
  theme_bw() +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="Replicate")) +
  scale_fill_brewer(palette = "Dark2") +
  facet_wrap(.~type, scales = "free_x", nrow = 1) ->
  pi.stats.box

ggsave(pi.stats.box, file = "pi.stats.box.pdf", w = 8, h = 3)

#####
pi.d.stats %>%
  group_by(samp) %>%
  summarize(mC = mean(read_depth))

#plot coverage
pi.d.stats %>%
  ggplot(aes(
    x=window-1,
    y=read_depth,
    color = samp,
  )) +
  ggtitle("Mapping Coverage") +
  ylab("Coverage X") +
  xlab("Genome window (Mbp)") +
  geom_line(linewidth = 0.8) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_brewer(palette = "Dark2") +
  facet_grid(~chr, scales = "free_x") ->
  depth.stats

ggsave(depth.stats, file = "depth.stats.pdf", w = 9, h = 3)

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
  facet_wrap(.~pop.group, scales = "free_x", ncol = 2) ->
  pi.stats.box

ggsave(pi.stats.box, file = "pi.stats.box.pdf", w = 4.5, h = 5)

###
pi.d.stats %>%
  group_by(pop.group, type) %>%
  summarise(mean.pi = mean(Pi)*100)


pi.d.stats %>%
  group_by( bio.rep) %>%
  summarise(mean.pi = mean(Pi)*100)

pi.d.stats %<>% mutate(parental_t = bio.rep==0)

lm(Pi~type, data =pi.d.stats) %>% summary
lm(Pi~type, data =pi.d.stats) %>% aov %>% summary  

lm(Pi~bio.rep*chr, data =pi.d.stats) %>% summary
lm(Pi~bio.rep*chr, data =pi.d.stats) %>% aov %>% TukeyHSD  
