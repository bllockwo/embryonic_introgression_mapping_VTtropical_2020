### Make thesis figures

library(tidyverse)
library(magrittr)
library(data.table)
library(reshape2)
library(foreach)
library(forcats)
library(lme4)

###load data
dat  <- fread("/netfiles/nunezlab/D_melanogaster_resources/Datasets/Embryonic_Thermal_Tolerance/Final_data.txt")

#dat %<>%
#  filter(Collector == "EB")
### Plot 1: We will plot all the lines in a single plot
dat.wBinom <-
  foreach(i=1:dim(dat)[1],
          .combine = "rbind",
          .errorhandling = "remove")%do%{
            
            tmp <- dat[i,]
            
            test <- prop.test(x=tmp$Hatched, n=tmp$Sample_size, 
                             p =0.5)
            
            data.frame(tmp,
                       p.val=test$p.value,
                       uci=test$conf.int[1],
                       lci=test$conf.int[2]
                       )
          }

### Global test --- ANOVA
model <- glm(cbind(Hatched, (Sample_size-Hatched)) ~ Genotype*SP70_allele, 
    data = filter(dat.wBinom, Wolbachia_Status == "No"), family= binomial)

anova(model)

t.test(Proportion~SP70_allele, data = filter(dat.wBinom, Wolbachia_Status == "No"))
###########


dat.wBinom %>%
  filter(Wolbachia_Status == "No") %>%
  ggplot(aes(
    x=fct_reorder(Line, Proportion),
    y=Proportion,
    ymin=uci,
    ymax=lci,
    #shape = Genotype,
  )) + geom_errorbar(width = 0.3) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.5) +
  theme_bw() +
  theme(axis.text.x = 
          element_text(angle = 90, 
                       vjust = 0.5, 
                       hjust=1)) +
  facet_grid(.~SP70_allele, scales = "free_x") +
  guides(color=guide_legend(title="Inversion")) +
  xlab("Line") + ylab("Proportion Hatched") ->
  data.plot
ggsave(data.plot, file = "data.plot.pdf", h = 3.5, w = 4)

dat.wBinom %>%
  filter(Wolbachia_Status == "No") %>%
  ggplot(aes(
    x=SP70_allele,
    color=SP70_allele,
    y=Proportion,
  )) + theme_bw() +
 geom_boxplot() ->
  data.plot.box
ggsave(data.plot.box, file = "data.plot.box.pdf", h = 3.5, w = 4)
  
