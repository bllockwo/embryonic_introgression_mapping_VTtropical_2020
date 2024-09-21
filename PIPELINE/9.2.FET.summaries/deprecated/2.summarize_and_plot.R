##### Summarize and aggregate the FET P-vals sim vs real.

library(tidyverse)
library(magrittr)
library(foreach)
library(vroom)
library(foreach)
library(data.table)
library(gmodels)


reco_windows_f="/gpfs2/scratch/jcnunez/fst_brent/simulations/reco_wins.master.txt"
reco_windows = fread(reco_windows_f)

#### Binning analysis -- deprecated -- stays as a commented out code.

grid <- 0.01
my_seq = data.table(min_p=seq(from=0, to=1-grid, by=grid), max_p=seq(from=grid, to=1, by=grid))
setkey(tmp, chr)

####
real_dat <-"/netfiles02/lockwood_lab/IntrogressionProject/FET_output/Fisher_test_results_All_crosses.txt"
real_FET <- fread(real_dat)
names(real_FET)[1] = "chr"

pops <- c("skf1_f.test_pval", 
         "skf2_f.test_pval", 
         "skf3_f.test_pval",
         "vt8f1_f.test_pval", 
         "vt8f2_f.test_pval", 
         "vt8f3_f.test_pval")

real_enrich <-
  foreach(pop.i=pops, .combine="rbind")%do%{
    
    tmp <- real_FET %>% dplyr::select(chr,POS,all_of(pop.i))  
    setDT(tmp)
    setnames(tmp, 2, "real_pos")
    setnames(tmp, 3, "p_fet")
    
    L = dim(tmp)[1] 
    tmp %>%
      arrange(p_fet) %>% 
      as.data.frame() %>%
      mutate(rank = seq(from=1,  to = L, by = 1)) %>% 
      mutate(rank_norm = rank/L) %>% 
      arrange(as.numeric(real_pos))->
      tmp.rnf
    
    message(pop.i)
    
    foreach(k=1:dim(reco_windows)[1], .combine = "rbind")%do%{
      
      chr.i=reco_windows$V1[k]
      sta =reco_windows$V3[k]
      end =reco_windows$V4[k]
      
      tmp.rnf %>%
        filter(chr == chr.i) %>%
        filter(real_pos >= sta & real_pos <= end ) ->
        tmp.slice
      
      p.i=0.01
      tmp.slice %>%
        summarise(N = sum(rank_norm < p.i),
                  Lw = n()) %>%
        mutate(p.binom = c(binom.test(N, Lw, p.i)$p.value)) %>%
        mutate(type=pop.i,
               id= k,
               chr=chr.i,
               win.sta=sta,
               win.edn=end
        ) -> o
      return(o)
      
    }
    
##foreach(chr.i=unique(tmp$chr), .combine="rbind")%do%{
##  message(paste(chr.i, pop.i, sep=" / "))
##  
##  tmp <- real_FET %>% dplyr::select(chr, all_of(pop.i))
##  setDT(tmp)
##  setnames(tmp, 2, "p_fet")
##  
##  new_tmp <- tmp[chr==chr.i][my_seq, .(N = .N), on = .(p_fet > min_p, p_fet <= max_p), by = .EACHI]
##  new_tmp[,chr:=chr.i]
##  new_tmp[,pop:=pop.i]
##  setnames(new_tmp, 1:2, c("p.low", "p.high"))
##  
##  tot <- sum(new_tmp$N)
##  sig <- new_tmp$N[which(new_tmp$p.high == grid)]
##  
##  data.frame(chr = chr.i,
##             pop = pop.i,
##             type = "real",
##             tot,
##             sig,
##             prop = sig/tot*100
##             ) -> o
##  
##  return(o)
##
##  }
    
    }

save(real_enrich, file = "real_enrich.FET.ag.win.Rdata")


#### window analysis
###


FETfils = system("ls /gpfs2/scratch/jcnunez/fst_brent/simulations/FET_SiMS_Ps", 
                 intern = T)

sims_win_ag = 
foreach(i=FETfils,
        .combine = "rbind"
        )%do%{
 tmp=get(load(
            paste("/gpfs2/scratch/jcnunez/fst_brent/simulations/FET_SiMS_Ps/", i, sep = "")))
 setDT(tmp)
 
 L = dim(tmp)[1] 
 
 tmp %>%
   arrange(p_fet) %>% 
   as.data.frame() %>%
   mutate(rank = seq(from=1,  to = L, by = 1)) %>% 
   mutate(rank_norm = rank/L) %>% 
   arrange(as.numeric(real_pos))->
   tmp.rnf
 
 message(i)
 
 foreach(k=1:dim(reco_windows)[1], .combine = "rbind")%do%{
   
   chr.i=reco_windows$V1[k]
   sta =reco_windows$V3[k]
   end =reco_windows$V4[k]
   
   tmp.rnf %>%
     filter(chr == chr.i) %>%
     filter(real_pos >= sta & real_pos <= end ) ->
     tmp.slice
   
   p.i=0.01
   tmp.slice %>%
     summarise(N = sum(rank_norm < p.i),
               Lw = n()) %>%
     mutate(p.binom = c(binom.test(N, Lw, p.i)$p.value)) %>%
     mutate(type="simulation",
             id= k,
            chr=chr.i,
            win.sta=sta,
            win.edn=end
            ) -> o
   return(o)
   
 }
 
 #foreach(chr.i=unique(tmp$chr), .combine="rbind")%do%{
 # message(paste(i, chr.i, sep=" / "))
 #  
 #  new_tmp <- tmp[chr==chr.i][my_seq, .(N = .N), on = .(p_fet > min_p, p_fet <= max_p), by = .EACHI]
 #  new_tmp[,chr:=chr.i]
 #  #new_tmp[,inv:=inv.i]
 #  #new_tmp[,perm:=perm.i]
 #  #new_tmp[,perm_method:="new"]
 #  setnames(new_tmp, 1:2, c("p.low", "p.high"))
 #  
 #  tot <- sum(new_tmp$N)
 #  sig <- new_tmp$N[which(new_tmp$p.high == grid)]
 #  
 #  data.frame(chr = chr.i,
 #             pop = paste("sim", i, sep = "_"),
 #             type = "sim",
 #             tot,
 #             sig,
 #             prop = sig/tot*100
 #  ) -> o
 #  
 #  return(o)
 #  
 #}
          
}

save(sims_win_ag, file = "sims.FET.ag.win.Rdata")

######
sims_win_ag %>%
  group_by(chr, win.sta, win.edn) %>%
  summarise(u95 = quantile(p.binom, 0.95)) ->
  sims_win_ag.u95


######
load("real_enrich.FET.ag.win.Rdata")
load("sims.FET.ag.win.Rdata")

real_enrich %>% head
sims_win_ag.u95 %>% head

ggplot() +
  geom_vline(xintercept = 19115753/1e6) +
  geom_vline(xintercept = 20615753/1e6) +
  geom_line(
    data=real_enrich,
    aes(
      x=( win.sta + win.edn)/2e6,
      y=-log10(p.binom),
      color=type
    )) +
  geom_line(
    data=sims_win_ag.u95,
    aes(
      x=( win.sta + win.edn)/2e6,
      y=-log10(u95)
    ), color = "grey", linewidth = 1.9) + 
  scale_color_brewer(palette = "Set3") +
  theme_minimal() +  
  facet_grid(type~chr, scales = "free_x") -> 
  FET_enrich

ggsave(FET_enrich, file = "FET_enrich.pdf")
