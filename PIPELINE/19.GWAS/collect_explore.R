### libraries
library(data.table)
library(tidyverse)
library(foreach)


autocor = function(x) {
  abs(cor(x[-1], x[-length(x)]))
}
thresUnif = function(L, cor, xi, alpha = 0.05) {
  a.poly.coef = rbind(c(-5.5, 2.47, 2.04, 0.22), c(6.76, 
                                                   -4.16, -5.76, -4.08), c(-5.66, -1.82, 1.04, 1.16), 
                      c(-2.51, -4.58, -6.95, -9.16))
  b.poly.coef = rbind(c(-1.22, 0.37, 2.55, 3.45), c(3.17, 
                                                    2.14, -0.02, -0.98), c(-1.99, -2.35, -2.31, -2.33))
  if (sum(xi != 1:4) == 4) {
    print("xi must be equal to 1, 2, 3 or 4")
    thres = NULL
  }
  else {
    tmp.cor2 = cor^2
    a = log(L) + a.poly.coef[1, xi] * (cor^3) + a.poly.coef[2, 
                                                            xi] * (cor^2) + a.poly.coef[3, xi] * cor + a.poly.coef[4, 
                                                                                                                   xi]
    b = b.poly.coef[1, xi] * (cor^2) + b.poly.coef[2, 
                                                   xi] * cor + b.poly.coef[3, xi]
    thres = (log(-log(1 - alpha)) - a)/b
  }
  return(thres)
}

files = 
system("ls /gpfs2/scratch/jcnunez/fst_brent/GWAS_EB_OK/GWAS_out/batch.*.txt", 
       intern = T)

o = foreach(i=files,
            .combine = "rbind")%do%{
    tmp <- fread(i)
            }

o %>% filter(pos == 15501637)

o %>%
  group_by(CHR) %>%
  arrange(POS) -> o

o %>%
  filter(CHR == "X") %>% filter(POS == 15501637 )
####
o %>%
ggplot(aes(
  x=POS,
  y=-log10(PVAL))) + geom_point() +
  facet_grid(~CHR) ->
  gwas
ggsave(gwas, file = "gwas.png", w=9, h=4)

###
pstmp=as.data.frame(o[,2:3])
pi=o$AF
pval=-log10(o$PVAL)
####
png("lindley.png")
gwas.lcos=compute.local.scores(pstmp,
                               snp.pi=pi,
                               snp.pvalue = pval,
                               xi=1,manplot = T,main="xi=2")
dev.off()

