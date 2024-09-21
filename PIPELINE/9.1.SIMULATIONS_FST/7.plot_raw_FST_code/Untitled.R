ST1 <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/TABLES_data/TabS1.metadat.txt")
ST1inv <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/TABLES_data/TabS1.inversions_frequency.txt")

full_join(ST1, ST1inv) -> newTs1

write.table(newTs1, file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/DESTv2_data_paper/TABLES_data/TabS1.updJun25.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")

eff_cov = function(C, N){
  
  top = N*C
  bottom = N+C-1
  
  return(top/bottom)
}

newTs1 %>% filter(set != "dgn") %>% .$Cov %>% min(na.rm = T)
newTs1 %>% filter(set != "dgn") %>% .$Cov %>% max(na.rm = T)

Untitled