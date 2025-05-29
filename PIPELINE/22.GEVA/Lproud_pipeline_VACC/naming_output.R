library(tidyverse)
library(magrittr)
library(foreach)
library(data.table)


###

# Load the data
fil <- "VA.cm_GEVA_complete/VA.cm_GEVA_complete.txt"
datg <- fread(fil)
names(datg) = c("id","position", "MarkerID", "Clock", "Filtered", "N_concordant", "N_discordant", "PostMean", "PostMode", "PostMedian", "array_ID")

write.table(datg, file = "VA.cm_GEVA_named", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

