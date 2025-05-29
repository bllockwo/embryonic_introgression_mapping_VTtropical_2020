# Load necessary libraries
library(tidyverse)
library(magrittr)
library(foreach)
library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)

### By Luke Proud, May 13, 2025


# Load the data
file = "VA.cm_GEVA_named"
datg <- fread(file)


datg.mutated = datg %>% separate(id, into = c("chr", "array_start", "array_end"))
save(datg.mutated, file = "datg.mutated.Rdata")

# Age by chromosome boxplot
age_by_chromosome = ggplot(datg.mutated, aes(x=chr, y=PostMean)) + 
  geom_boxplot() +
  theme_classic() +
  labs(title = "Age of mutations within chromosomes", x="Chromosome", y="Age")
ggsave("allele_age_by_chromosome.png", plot = age_by_chromosome, width = 8, height = 6)




#finding the ages of interesting SNPs
get_allele_ages <- function(SNP, data) {
  # Filter the dataset for chromosome 2R and the specified positions
  data %>%
    filter(chr == '2R', position %in% SNP) %>%
    select(Position = position, Age = PostMean) %>%
    # Return empty data frame if no matches
    {if (nrow(.) == 0) {
      message("No alleles found at the specified positions in chromosome 2R.")
      data.frame(Position = numeric(0), Age = numeric(0))
    } else {
      .
    }}
}

# Get allele ages for specified SNPs
allele_positions <- c(20551633, 20550549, 20551182, 20551519, 20550323, 
                      20552286, 20551486, 20551129, 20550070, 20550046, 20550403)

allele_ages <- get_allele_ages(allele_positions, datg.mutated)

# Write the allele ages to a CSV file
write.csv(allele_ages, "ages_for_alleles_of_interest_2R.csv", row.names = FALSE)

# Plotting age density with allele points
age_density_2R <- datg.mutated %>%
  filter(chr == "2R") %>%
  ggplot(aes(x = PostMean)) +
  geom_density(alpha = 0.5, fill = "skyblue") +
  # Adding the red points for other SNPs
  geom_point(data = filter(allele_ages, Position != 20551633), aes(x = Age, y = 0), 
             color = "red", size = 3) +
             # Adding the blue point for the special SNP (20551633)
  geom_point(data = filter(allele_ages, Position == 20551633), aes(x = Age, y = 0), 
             color = "blue", size = 3) +
  theme_classic() +
  labs(title = "Age Density Plot for Chromosome 2R", x = "Age", y = "Density")

# Save the plot
ggsave("age_density_2R_7.png", plot = age_density_2R, width = 8, height = 6)





#estimating mean allele age per chromosome
estimate_age_stats <- function(data) {
  # Group data by chromosome and calculate mean and standard deviation
  age_stats <- data %>%
    group_by(chr) %>%
    summarize(
      Chromosome = unique(chr),
      Mean_Age = mean(PostMean, na.rm = TRUE),
      SD_Age = sd(PostMean, na.rm = TRUE)
    ) %>%
    select(Chromosome, Mean_Age, SD_Age)  # Select and order columns
  return(age_stats)
}

age_stats <- estimate_age_stats(datg.mutated)

write.csv(age_stats, "age_stats_by_chromosome.csv", row.names = FALSE)