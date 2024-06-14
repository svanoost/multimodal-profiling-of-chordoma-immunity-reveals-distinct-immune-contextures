## Script for calculating the distribution of HLA class I among the identified T cell groups
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)

# Set working directory
master.location <- setwd(master.location)

# Read the HLA class I expression patterns
HLA <- read.delim("./input_files/survival_data_chordomas.tsv")

# Count the different expression patterns per T cell group
HLA_dist <- HLA %>% 
  group_by(IF.group) %>%
  dplyr::count(HLA)

# Remove all samples either without the necessary data
HLA_dist <- HLA_dist[!is.na(HLA_dist$IF.group),]
HLA_dist <- HLA_dist[!is.na(HLA_dist$HLA),]

# Calculate the percentage of the distinct patterns per group
HLA_dist$perc <- 0

# Calculate the total number in the T cell high group
sum_high <- sum(HLA_dist[HLA_dist$IF.group == "High", "n"])

# Calculate the percentage with "defective" HLA class I
HLA_dist[HLA_dist$IF.group == "High" & HLA_dist$HLA == "Defective", "perc"] <- 
  (HLA_dist[HLA_dist$IF.group == "High" & HLA_dist$HLA == "Defective", "n"]*100) / sum_high

# Calculate the percentage with "positive" HLA class I
HLA_dist[HLA_dist$IF.group == "High" & HLA_dist$HLA == "Positive", "perc"] <- 
  (HLA_dist[HLA_dist$IF.group == "High" & HLA_dist$HLA == "Positive", "n"]*100) / sum_high

# Calculate the total number in the T cell low group
sum_low <- sum(HLA_dist[HLA_dist$IF.group == "Low", "n"])

# Calculate the percentage with "defective" HLA class I
HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Defective", "perc"] <- 
  (HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Defective", "n"]*100) / sum_low

# Calculate the percentage with "positive" HLA class I
HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Positive", "perc"] <- 
  (HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Positive", "n"]*100) / sum_low

# Calculate the percentage with "weak" HLA class I
HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Weak", "perc"] <- 
  (HLA_dist[HLA_dist$IF.group == "Low" & HLA_dist$HLA == "Weak", "n"]*100) / sum_low

#### Save the HLA class I expression pattern distribution for Figure 5B ####
save(HLA_dist, file = "./analysis_files/HLA_distributions.RData")
