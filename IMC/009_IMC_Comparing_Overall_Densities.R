## Script for calculating the overall cell densities of the T cell high and the T cell low group from the imaging mass cytometry data
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

# Load cell counts from all phenotypes
load(file = "./analysis_files/IMC_all_phenotype_counts.RData")

# Create vectors to group stromal cells and tumour cells
stroma <- c("PDPN+ Stromal cells", "Stromal cells", "CD39+ Stromal cells", 
            "Ki-67+ Stromal cells", "Vessels", "CD39+ Vessels")
tumour <- c("Tumour cells", "CD56+ Tumour cells", "CD57+ Tumour cells")

# Create a column for the grouping of cell types
melt_counts$cell_group <- "Immune"

# Assign the correct group to Stromal cells and Tumour cells
melt_counts[melt_counts$phenotype %in% stroma, "cell_group"] <- "Stroma"
melt_counts[melt_counts$phenotype %in% tumour, "cell_group"] <- "Tumour"

# Calculate the overall cell density for the 3 major cell groups
melt_overall <- melt_counts %>% 
  select(-phenotype) %>%
  group_by(L_ID, cell_group, group) %>%
  mutate(L_ID = L_ID, group = group, cell_group = cell_group, cell.count.sum = sum(cell.count)) %>%
  select(-cell.count) %>%
  unique()

#### Save overall cell densities for visualisation in Supplementary Figure 5C ####
save(melt_overall, file = "./analysis_files/IMC_overall_density_results.RData")

#### Perform student's t-test to compare T cell high with T cell low group ####
# Create an empty vector and data frame for statistical testing
test <- NA
results <- data.frame(phenotype = 0, p.val = 0, t.test = 0, mean_high = 0, mean_low = 0, median_high = 0, median_low = 0)

# Create a vector with the groups
cell_groups <- unique(melt_overall$cell_group)

# Loop through the major cell groups
for(i in 1:length(cell_groups)){
  test <- t.test(melt_overall[melt_overall$group == "High" & melt_overall$cell_group == cell_groups[i], "cell.count.sum"], 
                 melt_overall[melt_overall$group == "Low" & melt_overall$cell_group == cell_groups[i], "cell.count.sum"])
  results[i, "phenotype"] <- cell_groups[i]
  results[i, "p.val"] <- test$p.value
  results[i, "t.test"] <- test$statistic
  results[i, "mean_high"] <- mean(melt_overall[melt_overall$group == "High" & 
                                                 melt_overall$cell_group == cell_groups[i], 
                                               "cell.count.sum"]$cell.count.sum)
  results[i, "mean_low"] <- mean(melt_overall[melt_overall$group == "Low" & 
                                                melt_overall$cell_group == cell_groups[i], 
                                              "cell.count.sum"]$cell.count.sum)
  results[i, "median_high"] <- median(melt_overall[melt_overall$group == "High" & 
                                                     melt_overall$cell_group == cell_groups[i], 
                                                   "cell.count.sum"]$cell.count.sum)
  results[i, "median_low"] <- median(melt_overall[melt_overall$group == "Low" & 
                                                    melt_overall$cell_group == cell_groups[i], 
                                                  "cell.count.sum"]$cell.count.sum)
}

#### Save student's t-test results ####
write.table(results, "./output_files/IMC_overall_density_comparison_results.tsv")
