## Script for analysing the major T cell groups from the imaging mass cytometry data 
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

# Load the cell phenotype data & sample annotation files for primary tumours only 
load(file = "./analysis_files/IMC_Data_Primary_Tumours.RData")

# Filter samples with an ICR score & select T cell phenotypes
df <- as.data.frame(mat)
df$phenotype <- row.names(df)

# Melt the data frame for ggplot visualisation
melt_counts <- melt(df, value.name = "cell.count", variable.name = "L_ID")

# Join the phenotype counts with the sample annotation
melt_counts <- left_join(melt_counts, anno_samples, by = "L_ID")

#### Save all cell counts data frame ####
save(melt_counts, file = "./analysis_files/IMC_all_phenotype_counts.RData")

# Filter for only T cell phenotypes
melt_Tcells <- melt_counts[grepl("T cells", melt_counts$phenotype),]

#### Save T cell phenotypes for visualisation in Figure 2D ####
save(melt_Tcells, file = "./analysis_files/IMC_Tcell_counts.RData")

#### Student's t-test comparing phenotypes between T cell high and T cell low ####
# Remove T cell phenotypes
melt_pheno <- melt_counts[!grepl("T cells", melt_counts$phenotype),]

# Create a vector with the phenotypes
my_phenos <- unique(melt_pheno$phenotype)

# Create an empty vector and data frame for statistical testing
test <- NA
results <- data.frame(phenotype = 0, p.val = 0, t.test = 0, mean_high = 0, mean_low = 0, median_high = 0, median_low = 0)

# Loop through the phenotypes and perform student's t-test
for(i in 1:length(my_phenos)){
  test <- t.test(melt_pheno[melt_pheno$group == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"], 
                 melt_pheno[melt_pheno$group == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "phenotype"] <- my_phenos[i]
  results[i, "p.val"] <- test$p.value
  results[i, "t.test"] <- test$statistic
  results[i, "mean_high"] <- mean(melt_pheno[melt_pheno$group == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "mean_low"] <- mean(melt_pheno[melt_pheno$group == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "median_high"] <- median(melt_pheno[melt_pheno$group == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "median_low"] <- median(melt_pheno[melt_pheno$group == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
}

# Correct for multiple testing with the false discovery rate / Benjamini-Hochberg
results$FDR <- p.adjust(results$p.val, method = "fdr", n = nrow(results))

# Filter only significant phenotypes for visualisation
significant_phenos <- results[results$FDR < 0.05, "phenotype"]
melt_pheno <- melt_pheno[melt_pheno$phenotype %in% significant_phenos,]

#### Save the significant phenotype cell counts for visualisation in Figure 3A & Supplementary Figure 5A ####
save(melt_pheno, file = "./analysis_files/IMC_Tcell_group_comparison_counts.RData")

#### Save results from the student's t-test with FDR correction ####
write.table(results, "./output_files/IMC_Tcell_groups_FDR_results.tsv", sep = "\t", row.names = FALSE)
