## Script for Calculating the percentage of Ki-67+ T cells per phenotype
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

# Load the T cell phenotypes from the immunofluorescence data, with separate Ki-67+ subtypes
IF.counts <- read.delim("./input_files/mean_Tcell_counts_IF_Ki67_separate.tsv")
row.names(IF.counts) <- IF.counts$phenotype

# Transpose the data frame and rename the columns before calculating the percentages
IF.counts <- as.data.frame(t(IF.counts[, -1]))
colnames(IF.counts)[c(2, 5, 6)] <- c("Ki67_CD4_Tcells", 
                                     "Ki67_CD8_Tcells",
                                     "Ki67_Tregs")

# Calculate the total number of T cells per phenotype
IF.counts[, "total_cd8"] <- IF.counts[, "CD8_Tcells"]+IF.counts[, "Ki67_CD8_Tcells"]
IF.counts[, "total_cd4"] <- IF.counts[, "CD4_Tcells"]+IF.counts[, "Ki67_CD4_Tcells"]

# Calculate the percentages of Ki-67+ T cells
IF.counts[, "perc_cd8"] <- IF.counts[, "Ki67_CD8_Tcells"]*100/IF.counts[, "total_cd8"]
IF.counts[, "perc_cd4"] <- IF.counts[, "Ki67_CD4_Tcells"]*100/IF.counts[, "total_cd4"]

# Create a data frame for to plot the percentages with ggplot
my_melt <- melt(IF.counts[, c("perc_cd4", "perc_cd8")])

#### Supplementary Figure 6B ####
# Plot the percentages of Ki67+ T cells in boxplots
ggplot()+
  geom_boxplot(data = my_melt, aes(x = variable, y = value))+
  theme_bw()+
  scale_x_discrete(labels = c("CD4+", "CD8+"))+
  xlab("T cell phenotype")+ylab("Percentage of Ki-67+ cells")+
  theme(axis.text = element_text(color = "black"),
        panel.grid = element_blank())

# Compare the percentages of Ki-67+CD4+ T cells with those of Ki-67+CD8+ T cells
t.test(IF.counts$perc_cd4, IF.counts$perc_cd8)



