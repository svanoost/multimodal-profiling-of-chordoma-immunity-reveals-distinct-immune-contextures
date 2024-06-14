## Script for correlating the CD4+ T cell infiltration with the CD8+ T cell infiltration in the combined cohort. 
## Also includes the calculation of the CD4/CD8 ratio
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)

# Set working directory
master.location <- setwd(master.location)

# Read T cell count data from the immunofluorescence data
IF_counts <- read.delim("./input_files/mean_Tcell_counts_IF.tsv")
row.names(IF_counts) <- IF_counts$phenotype

# Transpose the data frame for visualisation and ratio calculation
IF_counts <- as.data.frame(t(IF_counts[,-1]))
colnames(IF_counts) <- c("CD4", "CD8", "Tregs", "KI67")

#### Supplementary Figure 6A ####
# Plot the correlation between CD4+ and CD8+ T cells in the combined cohort
ggscatter(IF_counts, x = "CD8", y = "CD4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CD8+ T cells", ylab = "CD4+ T cells")

# Calculate the CD4/CD8 ratio
IF_counts$ratio <- IF_counts$CD4 / IF_counts$CD8
IF_counts$type <- "Chordoma"

#### Supplementary Figure 6C ####
# Plot the CD4/CD8 ratio for the combined cohort
ggplot(data = IF_counts, aes(x = type, y= ratio))+
  geom_boxplot()+
  theme_bw()+
  xlab(NULL)+ylab("CD4/CD8 ratio")+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))
