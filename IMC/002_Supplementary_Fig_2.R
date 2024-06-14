## Script for visualisation of the imaging mass cytometry data overview
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Load combined heatmap
load(file = "./analysis_files/Supplementary_Fig_2_IMC_Heatmap.RData")

#### Supplementary Figure 2 ####
# Plot the combined heatmaps
draw(ht_list, 
     row_title = "Samples", 
     row_title_gp = gpar(font = 1, fontsize = 7), 
     column_title = "Phenotypes",
     column_title_gp = gpar(font = 1, fontsize = 7),
     ht_gap = unit(c(0.8), "cm"),
     adjust_annotation_extension = FALSE)
