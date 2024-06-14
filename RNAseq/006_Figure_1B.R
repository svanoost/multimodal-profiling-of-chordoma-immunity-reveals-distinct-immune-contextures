### Script for visualisation of the ICR score of chordomas
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Load normalised data and design file
load("./analysis_files/DESeq2_Normalised_Data_Chordomas.RData")

# Create top annotation for the heatmap
ha1 <- HeatmapAnnotation(df = design_chord[, c("Diagnosis", "ICR.score")], 
                         col = list(Diagnosis = c("Chordoma" = "brown1"),
                                    ICR.score = c("High" = "#E7B800", "Low" = "#2E9FDF")))

#### Figure 1B ####
# Plot the ICR heatmap
Heatmap(Density.Z,
        name = "Expression (Z-Score)",
        column_title = "Samples",
        row_title = "Immune-related genes",
        show_column_names = FALSE,
        top_annotation = ha1,
        row_labels = ICR_genes$Gene.name, # use the HGNC symbols as labels
        row_names_gp = gpar(fontface = "italic"),
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
