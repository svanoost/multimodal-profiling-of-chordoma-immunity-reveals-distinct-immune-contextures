## Script for visualisation of the T cell groups from the imaging mass cytometry data  
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Load the cell phenotype data & sample annotation files for primary tumours only 
load(file = "./analysis_files/IMC_Data_Primary_Tumours.RData")

# Create list with annotation colours
# The used annotations in the scripts and input data are slightly different from the final figures e.g. location = Anatomical location
anno_colours <- list(location = c("Clivus" = "#4DAF4A", "Mobile Spine" = "#984EA3", "Sacrum" = "#FF7F00"), # Clivus = Skull base; Sacrum = Sacrococcygeal
                     radiotherapy = c("Neoadjuvant" = "black"),
                     ICR.score = c("High" = "#E7B800", "Low" = "#2E9FDF"),
                     Status = c("Alive" = "#66A61E", "Dead" = "#E7298A"),
                     group = c("High" = "#DC0018", "Low" = "#2D69A9")) # group = T cell infiltration 

# Keep only the T cell phenotypes & align the order with the annotation file
mat1 <- mat[grepl("T cells", row.names(mat)), anno_samples$L_ID]

# Calculate the Z-Score per phenotype
Density.Z = mat1
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat1[j,]-mean(mat1[j,]))/sd(mat1[j,]) # z-score the matrix
}

# Make sure the data frame has the same order as the annotation file
all(colnames(Density.Z) == anno_samples$L_ID)

# Create top annotation for heatmap
ha1 <- HeatmapAnnotation(df = anno_samples[, c("location", "radiotherapy", "ICR.score", "Status", "group")], 
                         col = anno_colours,
                         na_col = "white",
                         gap = unit(c(1, 1, 1, 1), "mm"))

#### Figure 2B ####
# Plot heatmap with T cell phenotypes in 32 primary tumours
Heatmap(Density.Z,
        name = "Density (Z-Score)",
        column_title = "Samples",
        row_title = "T cell phenotypes",
        show_column_names = FALSE,
        top_annotation = ha1,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
