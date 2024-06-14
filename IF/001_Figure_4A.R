## Script for visualisation of the T cell groups from the immunofluorescence data  
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Read T cell count data from the immunofluorescence data & survival data
IF_counts <- read.delim("./input_files/mean_Tcell_counts_IF.tsv")
row.names(IF_counts) <- IF_counts$phenotype

anno_samples <- read.delim("./input_files/annotation_full_cohort_chordoma.tsv")
row.names(anno_samples) <- anno_samples$L_ID

# Create a list with the annotation colours
# The used annotations in the scripts and input data are slightly different from the final figures e.g. location = Anatomical location
anno_colours <- list(location = c("Clivus" = "#4DAF4A", "ExtraAxial" = "#FFFF33", 
                                  "Mobile Spine" = "#984EA3", "Sacrum" = "#FF7F00"), # Clivus = Skull base; Sacrum = Sacrococcygeal
                     radiotherapy = c("Neoadjuvant" = "black"),
                     IF.group = c("High" = "#DC0018", "Low" = "#2D69A9")) # group = T cell infiltration

# Remove the phenotype column
mat1 <- as.matrix(IF_counts[, -1])

# Calculate the Z-Score per T cell phenotype
Density.Z = mat1
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat1[j,]-mean(mat1[j,]))/sd(mat1[j,]) # z-score the matrix
}

# Make sure the data frame has the same order as the annotation file
anno_samples <- anno_samples[colnames(Density.Z),]
all(colnames(Density.Z) == anno_samples$L_ID)

# Create top annotation for heatmap
ha1 <- HeatmapAnnotation(df = anno_samples[, c("location", "radiotherapy", "IF.group")], 
                         col = anno_colours,
                         na_col = "white",
                         gap = unit(c(1, 1), "mm"))

#### Figure 4A ####
# Plot heatmap with T cell phenotypes for entire cohort
Heatmap(Density.Z,
        name = "Density (Z-Score)",
        column_title = "Samples",
        row_title = "T cell phenotypes",
        top_annotation = ha1,
        show_column_names = FALSE,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
