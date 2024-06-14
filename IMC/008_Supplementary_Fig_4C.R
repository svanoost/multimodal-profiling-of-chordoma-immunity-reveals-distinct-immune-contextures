## Script for visualsing the microenvironment during the course of disease from the imaging mass cytometry data
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)

# Set working directory
master.location <- setwd(master.location)

# Read the cell phenotype data
pheno_counts <- read.delim("./input_files/mean_phenotype_counts_chordoma.tsv")
row.names(pheno_counts) <- pheno_counts$phenotype

# Read the annotation files for the cell types and for the samples
anno_samples <- read.delim("./input_files/annotation_samples_chordoma.tsv")
row.names(anno_samples) <- anno_samples$L_ID

# Keep only the patients with multiple samples in the annotation file
anno_samples <- anno_samples[!is.na(anno_samples$Pat_ID),]

# Create list with annotation colours
# The used annotations in the scripts and input data are slightly different from the final figures e.g. location = Anatomical location
anno_colours <- list(location = c("Clivus" = "#4DAF4A", "Mobile Spine" = "#984EA3", "Sacrum" = "#FF7F00"), # Clivus = Skull base; Sacrum = Sacrococcygeal
                     radiotherapy = c("Neoadjuvant" = "black"),
                     PRM = c("Primary" = "#FBB4AE", "Recurrence" = "#B3CDE3", "Metastasis" = "#CCEBC5"), # PRM = Sample type
                     Pat_ID = c("P1" = "#05288C", "P2" = "#799303", "P3" = "#CE06D4", 
                                "P4" = "#067804", "P5" = "#F9D40F", "P6" = "#D84B02"),
                     group = c("High" = "#DC0018", "Low" = "#2D69A9")) # group = T cell infiltration

# Create a vector to remove Tumour cells for visualisation
tumour <- c("Tumour cells", "CD56+ Tumour cells", "CD57+ Tumour cells")

# Remove Tumour cells from the data frame & align the order with the annotation file
mat1 <- as.matrix(pheno_counts[!row.names(pheno_counts) %in% tumour, anno_samples$L_ID])

# Calculate the Z-Score per phenotype
Density.Z = mat1
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat1[j,]-mean(mat1[j,]))/sd(mat1[j,]) # z-score the matrix
}

# Check whether the order of the data frame is the same as the annotation file
all(colnames(Density.Z) == anno_samples$L_ID)

# Create a factor of the sample type for visualisation
anno_samples$PRM <- factor(anno_samples$PRM, 
                           levels = c("Primary", "Recurrence", "Metastasis"))

#### Supplementary Figure 4 ####
# Create top annotation for the heatmap
ha1 <- HeatmapAnnotation(df = anno_samples[, c("Pat_ID", "PRM", "location", "radiotherapy", "group")], 
                         col = anno_colours,
                         na_col = "white",
                         gap = unit(c(1, 1, 1, 1), "mm"))

# Plot the heatmap for all patient with multiple samples
Heatmap(Density.Z,
        name = "Density (Z-Score)",
        column_title = "Samples",
        row_title = "Phenotypes",
        show_column_names = FALSE,
        top_annotation = ha1,
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))
