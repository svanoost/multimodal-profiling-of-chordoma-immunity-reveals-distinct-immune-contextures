## Script for generating a combined heatmap of the imaging mass cytometry data
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

# Set working directory
master.location <- setwd(master.location)

# Read the cell phenotype data
pheno_counts <- read.delim("./input_files/mean_phenotype_counts_chordoma.tsv")
row.names(pheno_counts) <- pheno_counts$phenotype

# Read the annotation files for the cell types and for the samples
anno_samples <- read.delim("./input_files/annotation_samples_chordoma.tsv")
row.names(anno_samples) <- anno_samples$L_ID

anno_pheno <- read.delim("./input_files/annotation_phenotypes_chordoma.tsv")
row.names(anno_pheno) <- anno_pheno$phenotype

# Keep only the primary tumours in the annotation file
anno_samples <- anno_samples[!is.na(anno_samples$group),]

# Create list with annotation colours
# The used annotations in the scripts and input data are slightly different from the final figures e.g. location = Anatomical location
anno_colours <- list(location = c("Clivus" = "#4DAF4A", "Mobile Spine" = "#984EA3", "Sacrum" = "#FF7F00"), # Clivus = Skull base; Sacrum = Sacrococcygeal
                     radiotherapy = c("Neoadjuvant" = "black"),
                     ICR.score = c("High" = "#E7B800", "Low" = "#2E9FDF"))

cell.type.colors <- brewer.pal(9, "Set3")

# Keep only primary tumours in the phenotype table
mat <- as.matrix(pheno_counts[, anno_samples$L_ID])

# Check whether sample IDs correspond between the phenotype table and sample annotation file
all(colnames(mat) == anno_samples$L_ID)

# Calculate the Z-Score per phenotype
Density.Z = mat
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat[j,]-mean(mat[j,]))/sd(mat[j,]) # z-score the matrix
}

# OPTIONAL
# Cluster all primary tumours unsupervised to visualise the row order for the overview heatmap (already given in ./input_files)
Heatmap(t(Density.Z),
        name = "Density (Z-Score)",
        row_title = "Samples",
        column_title = "Phenotypes",
        clustering_method_rows = "ward.D",
        clustering_method_columns = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# Read row order from file
row_order_heatmap <- read.delim("./input_files/row_order_heatmap_chordoma.tsv")

# assign the correct row order to the Z-Score matrix and sample annotation file 
Density.Z <- Density.Z[, row_order_heatmap$L_ID]
anno_samples <- anno_samples[row_order_heatmap$L_ID,]

#### Create separate heatmaps for each cell type ####
# The code from line 67:359 will be very repetitve to create separate heatmaps for each cell types.
# The first (T cells) heatmap will be annotated with descriptions. 

#### T cells heatmap ####
# Filter samples based on the cell type
hmap_1 <- anno_pheno %>% filter(Cell.Type == "T")
hmap_1 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_1),] 

# Create an annotation block for the heatmap including box plots
ha_1 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[1]),
                                                   labels = "T cells", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_1),]), 
                                                                     gp = gpar(fill = cell.type.colors[1])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

# Create the row annotation for the samples. Only necessary to run once
hb_1 <- rowAnnotation(df = anno_samples[, c("location", "radiotherapy", "ICR.score")], 
                      col = anno_colours,
                      na_col = "white",
                      gap = unit(c(1, 1), "mm"))

# Create the heatmap for this cell type
h1 <- Heatmap(t(hmap_1), # transpose the heatmap so the phenotype boxplots will be on top
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_1,
              row_gap = unit(2, "mm"),
              row_title = " ",
              left_annotation = hb_1,
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

# Display the heatmap
h1

#### Lymphocyte heatmap ####
hmap_2 <- anno_pheno %>% filter(Cell.Type == "Lym")
hmap_2 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_2),] 

ha_2 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[2]),
                                                   labels = "Lym", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_2),]), 
                                                                     gp = gpar(fill = cell.type.colors[2])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h2 <- Heatmap(t(hmap_2),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_2,
              row_gap = unit(2, "mm"),row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h2

#### Monocytes/Dendritic cells heatmap ####
hmap_3 <- anno_pheno %>% filter(Cell.Type == "Mo")
hmap_3 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_3),] 

ha_3 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[3]),
                                                   labels = "MonoDCs", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_3),]), 
                                                                     gp = gpar(fill = cell.type.colors[3])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h3 <- Heatmap(t(hmap_3),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_3,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h3

#### Macrophages heatmap ####
hmap_4 <- anno_pheno %>% filter(Cell.Type == "Mac")
hmap_4 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_4),] 

ha_4 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[4]),
                                                   labels = "Macro", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_4),]), 
                                                                     gp = gpar(fill = cell.type.colors[4])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h4 <- Heatmap(t(hmap_4),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_4,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h4

#### Granulocytes heatmap ####
hmap_5 <- anno_pheno %>% filter(Cell.Type == "Gra")
hmap_5 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_5),] 

ha_5 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[5]),
                                                   labels = "Gra", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_5),]), 
                                                                     gp = gpar(fill = cell.type.colors[5])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h5 <- Heatmap(t(hmap_5),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_5,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h5

#### Other immune cells heatmap ####
hmap_6 <- anno_pheno %>% filter(Cell.Type == "Oth")
hmap_6 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_6),]

ha_6 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[6]),
                                                   labels = "Other", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_6),]), 
                                                                     gp = gpar(fill = cell.type.colors[6])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h6 <- Heatmap(t(hmap_6),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_6,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h6

#### Stromal cells heatmap ####
hmap_7 <- anno_pheno %>% filter(Cell.Type == "Stro")
hmap_7 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_7),]

ha_7 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[7]),
                                                   labels = "Stroma", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_7),]), 
                                                                     gp = gpar(fill = cell.type.colors[7])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h7 <- Heatmap(t(hmap_7),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_7,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h7

#### Vessels heatmap ####
hmap_8 <- anno_pheno %>% filter(Cell.Type == "Ves")
hmap_8 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_8),]

ha_8 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[8]),
                                                   labels = "Vessels", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_8),]), 
                                                                     gp = gpar(fill = cell.type.colors[8])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h8 <- Heatmap(t(hmap_8),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_8,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

h8

#### Tumour cells heatmap ####
hmap_9 <- anno_pheno %>% filter(Cell.Type == "Tum")
hmap_9 <- Density.Z[row.names(Density.Z) %in% row.names(hmap_9),] 

ha_9 <- HeatmapAnnotation(`Cell type` = anno_block(gp = gpar(fill = cell.type.colors[9]),
                                                   labels = "Tumour", 
                                                   height = unit(0.3, "cm"),
                                                   labels_gp = gpar(col = "black", fontsize = 7)),
                          `Mean cell densities / mm2` = anno_boxplot(t(mat[row.names(mat) %in% row.names(hmap_9),]), 
                                                                     gp = gpar(fill = cell.type.colors[9])),
                          annotation_label = "",
                          height = unit(6, "cm"),
                          gap = unit(1, "mm"))

h9 <- Heatmap(t(hmap_9),
              name = "Density (Z-Score)",
              border = TRUE,
              show_row_names = FALSE,
              column_names_gp = gpar(fontsize = 7),
              top_annotation = ha_9,
              row_gap = unit(2, "mm"),
              row_title = "",
              column_title = "",
              clustering_method_columns = "ward.D",
              clustering_method_rows = "ward.D",
              column_names_rot = -45,
              show_column_dend = FALSE,
              cluster_rows = FALSE,
              show_row_dend = FALSE,
              col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

#### Combine the heatmaps ####
ht_list <- h1 + h2 + h3 +h4 + h5 + h6 + h7 + h8 + h9

#### Save the combined heatmap & the phenotype table and sample annotation file for primary tumours only ####
save(ht_list, file = "./analysis_files/Supplementary_Fig_2_IMC_Heatmap.RData")
save(mat, anno_samples, file = "./analysis_files/IMC_Data_Primary_Tumours.RData")
