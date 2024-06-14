### Script for visualisation of variable gene expression results including all studied sarcoma subtypes 
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ComplexHeatmap)
library(circlize)
library(ggplot2)

# Set working directory
master.location <- setwd(master.location)

# Load design file, ICR score and variable genes Z-Scores
load(file = "./analysis_files/HGNC_and_Design_Sarcomas.RData")
load(file ="analysis_files/ICRscore_Sarcomas.RData")
load(file ="analysis_files/ZScore_VariableGenes_Sarcomas.RData")

#### Variable Gene expression across sarcomas ####
# Assign colours for annotation of heatmaps
anno_colours <- list(Diagnosis = c("Chordoma" = "brown1", "Chondrosarcoma" = "darkolivegreen3", "Myxofibrosarcoma" = "darkorchid2",
                                   "Osteosarcoma" = "cornflowerblue", "USTS" = "darkgoldenrod3"))

# Create an annotation data frame with only the diagnosis as variable
design_sarcs_annotation <- data.frame(Diagnosis = design_sarcs$Diagnosis, row.names = design_sarcs$L_ID)

# Check whether the order of the data frame corresponds to the annotation file
all(colnames(Density.Z) == row.names(design_sarcs_annotation))

#### Figure 1A ####
# Create top annotation for heatmap
ha1 <- HeatmapAnnotation(df = design_sarcs_annotation, col = anno_colours)

# Plot the heatmap with the most variable genes across sarcomas
Heatmap(Density.Z,
        name = "Expression (Z-Score)",
        column_title = "Samples",
        row_title = "Top 200 most vairably expressed genes",
        show_column_names = FALSE,
        top_annotation = ha1,
        row_split = 2, # By splitting the rows, the chordoma-associated cluster moves to the top op the heatmap
        row_gap = unit(0, "mm"),
        row_names_gp = gpar(fontsize = 5),
        clustering_method_columns = "ward.D",
        clustering_method_rows = "ward.D",
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")))

#### ICR score across sarcomas ####
#### Figure 1C ####
# Plot box plots with ICR score per sarcoma subtype
ggplot(data = ICR.score.df)+
  geom_boxplot(aes(x = Diagnosis, y = ICR.score, fill = Diagnosis))+
  xlab(NULL)+
  ylab("ICR score")+
  scale_fill_manual(name = "Diagnosis", 
                     breaks = c("Chondrosarcoma", "Chordoma", "Myxofibrosarcoma", "Osteosarcoma", "USTS"), 
                     values = c("darkolivegreen3", "brown1", "darkorchid2", 'cornflowerblue', "darkgoldenrod3"))+
  theme_bw()+
  theme(axis.text = element_text(color = "black"),
        legend.position = "none",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
