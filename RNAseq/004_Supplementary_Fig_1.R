## Script for visualisation of DGE results including all studied sarcoma subtypes 
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(EnhancedVolcano)
library(ggplot2)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

# Load DGE results, data frame with HGNC symbols & design file
load(file = "./analysis_files/DGE_results_Sarcomas.RData")
load(file = "./analysis_files/HGNC_and_Design_Sarcomas.RData")

#### Visualisation of differentially expressed immune-related genes ####
DEG.df$annot <- NA

# All identified immune-related DEGs in descending order based on the log2 fold change
immune_genes <- c("CXCL1", "CXCL6", "ITK", "IL5RA", "IL7", "CCR3", "HLA-DQA2", "IL18", "CD40LG", "C4A", "LAMP3", "CD1E", "CD226", 
                  "MARCO", "LYVE1", "IL17D", "IL11", "IL36RN", "VTCN1", "IL17B")

# Create an annotation column for all immune-related DEGs
DEG.df[DEG.df$hgnc_symbol %in% immune_genes, "annot"] <- DEG.df[DEG.df$hgnc_symbol %in% immune_genes, "hgnc_symbol"]

#### Supplementary Figure 1A ####
# Plot enhanced volcano plot with immune-related DEGs
EnhancedVolcano(DEG.df,
                lab = DEG.df$annot,
                x = 'log2FoldChange',
                y = 'padj',
                xlab = "Log2 fold change",
                ylab = "-Log10 adjusted P",
                legendLabels=c('NS','Log2 FC','adjusted P',
                               'adjusted P & Log2 FC'),
                legendLabSize = 7,
                axisLabSize = 7,
                title = NULL,
                subtitle = NULL,
                pCutoff = 0.002453282,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 1.8,
                xlim = c(-10, 10), # for visualisation of the highlighted genes, the axis are manually set
                ylim = c(0, 100), # for visualisation of the highlighted genes, the axis are manually set
                drawConnectors = TRUE,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                widthConnectors = 0.5)


#### Compare the normalised gene expression for the selected genes between chordomas and other sarcomas ####
df.degs <- df[row.names(df) %in% immune_genes, 1:59]
df.degs <- as.data.frame(t(df.degs))
df.degs$L_ID <- row.names(df.degs)

# Join DEG data frame with design file
df.degs <- left_join(df.degs, design_sarcs[, c("L_ID", "Diagnosis", "vs")])

# Create factor to order the sarcoma subtypes from Chordoma to other sarcomas (Bone sarcomas, Soft tissue sarcomas)
df.degs$Diagnosis <- factor(df.degs$Diagnosis, 
                            levels = c("Chordoma", "Chondrosarcoma", "Osteosarcoma", 'Myxofibrosarcoma', "USTS"))

# melt the DEG data frame for visualisation with box plots
melt_degs <- melt(df.degs, value.name = "value")

# Immune-related DEGs enriched in other sarcomas
sarcoma_only_genes <- c("MARCO", "LYVE1", "IL17D", "IL11", "IL36RN", "VTCN1", "IL17B")

# Create a column for annotation in which group the genes are enriched
melt_degs$DEG <- "Chordoma"

# Assign the group other sarcoma to the correct genes
melt_degs[melt_degs$variable %in% sarcoma_only_genes, "DEG"] <- "Other sarcoma"

# Create a factor to order the genes based on the log2 fold change in descending order
melt_degs$variable <- factor(melt_degs$variable, levels = immune_genes)

#### Supplementary Figure 1B ####
# Plot box plots with immune-related DEGs
ggplot(data = melt_degs, aes(y = value, x = variable))+
  geom_boxplot(aes(fill = Diagnosis))+
  xlab("Differentially expressed immune-related genes")+
  ylab("Log2 normalised gene expression")+
  scale_fill_manual(name = "Diagnosis", 
                     breaks = c("Chordoma", "Chondrosarcoma", "Osteosarcoma", "Myxofibrosarcoma", "USTS"), 
                     values = c("brown1", "darkolivegreen3", 'cornflowerblue', "darkorchid2", "darkgoldenrod3"))+
  facet_wrap(~DEG+variable, scale = "free", ncol = 5)+
  theme_bw()+
  theme(axis.text.x = element_text(face = "italic"),
      axis.text = element_text(colour = "black"),
      panel.grid = element_blank(),
      legend.position = "top")
