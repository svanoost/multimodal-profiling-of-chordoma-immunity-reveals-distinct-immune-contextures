## Script for visualisation of the T cell groups from the imaging mass cytometry data  
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)

# Set working directory
master.location <- setwd(master.location)

# Load IMC T cell count data & significant phenotype count data
load(file = "./analysis_files/IMC_Tcell_counts.RData")
load(file = "./analysis_files/IMC_Tcell_group_comparison_counts.RData")

#### Figure 2D ####
# Transform the phenotypes into a factor for visualisation
melt_Tcells$phenotype <- factor(melt_Tcells$phenotype, 
                               levels = c("T cells", 'CD4+ T cells', "CD8+ T cells", 
                                          "Ki-67+ T cells", "Regulatory T cells"))

# Plot the box plots comparing T cell phenotypes between both groups
ggplot(data = melt_Tcells, aes(x = phenotype, y = cell.count, fill = group))+
  geom_boxplot()+
  xlab("T cell phenotypes")+
  ylab("Mean cell density / mm2")+
  scale_x_discrete(labels = c("T cells" = "T cells",
                              "CD4+ T cells" = "CD4+",
                              "CD8+ T cells" = "CD8+",
                              "Ki-67+ T cells" = "Ki-67+",
                              "Regulatory T cells" = "Tregs"))+
  scale_fill_manual(name = "T cell infiltration", 
                    breaks = c("High", "Low"), 
                    values = c("#DC0018", "#2D69A9"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"))


#### Figure 3A ####
# Create a vector with the phenotypes to plot in Figure 3A
fig3a <- c("Dendritic cells", "Immune cell aggregates")

# Plot box plots for dendritic cells and immune cell aggregates
ggplot(data = melt_pheno[melt_pheno$phenotype %in% fig3a,], 
       aes(x = phenotype, y = cell.count, fill = group))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("Mean cell density / mm2")+
  scale_fill_manual(name = NULL, 
                    breaks = c("High", "Low"), 
                    labels = c("T cell high", "T cell low"),
                    values = c("#DC0018", "#2D69A9"))+
  facet_wrap(~phenotype, scale = "free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank())
