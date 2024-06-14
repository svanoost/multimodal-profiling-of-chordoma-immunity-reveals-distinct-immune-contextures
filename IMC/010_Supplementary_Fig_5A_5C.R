## Script for visualising the differences between the T cell high and the T cell low group from the imaging mass cytometry data
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)

# Set working directory
master.location <- setwd(master.location)

# Load significant phenotypes data and overall cell density data
load(file = "./analysis_files/IMC_Tcell_group_comparison_counts.RData")
load(file = "./analysis_files/IMC_overall_density_results.RData")

#### Supplementary Figure 5A ####
# Create a vector with the phenotypes from Figure 3A
fig3a <- c("Dendritic cells", "Immune cell aggregates")

# Plot box plots for the remaining significant phenotypes by excluding the ones from Figure 3A
ggplot(data = melt_pheno[!melt_pheno$phenotype %in% fig3a,], 
       aes(x = phenotype, y = cell.count, fill = group))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("Mean cell density / mm2")+
  scale_fill_manual(name = "T cell infiltration", 
                    breaks = c("High", "Low"),
                    values = c("#DC0018", "#2D69A9"))+
  facet_wrap(~phenotype, scale = "free",
             ncol = 2)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        strip.background = element_blank(),
        strip.text = element_blank())

#### Supplementary Figure 5C ####
# Plot box plots with overall cell densities, comparing T cell high with T cell low
ggplot(data = melt_overall, aes(x = cell_group, y = cell.count.sum, fill = group))+
  geom_boxplot()+
  xlab(NULL)+
  ylab("Mean cell density / mm2")+
  scale_x_discrete(labels = c("Immune" = "Immune cells",
                              "Stroma" = "Stromal cells",
                              "Tumour" = "Tumour cells"))+
  scale_fill_manual(name = "T cell infiltration", 
                    breaks = c("High", "Low"),
                    values = c("#DC0018", "#2D69A9"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(color="black"),
        legend.position = "none")
