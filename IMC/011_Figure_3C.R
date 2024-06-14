## Script for visualisation of the interaction analysis, comparing the T cell groups from the imaging mass cytometry data  
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)

# Read results from interaction analysis with ImaCyte & the order for the plot
df_comp <- read.delim("./input_files/interaction_analysis_Tcell_groups.tsv")
phenotypes_ordered <- read.delim("./input_files/axis_order_interaction_analysis.tsv")

#### Figure 3C ####
# Plot the tile plot comparing significant interactions between both T cell groups
ggplot(data = df_comp, 
       aes(x = Pheno_interest, y = Pheno_neighbour, fill = Group))+
  scale_x_discrete(limits=phenotypes_ordered$phenotypes)+ 
  scale_y_discrete(limits=phenotypes_ordered$phenotypes)+ 
  geom_tile(color = "white", lwd = 1, linetype = 1)+ 
  scale_fill_manual(name = "Interaction enriched in", 
                    breaks = c("High", "Both", "Low"), 
                    values = c("#DC0018", "grey", "#2D69A9"),
                    na.value = "white")+
  xlab("Phenotype of interest")+
  ylab("Phenotype in neighbourhood")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.text = element_text(color = "black"),
        panel.grid = element_blank())
