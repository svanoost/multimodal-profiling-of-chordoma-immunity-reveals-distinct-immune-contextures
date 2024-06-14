## Script for calculating the distribution of HLA class I among alternative groups based on CD8+ T cell density
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)

# Set working directory
master.location <- setwd(master.location)

# Read the HLA class I expression patterns
HLA <- read.delim("./input_files/survival_data_chordomas.tsv")

# Read T cell count data from the immunofluorescence data
IF_counts <- read.delim("./input_files/mean_Tcell_counts_IF.tsv")
row.names(IF_counts) <- IF_counts$phenotype

# Transpose the data frame for caluclating the median CD8+ T cell density
IF_counts <- as.data.frame(t(IF_counts[,-1]))
colnames(IF_counts) <- c("CD4", "CD8", "Tregs", "KI67")

# Create a column with the CD8 high and low separation based on the median CD8+ T cell count
IF_counts$group <- "High"
IF_counts[IF_counts$CD8 < median(IF_counts$CD8), "group"] <- "Low"

# Merge the HLA class I expression patterns with the CD8+ T cell groups
IF_counts$L_ID <- row.names(IF_counts)
HLA <- left_join(HLA, IF_counts[, c("L_ID", "group")])

# Count the different expression patterns per CD8+ T cell group
HLA_dist <- HLA %>% 
  group_by(group) %>%
  dplyr::count(HLA)

# Remove all samples either without the necessary data
HLA_dist <- HLA_dist[!is.na(HLA_dist$group),]
HLA_dist <- HLA_dist[!is.na(HLA_dist$HLA),]

# Calculate the percentage of the distinct patterns per group
HLA_dist$perc <- 0

# Calculate the total number in the CD8+ high group
sum_high <- sum(HLA_dist[HLA_dist$group == "High", "n"])

# Calculate the percentage with "defective" HLA class I
HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Defective", "perc"] <- 
  (HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Defective", "n"]*100) / sum_high

# Calculate the percentage with "positive" HLA class I
HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Positive", "perc"] <- 
  (HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Positive", "n"]*100) / sum_high

# Calculate the percentage with "positive" HLA class I
HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Weak", "perc"] <- 
  (HLA_dist[HLA_dist$group == "High" & HLA_dist$HLA == "Weak", "n"]*100) / sum_high

# Calculate the total number in the CD8+ low group
sum_low <- sum(HLA_dist[HLA_dist$group == "Low", "n"])

# Calculate the percentage with "defective" HLA class I
HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Defective", "perc"] <- 
  (HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Defective", "n"]*100) / sum_low

# Calculate the percentage with "positive" HLA class I
HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Positive", "perc"] <- 
  (HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Positive", "n"]*100) / sum_low

# Calculate the percentage with "weak" HLA class I
HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Weak", "perc"] <- 
  (HLA_dist[HLA_dist$group == "Low" & HLA_dist$HLA == "Weak", "n"]*100) / sum_low

#### Supplementary Figure 6D ####
# Plot the stacked bar plot with HLA class I expression per CD8+ T cell group
ggplot()+
  geom_bar(data = HLA_dist, aes(x = group, y = perc, fill = HLA), color = "black", stat = "identity")+
  #coord_flip()+
  ylab("Cumulative percentage of HLA class I status")+
  xlab(NULL)+
  theme_bw()+
  scale_x_discrete(labels = c("High" = "CD8 high",
                              "Low" = "CD8 low"))+
  scale_fill_manual(breaks = c("Positive", "Weak", "Defective"),
                    name = "HLA class I expression",
                    values = c("#FC8D59", "#FFFFBF", "#91BFDB"))+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
