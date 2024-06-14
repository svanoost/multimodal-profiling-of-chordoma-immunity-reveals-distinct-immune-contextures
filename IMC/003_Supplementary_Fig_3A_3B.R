## Script for visualisation of the T cell phenotypes from the imaging mass cytometry data in relation to the ICR score   
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)

# Set working directory
master.location <- setwd(master.location)

# Load the cell phenotype data & sample annotation files for primary tumours only 
load(file = "./analysis_files/IMC_Data_Primary_Tumours.RData")

# Filter samples with an ICR score & select T cell phenotypes
df <- as.data.frame(mat[grepl("T cells", row.names(mat)), 
                        colnames(mat) %in% anno_samples[!is.na(anno_samples$ICR),"L_ID"]])
df$phenotype <- row.names(df)

# Melt the data frame for ggplot visualisation
melt_pheno <- melt(df, value.name = "cell.count", variable.name = "L_ID")

# Join the phenotype counts with the sample annotation
melt_pheno <- left_join(melt_pheno, anno_samples, by = "L_ID")

#### Student's t-test comparing T cell phenotypes between ICR-high and ICR-low ####
# Create a vector with the T cell phenotypes
my_phenos <- unique(melt_pheno$phenotype)

# Create an empty vector and data frame for statistical testing
test <- NA
results <- data.frame(phenotype = 0, p.val = 0, t.test = 0, mean_high = 0, mean_low = 0, median_high = 0, median_low = 0)

# Loop through the T cell phenotypes and perform student's t-test
for(i in 1:length(my_phenos)){
  test <- t.test(melt_pheno[melt_pheno$ICR == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"], 
                 melt_pheno[melt_pheno$ICR == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "phenotype"] <- my_phenos[i]
  results[i, "p.val"] <- test$p.value
  results[i, "t.test"] <- test$statistic
  results[i, "mean_high"] <- mean(melt_pheno[melt_pheno$ICR == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "mean_low"] <- mean(melt_pheno[melt_pheno$ICR == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "median_high"] <- median(melt_pheno[melt_pheno$ICR == "High" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
  results[i, "median_low"] <- median(melt_pheno[melt_pheno$ICR == "Low" & melt_pheno$phenotype == my_phenos[i], "cell.count"])
}

#### Save results from the student's t-test ####
write.table(results, "./output_files/IMC_ICR_t_test_results.tsv", sep = "\t", row.names = FALSE)

#### Supplementary Figure 3A ####
# Make the phenotypes into a factor for visualisation
melt_pheno$phenotype <- factor(melt_pheno$phenotype, 
                               levels = c("T cells", 'CD4+ T cells', "CD8+ T cells", 
                                          "Ki-67+ T cells", "Regulatory T cells"))

# Plot the box plots comparing ICR-high with ICR-low samples
ggplot(data = melt_pheno, aes(x = ICR.score, y = cell.count))+
  geom_boxplot(aes(fill = ICR.score), outlier.shape = NA)+
  geom_jitter(width = 0.2)+
  xlab("T cell phenotypes")+
  ylab("Mean cell density / mm2")+
  scale_x_discrete(labels = NULL)+
  scale_fill_manual(name = "ICR score", 
                    breaks = c("High", "Low"), 
                    values = c("#E7B800", "#2E9FDF"))+
  facet_wrap(~phenotype, scale = "free",
             strip.position = "bottom")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text = element_text(colour = "black"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank())

#### Supplementary Figure 3B ####
# Create a data frame with CD4+ and CD8+ T cells for a correlation plot
cor_df <- as.data.frame(t(mat[c("CD4+ T cells", "CD8+ T cells"),]))
colnames(cor_df) <- c("CD4", "CD8")

# Plot the correlation between CD4+ and CD8+ T cells
ggscatter(cor_df, x = "CD8", y = "CD4", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "CD8+ T cells", ylab = "CD4+ T cells")
