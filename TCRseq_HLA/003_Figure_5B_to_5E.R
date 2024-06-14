## Script for calculating the Shannon diversity & unique number of clones from the TCR sequencing data ()
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(ggplot2)

# Set working directory
master.location <- setwd(master.location)

# Note: L6951 is used for patient P5 and L2038 is used for patient P6.
# L6951 is the primary tumour after radiotherapy
# L2038 is the second recurrence

# Read TCRseq results and HLA class I expression patterns
TCR.df <- read.delim("./input_files/TCRseq_measurements.tsv")
load(file = "./analysis_files/HLA_distributions.RData")

#### Figure 5B ####
# Plot the stacked bar plot with HLA class I expression per T cell group
ggplot()+
  geom_bar(data = HLA_dist, aes(x = IF.group, y = perc, fill = HLA), color = "black", stat = "identity")+
  #coord_flip()+
  ylab("Cumulative percentage of HLA class I status")+
  xlab(NULL)+
  theme_bw()+
  scale_x_discrete(labels = c("High" = "T cell high",
                              "Low" = "T cell low"))+
  scale_fill_manual(breaks = c("Positive", "Weak", "Defective"),
                    name = "HLA class I expression",
                    values = c("#FC8D59", "#FFFFBF", "#91BFDB"))+
  theme(legend.position = "none",
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

#### TCR sequencing results and relation to HLA class I expression ####
# Keep only the first lesion of each patient
TCR.df <- TCR.df[!is.na(TCR.df$IF.group),]

#### Perform student's t-test to compare T cell high with T cell low group ####
# Create an empty vector and data frame for statistical testing
test <- NA
results <- data.frame(variable = 0, p.val = 0, t.test = 0, mean_high = 0, mean_low = 0, median_high = 0, median_low = 0)

# Create a vector with the variables
variables <- c("Clones", "Shannon_diversity")

i=1
# Loop through both measurements
for(i in 1:length(variables)){
  test <- t.test(TCR.df[TCR.df$IF.group == "High", variables[i]], 
                 TCR.df[TCR.df$IF.group == "Low", variables[i]])
  results[i, "variable"] <- variables[i]
  results[i, "p.val"] <- test$p.value
  results[i, "t.test"] <- test$statistic
  results[i, "mean_high"] <- mean(TCR.df[TCR.df$IF.group == "High", 
                                         variables[i]])
  results[i, "mean_low"] <- mean(TCR.df[TCR.df$IF.group == "Low", 
                                        variables[i]])
  results[i, "median_high"] <- median(TCR.df[TCR.df$IF.group == "High", 
                                             variables[i]])
  results[i, "median_low"] <- median(TCR.df[TCR.df$IF.group == "Low", 
                                            variables[i]])
}

#### Save student's t-test results ####
write.table(results, "./output_files/TCRseq_t_test_results.tsv")

#### Figure 5C ####
# Plot the Shannon diversity for both T cell groups
ggplot(data = TCR.df, aes(x = IF.group, y = Shannon_diversity, fill = IF.group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  xlab(NULL)+
  ylab("Shannon diversity")+
  geom_hline(yintercept = 2.5, linetype = "dashed")+
  scale_x_discrete(labels = c("High" = "T cell high",
                              "Low" = "T cell low"))+
  scale_fill_manual(name = NULL, 
                    breaks = c("High", "Low"),
                    values = c("#DC0018", "#2D69A9"))+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

#### Figure 5D ####
# Plot the number of unique TCR clones for both T cell groups
ggplot(data = TCR.df, aes(x = IF.group, y = Clones, fill = IF.group))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(size = 2)+
  xlab(NULL)+
  ylab("Number of unique TCR clones")+
  
  scale_x_discrete(labels = c("High" = "T cell high",
                              "Low" = "T cell low"))+
  scale_fill_manual(name = NULL, 
                    breaks = c("High", "Low"),
                    values = c("#DC0018", "#2D69A9"))+
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")

#### Figure 5E ####
# Plot the box plots with the Shannon diversity per HLA class I expression pattern group
ggplot(data = TCR.df, aes(x = HLA_Class_I, y = Shannon_diversity, fill = HLA_Class_I))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter()+
  xlab(NULL)+
  ylab("Shannon diversity")+
  geom_hline(yintercept = 2.5, linetype = "dashed")+
  scale_x_discrete(limits = c("Positive", "Defective", "Weak"))+
  scale_fill_manual(breaks = c("Positive", "Defective", "Weak"),
                    name = "HLA class I expression",
                    values = c("#FC8D59", "#FFFFBF", "#91BFDB"))+ 
  theme_bw()+
  theme(axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        legend.position = "none")
