### Script for visualisation of the ICR score in relation to survival in chordomas
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survminer)
library(survival)

# Set working directory
master.location <- setwd(master.location)

# Read survival data
ICR_survival <- read.delim("./input_files/survival_data_chordomas.tsv")

# Calculate the disease-specific survival for the chordomas with ICR score, comparing ICR-high with ICR-low
fit <- survfit(Surv(DSS, status_DSS) ~ ICR.score, data = ICR_survival)

# Summarise the survival data
summary(fit)
fit

#### Figure 1D ####
# Plot the Kaplan-Meier curves for the ICR-high and ICR-low chordomas
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("ICR high", "ICR low"), # label the groups
           font.legend = c(12),
           risk.table.font = c(5),
           legend.title = "",
           xlim = c(0, 150), # use 150 months as max
           break.x.by = 30, # use 30 months as breaks
           palette = c("#E7B800", "#2E9FDF")) +
  xlab(label = "Time (months)")
