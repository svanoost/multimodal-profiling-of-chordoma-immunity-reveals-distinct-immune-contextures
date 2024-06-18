## Script for visualisation of the T cell groups  from the immunofluorescence data in relation to survival
# R version 4.0.2

# Note: follow up data is available upon reasonable request

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survminer)
library(survival)

# Set working directory
master.location <- setwd(master.location)

# # Read survival data
Tcell_survival <- read.delim("./input_files/survival_data_chordomas.tsv")

#### Disease-Specific Survival ####
# Calculate the disease-specific survival for both T cell groups
fit <- survfit(Surv(DSS, status_DSS) ~ IF.group, data = Tcell_survival)

# Summarise the survival data
summary(fit)
fit

#### Figure 4B ####
# Plot the Kaplan-Meier curves for the T cell high and T cell low chordomas
ggsurvplot(fit,
           pval = TRUE, 
           conf.int = FALSE,
           risk.table = TRUE,
           risk.table.y.text = FALSE,
           legend.labs = c("T cell high", "T cell low"), # label the groups
           font.legend = c(12),
           risk.table.font = c(5),
           legend.title = "",
           xlim = c(0, 150), # use 150 months as max
           break.x.by = 30, # use 30 months as breaks
           palette = c("#DC0018", "#2D69A9"))+
  xlab(label = "Time (months)")
