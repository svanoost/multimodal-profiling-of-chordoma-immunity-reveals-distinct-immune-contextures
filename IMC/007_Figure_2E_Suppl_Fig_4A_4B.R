### Script for visualisation of the T cell groups in relation to survival in chordomas
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survminer)
library(survival)

# Read survival data
Tcell_survival <- read.delim("./input_files/survival_data_chordomas.tsv")

# Calculate the disease-specific survival for the chordomas with ICR score, comparing ICR-high with ICR-low
fit <- survfit(Surv(DSS, status_DSS) ~ IMC.group, data = Tcell_survival)

# Summarise the survival data
summary(fit)
fit

#### Figure 2E ####
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

#### Supplementary Figure 4A ####
# Calculate the disease-specific survival for the chordomas with ICR score, comparing ICR-high with ICR-low
fit <- survfit(Surv(RFS, status_RFS) ~ IMC.group, data = Tcell_survival)

# Summarise the survival data
summary(fit)
fit

#### Figure 2E ####
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

#### Supplementary Figure 4B ####
# Calculate the disease-specific survival for the chordomas with ICR score, comparing ICR-high with ICR-low
fit <- survfit(Surv(MFS, status_MFS) ~ IMC.group, data = Tcell_survival)

# Summarise the survival data
summary(fit)
fit

#### Figure 2E ####
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
