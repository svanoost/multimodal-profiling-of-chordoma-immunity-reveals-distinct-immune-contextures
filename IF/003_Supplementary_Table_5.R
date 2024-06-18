## Script for the prognostic analysis for the T cell phenotypes from the immunofluorescence data
# R version 4.0.2

# Note: follow up data is available upon reasonable request

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(survminer)
library(survival)

# Set working directory
master.location <- setwd(master.location)

# Read survival data
survival <- read.delim("./input_files/survival_data_chordomas.tsv")

# Load the T cell phenotypes from the immunofluorescence data, with separate Ki-67+ subtypes
IF.counts <- read.delim("./input_files/mean_Tcell_counts_IF_Ki67_separate.tsv")
row.names(IF.counts) <- IF.counts$phenotype

# Merge the T cell phenotype counts with the survival data
IF.counts <- as.data.frame(t(IF.counts[, -1]))
IF.counts$L_ID <- row.names(IF.counts)
colnames(IF.counts)[c(2, 5, 6)] <- c("Ki67_CD4_Tcells", 
                                               "Ki67_CD8_Tcells",
                                               "Ki67_Tregs")
survival <- left_join(survival, IF.counts)

#### Univariate cox proportional hazard models for the disease-specific survival ####
# Create a vector with all covariates
covariates_DSS <- c("CD4_Tcells", "Ki67_CD4_Tcells", "CD8_Tcells", "Ki67_CD8_Tcells", 
                    "Tregs", "Ki67_Tregs")

# Create a formula to loop through all selected covariates
univ_formulas_DSS <- sapply(covariates_DSS,
                            function(x) as.formula(paste('Surv(DSS, status_DSS)~', x)))

# Loop through all cox proportional hazard models
univ_models_DSS <- lapply( univ_formulas_DSS, function(x){coxph(x, data = survival)})

# Adjust the names
names(univ_models_DSS) <- c("CD4+ T cells", "Ki-67+CD4+ T cells", "CD8+ T cells",
                            "Ki-67+CD8+ T cells", "Regulatory T cells", "Ki-67+ Regulatory T cells")

# Summarise the results from the cox proportinal hazard models
univ_results_DSS <- lapply(univ_models_DSS,
                           function(x){ 
                             x <- summary(x)
                           })

#### Save the results from the univariate analysis of the disease-specific survival ####
save(univ_results_DSS, file = "./output_files/coxph_results_univariate_DSS_Tcells.RData")