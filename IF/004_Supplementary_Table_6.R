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
survival <- read.delim("./input_files/survival_data_chordomas.tsv")

# Exclude the following patients for specific analyses
# disease recurrence: margin == no surgery (n = 11)
survival[survival$margin == "no surgery", "status_RFS"] <- NA

# surgical resection margin: no surgery (n = 11)
survival[survival$margin == "no surgery", "margin"] <- NA

# anatomical location: ExtraAxial (n = 1)
survival[survival$location == "ExtraAxial", "location"] <- NA

# timing of radiotherapy: Sandwich (n = 2) & Monotherapy (n = 13)
survival[!grepl(pattern = "juvant", survival$RT_timing), "RT_timing"] <- NA

#### Univariate cox proportional hazard models for the disease-specific survival ####
# Create a vector with all covariates
covariates_DSS <- c("age_at_diagnosis", "max_size", "location", "treatment", 
                    "RT_timing", "margin", "IF.group", "status_RFS")

# Create a formula to loop through all selected covariates
univ_formulas_DSS <- sapply(covariates_DSS,
                        function(x) as.formula(paste('Surv(DSS, status_DSS)~', x)))

# Loop through all cox proportional hazard models
univ_models_DSS <- lapply( univ_formulas_DSS, function(x){coxph(x, data = survival)})

# Adjust the names
names(univ_models_DSS) <- c("Age at diagnosis", "max tumour size (cm)", "Anatomical location (Skull base)",
                        "Treatment (RT alone)", "Timing of RT (Adjuvant)", "Surgical margin (R0)", 
                        "T cell infiltration (High)", "Disease recurrence (No)")

# Summarise the results from the cox proportinal hazard models
univ_results_DSS <- lapply(univ_models_DSS,
                       function(x){ 
                         x <- summary(x)
                       })

#### Save the results from the univariate analysis of the disease-specific survival ####
save(univ_results_DSS, file = "./output_files/coxph_results_univariate_DSS.RData")

#### Cox proportional hazard models for the recurrence-free survival ####
# Create a vector with all covariates
covariates_RFS <- c("age_at_diagnosis", "max_size", "location", "treatment", 
                    "RT_timing", "margin", "IF.group")

# Create a formula to loop through all selected covariates
univ_formulas_RFS <- sapply(covariates_RFS,
                            function(x) as.formula(paste('Surv(RFS, status_RFS)~', x)))

# Loop through all cox proportional hazard models
univ_models_RFS <- lapply( univ_formulas_RFS, function(x){coxph(x, data = survival[survival$treatment != "RT",])})

# Adjust the names
names(univ_models_RFS) <- c("Age at diagnosis", "max tumour size (cm)", "Anatomical location (Skull base)",
                            "Treatment (Surgery)", "Timing of RT (Adjuvant)", "Surgical margin (R0)", 
                            "T cell infiltration (High)")

# Summarise the results from the cox proportinal hazard models
univ_results_RFS <- lapply(univ_models_RFS,
                           function(x){ 
                             x <- summary(x)
                           })

# Anatomical location & surgical margin are significantly associated with recurrence-free survival
# Therefore, these variables are taken along for multivariate analysis

multi_models_RFS <- coxph(Surv(RFS, status_RFS) ~ location + margin, data = survival[survival$treatment != "RT",])
multi_results_RFS <- summary(multi_models_RFS)

#### Save the results from the univariate analysis of the recurrence-free survival ####
save(univ_results_RFS, multi_results_RFS, file = "./output_files/coxph_results_RFS.RData")

#### Univariate cox proportional hazard models for the metastasis-free survival ####
# Create a vector with all covariates
covariates_MFS <- c("age_at_diagnosis", "max_size", "location", "treatment", 
                    "RT_timing", "margin", "IF.group", "status_RFS")

# Create a formula to loop through all selected covariates
univ_formulas_MFS <- sapply(covariates_MFS,
                            function(x) as.formula(paste('Surv(MFS, status_MFS)~', x)))

# Loop through all cox proportional hazard models
univ_models_MFS <- lapply( univ_formulas_MFS, function(x){coxph(x, data = survival)})

# Adjust the names
names(univ_models_MFS) <- c("Age at diagnosis", "max tumour size (cm)", "Anatomical location (Skull base)",
                            "Treatment (RT alone)", "Timing of RT (Adjuvant)", "Surgical margin (R0)", 
                            "T cell infiltration (High)", "Disease recurrence (No)")

# Summarise the results from the cox proportinal hazard models
univ_results_MFS <- lapply(univ_models_MFS,
                           function(x){ 
                             x <- summary(x)
                           })

#### Save the results from the univariate analysis of the metastasis-free survival ####
save(univ_results_MFS, file = "./output_files/coxph_results_univariate_MFS.RData")
