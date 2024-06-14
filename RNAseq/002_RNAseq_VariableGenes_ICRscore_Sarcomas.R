## Script for processing of RNA-seq data including all studied sarcoma subtypes 
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(biomaRt)
library(dplyr)

# Set working directory
master.location <- setwd(master.location)

# Load normalised data, design file & DGE results
load("./analysis_files/DESeq2_Normalised_Data_Sarcomas.RData")

#### Annotate ensembl IDs with HGNC symbols using biomaRt ####
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

filters = listFilters(ensembl)

attributes = listAttributes(ensembl)

gene_names_biomart <- getBM(attributes=c("ensembl_gene_id", 
                                         "chromosome_name",
                                         "hgnc_symbol"), 
                            filters = "ensembl_gene_id", 
                            values = row.names(df.log), 
                            mart = ensembl) # use normalised data as input for biomaRt


# Create column with ensembl IDs
df.log$ensembl_gene_id <- row.names(df.log)

# Join the normalised data with the gene names
df <- left_join(df.log, gene_names_biomart, by = "ensembl_gene_id")

# Filter genes from chromosome Y, as these will be highly variable across the cohort
df <- df[!df$chromosome_name == "Y",]

# Remove NAs 
df <- df[!is.na(df$hgnc_symbol),]

# Assign ensembl IDs to genes without HGNC symbols
df[df$hgnc_symbol == "", "hgnc_symbol"] <- df[df$hgnc_symbol == "", "ensembl_gene_id"]

# Assign ensembl IDs to duplicated HGNC symbols
df$dups <- duplicated(df$hgnc_symbol)
df[df$dups == TRUE, "hgnc_symbol"] <- df[df$dups == TRUE, "ensembl_gene_id"]

row.names(df) <- df$hgnc_symbol

#### Calculate the variance per gene among all samples ####
df$var_by_row <- apply(df[, design_sarcs$L_ID], 1, var) # find variance of each gene for 59 samples
df <- df[order(-df$var_by_row),] # sort by variance

# Subset the most variable genes
mat1 = as.matrix(df[1:200, design_sarcs$L_ID])

# Calculate Z-Score per gene
Density.Z = mat1
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat1[j,]-mean(mat1[j,]))/sd(mat1[j,]) # z-score the matrix
}

#### Save Z-Scores for most variable genes ####
save(Density.Z, file ="analysis_files/ZScore_VariableGenes_Sarcomas.RData")

#### ICR score ####
ICR_genes <- read.table("./input_files/ICR_genes.tsv", sep = "\t", header = TRUE)

# Filter normalised data for ICR genes
df.ICR <- df[row.names(df) %in% ICR_genes$Gene.name, design_sarcs$L_ID]

# Calculate ICR score per sample
ICR.score <- c(0,0)
for(i in 1:ncol(df.ICR)){
  ICR.score[i] <- mean(as.numeric(df.ICR[, i]))
}

# Create ICR score data frame
ICR.score.df <- data.frame(L_ID = colnames(df.ICR), ICR.score = ICR.score)

# Join ICR score with design file
ICR.score.df <- left_join(ICR.score.df, 
                          design_sarcs[, c("L_ID", "Diagnosis")], 
                          by = "L_ID")

# Create a factor for the different subtypes, ordered from ICR-high to ICR-low
ICR.score.df$Diagnosis <- factor(ICR.score.df$Diagnosis, 
                                 levels = c("USTS", "Chordoma", "Myxofibrosarcoma", 
                                            "Osteosarcoma", "Chondrosarcoma"))

#### Save ICR scores ####
save(ICR.score.df, file ="analysis_files/ICRscore_Sarcomas.RData")

#### Further processing of differential gene expression results ####
DEG.df <- as.data.frame(res)

# Arrange genes based on significance
DEG.df <- DEG.df %>% arrange(padj)
DEG.df$ensembl_gene_id <- row.names(DEG.df)

# Join DGE results with HGNC symbols
DEG.df <- left_join(DEG.df, df[, c("ensembl_gene_id", "hgnc_symbol")])
DEG.df[is.na(DEG.df$hgnc_symbol), "hgnc_symbol"] <- DEG.df[is.na(DEG.df$hgnc_symbol), "ensembl_gene_id"]

#### Save DGE results ####
save(DEG.df, file = "./analysis_files/DGE_results_Sarcomas.RData")

#### save data frame with HGNC symbols & design file ####
save(df, design_sarcs, file = "./analysis_files/HGNC_and_Design_Sarcomas.RData")
