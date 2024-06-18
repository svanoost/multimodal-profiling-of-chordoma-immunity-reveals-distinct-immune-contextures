## Script for normalisation & processing of RNA-seq data of chordomas
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(DESeq2)

# Set working directory
master.location <- setwd(master.location)

# Read count table data
gene_expr <- read.delim("./input_files/raw_counts_chordoma.tsv") # This data can be found on GEO
row.names(gene_expr) <- gene_expr$feature

# Read design file
design_chord <- read.delim("./input_files/design_file_sarcomas.tsv")
row.names(design_chord) <- design_chord$L_ID
design_chord <- design_chord[design_chord$Diagnosis == "Chordoma",]

# Add an ICR score column
design_chord$ICR.score <- "High"
ICR.low <- c("L3382", "L6473", "L6617")
design_chord[design_chord$L_ID %in% ICR.low, "ICR.score"] <- "Low"
rm(ICR.low)

# Check whether sample IDs correspond between the count table and design file
gene_expr <- gene_expr[, design_chord$L_ID]
all(colnames(gene_expr) == design_chord$L_ID)

# Filter genes 
length(which(apply(gene_expr,1, function(x){any(x > 20)}))) # keep 36716
idx <- which(apply(gene_expr,1, function(x){any(x > 20)}))
gene_expr_filtered <- gene_expr[idx,]

#### Create DEseq2 object, normalise and filter the count table ####
dds <- DESeqDataSetFromMatrix(countData = gene_expr_filtered, colData = design_chord, design = ~ICR.score)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

counttable.norm.dds <- counts(dds, normalized=TRUE)
df.log <- as.data.frame(log2(counttable.norm.dds+1))

#### ICR score among chordomas ####
# Load ICR genes
ICR_genes <- read.table("./input_files/ICR_genes.tsv", sep = "\t", header = TRUE)

# Filter the normalised data for the ICR genes by using the ensembl gene ID
mat1 <- as.matrix(df.log[row.names(df.log) %in% ICR_genes$Gene.stable.ID, ])

# Calculate the Z-Score for the ICR genes
Density.Z = mat1
for(j in 1: nrow(Density.Z))  {
  Density.Z[j,] = (mat1[j,]-mean(mat1[j,]))/sd(mat1[j,]) # z-score the matrix
}

#### Save the Z-Score of the ICR genes together with the design file ####
save(Density.Z, design_chord, ICR_genes, file = "./analysis_files/DESeq2_Normalised_Data_Chordomas.RData")
