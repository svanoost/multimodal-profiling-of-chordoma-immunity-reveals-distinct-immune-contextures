## Script for normalisation & processing of RNA-seq data including all studied sarcoma subtypes 
# R version 4.0.2

#### Set up environment ####
rm(list = ls())

# List of required libraries
library(DESeq2)

# Set working directory
master.location <- setwd(master.location)

# Read count table data
gene_expr <- read.delim("./input_files/raw_counts_sarcomas.tsv") # osteosarcoma & chordoma data can be found on GEO

# Read design file
design_sarcs <- read.delim("./input_files/design_file_sarcomas.tsv")
row.names(design_sarcs) <- design_sarcs$L_ID

# Check whether sample IDs correspond between the count table and design file
all(colnames(gene_expr) == design_sarcs$L_ID)

# Filter genes 
length(which(apply(gene_expr,1, function(x){any(x > 20)})))
idx <- which(apply(gene_expr,1, function(x){any(x > 20)}))
gene_expr_filtered <- gene_expr[idx,]

#### Create DEseq2 object, normalise and filter the count table ####
dds <- DESeqDataSetFromMatrix(countData = gene_expr_filtered, colData = design_sarcs, design = ~vs)
dds <- estimateSizeFactors(dds)
idx <- rowSums(counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)

counttable.norm.dds <- counts(dds, normalized=TRUE)
df.log <- as.data.frame(log2(counttable.norm.dds+1))

#### Differential gene expression analysis ####
res <- results(dds, contrast=c("vs", "Chordoma","Others"))

#### Save the normalised data together with the design file & DGE results ####
save(df.log, design_sarcs, res, file = "./analysis_files/DESeq2_Normalised_Data_Sarcomas.RData")
