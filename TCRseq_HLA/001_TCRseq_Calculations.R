## Script for calculating the Shannon diversity & unique number of clones from the TCR sequencing data ()
# R version 4.0.2

# Note: raw data is available upon reasonable request

#### Set up environment ####
rm(list = ls())

# Set working directory
master.location <- setwd(master.location)

# Create a vector with path to the TCR clones summaries
clone_sum <- "./clone_summaries"

# Create a vector with all the clone summary files
files <- dir(clone_sum)

# create an empty list to load all the separate clone summaries
TCR_list <- list()

# Loop through all clone summaries and remove clones with less then 10 Plus AND Minus counts
for(i in 1:length(files)){
  TCR_list[[i]] <- read.csv(paste0(clone_sum, "/", files[i]))
  TCR_list[[i]] <- TCR_list[[i]] %>% 
    filter(Plus.Counts >= 10 & Minus.Counts >= 10)
  TCR_list[[i]]$updated_Frequency <- 0
}

# Loop through the clone summary file names to retrieve the sample IDs
L_ID <- c(0,0)
for(i in 1:length(files)){
  L_ID[i] <- substr(files[i], start = 14, stop = 18)
}

# Assigne the sample IDs to the list with the clone summaries and create a new data frame
names(TCR_list) <- L_ID
TCR.df <- as.data.frame(L_ID)

# Create columns for the number of unique clones and the Shannon diversity
TCR.df$Clones <- 0
TCR.df$Shannon_diversity <- 0

# Loop through the list of clone summaries and calculate the Shannon diversity
# Save the results in the data frame
for(j in 1:length(TCR_list)){
  for(i in 1:nrow(TCR_list[[j]])){
    TCR_list[[j]][i, "updated_Frequency"] <- TCR_list[[j]][i, "updated_Frequency"]/sum(TCR_list[[j]]$Total.Counts)
  }
  TCR_list[[j]]$ln <- TCR_list[[j]]$updated_Frequency*log(TCR_list[[j]]$updated_Frequency)
  TCR.df[j, c("Clones", "Shannon_diversity")] <- c(nrow(TCR_list[[j]]), sum(TCR_list[[j]]$ln)*-1)
}

#### Save the results from the TCR sequencing ####
write.table(TCR.df, file = "./input_files/TCRseq_measurements.tsv", sep = "\t")
