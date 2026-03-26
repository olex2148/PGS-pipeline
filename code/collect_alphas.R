library(data.table)
library(bigsnpr)
library(dplyr)



# List all files in the directory
dir <- "export/gwas_catalog_mle_mle2/raw_models/"
model_files <- list.files(dir, pattern = "raw_models.rds")

# Making the file for alphas:
# Initiate df, loop through each .rds file, do some processing of auto chains, and taking the mean of the alphas from the passed chains

alpha_df <- data.frame(accessionID = character(), alpha = numeric())

for (i in 1:length(model_files)) {
  # Read the model files
  raw_model <- readRDS(paste0(dir, model_files[i]))
  
  # Get the accession ID
  base_name <- gsub("\\.[^.]*$", "", gsub("__raw_models\\.rds$", "", model_files[i])) # Removing both raw_models.rds suffix as well as any .tbl/.txt/.tsv suffix before that
  print(base_name)
  
  # Do some processing on auto chains
  multi_auto <- raw_model$ldpred2 
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
  alpha_auto <- mean(sapply(multi_auto[keep], function(auto) auto$alpha_est))
  
  
  # cbinding accession ID and alpha
  alpha_row <- data.frame(cbind(base_name, round(alpha_auto, 4))) %>%
    rename(accessionID = base_name,
           alpha = V2)
  
  alpha_df <- rbind(alpha_df, alpha_row)
}

write.table(alpha_df,
            file = "export/gwas_catalog_alphas.txt",
            sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)

# Read in to make sure it looks as intended
finished_table <- fread("export/BOR_all_auto.txt")

