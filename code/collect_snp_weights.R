library(data.table)
library(bigsnpr)
library(dplyr)



# List all files in the directory
model_dir <- "export/gwas_catalog_mle/raw_models/"
munged_dir <- "export/gwas_catalog_mle/munged_sumstats/"

model_files <- list.files(model_dir, pattern = "raw_models.rds")
munged_sumstats <- list.files(munged_dir, pattern = "munged.rds") #for matching variants later, in model only weights sorted

# HapMap3+ - iPSYCH overlap
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = "data/hapmap3+", fname = "map_hm3_plus.rds")) %>% 
  mutate(beta = 0)

dosage <- snp_attach("~/") 
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Creating df with the overlap between the two datasets
info_snp <- snp_match(info, dosage$map)

# hapmap3plus_ipsych is the df we'll left_join all snp weights onto

# Making the file for auto betas:
# Initiate df, loop through each .rds file, do some processing of auto chains, and left join with hapmap3plus_ipsych

hapmap3plus_ipsych <- info_snp %>% select(rsid) 
for (i in 1:length(model_files)) {
  # Read the both model and munged sumstats files
  raw_model <- readRDS(paste0(model_dir, model_files[i]))
  df_beta = readRDS(paste0(munged_dir, munged_sumstats[i]))
  
  base_name <- gsub("\\.[^.]*$", "", gsub("_raw_models\\.rds$", "", model_files[i])) # Removing both raw_models.rds suffix as well as any .tbl/.txt/.tsv suffix before that
  print(base_name)
  
  # Do some processing on auto chains
  multi_auto <- raw_model$ldpred2 
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
  beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
  
  # Finding indices in iPSYCH of variants used in the current PGS model
  df_beta <- df_beta[, c("chr", "pos", "a0", "a1")]; df_beta$beta <- 1
  map_pgs <- snp_match(df_beta, dosage$map)
  
  idx_in_ipsych <- map_pgs[["_NUM_ID_"]]
  idx_in_pgs <- map_pgs[["_NUM_ID_.ss"]] # Since the above snp_match can remove a few variants, we also need to specify which indices to keep from the model
  
  # cbinding rsid with snp weights and then making it a dataframe
  snp_weights <- data.frame(cbind(dosage$map$SNP[idx_in_ipsych], round(beta_auto[idx_in_pgs], 7))) %>% # Lesser decimals yields a lot of zeroes
    rename(rsid = X1,
           !!base_name := X2)
  
  # Left join the snp weights with snp_list
  hapmap3plus_ipsych <- hapmap3plus_ipsych %>% 
    left_join(snp_weights, by = "rsid")
}

write.table(hapmap3plus_ipsych,
            file = "export/BOR_all_auto.txt",
            sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)

# Read in to make sure it looks as intended
finished_table <- fread("export/BOR_all_auto.txt")


# Doing the same for lassosum

# Loop through each .rds file, do some processing of auto chains, and left join with snp_list
hapmap3plus_ipsych_lasso <- info_snp %>% select(rsid, chr, pos, a0, a1) # Do we only need rsid?
for (i in 1:length(model_files)) {
  # Read the both model and munged sumstats files
  raw_model <- readRDS(paste0("results/alisha/daner_bor15_all_noDanish_24b/", model_files[i]))
  df_beta = readRDS(paste0("steps/munged_sumstats/alisha/", munged_sumstats[i]))
  
  base_name <- gsub("\\.[^.]*$", "", gsub("_raw_models\\.rds$", "", model_files[i]))
  print(base_name)
  
  # Renaming the n column to n_eff to generalise code
  if("n" %in% colnames(df_beta) & !"n_eff" %in% colnames(df_beta)){df_beta$n_eff = df_beta$n}
  
  # First creating ld cor
  tmp <- tempfile(tmpdir = "steps/corr", pattern = model_files[i])
  cat("Reading in LD blocks for chromosome ")
  for (chr in 1:22) {
    
    cat(chr, "..", sep = "")
    
    ## indices in df_beta
    ind.chr <- which(df_beta$chr == chr)
    ## indices in info
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    ## indices in corr_chr
    ind.chr3 <- match(ind.chr2, which(info$chr == chr))
    
    corr_chr <- readRDS(paste0("", chr, ".rds"))[ind.chr3, ind.chr3]
    
    if (chr == 1) {
      corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  cat("\n")
  
  # Pseudo-validation
  beta_lassosum <- raw_model$lassosum 
  params <- attr(beta_lassosum, "grid_param")
  scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
  beta_hat <- df_beta$beta / scale
  
  if("p" %in% colnames(df_beta)) {
    df_beta$p <- as.numeric(df_beta$p)
  }
  
  fdr <- fdrtool::fdrtool(df_beta$p, statistic = "pvalue", plot = FALSE)
  beta_hat_shrunk <- beta_hat * (1 - fdr$lfdr)
  
  params$auto_score <- apply(beta_lassosum, 2, function(beta) {
    cat(".")
    beta <- beta / scale
    bRb <- crossprod(beta, bigsparser::sp_prodVec(corr, beta))
    crossprod(beta, beta_hat_shrunk) / sqrt(bRb)
  })
  
  # Choosing best lassosum model
  best_lassosum <- params %>%
    mutate(id = row_number()) %>%
    arrange(desc(auto_score)) %>%
    slice(1) %>%
    pull(id) %>% 
    beta_lassosum[, .]
  
  # Finding indices in iPSYCH of variants used in the current PGS model
  df_beta <- df_beta[, c("chr", "pos", "a0", "a1")]; df_beta$beta <- 1
  map_pgs <- snp_match(df_beta, dosage$map)
  
  idx_in_ipsych <- map_pgs[["_NUM_ID_"]]
  idx_in_pgs <- map_pgs[["_NUM_ID_.ss"]]
  
  snp_weights <- data.frame(cbind(dosage$map$SNP[idx_in_ipsych], round(best_lassosum[idx_in_pgs], 5))) %>% # Lesser decimals yields a lot of zeroes
    rename(rsid = X1,
           !!base_name := X2)
  
  # Left join the snp weights with snp_list
  hapmap3plus_ipsych_lasso <- hapmap3plus_ipsych_lasso %>% 
    left_join(snp_weights, by = "rsid")
}

write.table(hapmap3plus_ipsych_lasso,
            file = "export/BOR_all_lassosum.txt",
            sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)

