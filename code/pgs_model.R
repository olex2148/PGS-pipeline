#' LDpred2-auto and lassosum2 on parsed sumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date
#' 
#' @description This script takes two arguments: An input file of parsed sumstats
#' (output of parser.R) as well as a base name for the output  files - as there are several
#' (models, scores, figures, foelgefil)
#' 
#' @Todo
#'    - Make general for continuous traits

suppressPackageStartupMessages({
  library(bigsnpr)
  library(dplyr)
  library(ggplot2)
  library(testit)
  library(openxlsx)
  library(data.table)
  library(rjson)
})

# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
munged_sumstats <- args[1]
base_path <- args[2]
foelgefil <- args[3]

base_name <- gsub(".rds", "", basename(base_path))

paths <- fromJSON(file = "data/paths.json")

# Loading data -------------------------------------------------------------------------------------------------
# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = paths$hapmap_path, fname = "map_hm3_plus.rds"))

# Reading in iPSYCH data
dosage <- snp_attach(paths$dosage_path) 
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Reading in parsed sumstats
# df_beta = readRDS(paths$test_parsed) # For testing
df_beta = readRDS(munged_sumstats)

# Running LDSC -------------------------------------------------------------------------------------------------------
cat("Running LDSC \n")

ld_size <- nrow(info)
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = ld_size,
                               chi2 = (beta / beta_se)^2, 
                               sample_size = n_eff,
                               ncores = nb_cores()))
h2_init <- ldsc[["h2"]]
cat("LDSC-estimated heritability on the observed scale:", h2_init, "\n")

# Reading in LD blocks -----------------------------------------------------------------------------------------------
corr_dir <- paste0("steps/corr/", base_name, ".rds")

corr <- runonce::save_run({
  for (chr in 1:22) {
    
    cat(chr, "..", sep = "")
    
    ## indices in df_beta
    ind.chr <- which(df_beta$chr == chr)
    
    ## indices in info
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    
    ## indices in corr_chr
    ind.chr3 <- match(ind.chr2, which(info$chr == chr))
    
    corr_chr <- readRDS(paste0(paths$ld_blocks_path, chr, ".rds"))[ind.chr3, ind.chr3]
    
    if (chr == 1) {
      corr <- as_SFBM(corr_chr, corr_dir, compact = TRUE)
    } else {
      corr$add_columns(corr_chr, nrow(corr))
    }
  }
  corr
}, file = corr_dir)

# LDpred2-auto ----------------------------------------------------------------------------------------------------
set.seed(72)
coef_shrink <- 0.9

repeat { 
  cat("Running LDpred2-auto with shrinkage coefficient", coef_shrink, "\n")
  
  multi_auto <- snp_ldpred2_auto(
    corr, df_beta, h2_init = h2_init,
    vec_p_init = seq_log(1e-4, 0.2, length.out = 50), burn_in = 500, num_iter = 500,
    use_MLE = FALSE, # for power/convergence issues, alpha
    report_step = 20, ncores = nb_cores(), allow_jump_sign = FALSE, shrink_corr = coef_shrink)
  
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
  
  perc_kept <- length(keep)/50
  cat(length(keep), "chains passed QC \n")
  
  if(perc_kept >= 0.2) break                       # At least 10 chains should pass QC
  coef_shrink <- coef_shrink - 0.1
  if(coef_shrink < 0.4) break                     # We won't allow a shrinkage coef smaller than 0.4
  cat("Rerunning \n")
}

# Average of kept chains
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Inspecting the first chain to pass QC (arbitrary)
auto <- multi_auto[keep][[1]]

q <- plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) + 
    theme_bigstatsr() + 
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)

ggsave(paste0(base_path, "_1st_kept_chain.jpeg"), q)

# Saving auto parameters -------------------------------------------------------------------------------------------
cat("h2 quantiles \n")
all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
(quant_h2 <- quantile(all_h2, c(0.5, 0.025, 0.975)))

cat("p quantiles \n")
all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
(quant_p <- quantile(all_p, c(0.5, 0.025, 0.975)))

cat("Alpha quantiles \n")
all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
quantile(all_alpha, c(0.5, 0.025, 0.975), na.rm = TRUE)

cat("r2 quantiles \n")
bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
  b1 <- bsamp[[ic]]
  
  Rb1 <- apply(b1, 2, function(x)
    coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
  b2 <- do.call("cbind", bsamp[-ic])
  
  b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
}))
(quant_r2 <- quantile(all_r2, c(0.5, 0.025, 0.975)))

saveRDS(list(r2 = all_r2, h2 = all_h2, alpha = all_alpha, p = all_p),
        paste0(base_path, "_auto_parameters.rds"))

# Lassosum ---------------------------------------------------------------------------------------------------------
cat("Running lassosum \n")
beta_lassosum <- snp_lassosum2(corr, df_beta, ncores = nb_cores())

# Pseudo-validation
params <- attr(beta_lassosum, "grid_param")
scale <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))
beta_hat <- df_beta$beta / scale

pval <- df_beta$p
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

# Saving the raw models                 
saveRDS(list(ldsc = ldsc, ldpred2 = multi_auto, lassosum = beta_lassosum),
        paste0(base_path, "_raw_models.rds"))

# Adding metrics to foelgefil
foelgefil_df <- read.xlsx(xlsxFile = foelgefil)

foelgefil_df$h2_init <- h2_init
foelgefil_df$h2_auto <- quant_h2[1]
foelgefil_df$h2_2.5 <- quant_h2[2]
foelgefil_df$h2_97.5 <- quant_h2[3]
foelgefil_df$p2_auto <- quant_p[1]
foelgefil_df$p_2.5 <- quant_p[2]
foelgefil_df$p_97.5 <- quant_p[3]
  
write.xlsx(foelgefil_df,
           file = foelgefil,
           rownames = FALSE)
  
# Predicting in iPSYCH ------------------------------------------------------------------------------------------------------

# Reading in PCs and info for covariates
pcs <- readRDS(paths$pcs_path)
meta <- fread(paths$meta_path)

# Computing sex and age
covariates_df <- dosage$fam %>% 
  left_join(meta[, c("fdato", "gender", "pid")], by = c("family.ID" = "pid")) %>% 
  mutate(sex = ifelse(gender == "F", 1, 0),
         # Time diff in years between present date and fdate
         age = lubridate::time_length(
           difftime(
             as.Date("01/01/2023", format = "%d/%m/%Y"), 
             as.Date(fdato, format = "%d/%m/%Y")), 
           "years")) %>% 
  select(-c(paternal.ID, maternal.ID, affection, gender))

# Checking order is preserved
# identical(dosage$fam$sample.ID, covariates_df$sample.ID)

cov <- cbind(covariates_df$sex, covariates_df$age, covariates_df$is_2012, pcs)

G <- dosage$genotypes

map_pgs <- df_beta[1:4]; map_pgs$beta <- 1
map_pgs2 <- snp_match(map_pgs, dosage$map)

# Compute scores for all individuals in iPSYCH
pred_auto <- big_prodVec(G,
                         beta_auto, # Model
                         ind.col = map_pgs2[["_NUM_ID_"]], # Indices in G of snps used in auto
                         ncores = nb_cores())

pred_lassosum <- big_prodVec(G,
                             best_lassosum,
                             ind.col = map_pgs2[["_NUM_ID_"]], # Indices in G of snps used in auto
                             ncores = nb_cores())

# cbind with family and sample ID and save
scores <- as.data.frame(cbind(covariates_df$family.ID, 
                              covariates_df$sample.ID, 
                              as.numeric(pred_auto),
                              as.numeric(pred_lassosum)))
colnames(scores) <- c("family.ID", "sample.ID", "ldpred2_pgs", "lassosum_pgs")
saveRDS(scores, paste0(base_path, "_scores.rds"))

cat("\n Finished computing scores. Models, auto model parameters, foelgefil, and auto + lassosum scores were saved in 4 distinct files in", base_path, "/")





