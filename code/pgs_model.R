#' LDpred2-auto and lassosum2 on parsed sumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date
#' 
#' @description This script takes two arguments: An input file of parsed sumstats
#' (output of parser.R) as well as a base name for the output  files - as there are several
#' (models, scores, figures, følgefil)
#' 
#' @Todo
#'    - Make general for continuous traits as well
#'    - Move n_eff check to parser

suppressPackageStartupMessages({
  library(bigsnpr)
  library(dplyr)
  library(ggplot2)
  library(testit)
  library(openxlsx)
  library(data.table)
})

# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
parsed_sumstats <- args[1]
base_path <- args[2]
base_name <- sapply(strsplit(base_path, split='/', fixed=TRUE), function(x) (x[3])) # Removing the path in front of base name

source("code/aux_scripts/input_paths.R")
set.seed(72)

# Loading data -------------------------------------------------------------------------------------------------

# Reading in parsed sumstats
sumstats = readRDS(parsed_sumstats)
# sumstats = readRDS("steps/parsed_sumstats/test_adhd.rds") # Example

# Reading in iPSYCH data
dosage <-  snp_attach(dosage_path)
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Reading in HapMap3+
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = hapmap_path, fname = "map_hm3_plus.rds"))

# Finding Hapmap and iPSYCH overlap with sumstats --------------------------------------------------------------

info_ipsych_external <- snp_match(
  sumstats, 
  dosage$map)                                                                                      # iPSYCH-external

info_snp <- snp_match(
  info_ipsych_external[, c("pos","chr","a0","a1","beta","beta_se","freq","p","n_eff", "info")], 
  info)                                                                                           # Hapmap overlap with remaining

# QC -----------------------------------------------------------------------------------------------------------

# First checking that the effective population size is in the parsed sumstats
assert("No effective population size in parsed sumstats",
       !all(is.na(info_snp$n_eff)))

# Filtering on beta_se and n_eff which should be present in all sumstats
info_snp2 <- info_snp %>% 
  filter(beta_se > 0 & n_eff > (0.7 * max(n_eff)))

# Other QC steps that cannot be performed in all sumstats (if they don't contain freq and/or info)
is_bad <- NA
sd_af <- NA

if( !all(is.na(info_snp2$freq)) ){         # If freq exists
  
  sd_af <- with(info_snp2, sqrt(2 * freq * (1 - freq)))
  sd_ss <- with(info_snp2, 2 / sqrt(n_eff * beta_se^2 + beta^2))
  sd_ss2 <- sd_ss / quantile(sd_ss, 0.999) * sqrt(0.5) 
  
  is_bad <- 
    sd_ss2 < (0.5 * sd_af) |
    sd_ss2 > (sd_af + 0.1) |
    sd_ss2 < 0.05 |
    sd_af < 0.05
  
} else {
  cat("No allele frequencies available in summary statistics. QC step not performed. \n")
}

if( !all(is.na(info_snp2$info)) ){       # If info exists
  is_bad <- is_bad | 
    info_snp2$info < 0.7
} else {
  cat("No INFO scores available in summary statistics. QC step not performed. \n")
}

cat("Number of variants to be filtered out:", sum(is_bad) + (nrow(info_snp) - nrow(info_snp2)), "\n")

if( all(is.na(is_bad)) ){               # if neither info nor freq exist
  df_beta <- info_snp2
} else {
  df_beta <- info_snp2[!is_bad, ]       # If one of them does
}

# Saving QC in plot if possible
if( length(sd_af) > 1 ){
  
  p <- ggplot(slice_sample(data.frame(sd_af, sd_ss2, is_bad), n = 50e4)) +
    geom_point(aes(sd_af, sd_ss2, color = is_bad), alpha = 0.5) +
    theme_bigstatsr(0.9) + 
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red", linewidth = 1.5) +
    labs(x = "Standard deviations in the reference set",
         y = "Standard deviations derived from the summary statistics",
         color = "To remove?")
  
  ggsave(paste0(base_path, "_QC.jpeg"), p)
}

# Running LDSC -------------------------------------------------------------------------------------------------------
cat("Running LDSC \n")
ldsc <- with(df_beta, snp_ldsc(ld, ld_size = length(ld),
                               chi2 = (beta / beta_se)^2, 
                               sample_size = n_eff,
                               ncores = nb_cores()))
h2_init <- ldsc[["h2"]]
cat("LDSC-estimated heritability on the observed scale:", h2_init, "\n")

# Reading in LD blocks -----------------------------------------------------------------------------------------------
tmp <- tempfile(tmpdir = "steps/corr")

cat("Reading in LD blocks for chromosome ")
for (chr in 1:22) {
  
  cat(chr, "..", sep = "")
  
  ## indices in df_beta
  ind.chr <- which(df_beta$chr == chr)
  ## indices in map_ldref
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in corr_chr
  ind.chr3 <- match(ind.chr2, which(info$chr == chr))
  
  corr_chr <- readRDS(paste0(ld_blocks, chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}
cat("\n")

# LDpred2-auto ------------------------------------------------------------------------------------------------------
coef_shrink <- 0.8

cat("Running LDpred2-auto with shrinkage coefficient", coef_shrink, "\n")
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_init,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 50), burn_in = 500, num_iter = 500,
  use_MLE = FALSE, # for power/convergence issues, alpha
  report_step = 20, ncores = nb_cores(), allow_jump_sign = FALSE, shrink_corr = coef_shrink)

cat("Running lassosum \n")
lasso <- snp_lassosum2(corr, df_beta, ncores = nb_cores())

saveRDS(list(ldsc = ldsc, ldpred2 = multi_auto, lassosum = lasso),
        paste0(base_path, "_raw_models.rds"))

# QC on chains from auto --------------------------------------------------------------------------------------------

range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- (range > (0.95 * quantile(range, 0.95)))

# Making sure some chains were kept
assert("No chains passed QC", sum(keep) != 0)

cat(sum(keep), "chains passed QC \n")

# Average of kept chains
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

# Inspecting the first chain to pass QC (arbitrary)
auto <- multi_auto[keep][[1]]

q <- plot_grid(
  qplot(y = auto$path_p_est) + 
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() + 
    labs(y = "p") + 
    qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
ggsave(paste0(base_path, "_1st_kept_chain.jpeg"), q)

# Saving parameters ------------------------------------------------------------------------------------------------

all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
quantile(all_h2, c(0.5, 0.025, 0.975))

all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
quantile(all_p, c(0.5, 0.025, 0.975))

# MLE
all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
quantile(all_alpha, c(0.5, 0.025, 0.975))

bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
  b1 <- bsamp[[ic]]
  
  Rb1 <- apply(b1, 2, function(x)
    coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
  b2 <- do.call("cbind", bsamp[-ic])
  
  b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
}))

quantile(all_r2, c(0.5, 0.025, 0.975))

saveRDS(list(r2 = all_r2, h2 = all_h2, alpha = all_alpha, p = all_p),
        paste0(base_path, "auto_parameters.rds"))

# Getting mean h2 and p from chains and saving følgefil
h2_mean <- mean(all_h2)
h2_se <- sd(all_h2) / sqrt(sum(keep)) # sd / sqrt(n)

p_mean <- mean(all_p)
p_se <- sd(all_p) / sqrt(sum(keep)) # sd / sqrt(n)

foelgefil <- data.frame(
  id = base_name,
  restrictions = NA,
  reported_trait = NA,
  pubmed_id = NA,
  first_author = NA,
  journal = NA,
  title = NA,
  publication_date = NA,
  N = NA,
  Ncase = NA,
  Ncontrol = NA,
  se = "beta",
  M_or = nrow(sumstats), # Variants in the original sumstats
  M_m = nrow(info_snp),  # Overlap with HapMap and iPSYCH
  M_qc = nrow(df_beta),  # Variants after QC
  m_h2 = h2_mean,        # Ldpred2 h2 mean from kept chains
  se_he = h2_se,         # h2 se
  m_p = p_mean,          # LDpred2 polygenicity mean from kept chains
  se_p = p_se,           # p se
  h2_init = h2_init      # LDSC h2
  
)

write.xlsx(foelgefil, 
           file = paste0("results/følgefiler/", base_name, "_følgefil.xlsx"), 
           rownames = FALSE)

# Predicting in iPSYCH ------------------------------------------------------------------------------------------------------

# Reading in PCs and info for covariates
pcs <- readRDS(pcs_path)
meta <- fread(meta_path)

# Computing sex and age
covariates_df <- dosage$fam %>% 
  left_join(meta[, c("fdato", "gender", "pid")], by = c("family.ID" = "pid")) %>% 
  mutate(sex = ifelse(gender == "F", 1, 0),
         # Time diff in years between present date and fdate
         age = lubridate::time_length(
           difftime(
             as.Date(Sys.Date(), format = "%d/%m/%Y"), 
             as.Date(fdato, format = "%d/%m/%Y")), 
           "years")) %>% 
  select(-c(paternal.ID, maternal.ID, affection, gender))

# Checking order is preserved
# identical(dosage$fam$sample.ID, covariates_df$sample.ID)

cov <- cbind(covariates_df$sex, covariates_df$age, covariates_df$is_2012, pcs)

G <- dosage$genotypes

# Finding indices of variants in G which is used in models
ipsych_sumstats_index <- snp_match(
  df_beta[, c("pos","chr","a0","a1","beta","beta_se","freq","p","n_eff", "info")], # Only some of the cols, to avoid duplicate NUM_SS
  dosage$map)
  
# Compute scores for all individuals in iPSYCH
pred_auto <- big_prodVec(G,
                         beta_auto, # Model
                         ind.col = ipsych_sumstats_overlap[["_NUM_ID_"]], # Indices in G of snps used in auto
                         ncores = nb_cores())

# cbind with family and sample ID and save
scores <- as.data.frame(cbind(covariates_df$family.ID, covariates_df$sample.ID, as.numeric(pred_auto)))
colnames(scores) <- c("family.ID", "sample.ID", "ldpred2_pgs")
saveRDS(scores, scores_out)

cat("Finished computing scores. Models, auto model paramters, følgefil, and auto scores were saved in 4 distinct files in", base_path, "/")

  



