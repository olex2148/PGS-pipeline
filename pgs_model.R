library(bigsnpr)
library(dplyr)
library(ggplot2)

# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
parsed_sumstats <- args[1]
auto_model_out <- args[2]
scores_out <- args[3]
parameters_out <- args[4]

setwd("~/NCRR-PRS/faststorage/osh/PGS/updates_and_adds/")
set.seed(72)

# Reading in ext sumstats
# df = readRDS(parsed_sumstats)
sumstats = readRDS("steps/parsed_sumstats/test_adhd.rds") #Example

# Reading in iPSYCH
dosage <-  snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/dosage_ipsych2015.rds")
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Reading in HapMap3+
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = "steaps/hapmap3+", fname = "map_hm3_plus.rds"))

# Finding Hapmap and iPSYCH overlap with sumstats
info_ipsych_external <- snp_match(sumstats, dosage$map) # iPSYCH-external
info_snp <- snp_match(info_ipsych_external[, c("pos","chr","a0","a1","beta","beta_se","freq","p","n_eff", "info")], info) #hapmap overlap with remaining

# QC
# TODO: Add check for the presence of the info and freq column

# Removing SEs below zero to avoid errors in sd_ss
info_snp2 <- info_snp %>% 
  filter(beta_se > 0)

sd_af <- with(info_snp2, sqrt(2 * freq * (1 - freq)))
sd_ss <- with(info_snp2, 2 / sqrt(n_eff * beta_se^2 + beta^2))
sd_ss2 <- sd_ss / quantile(sd_ss, 0.999) * sqrt(0.5) 
is_bad <- 
  sd_ss2 < (0.5 * sd_af) |
  sd_ss2 > (sd_af + 0.1) |
  sd_ss2 < 0.05 |
  sd_af < 0.05 |
  info_snp2$n_eff < (0.7 * max(info_snp2$n_eff))|
  info_snp2$info < 0.7

table(is_bad)

df_beta <- info_snp2[!is_bad, ]

#TODO: Make png non-hard coded
png(filename = "")
ggplot(slice_sample(data.frame(sd_af, ad_ss2, is_bad), n = 50e4)) +
  geom_point(aes(sd_af, sd_ss2, color = is_bad), alpha = 0.5) +
  theme_bigstatsr(0.9) + 
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red", size = 1.5) +
  labs(x = "Standard deviations in the reference set",
       y = "Standard deviations derived from the summary statistics",
       color = "To remove?")
dev.off()

ldsc <- with(df_beta, snp_ldsc(ld, ld_size = length(ld),
                               chi2 = (beta / beta_se)^2, #qnorm(1-p/2)^2,
                               sample_size = n_eff,
                               ncores = nb_cores()))
h2_est <- ldsc[["h2"]]
h2_liab <- h2_est * coef_to_liab(K_pop = 0.05) #TODO: Remove - k_pop is specific, just for testing

# LDpred2-auto
tmp <- tempfile(tmpdir = "steps/corr")
for (chr in 1:22) {
  
  cat(chr, "..", sep = "")
  
  ## indices in df_beta
  ind.chr <- which(df_beta$chr == chr)
  ## indices in map_ldref
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in corr_chr
  ind.chr3 <- match(ind.chr2, which(info$chr == chr))
  
  corr_chr <- readRDS(paste0("~/NCRR-PRS/faststorage/florian/ldpred2-inference/data/corr_hm3_plus/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}

coef_shrink <- 0.8

multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 50), burn_in = 500, num_iter = 500,
  use_MLE = FALSE, # for power/convergence issues, alpha
  report_step = 20, ncores = nb_cores(), allow_jump_sign = FALSE, shrink_corr = coef_shrink)

saveRDS(list(ldsc = ldsc, ldpred2 = multi_auto),
        auto_model_out)
saveRDS(list(ldsc = ldsc, ldpred2 = multi_auto),
        "results/models/adhd_test.rds")

###### Predicting in iPSYCH
# Reading in PCs and some iPSYCH meta info for covariates
pcs <- readRDS("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/PC.rds")
meta <- bigreadr::fread2("~/Register/2019_06/csv/ipsych2015design_v2.csv")

# Compute sex and age
# TODO: Remove adhd2015I, only used for testing
df <- dosage$fam %>% 
  left_join(meta[, c("fdato", "gender", "pid", "adhd2015I")], by = c("family.ID" = "pid")) %>% 
  mutate(sex = ifelse(gender == "F", 1, 0),
         # Time diff in years between present date and fdate
         age = lubridate::time_length(difftime(as.Date(Sys.Date(), format = "%d/%m/%Y"), as.Date(fdato, format = "%d/%m/%Y")), "years"))

# Checking order is preserved
identical(dosage$fam$sample.ID, df$sample.ID)

cov <- cbind(df$sex, df$age, df$is_2012, pcs)

G <- dosage$genotypes

# Finding indices of variants in G which is used in model
ipsych_sumstats_index <- snp_match(
  df_beta[, c("pos","chr","a0","a1","beta","beta_se","freq","p","n_eff", "info")], # Only some of the cols, to avoid duplicate NUM_SS
  dosage$map)

# QC on chains from auto
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- (range > (0.95 * quantile(range, 0.95)))
sum(keep)

# Average of kept chains
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))

#Inspect a chain that has passed QC
#TODO: Make universal, and save png of chains
auto <- multi_auto[[8]]
#TODO: make png name non-hard coded
png(filename = "")
plot_grid(
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
dev.off()
  
# Compute scores for all individuals in iPSYCH
pred_auto <- big_prodVec(G,
                         beta_auto, # Model
                         ind.col = ipsych_sumstats_overlap[["_NUM_ID_"]], # Indices in G of snps used in auto
                         ncores = nb_cores())
#cbind with family and sample ID and save

scores <- as.data.frame(cbind(df$family.ID, df$sample.ID, as.numeric(pred_auto)))
colnames(scores) <- c("family.ID", "sample.ID", "PGS")
saveRDS(scores, scores_out)

# Save parameters
all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
quantile(all_h2, c(0.5, 0.025, 0.975))

all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
quantile(all_p, c(0.5, 0.025, 0.975))

#MLE
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

#TODO: import code to make metrics file, save not in rds
saveRDS(list(r2 = all_r2, h2 = all_h2, alpha = all_alpha, p = all_p),
        parameters)
  



