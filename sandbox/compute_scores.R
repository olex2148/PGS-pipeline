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
                         ind.col = ipsych_sumstats_index[["_NUM_ID_"]], # Indices in G of snps used in auto
                         ncores = nb_cores())

pred_lassosum <- big_prodVec(G,
                             best_lassosum,
                             ind.col = ipsych_sumstats_index[["_NUM_ID_"]], # Indices in G of snps used in auto
                             ncores = nb_cores())

# TODO: Add scaling?
# auto_scaled <- (pred_auto - mean(pred_auto)) / sd(pred_auto)
# lassosum_scaled <- (pred_lassosum - mean(pred_lassosum)) / sd(pred_lassosum)

# cbind with family and sample ID and save
scores <- as.data.frame(cbind(covariates_df$family.ID, 
                              covariates_df$sample.ID, 
                              as.numeric(pred_auto),
                              as.numeric(pred_lassosum)))
colnames(scores) <- c("family.ID", "sample.ID", "ldpred2_pgs", "lassosum_pgs")
saveRDS(scores, paste0(base_path, "_scores.rds"))

cat("Finished computing scores. Models, auto model parameters, følgefil, and auto + lassosum scores were saved in 4 distinct files in", base_path, "/")


