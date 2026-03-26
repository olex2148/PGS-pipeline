library(bigsnpr)
library(dplyr)
library(data.table)

# fam object
dosage <-  snp_attach("~")

# Scores
prs_one <- readRDS("/")
prs_two <- readRDS("")
prs_three <- readRDS("")$prs

# Munged sumstats - for computing the weights in the next step
sumstats_prs_one <- readRDS("")
sumstats_prs_two <- readRDS("")

## PRS 1
n_cas_prs_one <- 18000
n_con_prs_one <- 42000

## PRS 2
n_cas_prs_two <- 15000
n_con_prs_two <- 35000

n_eff_prs_one <- 4 / (1 / n_cas_prs_one + 1 / n_con_prs_one)
n_eff_prs_two <- 4 / (1 / n_cas_prs_two + 1 / n_con_prs_two)

sumstats_prs_one$n_eff <- n_eff_prs_one
sumstats_prs_two$n_eff <- n_eff_prs_two

# And their respective weights (sqrt(neff))
w_prs_one <- mean(sqrt(sumstats_prs_one$n_eff)) # 103
w_prs_two <- mean(sqrt(sumstats_prs_two$n_eff)) # 177
w_prs_three <- sqrt(49200) # 209

# Centering and scaling - could also be with (scale() function)
prs_one_scaled <- (prs_one$ldpred2_pgs - mean(prs_one$ldpred2_pgs))/sd(prs_one$ldpred2_pgs)
prs_two_scaled <- (prs_two$ldpred2_pgs - mean(prs_two$ldpred2_pgs))/sd(prs_two$ldpred2_pgs)
prs_three_scaled <- (prs_three - mean(prs_three))/sd(prs_three)


# meta all
meta_int_decode_pgc <- prs_one_scaled * w_prs_one + prs_two_scaled * w_prs_two + prs_three_scaled * w_prs_three

# meta pgc, decode
meta_decode_pgc <- prs_one_scaled * w_prs_one + prs_two_scaled * w_prs_two


# Save
all_scores <- data.frame(
  IID = dosage$fam$sample.ID,
  FID = dosage$fam$family.ID,
  pgc = prs_one$ldpred2_pgs,
  decode = prs_two$ldpred2_pgs,
  pgc_decode = meta_decode_pgc,
  ipsych_5_cv = prs_three,
  ipsych_pgc_decode = meta_int_decode_pgc
  
)

write.csv(all_scores, "", 
          quote = FALSE, row.names = FALSE)

