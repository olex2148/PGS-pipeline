library(dplyr)
library(Dict)

# Command line arguments for this script
# args = commandArgs(trailingOnly=True)
#df = read.table(args[1], header = TRUE) #read.table slow, revisit fread, arrow not useful
# outputs = args[2]

setwd("~/NCRR-PRS/faststorage/osh/PGS/updates_and_adds/")
df = read.table("~/NCRR-PRS/faststorage/ExtSumStats/ADHD/daner_ADHD_meta_PGCeur_v1_filtered.meta", header = TRUE) # For testing
head(df)

#Renaming columns by fist creating a dict with the new and pissble original colnames
idkey <- Dict$new("a0" = c("a1","allele_1","allele1", "eff_allele", "effect_allele", "ea", "testallele"),
                  "a1" = c("a2", "allele_2", "allele2", "alt_allele", "alternative_allele", "other_allele", "nea", "oa"),
                  "chr" = c("chromosome", "chr_id", "chrom"),
                  "pos" = c("bp", "position", "bp_pos","base_pair_location", "posgrch37"),
                  "freq" = c("freq1","freq1.hapmap","frq","effect_allele_frequency","frq_a1","eaf"),
                  "p" = c("pval","p.value","p_val","p_value","p-value"),
                  "beta" = c("effect","effectsize","b","weight","est"),
                  "beta_se" = c("se","stderr","standard_error"),
                  "n_eff" = c("effective_n","neff","n_eff"),
                  "info" = c("impinfo","infoscore","info_score"),
                  "or" = c("oddsratio","odds_ratio"),
                  "n_cases" = c("ncas","total_cases","total_ncase","nca"),
                  "n_controls" = c("ncon", "total_ncontrol", "nco"),
                  "half_neff" = c("neffdiv2","neff_half"),
                  "freq_cases" = c("fcas", "eaf_cases"),
                  "freq_controls" = c("fcon", "eaf_controls"),
                  "n_total" = c("totaln", "totalsamplesize", "n", "total_n"),
                  "z" = c("zscore","z_score","z-score"))

colnames(df) <- tolower(colnames(df)) #reduce number of possible versions, lower case letters

# Renaming for all keys in the above dict: Looping over keys in dict, to see if any of their vals are present in colnames
for (key in idkey$keys) {
  colnames(df[colnames(df) %in% idkey[key]] <- key)
  
}

#TODO: check if continous, e.g. BMI, in which case just use N
# check if df has OR, freq, and effective N
if("or" %in% colnames(df)){
  df$beta = log(df$or);
  df$beta_se = df(df$beta/qnorm(1-df$p/2)) #beta/z
}

#TODO: ad if only total N
if(!"n_eff" %in% colnames(df)){
  if("n_cases" %in% colnames(df)){
    df$n_eff = 4/(1/df$n_cases + 1/df$n_controls)
  } else {
      if("half_neff" %in% colnames(df)){
        df$n_eff = df$half_neff * 2
      } else {
        df$n_eff = NA
    }
  }
}

if(!"freq" %in% colnames(df)){
  if("freq_cases" %in% colnames(df)) {
    df$freq = (dffreq_cases + df$freq_controls)/2
  } else {
    df$freq = NA
  } 
}

if(!"info" %in% colnames(df)){
  df$info = NA #So it's created for the command below if it's not already in sumstats
}

sumstats <- 
  select(chr, pos, a0, a1, beta, beta_se, p, n_eff, freq, info)

# Checking the renaming
head(sumstats)

# Saving the parserd DF in the outputfile
# saveRDS(sumstats, output)
saveRDS(sumstats, "steps/parsed_sumstats/test_adhd.rds") # for testing
