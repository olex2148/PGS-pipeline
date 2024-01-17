#' Generalized script for parsing many different summary statistics
#' 
#' @author Ole Sahlholdt Hansen
#' @date
#' 
#' @description This script takes two arguments: An input file of sumstats as well as an output name 
#' 
#' @Todo
#'   - Update dict to contain names from MungeSumstats too
#'   - Check if continuous

suppressPackageStartupMessages({
  library(dplyr)
  library(Dict)
  library(data.table)
})

# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
sumstats = fread(args[1])
output = args[2]

# sumstats = read.table("", header = TRUE) # For testing
head(sumstats)

#Renaming columns by fist creating a dict with the new and pissble original colnames
idkey <- Dict$new("a0" = c("a1","allele_1","allele1", "eff_allele", "effect_allele", "ea", "reference_allele", "testallele"),
                  "a1" = c("a2", "allele_2", "allele2", "alt_allele", "alternative_allele", "other_allele", "nea", "oa"),
                  "chr" = c("chromosome", "chr_id", "chrom"),
                  "pos" = c("bp", "position", "bp_pos","base_pair_location", "posgrch37"),
                  "freq" = c("freq1","freq1.hapmap", "frq", "effect_allele_frequency", "frq_a1", "eaf"),
                  "p" = c("pval", "p.value", "p_val", "p_value", "p-value"),
                  "beta" = c("effect","effectsize", "b", "est"),
                  "beta_se" = c("se","stderr", "standard_error"),
                  "n_eff" = c("effective_n", "neff"),
                  "info" = c("impinfo", "infoscore", "info_score"),
                  "or" = c("oddsratio", "odds_ratio"),
                  "n_cases" = c("ncas", "total_cases", "total_ncase", "nca"),
                  "n_controls" = c("ncon", "total_ncontrol", "nco"),
                  "half_neff" = c("neffdiv2", "neff_half"),
                  "freq_cases" = c("fcas", "eaf_cases"),
                  "freq_controls" = c("fcon", "eaf_controls"),
                  "n_total" = c("totaln", "totalsamplesize", "n", "total_n"),
                  "z" = c("zscore","z_score","z-score"))

colnames(sumstats) <- tolower(colnames(sumstats)) # To reduce number of possible versions

# Renaming for all keys in the above dict: Looping over keys in dict, to see if any of their vals are present in colnames
for (key in idkey$keys) {
  colnames(sumstats[colnames(sumstats) %in% idkey[key]] <- key)
  
}

# Check if frq columns are on the form frq_a_X and frq_u_X
colnames(sumstats)[grep("^fr?q_a_", colnames(sumstats))] <- "freq_cases"
colnames(sumstats)[grep("^fr?q_u_", colnames(sumstats))] <- "freq_controls"

# check if sumstats has OR, freq, and effective N
if("or" %in% colnames(sumstats)){
  sumstats$beta = log(sumstats$or);
  sumstats$beta_se = sumstats$beta/qnorm(1-sumstats$p / 2) # beta/z
}

if(!"n_eff" %in% colnames(sumstats)){
  if("n_cases" %in% colnames(sumstats)){
    sumstats$n_eff = 4/(1/sumstats$n_cases + 1/sumstats$n_controls)
  } else {
      if("half_neff" %in% colnames(sumstats)){
        sumstats$n_eff = sumstats$half_neff * 2
      } else {
        sumstats$n_eff = NA
    }
  }
}

if(!"freq" %in% colnames(sumstats)){
  if("freq_cases" %in% colnames(sumstats)) {
    sumstats$freq = (sumstatsfreq_cases + sumstats$freq_controls)/2
  } else {
    sumstats$freq = NA
  } 
}

if(!"info" %in% colnames(sumstats)){
  sumstats$info = NA #So it's created for the command below if it's not already in sumstats
}

parsed_sumstats <- sumstats %>% 
  select(chr, pos, a0, a1, beta, beta_se, p, n_eff, freq, info)

# Checking the renaming
head(parsed_sumstats)

# Saving the parsed sumstats in the outputfile
saveRDS(parsed_sumstats, output)
# saveRDS(parsed_sumstats, "steps/parsed_sumstats/test_adhd.rds") # for testing
