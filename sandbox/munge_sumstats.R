#' Generalized script for parsing many different summary statistics using MungeSumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date 02-02-2024
#' 
#' @description This script takes two arguments: An input file of sumstats as well as an output name 
#' 
#' @Todo
#'   - Prioritize Neff_half or 4/(1/cases + 1/controls)?


# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18") # Version needed for MungeSumstats
# BiocManager::install("MungeSumstats")
# 
# # Installing ref genomes
# BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37", "BSgenome.Hsapiens.1000genomes.hs37d5",
#                      "SNPlocs.Hsapiens.dbSNP144.GRCh38", "BSgenome.Hsapiens.NCBI.GRCh38")

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MungeSumstats)
  library(bigsnpr)
  library(testit)
  library(ggplot2)
})

# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
sumstats = read_sumstats(args[1])
base_name = args[2]

# Some names of output
output = paste0("steps/munged_sumstats/", base_name, "_munged.rds") # Final output of script
tmp = paste0("steps/tmp/", base_name, ".tsv.gz")                    # File for the result of format_sumstats, which is deleted 

source("code/aux/input_paths.R")

# sumstats = read_sumstats(test_path) # For testing
head(sumstats)

# Formatting sumstats using MungeSumstats ---------------------------------------------------------------------------------
# Inferring ref genome  -- list() because "get_genome_build" doesn't exist
ref_genome <- get_genome_builds(sumstats_list = list(ss1 = sumstats), dbSNP = 144, nThread = nb_cores())$ss1  # ~2 mins

reformatted <- format_sumstats(path=sumstats,                                           # ~8-10 mins
                               ref_genome=ref_genome, dbSNP = 144,                       # Detected ref genome
                               convert_ref_genome = "GRCh37",                            # Convert to HapMap3+ build if not already GRCh37
                               impute_beta = TRUE, impute_se = TRUE,                     # Convert OR to beta and se to beta_se
                               INFO_filter = 0.7,                                        # Filtering on INFO score - could be more strict (default 0.9)
                               nThread = nb_cores(),
                               mapping_file = sumstatsColHeaders,                        # Local mapping file
                               return_data = TRUE, return_format = "data.table",
                               save_path = tmp, force_new = TRUE) %>%                    
          rename(POS = BP, A0 = A1, A1 = A2, BETA_SE = SE) %>%                           # Renaming to fit LDpred2 format
          mutate(CHR = ifelse(CHR == "X", 23, ifelse(CHR == "Y", 24, as.numeric(CHR))))  # Converting X and Y chr to 23 and 24 

# Deleting the file that has been written to disc
file.remove(tmp)

# Finding Hapmap and iPSYCH overlap with sumstats ----------------------------------------------------------------------------

# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = hapmap_path, fname = "map_hm3_plus.rds"))

# Reading in iPSYCH data
dosage <- snp_attach(dosage_path) 
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Making colnames lower case for snp_match
colnames(reformatted) <- tolower(colnames(reformatted)) 

# Finding sumstats/iPSYCH overlap
snp_info <- snp_match(reformatted, dosage$map)                                                                                      

# Finding HapMap3+ overlap
in_test <- vctrs::vec_in(snp_info[, c("chr", "pos")], info[, c("chr", "pos")])
snp_info <- snp_info[in_test, ]                          

# Some manual checks -------------------------------------------------------------------------------------------------------

# Odds ratio -------------------------------------------------

# If reported effect size is odds ratio, MungeSumstats does not convert se, as long as SE exists
if("or" %in% colnames(snp_info)){
  snp_info$beta_se = snp_info$beta/qnorm(1-snp_info$p / 2) # beta/z
}

# Effective population size ----------------------------------

if(!"n_eff" %in% colnames(snp_info)){
  if("n_cas" %in% colnames(snp_info)){
    snp_info$n_eff = 4/(1/snp_info$n_cas + 1/snp_info$n_con)
  } else {
      if("neff_half" %in% colnames(snp_info)){
        snp_info$n_eff = snp_info$neff_half * 2
      }
  }
}

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "n_eff" %in% colnames(snp_info))

# Allele frequency ------------------------------------------

# Check if frq columns are on the form frq_a_X and frq_u_X
colnames(snp_info)[grep("^fr?q_a_", colnames(snp_info))] <- "frq_cas"
colnames(snp_info)[grep("^fr?q_u_", colnames(snp_info))] <- "frq_con"

if(!"frq" %in% colnames(snp_info)){
  if("frq_cas" %in% colnames(snp_info)) {
    snp_info$frq = (snp_info$frq_cas + snp_info$frq_con)/2
  }
}

# QC ------------------------------------------------------------------------------------------------------------------

# Beta, beta_se, and n_eff  ---------------------------------------
df_beta <- snp_info %>% 
  filter(beta != 0 & beta_se > 0 & n_eff > (0.7 * max(n_eff)))

# Allele frequency -----------------------------------------------
if("frq" %in% colnames(df_beta)){         # If freq exists
  
  sd_af <- with(df_beta, sqrt(2 * frq * (1 - frq)))
  sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
  sd_ss2 <- sd_ss / quantile(sd_ss, 0.999) * sqrt(0.5) 
  
  is_bad <- 
    sd_ss2 < (0.5 * sd_af) |
    sd_ss2 > (sd_af + 0.1) |
    sd_ss2 < 0.05 |
    sd_af < 0.05
  
  p <- ggplot(slice_sample(data.frame(sd_af, sd_ss2, is_bad), n = 50e4)) +
    geom_point(aes(sd_af, sd_ss2, color = is_bad), alpha = 0.5) +
    theme_bigstatsr(0.9) + 
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red", linewidth = 1.5) +
    labs(x = "Standard deviations in the reference set",
         y = "Standard deviations derived from the summary statistics",
         color = "To remove?")
  
  ggsave(paste0("results/", base_name, "/", base_name, "_QC.jpeg"), p)
  
  df_beta <- df_beta[!is_bad, ] 
  
} else {
  cat("No allele frequencies available in summary statistics. QC step not performed. \n")
}

# Making sure there are at least 500K variants in sumstats -----------------------------------------------------------------
assert("Less than 500K variants remaining in summary statistics following QC and Hapmap3+/iPSYCH overlap.", 
       nrow(df_beta) > 500000)
cat(nrow(df_beta), "variants remaining in munged sumstats.")

# Saving the parsed sumstats in the outputfile ------------------------------------------------------------------------------

head(df_beta)
saveRDS(df_beta, output)

# saveRDS(df_beta, "steps/munged_sumstats/test_adhd.rds") # for testing
