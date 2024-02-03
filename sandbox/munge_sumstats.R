#' Generalized script for parsing many different summary statistics using MungeSumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date 02-02-2024
#' 
#' @description This script takes two arguments: An input file of sumstats as well as an output name 
#' 
#' @Todo
#'   - Prioritize Neff_half or 4/(1/cases + 1/controls)?


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18") # Version needed for MungeSumstats
BiocManager::install("MungeSumstats")

# Installing ref genomes
BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37", "BSgenome.Hsapiens.1000genomes.hs37d5",
                     "SNPlocs.Hsapiens.dbSNP144.GRCh38", "BSgenome.Hsapiens.NCBI.GRCh38")

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
output = paste0("steps/munged_sumstats/", base_name, "_munged.rds")

# sumstats = read_sumstats("data/daner_bip_pgc3_nm_noukbiobank") # For testing
head(sumstats)

# Inferring ref genome
sumstats_list <- list(ss1 = sumstats) # Because "get_genome_build" doesn't exist
ref_genome <- get_genome_builds(sumstats_list = sumstats_list, dbSNP = 144)$ss1     # ~2 mins

# Formatting sumstats using MungeSumstats ---------------------------------------------------------------------------------
reformatted <- format_sumstats(path=sumstats,                                       # ~10 mins
                               ref_genome=ref_genome, dbSNP = 144,                       # Detected ref genome
                               convert_ref_genome = "GRCh37",                            # Convert to HapMap3+ build if not already GRCh37
                               impute_beta = TRUE, impute_se = TRUE,                     # Convert OR to beta and se to beta_se
                               INFO_filter = 0.7,                                        # Filtering on INFO score - could be more strict (default 0.9)
                               nThread = nb_cores(),
                               return_data = TRUE, return_format = "data.table") %>% 
          rename(POS = BP, A0 = A1, A1 = A2) %>%                                         # Renaming to fit LDpred2 format
          mutate(CHR = ifelse(CHR == "X", 23, ifelse(CHR == "Y", 24, as.numeric(CHR))))  # Converting X and Y chr to 23 and 24 

# Some manual checks -------------------------------------------------------------------------------------------------------

# Odds ratio -------------------------------------------------

# If reported effect size is odds ratio, MungeSumstats 
# does not convert se, as long as SE exists
if("OR" %in% colnames(reformatted)){
  reformatted$BETA_SE = reformatted$BETA/qnorm(1-reformatted$P / 2) # beta/z
}

# Effective population size ----------------------------------

if(!"N_EFF" %in% colnames(reformatted)){
  if("N_CAS" %in% colnames(reformatted)){
    reformatted$N_EFF = 4/(1/reformatted$N_CAS + 1/reformatted$N_CON)
  } else {
      if("NEFF_HALF" %in% colnames(reformatted)){
        reformatted$N_EFF = reformatted$NEFF_HALF * 2
      }
  }
}

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "N_EFF" %in% colnames(reformatted))

# Allele frequency ------------------------------------------

# Check if frq columns are on the form frq_a_X and frq_u_X
colnames(reformatted)[grep("^FR?Q_A_", colnames(reformatted))] <- "FRQ_CAS"
colnames(reformatted)[grep("^FR?Q_U_", colnames(reformatted))] <- "FRQ_CON"

if(!"FRQ" %in% colnames(reformatted)){
  if("FRQ_CAS" %in% colnames(reformatted)) {
    reformatted$FRQ = (reformatted$FRQ_CAS + reformatted$FRQ_CON)/2
  }
}

# QC ------------------------------------------------------------------------------------------------------------------

# Beta, beta_se, and n_eff  ---------------------------------------
parsed_sumstats <- reformatted %>% 
  filter(BETA != 0, BETA_SE > 0 & N_EFF > (0.7 * max(N_EFF)))

# Allele frequency -----------------------------------------------
if("FRQ" %in% colnames(parsed_sumstats)){         # If freq exists
  
  sd_af <- with(parsed_sumstats, sqrt(2 * FRQ * (1 - FRQ)))
  sd_ss <- with(parsed_sumstats, 2 / sqrt(N_EFF * BETA_SE^2 + BETA^2))
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
  
  parsed_sumstats <- parsed_sumstats[!is_bad, ] 
  
} else {
  cat("No allele frequencies available in summary statistics. QC step not performed. \n")
}


# Finding Hapmap and iPSYCH overlap with sumstats ----------------------------------------------------------------------------

# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = hapmap_path, fname = "map_hm3_plus.rds"))

# Making colnames lower case for snp_match
colnames(parsed_sumstats) <- tolower(colnames(parsed_sumstats)) 

# Finding sumstats/iPSYCH overlap
df_beta <- snp_match(parsed_sumstats, dosage$map)                                                                                      

# Finding HapMap3+ overlaps
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], info[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]                                                                                      

# Making sure there are at least 500K variants in sumstats
assert("Less than 500K variants remaining in summary statistics following QC and Hapmap3+/iPSYCH overlap.", 
       nrow(parsed_sumstats) > 500000)
cat(nrow(parsed_sumstats), "variants remaining in munged sumstats.")

# Saving the parsed sumstats in the outputfile ------------------------------------------------------------------------------

head(parsed_sumstats)
saveRDS(parsed_sumstats, output)

# saveRDS(parsed_sumstats, "steps/parsed_sumstats/test_adhd.rds") # for testing
