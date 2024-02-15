#' Generalized script for parsing many different summary statistics using MungeSumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date 02-02-2024
#' 
#' @description This script takes two arguments: An input file of sumstats as well as an output name 
#' 
#' @Todo

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.18") # Version needed for MungeSumstats
# BiocManager::install("MungeSumstats")
# 
# Installing ref genomes
# options(timeout=2000)
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38")
# BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")
# BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37")
# BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5")

# Libs ----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(MungeSumstats)
  library(bigsnpr)
  library(testit)
  library(ggplot2)
  library(rjson)
})

# Input and output -----------------------------------------------------------------
paths <- fromJSON(file = "data/paths.json")
load("data/sumstatsColHeaders.rda")
# sumstats = read_sumstats(paths$test_path) # For testing

# Command line arguments for this script
args <- commandArgs(trailingOnly = TRUE)
sumstats <- read_sumstats(args[1])
output_path <- args[2]
res_folder <- args[3]

# Removing path and suffix from input str
base_name <- gsub("_munged.rds", "", basename(output_path))

# Standardizing header
sumstats <- standardise_header(sumstats, mapping_file = sumstatsColHeaders, return_list = FALSE)

# Inferring reference genome and performing lift_over if necessary -------------------------------------------------------
ref_genome <- get_genome_builds(sumstats_list = list(ss1 = sumstats))$ss1

if(ref_genome != "GRCH37") {
  sumstats <- liftover(sumstats_dt = sumstats,
                       ref_genome = ref_genome,
                       convert_ref_genome = "GRCh37")
}

# Renaming to fit snp_match format and filtering away sex chromosomes ---------------------------------------------------
reformatted <- sumstats %>%
  rename(POS = BP, A0 = A1, A1 = A2, BETA_SE = SE) %>%
  filter(CHR %in% 1:22) %>%
  mutate(CHR = as.numeric(CHR))

# Some manual checks -----------------------------------------------------------------------------------------------------
colnames(reformatted) <- tolower(colnames(reformatted))

# Odds ratio -------------------------------------------------
# If reported effect size is odds ratio, MungeSumstats does not convert se, as long as SE exists
if("or" %in% colnames(reformatted)){
  reformatted$beta = with(reformatted, log(or))
  # TODO: Test if SE already there. If so, check that p-values match. Otherwise recompute beta_se.
  reformatted$beta_se = with(reformatted, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
}

# Effective population size ----------------------------------

if(!"n_eff" %in% colnames(reformatted)){
  reformatted$n_eff = if("neff_half" %in% colnames(reformatted)){
    with(reformatted, neff_half * 2)
  } else if("n_cas" %in% colnames(reformatted)){
    with(reformatted, 4/(1/n_cas + 1/n_con))
  } else {
    # TODO: What about continuous outcomes? Find N. Maybe read from some text file where you manually store this.
  }  
}

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "n_eff" %in% colnames(reformatted))

# Z score ----------------------------------------------------
if(!"beta" %in% colnames(reformatted) & "z" %in% colnames(reformatted)){
  # TODO: What I would do: get beta_se from n_eff and freq
  # 2 * frq (1 - frq) ~ 4 / (n_eff * beta_se^2)
  # get beta from p-val -> |z| and beta_se and sign(z)
    #reformatted$beta = with(reformatted, z / sqrt(2*frq*(1-frq)(n_eff + z^2)))  
    #reformatted$beta_se = with(reformatted, beta/qnorm(1-p/2)) # beta/z
  # TODO: what if `$frq` is missing? Then using the frequencies from `info` later
  # and reverse the freq if needed (cf. https://github.com/privefl/paper-infer/blob/main/code/prepare-sumstats/MDD.R#L43-L44)
}

# Allele frequency ------------------------------------------

# Check if frq columns are on the form frq_a_X and frq_u_X
colnames(reformatted)[grep("^fr?q_a_", colnames(reformatted))] <- "frq_cas"
colnames(reformatted)[grep("^fr?q_u_", colnames(reformatted))] <- "frq_con"

if(!"frq" %in% colnames(reformatted)){
  if("frq_cas" %in% colnames(reformatted)) {
    reformatted$frq = with(reformatted, (frq_cas + frq_con)/2)
  }
}

# Finding Hapmap overlap with sumstats -----------------------------------------------------------------------------------
# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = paths$hapmap_path, fname = "map_hm3_plus.rds"))

# Finding sumstats/HapMap3+ overlap
# TODO: can you find one example where this is needed
snp_info <- snp_match(reformatted, info, match.min.prop = 0.1)                                              

cat(nrow(snp_info), "variants in overlap with HapMap3+. \n")

# QC ------------------------------------------------------------------------------------------------------------------

# Beta, beta_se, and n_eff  ---------------------------------------
df_beta <- snp_info %>% 
  filter(beta != 0 & beta_se > 0 & n_eff > (0.7 * max(n_eff)))

# INFO score ----------------------------------------------------
if("info" %in% colnames(df_beta)) {
  df_beta <- df_beta %>%
    filter(info > 0.7)
}

# Allele frequency -----------------------------------------------
# TODO: always define this `$frq` -> can get it from `info` as before
if("frq" %in% colnames(df_beta)){         # If freq exists
  
  sd_af <- with(df_beta, sqrt(2 * frq * (1 - frq)))
  sd_ss <- with(df_beta, 2 / sqrt(n_eff * beta_se^2 + beta^2))
  sd_ss2 <- sd_ss / quantile(sd_ss, 0.999) * sqrt(0.5) 
  
  is_bad <- 
    sd_ss2 < (0.7 * sd_af) |
    sd_ss2 > (sd_af + 0.1) |
    sd_ss2 < 0.05 |
    sd_af < 0.05
  
  p <- ggplot(slice_sample(data.frame(sd_af, sd_ss2, is_bad), n = 100e3)) +
    geom_point(aes(sd_af, sd_ss2, color = is_bad), alpha = 0.5) +
    theme_bigstatsr(0.9) + 
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red", linewidth = 1.5) +
    labs(x = "Standard deviations derived from the allele frequencies",
         y = "Standard deviations derived from the summary statistics",
         color = "To remove?")
  
  ggsave(paste0(res_folder, "/", base_name, "_QC.jpeg"), p)
  
  df_beta <- df_beta[!is_bad, ] 
  
} else {
  cat("No allele frequencies available in summary statistics. QC step not performed. \n")
}
cat(nrow(df_beta), "variants remaining following QC. \n")

# Restricting to iPSYCH variants -----------------------------------------------------------------------------------------
# Reading in iPSYCH data
dosage <- snp_attach(paths$dosage_path) 
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Finding iPSYCH overlap
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], dosage$map[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]                          

# Making sure there are at least 60K variants in sumstats -----------------------------------------------------------------
# assert("Less than 60K variants remaining in summary statistics following QC and Hapmap3+/iPSYCH overlap.",
#        nrow(df_beta) > 60000)
cat(nrow(df_beta), "variants remaining in munged sumstats. \n")

# Saving the parsed sumstats in the outputfile ---------------------------------------------------------------------------
head(df_beta)
saveRDS(df_beta, output_path)

# saveRDS(df_beta, paths$test_parsed) # for testing

# Foelgefil --------------------------------------------------------------------------------------------------------------
foelgefil_df <- data.frame(
  ID = base_name,
  Restrictions = NA,
  Reported_Trait = NA,
  PubMed_ID = NA,
  First_Author = NA,
  Journal = NA,
  Title = NA,
  Publication_Date = NA,
  N = NA,
  Ncase = NA,
  Ncontrol = NA,
  se = "beta",
  M_or = nrow(sumstats),     # Variants in original sumstats
  M_m = nrow(snp_info),      # Overlap with HapMap3+
  M_qc = nrow(df_beta)       # Variants after QC and iPSYCH overlap
)

write.table(foelgefil_df,
           file = paste0(res_folder, "/", base_name, "_foelgefil.csv"),
           sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)
