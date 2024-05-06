#' Generalized script for parsing many different summary statistics using MungeSumstats
#' 
#' @author Ole Sahlholdt Hansen
#' @date 02-02-2024
#' 
#' @description This script takes two arguments: An input file of sumstats as well as an output name 
#' 
#' @Todo

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
  library(gwasrapidd)
  library(stringr)
})

# Input and output -------------------------------------------------------------------------------------------------------
paths <- fromJSON(file = "data/paths.json")
load(paths$col_header)
source(paths$get_n_function)
source(paths$create_model_info)
source(paths$or_to_beta_function)
source(paths$snp_match_names)
source(paths$z_to_beta_function)
source(paths$frq_col_function)
source(paths$n_col_function)
# sumstats = read_sumstats(paths$test_path) # For testing

# Command line arguments for this script
args <- commandArgs(trailingOnly = TRUE)
sumstats <- read_sumstats(args[1])
output_path <- args[2]
res_folder <- args[3]

# Outputting column names before munging
colnames(sumstats)

# Removing path and suffix from input str
base_name <- gsub("_munged.rds", "", basename(output_path))

# Accession ID to find in gwas catalog using gwasrapidd ------------------------------------------------------------------
accession_id <- str_match(args[1], "accession\\s*(.*?)\\s*_")[,2]

# If accession ID present, use to get info for model info
if(!is.na(accession_id)) {
  study_info <- get_studies(study_id = accession_id)
  num_inds <- get_n_cas_con(study_info@studies$initial_sample_size) # Get number of cases and controls or n from sample size string
  
} else {
  study_info <- NA
  num_inds <- list("n" = NA, "n_cas" = NA, "n_con" = NA)
}

model_info_df <- create_model_info(accession_id = accession_id) %>% 
  mutate(ID = base_name,
         M_Input = nrow(sumstats))

# Standardizing header --------------------------------------------------------------------------------------------------
# Discarding non-harmonised columns where harmonised version exists
if ("hm_pos" %in% colnames(sumstats)) {
  sumstats <- sumstats %>% select(-base_pair_location)
}
harmonised_exists <- colnames(sumstats)[which(paste0("hm_", colnames(sumstats)) %in% colnames(sumstats))]
sumstats <- sumstats %>% 
  select(-all_of(harmonized_exists))
colnames(sumstats) <- sub("^hm_", "", colnames(sumstats))

sumstats <- standardise_header(sumstats, mapping_file = sumstatsColHeaders, return_list = FALSE)
head(sumstats)

# Inferring reference genome and performing lift_over if necessary -------------------------------------------------------
# If SNP col is chr:pos
if("SNP" %in% colnames(sumstats) & !"RSID" %in% colnames(sumstats)){
  if(all(grepl(":", sumstats$SNP))) { 
    sumstats <- tidyr::separate(sumstats,                                        # SLOW
                                col = SNP, into = c("CHR", "BP"), sep = ":", 
                                convert = TRUE, extra = "drop")                  # If the col contains more than one : - e.g. chr:pos:a0:a1, :a0:a1 is dropped
  }
}

# if(all(c("SNP", "CHR", "BP") %in% colnames(sumstats))) {  # MungeSumstats needs these three cols to infer ref
#   ref_genome <- get_genome_builds(sumstats_list = list(ss1 = sumstats))$ss1
#   
#   if(ref_genome != "GRCH37") {
#     sumstats <- liftover(sumstats_dt = sumstats,
#                          ref_genome = ref_genome,
#                          convert_ref_genome = "GRCh37")
#   }
# }

# Some manual checks ---------------------------------------------------------------------------------------------------

# Small edits for snp_match format --
sumstats <- snp_match_format(sumstats)

# Allele frequency --
sumstats <- check_frq_col(sumstats, study_info, num_inds)

# Effective population size ---
check_n_col_list <- check_n_col(sumstats, num_inds, model_info_df)

sumstats <- check_n_col_list["sumstats"]$sumstats
model_info_df <- check_n_col_list["model_info"]$model_info

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "n_eff" %in% colnames(sumstats) | "n" %in% colnames(sumstats))

# Other effect size than beta -------------------------------------------------------------------------------------------------------------
# Odds ratio
sumstats <- or_to_beta(sumstats)

# Z score
sumstats <- z_to_beta(sumstats)

# Finding Hapmap overlap with sumstats -----------------------------------------------------------------------------------
# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = paths$hapmap_path, fname = "map_hm3_plus.rds"))

# Finding sumstats/HapMap3+ overlap
snp_info <- snp_match(sumstats, info, match.min.prop = 0.1, join_by_pos = FALSE) %>%
  select(-c(pos_hg18, pos_hg38))

cat(nrow(snp_info), "variants in overlap with HapMap3+. \n")

# If frq still not in sumstats, use af_UKBB 
if(!"frq" %in% colnames(snp_info)) { 
  snp_info$frq = snp_info$af_UKBB
}

# Saving in model info
model_info_df$M_HapMap <-  nrow(snp_info)

# QC -------------------------------------------------------------------------------------------------------------------------------------
snp_info$p <- as.numeric(snp_info$p)
if("n_eff" %in% colnames(snp_info)) {
  df_beta <- snp_info %>% 
    filter(beta != 0, beta_se > 0,
           n_eff > (0.7 * max(n_eff))) %>% 
    mutate(sd_af = sqrt(2 * frq * (1 - frq)),
           sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2),
           sd_ss2 = sd_ss/quantile(sd_ss, 0.999) * sqrt(0.5))
} else {
  df_beta <- snp_info %>% 
    filter(beta != 0, beta_se > 0) %>% 
    mutate(sd_af = sqrt(2 * frq * (1 - frq)),
           sd_ss = 2 / sqrt(n * beta_se^2 + beta^2),
           sd_ss2 = sd_ss/quantile(sd_ss, 0.999) * sqrt(0.5))
}

if("info" %in% colnames(df_beta)) {
  df_beta <- filter(df_beta, info > 0.7)
}

# TODO: Incorporate this step into QC
# df_beta$frq2 <- ifelse(df_beta$beta * sumstats$beta[df_beta$`_NUM_ID_.ss`] < 0,
#                           1 - df_beta$frq, df_beta$frq)
# diff <- with(df_beta, abs(af_UKBB - frq2))
  
df_beta$is_bad <- with(df_beta,
                      #  diff > 0.05 |
                       sd_ss2 < (0.7 * sd_af) | 
                       sd_ss2 > (sd_af + 0.1) |
                       sd_ss2 < 0.05 | 
                       sd_af < 0.05)

p <- qplot(sd_af, sd_ss2, color = is_bad, alpha = I(0.5),
           data = slice_sample(df_beta, n = 50e3)) +
           theme_bigstatsr() +
           coord_equal() +
           scale_color_viridis_d(direction = -1) +
           geom_abline(linetype = 2, color = "red") +
           labs(x = "Standard deviations derived from the allele frequencies",
                y = "Standard deviations derived from the summary statistics",
                color = "To remove?")

ggsave(paste0(res_folder, "/", base_name, "_QC.jpeg"), p)

df_beta <- df_beta %>% 
  filter(!is_bad) %>% 
  select(-c(is_bad, sd_af, sd_ss, sd_ss2, af_UKBB))
  

cat(nrow(df_beta), "variants remaining following QC. \n")

# Restricting to iPSYCH variants -----------------------------------------------------------------------------------------
# Reading in iPSYCH data
dosage <- snp_attach(paths$dosage_path) 
dosage$map <- dosage$map %>% 
  rename("chr" = "CHR", "pos" = "POS", "a0" = "a1", "a1" = "a2")

# Finding iPSYCH overlap
in_test <- vctrs::vec_in(df_beta[, c("chr", "pos")], dosage$map[, c("chr", "pos")])
df_beta <- df_beta[in_test, ]                          

# Making sure there is not no variants in sumstats ----------------------------------------------------------------------
assert("Less than 60K variants remaining in summary statistics following QC and Hapmap3+/iPSYCH overlap.",
       nrow(df_beta) > 10000)
cat(nrow(df_beta), "variants remaining after restricting to iPSYCH variants. \n")

# Saving in model info
model_info_df$M_QC <- nrow(df_beta)

# Saving the parsed sumstats in the outputfile ---------------------------------------------------------------------------
head(df_beta)
saveRDS(df_beta, output_path)

# saveRDS(df_beta, paths$test_parsed) # for testing

# Saving model info --------------------------------------------------------------------------------------------------------
write.table(model_info_df,
           file = paste0(res_folder, "/", base_name, "_model_info.csv"),
           sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)
