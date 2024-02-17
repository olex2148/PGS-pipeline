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
  library(gwasrapidd)
  library(stringr)
})

# Input and output -------------------------------------------------------------------------------------------------------
paths <- fromJSON(file = "data/paths.json")
load(paths$col_header)
source(paths$get_n_function)
# sumstats = read_sumstats(paths$test_path) # For testing

# Command line arguments for this script
args <- commandArgs(trailingOnly = TRUE)
sumstats <- read_sumstats(args[1])
output_path <- args[2]
res_folder <- args[3]

# Removing path and suffix from input str
base_name <- gsub("_munged.rds", "", basename(output_path))

# Accession ID to find in gwas catalog using gwasrapidd ------------------------------------------------------------------
accession_id <- str_match(args[1], "accession\\s*(.*?)\\s*_")[,2]

# Not all requested sumstats had accession ID
if(!is.na(accession_id)) {
  study_info <- get_studies(study_id = accession_id)
  
  num_inds <- get_n_cas_con(study_info@studies$initial_sample_size)
} else {
  study_info <- NA
}




# Standardizing header --------------------------------------------------------------------------------------------------
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
# If reported effect size is odds ratio
if("or" %in% colnames(reformatted)){
  reformatted$beta = with(reformatted, log(or))
  
  # If SE already there, check that derived and reported p-values match. Otherwise recompute beta_se.
  if("beta_se" %in% colnames(reformatted)) {
    z <- with(reformatted, beta/beta_se)
    derived_pval <- 2 * (1 - pnorm(abs(z)))
    pval_cor <- cor(derived_pval, reformatted$p) # Maybe other comparison method? 
    
    # If cor is poor, the reported SE is not SE of log(OR), but probably SE of OR. Therefore, compute SE of beta
    if(pval_cor < 0.9) { # Fitting threshold?
      reformatted$beta_se = with(reformatted, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
    }
  # if SE not there, estimate with beta/z
  } else if (!"beta_se" %in% colnames(reformatted)) {
    reformatted$beta_se = with(reformatted, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
  }
}

# Effective population size ----------------------------------
if(!"n_eff" %in% colnames(reformatted)){
  reformatted$n_eff = if("neff_half" %in% colnames(reformatted)){
    with(reformatted, neff_half * 2)
  } else if("n_cas" %in% colnames(reformatted)){
    with(reformatted, 4/(1/n_cas + 1/n_con))
  } else if(!is.na(study_info)){
    with(num_inds, 4/(1/n_cas + 1/n_con))
  }
}

# Total population size for continuous traits ----------------
if(!"n" %in% colnames(reformatted)) {
  reformatted$n = if(!is.na(study_info)){
    sum(study_info@ancestries$number_of_individuals)
  }
}

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "n_eff" %in% colnames(reformatted) | "n" %in% colnames(reformatted))

# Z score ----------------------------------------------------
# if(!"beta" %in% colnames(reformatted) & "z" %in% colnames(reformatted)){
  # TODO: What I would do: get beta_se from n_eff and freq
  # 2 * frq (1 - frq) ~ 4 / (n_eff * beta_se^2)
  # get beta from p-val -> |z| and beta_se and sign(z)
    #reformatted$beta = with(reformatted, z / sqrt(2*frq*(1-frq)(n_eff + z^2)))  
    #reformatted$beta_se = with(reformatted, beta/qnorm(1-p/2)) # beta/z
  # TODO: what if `$frq` is missing? Then using the frequencies from `info` later
  # and reverse the freq if needed (cf. https://github.com/privefl/paper-infer/blob/main/code/prepare-sumstats/MDD.R#L43-L44)
# }

# Allele frequency ------------------------------------------

# Check if frq columns are on the form frq_a_X and frq_u_X
colnames(reformatted)[grep("^fr?q_a_", colnames(reformatted))] <- "frq_cas"
colnames(reformatted)[grep("^fr?q_u_", colnames(reformatted))] <- "frq_con"

cases <- 
controls <- 

# Calc frq and weighting by number of cases and controls
if(!"frq" %in% colnames(reformatted)){
  if("frq_cas" %in% colnames(reformatted)) {
    reformatted$frq = with(reformatted, (freq1 * cases + freq2 * controls) / (cases + controls))
  }
}

# Finding Hapmap overlap with sumstats -----------------------------------------------------------------------------------
# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = paths$hapmap_path, fname = "map_hm3_plus.rds"))

# Finding sumstats/HapMap3+ overlap
# TODO: can you find one example where this is needed (match.min.prop)
snp_info <- snp_match(reformatted, info)                                              

cat(nrow(snp_info), "variants in overlap with HapMap3+. \n")


# If frq and info not in sumstats now, use frq and info scores from Hapmap ----------------------------------------------
if(!"frq" %in% colnames(df_beta)){         
  reformatted$frq = with(reformatted, (x + y)/2)
}

if("info" %in% colnames(df_beta)) {
  reformatted$frq = reformatted$something
}

# QC ------------------------------------------------------------------------------------------------------------------

df_beta <- snp_info %>% 
  # Beta, beta_se, n_eff, info  -----------------------------------
  filter(beta != 0,
         beta_se > 0,
         n_eff > (0.7 * max(n_eff)),
         info > 0.7
         ) %>% 
  # sd of af compared to sd of ss
  mutate(sd_af = sqrt(2 * frq * (1 - frq)),
         sd_ss = 2 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss2 = sd_ss/quantile(sd_ss, 0.999) * sqrt(0.5))

# TODO: Incorporate
# df_beta$freq2 <- ifelse(df_beta$beta * sumstats$beta[df_beta$`_NUM_ID_.ss`] < 0,
#                           1 - df_beta$freq, df_beta$freq)

df_beta$is_bad <- with(df_beta,
                       sd_ss2 < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                       sd_ss2 < 0.1 | sd_af < 0.05)

p <- qplot(sd_af, sd_ss2, color = is_bad, alpha = I(0.5),
           data = slice_sample(info_snp2, n = 50e3)) +
           theme_bigstatsr() +
           coord_equal() +
           scale_color_viridis_d(direction = -1) +
           geom_abline(linetype = 2, color = "red") +
           labs(x = "Standard deviations derived from the allele frequencies",
                y = "Standard deviations derived from the summary statistics",
                color = "To remove?")

ggsave(paste0(res_folder, "/", base_name, "_QC.jpeg"), p)

df_beta <- filter(df_beta, !is_bad)
  

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
if(!is.na(study_info)) {
  foelgefil_df <- data.frame(
    ID = base_name,
    Accession_ID = accession_id,
    Reported_Trait = study_info@studies$reported_trait,
    PubMed_ID = study_info@publications$pubmed_id,
    First_Author = study_info@publications$author_fullname,
    Journal = study_info@publications$publication,
    Title = study_info@publications$title,
    Publication_Date = study_info@publications$publication_date,
    N = sum(study_info@ancestries$number_of_individuals),
    Ncase = sum_cases,
    Ncontrol = sum_controls,
    Ancestral_Group = study_info@ancestral_groups$ancestral_group,
    M_input = nrow(sumstats),     # Variants in original sumstats
    M_hapmap = nrow(snp_info),    # Overlap with HapMap3+
    M_qc = nrow(df_beta)          # Variants after QC and iPSYCH overlap
  )
} else {
  foelgefil_df <- data.frame(
    ID = base_name,
    M_or = nrow(sumstats),     # Variants in original sumstats
    M_m = nrow(snp_info),      # Overlap with HapMap3+
    M_qc = nrow(df_beta)       # Variants after QC and iPSYCH overlap
  )
}
write.table(foelgefil_df,
           file = paste0(res_folder, "/", base_name, "_foelgefil.csv"),
           sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)
