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
  num_inds <- get_n_cas_con(study_info@studies$initial_sample_size) # Get number of cases and controls from sample size string
} else {
  study_info <- NA
}

# Standardizing header --------------------------------------------------------------------------------------------------
sumstats <- standardise_header(sumstats, mapping_file = sumstatsColHeaders, return_list = FALSE)
head(sumstats)

# Inferring reference genome and performing lift_over if necessary -------------------------------------------------------
# If SNP col is chr:pos
if(all(grepl(":", sumstats$SNP))) { 
  sumstats <- tidyr::separate(sumstats,                                        # SLOW
                              col = SNP, into = c("CHR", "BP"), sep = ":", 
                              convert = TRUE, extra = "drop")                  # If the col contains more than one : - e.g. chr:pos:a0:a1, :a0:a1 is dropped
}

if(all(c("SNP", "CHR", "BP") %in% colnames(sumstats))) {  # MungeSumstats needs these three cols to infer ref
  ref_genome <- get_genome_builds(sumstats_list = list(ss1 = sumstats))$ss1
  
  if(ref_genome != "GRCH37") {
    sumstats <- liftover(sumstats_dt = sumstats,
                         ref_genome = ref_genome,
                         convert_ref_genome = "GRCh37")
  }
}

# Renaming to fit snp_match format and filtering away sex chromosomes ---------------------------------------------------
colnames(sumstats) <- tolower(colnames(sumstats))

if(all(c("a1", "a2") %in% colnames(sumstats))){  sumstats <- rename(sumstats, a0 = a1, a1 = a2)}
if("bp" %in% colnames(sumstats)){                sumstats <- rename(sumstats, pos = bp)}
if("snp" %in% colnames(sumstats)){               sumstats <- rename(sumstats, rsid = snp)}
if("se" %in% colnames(sumstats)){                sumstats <- rename(sumstats, beta_se = se)}
if("chr" %in% colnames(sumstats)){               sumstats <- filter(sumstats, chr %in% 1:22) %>%  mutate(chr = as.numeric(chr))}

# Removing some redundant cols
if("direction" %in% colnames(sumstats)){         sumstats <- select(sumstats, !direction)}
if("ngt" %in% colnames(sumstats)){               sumstats <- select(sumstats, !ngt)}
if("hetisqt" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetisqt)}
if("hetdf" %in% colnames(sumstats)){             sumstats <- select(sumstats, !hetdf)}
if("hetpval" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetpval)}

# Odds ratio -------------------------------------------------------------------------------------------------------------
# If reported effect size is odds ratio
if("or" %in% colnames(sumstats)){
  sumstats$beta = with(sumstats, log(or))
  sumstats <- select(sumstats, !or)
  
  # If SE already there, check that derived and reported p-values match. Otherwise recompute beta_se.
  if("beta_se" %in% colnames(sumstats)) {
    z <- with(sumstats, beta/beta_se)
    derived_pval <- 2 * (1 - pnorm(abs(z)))
    pval_cor <- cor(derived_pval, sumstats$p) # Maybe other comparison method? 
    
    # If cor is poor, the reported SE is not SE of log(OR), but probably SE of OR. Therefore, compute SE of beta
    if(pval_cor < 0.9) { # Fitting threshold?
      sumstats$beta_se = with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
    }
  # if SE not there, estimate with beta/z
  } else if (!"beta_se" %in% colnames(sumstats)) {
    sumstats$beta_se = with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
  }
}

# Finding Hapmap overlap with sumstats -----------------------------------------------------------------------------------
# Reading in HapMap3+ 
info <- readRDS(runonce::download_file(
  "https://figshare.com/ndownloader/files/37802721",
  dir = paths$hapmap_path, fname = "map_hm3_plus.rds"))

# Making sure a0 and a1 are upper case
sumstats$a0 <- toupper(sumstats$a0)
sumstats$a1 <- toupper(sumstats$a1)

# Finding sumstats/HapMap3+ overlap
if(all(c("chr", "pos") %in% colnames(sumstats))) {
  snp_info <- snp_match(sumstats, info, match.min.prop = 0.1) %>%
    select(-c(pos_hg18, pos_hg38))
} else { # Some files have only SNP/rsid - attempt to snp_match without pos then
  snp_info <- snp_match(sumstats, info, match.min.prop = 0.1, join_by_pos = FALSE) %>%
    select(-c(pos_hg18, pos_hg38))
}


cat(nrow(snp_info), "variants in overlap with HapMap3+. \n")

# Allele frequency -------------------------------------------------------------------------------------------------------

# Check if frq columns are on the form frq_a_X and frq_u_X (PGC format)
frq_cas_col <- grep("^fr?q_a_", colnames(snp_info))
frq_con_col <- grep("^fr?q_u_", colnames(snp_info))

col_cas <- as.numeric(str_extract(colnames(snp_info)[frq_cas_col], "\\d+"))
col_con <- as.numeric(str_extract(colnames(snp_info)[frq_con_col], "\\d+"))

colnames(snp_info)[frq_cas_col] <- "frq_cas"
colnames(snp_info)[frq_con_col] <- "frq_con"

# Calc frq and weight by number of cases and controls
if(!"frq" %in% colnames(snp_info)){

  # Start by looking for frqs split between cases and controls
  if("frq_cas" %in% colnames(snp_info)) {

    # Calculating weighted mean of the two frequencies using n_cas n_con
    if("n_cas" %in% colnames(snp_info)) {
      snp_info$frq = with(snp_info, (frq_cas * n_cas + frq_con * n_con) / (n_cas + n_con))

    # If n_cas n_con not there, look up in gwas catalog using gwasrapidd
    } else if(!is.na(study_info)) {
      snp_info$frq = with(snp_info, (frq_cas * num_inds$n_cas + frq_con * num_inds$n_con) / (num_inds$n_cas + num_inds$n_con))

    # Otherwise use cas con from frq cols (PGC format)
    } else if(length(frq_cas_col) > 0) {
      snp_info$frq = with(snp_info, (frq_cas * col_cas + frq_con * col_con) / (col_cas + col_con))
    }
      
  # Removing frq_cas frq_con afterwards
  snp_info <- select(snp_info, -c(frq_cas, frq_con))

  # If frq_cas frq_con not in sumstats, use af_UKBB instead
  } else { 
    snp_info$frq = snp_info$af_UKBB
  }
}

# Effective population size ------------------------------------------------------------------------------------------------------
if(!"n_eff" %in% colnames(snp_info)){

  # Prioritizing estimation from neff_half
  if("neff_half" %in% colnames(snp_info)){
    snp_info$n_eff = with(snp_info, neff_half * 2)

    # Deleting col afterwards - with n_cas n_con if present
    snp_info <- select(snp_info, !neff_half)
    if("n_cas" %in% colnames(snp_info)){snp_info <- select(snp_info, -c(n_cas, n_con))}

  # If no neff_half, check for n_cas n_con  
  } else if("n_cas" %in% colnames(snp_info)){
    snp_info$n_eff = with(snp_info, 4/(1/n_cas + 1/n_con))

    # Then delete cols
    snp_info <- select(snp_info, -c(n_cas, n_con))

  # Otherwise, look up in gwas catalog using gwasrapidd  
  } else if(!is.na(study_info)){
    snp_info$n_eff = with(num_inds, 4/(1/n_cas + 1/n_con))

  # Last resort using numbers from frq_a_cas frq_u_con (PGC)
  } else if(length(frq_cas_col) > 0){
    snp_info$n_eff = 4/(1/col_cas + 1/col_con)
  }
}

# Total population size for continuous traits ------------------------------------------------------------------------------------------
if(!"n" %in% colnames(snp_info)) {
  snp_info$n = if(!is.na(study_info)){
    sum(study_info@ancestries$number_of_individuals)
  }
}

# Making sure its in the sumstats 
# - otherwise should be added manually
assert("No effective population size in parsed sumstats",
       "n_eff" %in% colnames(snp_info) | "n" %in% colnames(snp_info))

# Z score -------------------------------------------------------------------------------------------------------------------------------
# if(!"beta" %in% colnames(snp_info) & "z" %in% colnames(snp_info)){
# TODO: What I would do: get beta_se from n_eff and freq
# 2 * frq (1 - frq) ~ 4 / (n_eff * beta_se^2)
# get beta from p-val -> |z| and beta_se and sign(z)
#snp_info$beta = with(snp_info, z / sqrt(2*frq*(1-frq)(n_eff + z^2)))  
#snp_info$beta_se = with(snp_info, beta/qnorm(1-p/2)) # beta/z


# QC -------------------------------------------------------------------------------------------------------------------------------------
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

df_beta$frq2 <- ifelse(df_beta$beta * sumstats$beta[df_beta$`_NUM_ID_.ss`] < 0,
                          1 - df_beta$frq, df_beta$frq)
diff <- with(df_beta, abs(af_UKBB - frq2))
  
df_beta$is_bad <- with(df_beta,
                       diff > 0.05 |
                       sd_ss2 < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                       sd_ss2 < 0.1 | sd_af < 0.05)

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

df_beta <- df_beta %>% filter(!is_bad) %>% select(-c(sd_af, sd_ss, sd_ss2, af_UKBB))
  

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
    N_Cases = num_inds$n_cas,
    N_Controls = num_inds$n_con,
    Ancestral_Group = study_info@ancestral_groups$ancestral_group,
    M_Input = nrow(sumstats),     # Variants in original sumstats
    M_HapMap = nrow(snp_info),    # Overlap with HapMap3+
    M_QC = nrow(df_beta)          # Variants after QC and iPSYCH overlap
  )
} else {
  foelgefil_df <- data.frame(
    ID = base_name,
    Accession_ID = NA,
    Reported_Trait = NA,
    PubMed_ID = NA,
    First_Author = NA,
    Journal = NA,
    Title = NA,
    Publication_Date = NA,
    N = NA,
    N_Cases = NA,
    N_Controls = NA,
    Ancestral_Group = NA,
    M_Input = nrow(sumstats),     # Variants in original sumstats
    M_HapMap = nrow(snp_info),    # Overlap with HapMap3+
    M_QC = nrow(df_beta)          # Variants after QC and iPSYCH overlap
  )
}
write.table(foelgefil_df,
           file = paste0(res_folder, "/", base_name, "_foelgefil.csv"),
           sep = "\t", row.names = FALSE, append = FALSE, quote = FALSE)
