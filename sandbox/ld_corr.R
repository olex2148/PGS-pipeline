tmp <- tempfile(tmpdir = "data/corr", pattern = "hm3+_corr")

for (chr in 1:22) {
  
  cat(chr, "..", sep = "")
  
  corr_chr <- readRDS(paste0(ld_blocks, chr, ".rds"))
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}
corr$save()

# cat("Reading in LD blocks for chromosome ")
# for (chr in 1:22) {
#   
#   cat(chr, "..", sep = "")
#   
#   ## indices in df_beta
#   ind.chr <- which(df_beta$chr == chr)
#   ## indices in map_ldref
#   ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
#   ## indices in corr_chr
#   ind.chr3 <- match(ind.chr2, which(info$chr == chr))
#   
#   corr_chr <- readRDS(paste0(ld_blocks, chr, ".rds"))[ind.chr3, ind.chr3]
#   
#   if (chr == 1) {
#     corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
#   } else {
#     corr$add_columns(corr_chr, nrow(corr))
#   }
# }