#' Small checks to fit snp_match() format
#' 
#' @author Ole SH
#' 
#' @name snp_match_format
#' @param sumstats
#' 
#' @returns sumstats with names snp_match requires, uppercase alleles + removal of some redundant cols

snp_match_format <- function(sumstats){
  # Renaming
  colnames(sumstats) <- tolower(colnames(sumstats))
  
  if(all(c("a1", "a2") %in% colnames(sumstats))){                        sumstats <- rename(sumstats, a0 = a1, a1 = a2)}
  if("bp" %in% colnames(sumstats)){                                      sumstats <- rename(sumstats, pos = bp)}
  if("snp" %in% colnames(sumstats) & !"rsid" %in% colnames(sumstats)){   sumstats <- rename(sumstats, rsid = snp)}
  if("se" %in% colnames(sumstats)){                                      sumstats <- rename(sumstats, beta_se = se)}
  
  # Deleting sex chromosomes
  if("chr" %in% colnames(sumstats)){                                     sumstats <- sumstats %>% filter(chr %in% 1:22) %>%  mutate(chr = as.numeric(chr))}
  
  # Making sure a0 and a1 are upper case
  sumstats$a0 <- toupper(sumstats$a0)
  sumstats$a1 <- toupper(sumstats$a1)
  
  # Removing some redundant cols
  if("direction" %in% colnames(sumstats)){         sumstats <- select(sumstats, !direction)}
  if("ngt" %in% colnames(sumstats)){               sumstats <- select(sumstats, !ngt)}
  if("hetisqt" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetisqt)}
  if("hetdf" %in% colnames(sumstats)){             sumstats <- select(sumstats, !hetdf)}
  if("hetpval" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetpval)}
  
  return(sumstats)
}