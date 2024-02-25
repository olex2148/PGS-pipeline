#' Some renaming
#' 
#' @author Ole SH
#' 
#' @name snp_match_names
#' @param sumstats
#' 
#' @returns sumstats with names snp_match requires + removal of some redundant cols

snp_match_names <- function(sumstats){
  colnames(sumstats) <- tolower(colnames(sumstats))
  
  if(all(c("a1", "a2") %in% colnames(sumstats))){                        sumstats <- rename(sumstats, a0 = a1, a1 = a2)}
  if("bp" %in% colnames(sumstats)){                                      sumstats <- rename(sumstats, pos = bp)}
  if("snp" %in% colnames(sumstats) & !"rsid" %in% colnames(sumstats)){   sumstats <- rename(sumstats, rsid = snp)}
  if("se" %in% colnames(sumstats)){                                      sumstats <- rename(sumstats, beta_se = se)}
  
  # Removing some redundant cols
  if("direction" %in% colnames(sumstats)){         sumstats <- select(sumstats, !direction)}
  if("ngt" %in% colnames(sumstats)){               sumstats <- select(sumstats, !ngt)}
  if("hetisqt" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetisqt)}
  if("hetdf" %in% colnames(sumstats)){             sumstats <- select(sumstats, !hetdf)}
  if("hetpval" %in% colnames(sumstats)){           sumstats <- select(sumstats, !hetpval)}
  
  return(sumstats)
}