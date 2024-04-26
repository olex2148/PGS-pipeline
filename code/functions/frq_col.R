#' Making sure frq col exists
#' 
#' @author Ole SH
#' 
#' @name check_frq_col
#'
#' @param sumstats
#' @param study_info
#' @param num_inds
#' 
#' @returns sumstats with frq col

check_frq_col <- function(sumstats, study_info, num_inds){
  # Check if frq columns are on the form frq_a_X and frq_u_X (PGC format)
  frq_cas_col <- grep("^fr?q_a_", colnames(sumstats))
  frq_con_col <- grep("^fr?q_u_", colnames(sumstats))
  
  col_cas <- as.numeric(str_extract(colnames(sumstats)[frq_cas_col], "\\d+"))
  col_con <- as.numeric(str_extract(colnames(sumstats)[frq_con_col], "\\d+"))
  
  colnames(sumstats)[frq_cas_col] <- "frq_cas"
  colnames(sumstats)[frq_con_col] <- "frq_con"
  
  # Calc frq and weight by number of cases and controls
  if(!"frq" %in% colnames(sumstats)){
    
    # Start by looking for frqs split between cases and controls
    if("frq_cas" %in% colnames(sumstats)) {
      
      # Calculating weighted mean of the two frequencies using n_cas n_con
      if("n_cas" %in% colnames(sumstats)) {
        sumstats$frq = with(sumstats, (frq_cas * n_cas + frq_con * n_con) / (n_cas + n_con))
        
      # If n_cas n_con not there, look up in gwas catalog using gwasrapidd
      } else if(!is.na(study_info) & !is.na(num_inds$n_cas)) {
        sumstats$frq = with(sumstats, (frq_cas * num_inds$n_cas + frq_con * num_inds$n_con) / (num_inds$n_cas + num_inds$n_con))
        
        # Otherwise use cas con from frq cols (PGC format)
      } else if(length(frq_cas_col) > 0) {
        sumstats$frq = with(sumstats, (frq_cas * col_cas + frq_con * col_con) / (col_cas + col_con))
      }
      
      # Removing frq_cas frq_con afterwards
      sumstats <- select(sumstats, -c(frq_cas, frq_con))
    }
  }
  return(sumstats)
}