check_n_col <- function(sumstats, num_inds, foelgefil_df){
  # Check if frq columns are on the form frq_a_X and frq_u_X (PGC format)
  frq_cas_col <- grep("^fr?q_a_", colnames(sumstats))
  frq_con_col <- grep("^fr?q_u_", colnames(sumstats))
  
  col_cas <- as.numeric(str_extract(colnames(sumstats)[frq_cas_col], "\\d+"))
  col_con <- as.numeric(str_extract(colnames(sumstats)[frq_con_col], "\\d+"))
  
  colnames(sumstats)[frq_cas_col] <- "frq_cas"
  colnames(sumstats)[frq_con_col] <- "frq_con"

  if(!"n_eff" %in% colnames(sumstats)){

    # Prioritizing estimation from neff_half
    if("neff_half" %in% colnames(sumstats)){
      sumstats$n_eff = with(sumstats, neff_half * 2)
      
      # Deleting col afterwards - with n_cas n_con if present
      sumstats <- select(sumstats, !neff_half)
      if("n_cas" %in% colnames(sumstats)){
        
        # Saving the info before deleting
        foelgefil_df$N_Cases <- mean(sumstats$n_cas)
        foelgefil_df$N_Controls <- mean(sumstats$n_con)
        
        sumstats <- select(sumstats, -c(n_cas, n_con))
      }
      
      # If no neff_half, check for n_cas n_con  
    } else if("n_cas" %in% colnames(sumstats)){
      sumstats$n_eff = with(sumstats, 4/(1/n_cas + 1/n_con))
      
      # Saving the info before deleting
      foelgefil_df$N_Cases <- mean(sumstats$n_cas)
      foelgefil_df$N_Controls <- mean(sumstats$n_con)
      
      # Then delete cols
      sumstats <- select(sumstats, -c(n_cas, n_con))
      
      # Otherwise, look up in gwas catalog using gwasrapidd  
    } else if(!is.na(num_inds$n_cas)){
      sumstats$n_eff = with(num_inds, 4/(1/n_cas + 1/n_con))
      
      # Last resort using numbers from frq_a_cas frq_u_con (PGC)
    } else if(length(frq_cas_col) > 0){
      sumstats$n_eff = 4/(1/col_cas + 1/col_con)
      
      # Saving the info
      foelgefil_df$N_Cases <- col_cas
      foelgefil_df$N_Controls <- col_con
    }
      
    
    # Total population size for continuous traits 
    if(!"n" %in% colnames(sumstats) & !"n_eff" %in% colnames(sumstats)) {
      if(!is.na(num_inds$n)){
        sumstats$n =  num_inds$n
      }
    }
    
  }
  return(list("sumstats" = sumstats, "foelgefil" = foelgefil_df))
}