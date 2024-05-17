check_n_col <- function(sumstats, sample_size, model_info_df){
  col_cas <- as.numeric(str_extract(colnames(sumstats)[frq_cas], "\\d+"))
  col_con <- as.numeric(str_extract(colnames(sumstats)[frq_con], "\\d+"))

  if(!"n_eff" %in% colnames(sumstats)){

    # Prioritizing estimation from neff_half
    if("neff_half" %in% colnames(sumstats)){
      sumstats$n_eff = with(sumstats, neff_half * 2)
      
      # Saving the info before deleting later
      if(all(c("n_cas", "n_con") %in% colnames(sumstats))){
        model_info_df$frq_cas <- mean(sumstats$n_cas)
        model_info_df$frq_con <- mean(sumstats$n_con)
      }
      
      # If no neff_half, check for n_cas n_con  
    } else if(all(c("n_cas", "n_con") %in% colnames(sumstats))){
      sumstats$n_eff = with(sumstats, 4/(1/n_cas + 1/n_con))
      
      # Saving the info before deleting
      model_info_df$frq_cas <- mean(sumstats$n_cas)
      model_info_df$frq_con <- mean(sumstats$n_con)
      
      # Otherwise, look up in gwas catalog using gwasrapidd  
    } else if(!is.na(sample_size$n_eff)){
      sumstats$n_eff = sample_size$n_eff
      
      # Last resort using numbers from frq_a_cas frq_u_con (PGC)
    } else if(all(is.na(c(col_cas, col_con)))){
      sumstats$n_eff = 4/(1/col_cas + 1/col_con)
      
      # Saving the info
      model_info_df$frq_cas <- col_cas
      model_info_df$frq_con <- col_con
    }
      
    
    # Total population size for continuous traits 
    if(!"n" %in% colnames(sumstats) & !"n_eff" %in% colnames(sumstats)) {
      if(!is.na(sample_size$n)){
        sumstats$n =  sample_size$n
      }
    }
    
  }
  return(list("sumstats" = sumstats, "model_info" = model_info_df))
}