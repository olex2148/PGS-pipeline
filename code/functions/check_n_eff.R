check_n_eff <- function(sumstats, sample_size, model_info_df){

  if(!"n_eff" %in% colnames(sumstats)){

    # Prioritizing estimation from neff_half
    if("neff_half" %in% colnames(sumstats)){
      sumstats$n_eff = with(sumstats, neff_half * 2)
      
      # Saving the info before deleting later
      if(all(c("n_cas", "n_con") %in% colnames(sumstats))){
        model_info_df$frq_cas_median <- mean(sumstats$n_cas)
        model_info_df$frq_con_median <- mean(sumstats$n_con)
      }
      
      # If no neff_half, check for n_cas n_con  
    } else if(all(c("n_cas", "n_con") %in% colnames(sumstats))){
      sumstats$n_eff = with(sumstats, 4/(1/n_cas + 1/n_con))
      
      # Saving the info before deleting
      model_info_df$frq_cas_median <- mean(sumstats$n_cas)
      model_info_df$frq_con_median <- mean(sumstats$n_con)
      
      # Otherwise, look up in gwas catalog using gwasrapidd, or if added in check_frq  
    } else if(!is.na(sample_size$n_eff)){
      sumstats$n_eff = sample_size$n_eff
      
    } else if(!is.na(sample_size$n_cas)){
      sumstats$n_eff = with(sample_size, 4/(1/n_cas + 1/n_con))
      
      # Last resort using numbers from frq_a_cas frq_u_con (PGC)
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
