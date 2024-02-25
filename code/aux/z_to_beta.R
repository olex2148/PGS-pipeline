#' Z to beta function
#' 
#' From MungeSumstats R package
#' 

z_to_beta <- function(sumstats){
  if("z" %in% colnames(sumstats) & !"beta" %in% colnames(sumstats)) {
    # In case binary pheno, renaming n_eff to n to ease of following checks
    if("n_eff" %in% colnames(sumstats)){sumstats$n = sumstats$n_eff}
    
    if(all(c("z", "p", "n") %in% colnames(sumstats))){
      cat("Deriving BETA from Z, N, and P")
      sumstats$beta <- with(sumstats, z/sqrt(qchisq(p, n)))
      
    } else if (all(c("z", "n", "frq") %in% colnames(sumstats))){
      cat("Deriving BETA from Z, N, and FRQ")
      sumstats$beta <- with(sumstats, z/sqrt(2 * frq * (1 - frq) * (n + z^2)))
    }
    
    sumstats$beta_se <- with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
  }
  return(sumstats)
}