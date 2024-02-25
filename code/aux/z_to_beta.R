#' Z to beta function
#' 
#' From MungeSumstats R package
#' 

z_to_beta <- function(sumstats){
  
  if(all(c("z", "p", "n")) %in% colnames(sumstats)){
    cat("Deriving BETA from Z, N, and P")
    sumstats$beta <- with(sumstats, z/sqrt(qchisq(p, n)))

  } else if (all(c("z", "n", "frq")) %in% colnames(sumstats)){
    cat("Deriving BETA from Z, N, and FRQ")
    sumstats$beta <- with(sumstats, z/sqrt(2 * frq * (1 - frq) * (n + z^2)))
  }
}