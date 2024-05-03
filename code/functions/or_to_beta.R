#' From OR to beta and beta_se
#' 
#' @author Ole SH
#' 
#' @name or_to_beta
#' @param sumstats
#' 
#' @returns sumstats with beta and beta_se

or_to_beta <- function(sumstats){
  
  # Checking for OR
  if("or" %in% colnames(sumstats)){
    cat("Odds ratio found in sumstats, imputing beta.")
    sumstats <- sumstats[complete.cases(sumstats$or), ]
    sumstats$beta = with(sumstats, log(or))
    sumstats <- select(sumstats, !or) # Removing OR col
    
    # If SE already there, check that derived and reported p-values match. Otherwise recompute beta_se.
    if("beta_se" %in% colnames(sumstats)) {
      z <- with(sumstats, beta/beta_se)
      derived_pval <- 2 * (1 - pnorm(abs(z)))
      
      if (!any(is.na(derived_pval)) && !any(is.na(sumstats$p))) {
        pval_cor <- cor(derived_pval, sumstats$p) # Maybe other comparison method? 
      } else {
        pval_cor <- NA
      }
      
      # If cor is poor, the reported SE is not SE of log(OR), but probably SE of OR. Therefore, compute SE of beta
      if(is.na(pval_cor) || pval_cor < 0.9) { # Fitting threshold?
        sumstats$beta_se = with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
      }
      # if SE not there, estimate with beta/z
    } else if (!"beta_se" %in% colnames(sumstats)) {
      sumstats$beta_se = with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
    }
  }
  # If effect size is already beta, but beta_se is not there
  if(!"beta_se" %in% colnames(sumstats) & "beta" %in% colnames(sumstats)) {
    sumstats$beta_se = with(sumstats, abs(beta) / qnorm(pmax(p, .Machine$double.xmin) / 2, lower.tail = FALSE)) # beta/z
  }
  return(sumstats)
}