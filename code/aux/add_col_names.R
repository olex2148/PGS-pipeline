#' Summary Statistics Column Headers

#' # Make additions to sumstatsColHeaders using github version of MungeSumstats-
#' # shown is an example of adding columns for Standard Error (SE)

data("sumstatsColHeaders")
new_cols <- data.frame("Uncorrected"=c("NCAS",   "NCON", "FCAS",   "EAF_CASES", "FCON", "EAF_CONTROLS", "NEFFDIV2", "IMPINFO", "INFOSCORE", "DIRE", "ALT_ALLELE", "b", "RISK_ALLELE"),
                        "Corrected"= c("N_CAS", "N_CON", "FRQ_CAS", "FRQ_CAS", "FRQ_CON", "FRQ_CON",     "NEFF_HALF", "INFO",  "INFO",   "DIRECTION", "A1", "beta", "A2"))
sumstatsColHeaders <- rbind(sumstatsColHeaders,new_cols)

#Once additions are made, order & save the new mapping dataset
#now sort ordering -important for logic that
# uncorrected=corrected comes first

sumstatsColHeaders$ordering <-
    sumstatsColHeaders$Uncorrected==sumstatsColHeaders$Corrected
sumstatsColHeaders <-
    sumstatsColHeaders[order(sumstatsColHeaders$Corrected,
                             sumstatsColHeaders$ordering,decreasing = TRUE),]
rownames(sumstatsColHeaders)<-1:nrow(sumstatsColHeaders)
sumstatsColHeaders$ordering <- NULL
#' #manually move FRWQUENCY to above MAR - github issue 95
frequency <- sumstatsColHeaders[sumstatsColHeaders$Uncorrected=="FREQUENCY",]
maf <- sumstatsColHeaders[sumstatsColHeaders$Uncorrected=="MAF",]
if(as.integer(rownames(frequency))>as.integer(rownames(maf))){
  sumstatsColHeaders[as.integer(rownames(frequency)),] <- maf
  sumstatsColHeaders[as.integer(rownames(maf)),] <- frequency
}
#' usethis::use_data(sumstatsColHeaders,overwrite = TRUE, internal=TRUE)
save(sumstatsColHeaders,
      file="data/sumstatsColHeaders.rda")
