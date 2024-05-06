logs_path <- ".gwf/logs/"
filenames <- list.files(logs_path)

gwf_status_path <- "old logs/gwf status latest"
gwf_status <- readlines(gwf_status_path)

failed <- gwf_status[which(regexpr("failed", gwf_status) > 0)]
failed_ids <- regmatches(failed, regexpr("GCST\\d+", failed))

failed_munge <- failed[which(regexpr("munge", failed) > 0)]
failed_munge_ids <- regmatches(failed_munge, regexpr("GCST\\d+", failed_munge))

failed_ldpred2 <- failed[which(regexpr("ldpred2", failed) > 0)]
failed_ldpred2_ids <- regmatches(failed_ldpred2, regexpr("GCST\\d+", failed_ldpred2))

missing_beta <- c()
not_enough_variants <- c()
no_n <- c()
timedout <- c()

for (id in failed_munge_ids) {
  path <- paste0(logs_path, "munge_accession", id, "_.stderr")
  if (file.exists(path)) {
    log <- readLines(path)
    if (regexpr("Expected 'chr, rsid, a0, a1, beta'", log[length(log) -1]) > 0) {
      missing_beta <- c(missing_beta, id)
    }
    if (regexpr("Not enough variants have been matched", log[length(log) -1]) > 0) {
      not_enough_variants <- c(not_enough_variants, id)
    }
    if (regexpr("nrow\\(df_beta\\) > 10000 is not TRUE", log[length(log) -1]) > 0) {
      not_enough_variants <- c(not_enough_variants, id)
    }
    if (regexpr("No variant has been matched", log[length(log) -1]) > 0) {
      not_enough_variants <- c(not_enough_variants, id)
    }
    if (regexpr("colnames", log[length(log) -1]) > 0) {
      no_n <- c(no_n, id)
    }
  }
}

for (id in failed_ldpred2_ids) {
  path <- paste0(logs_path, "ldpred2_accession", id, "_.stderr")
  if (file.exists(path)) {
    log <- readLines(path)
    if (regexpr("DUE TO TIME LIMIT", log[length(log) -1]) > 0) {
      timedout <- c(timedout, id)
    }
  }
}

print(paste0(length(missing_beta), "/", length(failed_munge), "failed munging due to missing betas"))
print(paste0(length(not_enough_variants), "/", length(failed_munge), "failed munging due to insufficient variants"))
print(paste0(length(no_n), "/", length(failed_munge)), "failed munging due to no available sample size")
print(paste0(length(timedout), "/", length(failed_ldpred2), "failed modelling due to timeout"))

other_failures <- failed_ids[!(failed_ids %in% Reduce(union, list(missing_beta, no_n, not_enough_variants, timedout)))]
print(paste0(other_failures), "/", length(failed_ids))