
# Command line arguments for this script
args = commandArgs(trailingOnly = TRUE)
model_info <- args[1]
model_info_qc <- args[2]

#Read the csv
model_info <- read.delim(model_info, check.names = FALSE, stringsAsFactors = FALSE)


#QC FUNCTION


check_ldpred2_qc <- function(model_info) {

  qc <- model_info
  
  # 1. INITIAL HERITABILITY
  qc$h2_init_qc <- ifelse(
    is.na(qc$h2_init), "WARN_missing",
    ifelse(
      qc$h2_init <= 0 | qc$h2_init > 1,
      "FAIL",
      "PASS"
    )
  )

  # 2. HERITABILITY POWER / PRECISION
  qc$h2_power_qc <- ifelse(
    is.na(qc$h2_init_z), "WARN_missing",
    ifelse(
      qc$h2_init_z < 2, "WARN_low_power",
      ifelse(
        qc$h2_init_z < 4, "WARN_moderate_power",
        "PASS"
      )
    )
  )


  # 3. LDPRED2-AUTO HERITABILITY
  qc$h2_auto_qc <- ifelse(
    is.na(qc$h2_auto), "WARN_missing",
    ifelse(
      qc$h2_auto <= 0 | qc$h2_auto > 1,
      "FAIL",
      "PASS"
    )
  )

  # 4. NUMBER OF RETAINED CHAINS(out of 50)
  qc$chains_qc <- ifelse(
    is.na(qc$n_chains), "WARN_missing",
    ifelse(
      qc$n_chains < 25, "FAIL",
      ifelse(
        qc$n_chains < 30, "WARN",
        "PASS"
      )
    )
  )

  # 5. LD SHRINKAGE
  qc$shrinkage_qc <- ifelse(
    is.na(qc$shrinkage), "WARN_missing",
    ifelse(
      qc$shrinkage >= 0.8,
      "PASS",
      "WARN_strong_shrinkage"
    )
  )

  # 6. POLYGENICITY
  qc$p_auto_qc <- ifelse(
    is.na(qc$p_auto), "WARN_missing",
    ifelse(
      qc$p_auto <= 0 | qc$p_auto > 1,
      "FAIL",
      "PASS"
    )
  )

  # 7. LDSC RATIO


  if ("ldsc_ratio" %in% colnames(qc)) {

    qc$ldsc_ratio_qc <- ifelse(
      is.na(qc$ldsc_ratio), "WARN_missing",
      ifelse(
        qc$ldsc_ratio > 0.20, "FAIL",
        ifelse(
          qc$ldsc_ratio >= 0.10,
          "WARN", "PASS")
      )
    )

  } else {

    qc$ldsc_ratio_qc <- "WARN_missing"

  }


  # 8. OVERALL QC

  qc_checks <- c(
    "h2_init_qc",
    "h2_power_qc",
    "h2_auto_qc",
    "chains_qc",
    "shrinkage_qc",
    "p_auto_qc",
    "ldsc_ratio_qc"
  )


  qc$overall_qc <- apply(
    qc[, qc_checks, drop = FALSE],
    1,
    function(x) {
      if (any(grepl("^FAIL", x))) return("FAIL")
      if (any(grepl("^WARN", x))) return("WARN")
      return("PASS")
    }
  )


# REORDER QC COLUMNS
move_after <- function(data, column, after) {
  if (column %in% names(data) && after %in% names(data)) {
    cols <- names(data)
    cols <- cols[cols != column]
    pos <- match(after, cols)
    cols <- append(cols, column, after = pos)
    data <- data[, cols, drop = FALSE]
  }

  return(data)
}


qc <- move_after(qc, "h2_init_qc",    "h2_init")
qc <- move_after(qc, "h2_power_qc",   "h2_init_z")
qc <- move_after(qc, "h2_auto_qc",    "h2_auto")
qc <- move_after(qc, "chains_qc",     "n_chains")
qc <- move_after(qc, "shrinkage_qc",  "shrinkage")
qc <- move_after(qc, "p_auto_qc",     "p_auto")
qc <- move_after(qc, "ldsc_ratio_qc", "ldsc_ratio")

#Put overall last
qc <- qc[, c (setdiff(names(qc), "overall_qc"),
            "overall_qc"),
            drop = FALSE]
      return(qc)
}


#Add the QC columns
qc <- check_ldpred2_qc(model_info)

#Save

write.table( qc, file = model_info_qc,
    sep= "\t", row.names = FALSE,
    quote = FALSE, na = "NA")

