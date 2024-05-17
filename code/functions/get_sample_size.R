#' Split gwasrapidd initial_sample_size strings into cases and controls
#' 

get_sample_size <- function(initial_sample_size_str, replication_sample_size_str) {
  require(stringr)
  
  parse <- function(sample_size_str) {
    # Removing commas in the numbers
    sample_size_str <- gsub("(\\d),(\\d)", "\\1\\2", sample_size_str)
    
    # Split the string into parts
    parts <- unlist(strsplit(sample_size_str, ",", perl = TRUE))
    
    # Initialize
    sum_cases <- sum_controls <- sum_inds <- 0
    
    # Loop through the parts and sum cases and controls separately
    for (part in parts) {
      if(grepl("cases", part)) {
        sum_cases <- sum_cases + as.numeric(gsub("[^0-9]", "", part))
      } else if(grepl("controls", part)) {
        sum_controls <- sum_controls + as.numeric(gsub("[^0-9]", "", part))
      } else if(grepl("(individuals|men|women|males|females)", part)) {
        sum_inds <- sum_inds + as.numeric(gsub("[^0-9]", "", part))
      }
    }
    
    return(list("n" = sum_inds, "n_cas" = sum_cases, "n_con" = sum_controls))
  }
  
  initial <- parse(initial_sample_size_str)
  replication <- parse(replication_sample_size_str)
  
  n <- initial["n"] + replication["n"]
  n_cas <- initial["n_cas"] + replication["n_cas"]
  n_con <- initial["n_con"] + replication["n_con"]
  n_eff <- 4 / (1 / n_cas + 1 / n_con)
  n_bin <- n_cas + n_con
  
  res <- list("n" = n, "n_cas" = n_cas, "n_con" = n_con, "n_eff" = n_eff, "n_bin" = n_bin)
  return(lapply(res, function(elem) ifelse(elem == 0, NA, elem)))
}