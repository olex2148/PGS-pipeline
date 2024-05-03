#' Split gwasrapidd initial_sample_size strings into cases and controls
#' 

get_n_cas_con <- function(initial_sample_size_str) {
  require(stringr)
  
  # Removing commas in the numbers
  sample_size_str <- gsub("(\\d),(\\d)", "\\1\\2", initial_sample_size_str)
  
  # Split the string into parts
  parts <- unlist(strsplit(sample_size_str, "(,\\s*\\d)", perl = TRUE))
  
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
  n_cas <- ifelse(sum_cases == 0, NA, sum_cases)
  n_con <- ifelse(sum_controls == 0, NA, sum_controls)
  n <- ifelse(sum_inds == 0, NA, sum_inds)
  
  return(list("n" = n, "n_cas" = n_cas, "n_con" = n_con))

}