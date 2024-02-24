#' Split gwasrapidd initial_sample_size strings into cases and controls
#' 

get_n_cas_con <- function(initial_sample_size_str) {
  require(stringr)
  
  # Removing commas in the numbers
  sample_size_str <- gsub("(\\d),(\\d)", "\\1\\2", initial_sample_size_str)
  
  # Split the string into parts
  parts <- unlist(strsplit(sample_size_str, ", "))
  
  # If theres only one sentence in this variable, continuous trait and will return only n
  if(length(parts) == 1) {
    n <- as.numeric(gsub("[^0-9]", "", parts))
    
    return(list("n" = n, "n_cas" = NA, "n_con" = NA))  
           
  } else { 
  # Otherwise return ncas ncon
    
    # Initialize
    sum_cases <- 0
    sum_controls <- 0
    
    # Loop through the parts and sum cases and controls separately
    for (part in parts) {
      if (grepl("cases$", part)) {
        sum_cases <- sum_cases + as.numeric(gsub("[^0-9]", "", part))
      } else if (grepl("controls$", part)) {
        sum_controls <- sum_controls + as.numeric(gsub("[^0-9]", "", part))
      }
    }
    return(list("n" = NA, "n_cas" = sum_cases, "n_con" = sum_controls))
  }
}