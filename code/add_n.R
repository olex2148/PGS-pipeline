#' Script for adding n_cas and n_con to raw sumstats files

library(dplyr)
library(data.table)
library(MungeSumstats)

sumstats_path <- ""

sumstats <- read_sumstats(sumstats_path)

sumstats$n_cas <- 3000
sumstats$n_con <- 17000

write_sumstats(sumstats, sumstats_path)
