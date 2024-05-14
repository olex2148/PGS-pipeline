library(httr)
library(rvest)
library(progress)
library(dplyr)
library(stringr)
library(readr)
library(bigreadr)

paths <- fromJSON(file = "data/paths.json")
source(paths$get_n)

runonce::download_file("https://www.ebi.ac.uk/gwas/api/v2/summaryStatistics/studies/download",
                       "data/",
                       "summary_statistics_table_export.tsv",
                       overwrite = TRUE)
gwascatalog <- bigreadr::fread2("data/summary_statistics_table_export.tsv")
gwascatalog$publicationDate <- as.Date(gwascatalog$publicationDate)

extract_european_sample_size <- function(text) {
  if (is.na(text)) return(NA)
  number <- str_extract(text, "\\d+(?= European)")
  if (is.na(number)) return(NA)
  as.numeric(number)
}

with_NAs <- function(a, b, exp) {
  if (!(is.na(a) || is.na(b))) {
    exp(a, b)
  } else NA
}

gwascatalog$europeanSampleSize <- sapply(gwascatalog$discoverySampleAncestry, extract_european_sample_size)

filtered_gwascatalog <- gwascatalog %>%
  filter(publicationDate >= as.Date("2014-01-01")) %>%
  filter(!is.na(europeanSampleSize) & europeanSampleSize > 10000) %>%
  filter(associationCount > 0) %>%
  filter(!(genotypingTechnologies %in% 
             c("Exome genotyping array", "Exome genotyping array,Exome-wide sequencing")))

filtered_gwascatalog$n_cont <-
  filtered_gwascatalog$n_cas <-
  filtered_gwascatalog$n_con <-
  filtered_gwascatalog$n_total_bin <- 
  filtered_gwascatalog$n_eff_bin <- NA


for (i in seq(nrow(filtered_gwascatalog))) {
  parsed <- get_n(filtered_gwascatalog[i, ]["initialSampleDescription"], filtered_gwascatalog[i, ]["replicateSampleDescription"])
  n <- unlist(parsed["n"])
  n_cas <- unlist(parsed["n_cas"])
  n_con <- unlist(parsed["n_con"])
  filtered_gwascatalog[i, ]$n_cont <- n
  filtered_gwascatalog[i, ]$n_cas <- n_cas
  filtered_gwascatalog[i, ]$n_con <- n_con
  filtered_gwascatalog[i, ]$n_total_bin <- with_NAs(n_cas, n_con, function (a, b) a+b)
  filtered_gwascatalog[i, ]$n_eff_bin <- with_NAs(n_cas, n_con, function (a, b) 4 / (1 / b + 1 / b)) 
}

urls <- filtered_gwascatalog$summaryStatistics

fetch <- function(urls) {
  contents <- c()
  for (url in urls) {
    response <- httr::GET(url)
    text <- rvest::html_text(httr::content(response))
    contents <- append(contents, text)
  }
  return(data.frame(
    url = urls,
    content = contents
  ))
}

fetched <- fetch(paste0(urls, "/harmonised/"))
filenames <- "(\\d+-)?GCST\\d+\\S+\\.h\\.tsv\\.gz"
matches <- regexpr(filenames, fetched$content)
sumstats <- data.frame(
  url = fetched$url[which(matches > 0)],
  filename = regmatches(fetched$content, matches)
)
sumstats_metadata <- filtered_gwascatalog[which(matches > 0), ]

options(timeout=3600)
for (i in seq_along(sumstats$url[1:length(seq_along(sumstats$url))])) {
  url <- paste0(sumstats$url[i], sumstats$filename[i])
  dir <- file.path("<INSERT PATH HERE>")
  fname <- paste0("accession", (regmatches(sumstats$filename[i], regexpr("GCST\\d+", sumstats$filename[i]))), "_.gz")
  runonce::download_file(url, dir, fname, overwrite = FALSE)
  cat(paste0(i, "/", nrow(sumstats), "\n"))
}