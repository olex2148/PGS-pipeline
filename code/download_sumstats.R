library(httr)
library(rvest)
library(progress)
library(dplyr)
library(stringr)
library(readr)
library(bigreadr)
# 
# paths <- fromJSON(file = "data/paths.json")
# source(paths$get_n)
# source("code\\functions\\get_n.R")



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
  filter(!is.na(europeanSampleSize) & europeanSampleSize > 20000) %>% 
  filter(associationCount > 0) %>% 
  filter(!(genotypingTechnologies %in% 
             c("Exome genotyping array", "Exome genotyping array,Exome-wide sequencing")))

sample_size <- t(sapply(filtered_gwascatalog$initialSampleDescription, get_n))
rownames(sample_size) <- NULL
sample_size <- data.frame(sample_size)
filtered_gwascatalog <- filtered_gwascatalog %>% mutate(n_cont = apply(sample_size, 1, function(row) with_NAs(row[["n"]], 0, exp = function(a, b) a)),
                                                        n_cas = apply(sample_size, 1, function(row) with_NAs(row[["n_cas"]], 0, exp = function(a, b) a)),
                                                        n_con = apply(sample_size, 1, function(row) with_NAs(row[["n_con"]], 0, exp = function(a, b) a)),
                                                        n_cc = apply(sample_size, 1, function(row) with_NAs(row[["n_cas"]], row[["n_con"]], exp = function(a, b) a + b)),
                                                        n_eff = apply(sample_size, 1, function(row) with_NAs(row[["n_cas"]], row[["n_con"]], exp = function(a, b) 4 / (1 / a + 1 / b))))



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

options(timeout=3600)
for (i in seq_along(sumstats$url[1:length(seq_along(sumstats$url))])) {
  url <- paste0(sumstats$url[i], sumstats$filename[i])
  dir <- file.path("<INSERT PATH HERE>")
  fname <- paste0("accession", (regmatches(sumstats$filename[i], regexpr("GCST\\d+", sumstats$filename[i]))), "_.gz")
  runonce::download_file(url, dir, fname, overwrite = FALSE)
  cat(paste0(i, "/", nrow(sumstats), "\n"))
}