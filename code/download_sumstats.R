library(httr)
library(rvest)
library(progress)
library(dplyr)
library(stringr)
library(readr)
library(bigreadr)

gwascatalog <- bigreadr::fread2("summary_statistics_table_export.tsv")
gwascatalog$publicationDate <- as.Date(gwascatalog$publicationDate)

extract_european_sample_size <- function(text) {
  if (is.na(text)) return(NA)
  number <- str_extract(text, "\\d+(?= European)")
  if (is.na(number)) return(NA)
  as.numeric(number)
}

gwascatalog$europeanSampleSize <- sapply(gwascatalog$discoverySampleAncestry, extract_european_sample_size)

filtered_gwascatalog <- gwascatalog %>%
  filter(publicationDate >= as.Date("2014-01-01")) %>%
  filter(!is.na(europeanSampleSize) & europeanSampleSize > 20000) %>% 
  filter(associationCount > 0) %>% 
  filter(!(genotypingTechnologies %in% 
             c("Exome genotyping array", "Exome genotyping array,Exome-wide sequencing")))

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