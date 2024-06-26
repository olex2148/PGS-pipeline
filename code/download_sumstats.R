library(httr)
library(rvest)
library(progress)
library(dplyr)
library(stringr)
library(readr)
library(bigreadr)
library(rjson)

paths <- rjson::fromJSON(file = "data/paths.json")
source(paths$get_sample_size)

runonce::download_file("https://www.ebi.ac.uk/gwas/api/v2/summaryStatistics/studies/download",
                       "data/",
                       "summary_statistics_table_export.tsv",
                       overwrite = TRUE)
gwascatalog <- bigreadr::fread2("data/summary_statistics_table_export.tsv")
gwascatalog$publicationDate <- as.Date(gwascatalog$publicationDate)

extract_european_sample_size <- function(text) {
  if (is.na(text)) return(NA)
  number <- stringr::str_extract(text, "\\d+(?= European)")
  if (is.na(number)) return(NA)
  as.numeric(number)
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
  filtered_gwascatalog$n_eff <- 
  filtered_gwascatalog$n_bin <- NA


for (i in seq(nrow(filtered_gwascatalog))) {
  parsed <- get_sample_size(filtered_gwascatalog[i, ]["initialSampleDescription"], filtered_gwascatalog[i, ]["replicateSampleDescription"])
  filtered_gwascatalog[i, ]$n_cont <- parsed$n
  filtered_gwascatalog[i, ]$n_cas <- parsed$n_cas
  filtered_gwascatalog[i, ]$n_con <- parsed$n_con
  filtered_gwascatalog[i, ]$n_eff <- parsed$n_eff
  filtered_gwascatalog[i, ]$n_bin <- parsed$n_bin
}

filtered_gwascatalog <- filtered_gwascatalog %>% filter(n_eff > 10000 | n_cont > 10000)

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
saveRDS(sumstats_metadata, "data/sumstats_metadata.rds")

options(timeout=3600)
for (i in seq_along(sumstats$url[1:length(seq_along(sumstats$url))])) {
  url <- paste0(sumstats$url[i], sumstats$filename[i])
  dir <- file.path("data/sumstats/gwascatalog")
  fname <- paste0("accession", (regmatches(sumstats$filename[i], regexpr("GCST\\d+", sumstats$filename[i]))), "_.gz")
  runonce::download_file(url, dir, fname, overwrite = FALSE)
  cat(paste0(i, "/", nrow(sumstats), "\n"))
}