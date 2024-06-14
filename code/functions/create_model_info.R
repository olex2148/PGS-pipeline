
create_model_info <- function(accession_id){
  model_info_df <- data.frame(
    id = NA,
    accession_id = accession_id,
    reported_trait = NA,
    pubmed_id = NA,
    first_author = NA,
    journal = NA,
    title = NA,
    publication_date = NA,
    n = NA,
    n_cas = NA,
    n_con = NA,
    n_eff = NA,
    frq_cas_median = NA,
    frq_con_median = NA,
    m_input = NA,
    m_hapmap = NA,
    m_qc = NA          
  )
  
  # if using gwasrapidd
  # if(!is.na(accession_id)) {
  #   study_info <- get_studies(study_id = accession_id)
  #   sample_size <- get_sample_size(study_info@studies$initial_sample_size, study_info@studies$replication_sample_size)
  #   # Saving variables
  #   model_info_df <- model_info_df %>% 
  #     mutate(
  #       accession_id = accession_id,
  #       reported_trait = study_info@studies$reported_trait,
  #       pubmed_id = study_info@publications$pubmed_id,
  #       first_author = study_info@publications$author_fullname,
  #       journal = study_info@publications$publication,
  #       title = study_info@publications$title,
  #       publication_date = as.Date(study_info@publications$publication_date, origin = "1970-01-01"),
  #       n = ifelse(!is.na(sample_size$n), sample_size$n, sample_size$n_bin),
  #       n_cas = sample_size$n_cas,
  #       n_con = sample_size$n_con,
  #       n_eff = sample_size$n_eff
  #     )
  # }
  
  # if using locally stored gwascatalog metadata
  # metadata <- readRDS("data/sumstats_metadata.rds")
  # info <- metadata[metadata$accessionId == accession_id, ]
  # 
  # model_info_df <- model_info_df %>% 
  #       mutate(
  #         reported_trait = info$reportedTrait,
  #         pubmed_id = info$pubmedId,
  #         first_author = info$firstAuthor,
  #         journal = info$journal,
  #         title = info$title,
  #         publication_date = as.Date(info$publicationDate, origin = "1970-01-01"),
  #         n = ifelse(!is.na(info$n_cont), info$n_cont, info$n_bin),
  #         n_cas = info$n_cas,
  #         n_con = info$n_con,
  #         n_eff = info$n_eff
  #       )
  
  return(model_info_df)
}