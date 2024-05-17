
create_model_info <- function(accession_id){
  model_info_df <- data.frame(
    id = NA,
    accession_id = NA,
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
    frq_cas = NA,
    frq_con = NA,
    m_input = NA,
    m_hapmap = NA,
    m_qc = NA          
  )
  
  if(!is.na(accession_id)) {
    study_info <- get_studies(study_id = accession_id)
    sample_size <- get_sample_size(study_info@studies$initial_sample_size, study_info@studies$replication_sample_size)
    # Saving variables
    model_info_df <- model_info_df %>% 
      mutate(
        accession_id = accession_id,
        reported_trait = study_info@studies$reported_trait,
        pubmed_id = study_info@publications$pubmed_id,
        first_author = study_info@publications$author_fullname,
        journal = study_info@publications$publication,
        title = study_info@publications$title,
        publication_date = as.Date(study_info@publications$publication_date, origin = "1970-01-01"),
        n = ifelse(!is.na(sample_size$n), sample_size$n, sample_size$n_bin),
        n_cas = sample_size$n_cas,
        n_con = sample_size$n_con,
        n_eff = sample_size$n_eff
      )
  }
  
  return(model_info_df)
}