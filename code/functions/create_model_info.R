
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
    n_cases = NA,
    n_controls = NA,
    m_input = NA,
    m_hapmap = NA,
    m_qc = NA          
  )
  
  if(!is.na(accession_id)) {
    study_info <- get_studies(study_id = accession_id)
    num_inds <- get_n_cas_con(study_info@studies$initial_sample_size)
    # Saving variables
    model_info_df <- model_info_df %>% 
      mutate(
        accession_id = accession_id,
        reported_trait = study_info@studies$reported_trait,
        pubmed_id = study_info@publications$pubmed_id,
        first_author = study_info@publications$author_fullname,
        journal = study_info@publications$publication,
        title = study_info@publications$title,
        publication_date = study_info@publications$publication_date
      )
    
    if(!is.na(num_inds$n)){
      model_info_df <- model_info_df %>% 
        mutate(
          n = num_inds$n,
          n_cases = num_inds$n_cas,
          n_controls = num_inds$n_con
        )
    } else {
      model_info_df <- model_info_df %>% 
        mutate(
          n = num_inds$n_cas + num_inds$n_con,
          n_cases = num_inds$n_cas,
          n_controls = num_inds$n_con
        )
    }
  }
  
  return(model_info_df)
}