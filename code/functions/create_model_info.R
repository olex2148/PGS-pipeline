
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
  
  
  
  return(model_info_df)
}
