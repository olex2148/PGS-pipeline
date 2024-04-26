
create_foelgefil <- function(accession_id){
  foelgefil_df <- data.frame(
    ID = NA,
    Accession_ID = NA,
    Reported_Trait = NA,
    PubMed_ID = NA,
    First_Author = NA,
    Journal = NA,
    Title = NA,
    Publication_Date = NA,
    N = NA,
    N_Cases = NA,
    N_Controls = NA,
    M_Input = NA,
    M_HapMap = NA,
    M_QC = NA          
  )
  
  if(!is.na(accession_id)) {
    study_info <- get_studies(study_id = accession_id)
    num_inds <- get_n_cas_con(study_info@studies$initial_sample_size)
    # Saving variables
    foelgefil_df <- foelgefil_df %>% 
      mutate(
        Accession_ID = accession_id,
        Reported_Trait = study_info@studies$reported_trait,
        PubMed_ID = study_info@publications$pubmed_id,
        First_Author = study_info@publications$author_fullname,
        Journal = study_info@publications$publication,
        Title = study_info@publications$title,
        Publication_Date = study_info@publications$publication_date
      )
    
    if(!is.na(num_inds$n)){
      foelgefil_df <- foelgefil_df %>% 
        mutate(
          N = num_inds$n,
          N_Cases = num_inds$n_cas,
          N_Controls = num_inds$n_con
        )
    } else {
      foelgefil_df <- foelgefil_df %>% 
        mutate(
          N = num_inds$n_cas + num_inds$n_con,
          N_Cases = num_inds$n_cas,
          N_Controls = num_inds$n_con
        )
    }
  }
  
  return(foelgefil_df)
}