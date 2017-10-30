make_table <- function(interactions = get_interactions(), ontologies = get_ontologies()) {
  # interactions = get_interactions(); ontologies = get_ontologies()
  lookup <- make_lookup(interactions, ontologies)
  
  interactions <- interactions %>%
    dplyr::mutate(interaction = row_number()) %>%
    dplyr::select(interaction, starts_with("INTERACTOR")) %>%
    tidyr::gather("interactor_side", "interactor_id", starts_with("INTERACTOR")) %>%
    dplyr::left_join(lookup_interactions, "interactor_id") 
  
  
    dplyr::left_join(lookup_ontologies, "aliases")
}
