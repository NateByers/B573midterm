

make_table <- function(interactions = get_interactions(), ontologies = get_ontologies()) {
  # interactions = get_interactions(); ontologies = get_ontologies()
  lookup <- make_lookup(interactions, ontologies)
  
  interactions <- interactions %>%
    dplyr::mutate(interaction = row_number()) %>%
    dplyr::select(interaction, starts_with("OFFICIAL")) %>%
    tidyr::gather("interactor", "official_symbol", starts_with("OFFICIAL")) %>%
    dplyr::left_join(lookup, "official_symbol") 
  
  interactions
}
