make_table <- function(interactions = get_interactions(), ontologies = get_ontologies()) {
  # interactions = get_interactions(); ontologies = get_ontologies()
  lookup_interactions <- make_interaction_lookup(interactions)
  lookup_ontologies <- make_ontologies_lookup(ontologies)
  
  interactions <- interactions %>%
    dplyr::mutate(interaction = row_number()) %>%
    dplyr::select(interaction_id, starts_with("INTERACTOR")) %>%
    tidyr::gather("interactor_id", "alias", starts_with("INTERACTOR")) %>%
    dplyr::left_join(lookup_interactions, "interactor_id")
}
