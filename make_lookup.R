make_interaction_lookup <- function(interactions = get_interactions()) {
  # the %>% function is from the magrittr package, which is imported with the dplyr package

  lookup <- lapply(c("A", "B"), function(x, interactions) {
    # x = "A"
    
    interactions %>%
      dplyr::select(ends_with(x)) %>%
      dplyr::rename_all(funs(sub(paste0("_", x), "", .))) %>%
      dplyr::mutate(aliases = paste(OFFICIAL_SYMBOL, ALIASES_FOR, sep = "|")) %>%
      dplyr::select(OFFICIAL_SYMBOL, aliases) %>%
      dplyr::rename(official_symbol = OFFICIAL_SYMBOL) %>%
      tidyr::separate_rows(aliases, sep = "\\|") %>%
      dplyr::mutate(aliases = trimws(aliases, "both")) %>%
      dplyr::distinct() 
    
  }, interactions = interactions)
  
  lookup <- Reduce(rbind, lookup) %>%
    dplyr::distinct()
  
  lookup
}

make_ontologies_lookup <- function(ontologies = get_ontologies()) {
  
  lookup <- ontologies %>%
    dplyr::select(DB_Object_Synonym, DB_Object_Name) %>%
    dplyr::rename(aliases = DB_Object_Synonym, protein_function = DB_Object_Name) %>%
    tidyr::separate_rows(aliases, sep = "\\|") %>%
    dplyr::mutate(aliases = trimws(aliases, "both")) %>%
    dplyr::distinct() 
  
}

make_lookup <- function(interactions = get_interactions(), ontologies = get_ontologies()) {
  lookup_interactions <- make_interaction_lookup(interactions)
  lookup_ontologies <- make_ontologies_lookup(ontologies)
  
  lookup <- inner_join(lookup_interactions, lookup_ontologies, "aliases") %>%
    dplyr::select(-aliases) %>%
    dplyr::mutate(protein_function = gsub("-|,", " ", protein_function),
                  protein_function = trimws(protein_function, "both"),
                  protein_function = gsub("\\s{2,}", " ", protein_function)) %>%
    distinct() %>%
    arrange(protein_function)
  
  lookup
}

