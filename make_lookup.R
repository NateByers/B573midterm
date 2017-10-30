make_interaction_lookup <- function(interactions = get_interactions()) {
  # the %>% function is from the magrittr package, which is imported with the dplyr package

  lookup <- lapply(c("A", "B"), function(x, interactions) {
    # x = "A"
    
    interactions %>%
      dplyr::select(ends_with(x)) %>%
      dplyr::rename_all(funs(sub(paste0("_", x), "", .))) %>%
      dplyr::mutate(aliases = paste(OFFICIAL_SYMBOL, ALIASES_FOR, sep = "|")) %>%
      dplyr::select(INTERACTOR, aliases) %>%
      dplyr::rename(interactor_id = INTERACTOR) %>%
      tidyr::separate_rows(aliases, sep = "\\|") %>%
      dplyr::mutate(aliases = trimws(aliases, "both"),
                    aliases = tolower(aliases)) %>%
      dplyr::distinct() 
    
  }, interactions = interactions)
  
  lookup <- Reduce(rbind, lookup) %>%
    dplyr::distinct()
  
  lookup
}

make_ontologies_lookup <- function(ontologies = get_ontologies()) {
  
  lookup <- ontologies %>%
    dplyr::select(DB_Object_ID, DB_Object_Synonym, GO_ID) %>%
    dplyr::rename(gene_id = DB_Object_ID, aliases = DB_Object_Synonym, go_id = GO_ID) %>%
    tidyr::separate_rows(aliases, sep = "\\|") %>%
    dplyr::mutate(aliases = trimws(aliases, "both"),
                  aliases = tolower(aliases)) %>%
    dplyr::distinct() 
  
}

