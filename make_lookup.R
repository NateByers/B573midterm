make_interaction_lookup <- function(interactions = get_interactions()) {
  # the %>% function is from the magrittr package, which is imported with the dplyr package

  lookup <- lapply(c("A", "B"), function(x, interactions) {
    # x = "A"
    
    interactions %>%
      dplyr::select(ends_with(x)) %>%
      dplyr::rename_all(funs(sub(paste0("_", x), "", .))) %>%
      dplyr::mutate(ALIASES = paste(OFFICIAL_SYMBOL, ALIASES_FOR, sep = "|")) %>%
      dplyr::select(OFFICIAL_SYMBOL, ALIASES) %>%
      tidyr::separate_rows(ALIASES, sep = "\\|") %>%
      dplyr::distinct() 
    
  }, interactions = interactions)
  
  lookup <- Reduce(rbind, lookup) %>%
    dplyr::distinct()
  
  lookup
}

