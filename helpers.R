get_interactions <- function() {
  interactions <- read.table("BIOGRID-ORGANISM-Homo_sapiens-3.1.91.tab.txt",
                             stringsAsFactors = FALSE, skip = 35, header = TRUE, sep = "\t",
                             comment.char = "", quote = "") 
  
  for(i in names(interactions)) {
    if(class(interactions[[i]]) == "character") {
      interactions[[i]] <- tolower(interactions[[i]])
    }
  }
  
  interactions
}


get_ontologies <- function() {
  
  ontologies <- readr::read_delim("gene_association.goa_human.txt", "\t", quote = "",
                                  skip = 23) 
  
  for(i in names(ontologies)) {
    if(class(ontologies[[i]]) == "character") {
      ontologies[[i]] <- tolower(ontologies[[i]])
    }
  }
  
  # get the headers from the README is here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README
  
  headers <- readLines("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README")
  headers <-headers[(grep("^GAF2.1", headers) + 3):(grep("^GPAD1.1", headers) - 2)]
  headers <- sub("\\t\\d{1,2}\\s{1,7}", "", headers)
  headers <- gsub("\\(|\\)|:", "_", headers)
  headers <- gsub(" ", "", headers)
  
  names(ontologies) <- headers
  
  ontologies
}

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
    dplyr::select(DB_Object_Synonym, GO_ID) %>%
    dplyr::rename(aliases = DB_Object_Synonym, protein_function = GO_ID) %>%
    tidyr::separate_rows(aliases, sep = "\\|") %>%
    dplyr::mutate(aliases = trimws(aliases, "both")) %>%
    dplyr::distinct() 
}

make_lookup <- function(interactions = get_interactions(), ontologies = get_ontologies()) {
  lookup_interactions <- make_interaction_lookup(interactions)
  lookup_ontologies <- make_ontologies_lookup(ontologies)
  
  lookup <- inner_join(lookup_interactions, lookup_ontologies, "aliases") %>%
    dplyr::select(-aliases) %>%
    distinct() 
  
  lookup
}


get_protein_function <- function(interactions = get_interactions(), 
                                 lookup = make_lookup(), protein) {
  # protein <- "EXOSC4"
  protein <- tolower(protein)
  
  interactions <- interactions %>%
    dplyr::mutate(interaction = row_number()) %>%
    dplyr::select(interaction, starts_with("OFFICIAL")) %>%
    tidyr::gather("interactor", "official_symbol", starts_with("OFFICIAL"))
  
  protein_df <- interactions %>%
    dplyr::filter(official_symbol == protein)
  
  table <- interactions %>%
    dplyr::semi_join(protein_df, "interaction") %>%
    dplyr::filter(official_symbol != protein) %>%
    dplyr::inner_join(lookup, "official_symbol") %>%
    dplyr::group_by(protein_function) %>%
    dplyr::summarize(number_of_interactions = n()) %>%
    dplyr::filter(number_of_interactions > 1) %>%
    dplyr::arrange(desc(number_of_interactions)) %>%
    dplyr::mutate(official_symbol = toupper(protein),
                  protein_function = toupper(protein_function))
    
  table
}


return_function_text <- function(interactions, lookup, protein) {
  
  protein_df <- get_protein_function(interactions, lookup, protein) %>%
    dplyr::select(-official_symbol) %>%
    as.data.frame()
  
  
  cat(paste0("Protein Function(s) for ", protein,":\n"))
  
  if(nrow(protein_df)) {
     return(protein_df)
  } else {
    cat("Not enough information to assign a function")
  }
  
}

read_protein <- function(interactions, lookup) {
  
  protein <- readline(prompt = "Enter protein: ") 
  
  return_function_text(interactions, lookup, protein)
}

start_shiny <- function() {
  
  start <- readline(prompt = "Start shiny app [y/n]?")
  
  if(tolower(start) %in% c("y", "yes")) {
    save(interactions, lookup, file = "shiny_data.rda")
    if(!require(shiny)) {
      install.packages("shiny")
    }
    shiny::runApp()
  }
}




