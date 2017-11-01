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

make_table <- function(interactions = get_interactions(), ontologies = get_ontologies(),
                       lookup = make_lookup()) {
  
  interactions <- interactions %>%
    dplyr::mutate(interaction = row_number()) %>%
    dplyr::select(interaction, starts_with("OFFICIAL")) %>%
    tidyr::gather("interactor", "official_symbol", starts_with("OFFICIAL")) %>%
    dplyr::left_join(lookup, "official_symbol") %>%
    dplyr::arrange(interaction, interactor)
  
  interactions
}

get_protein_function <- function(table = make_table(), protein) {
  
  protein_df <- table %>%
    dplyr::filter(official_symbol == protein)
  
  table <- table %>%
    semi_join(protein_df, "interaction") %>%
    filter(official_symbol != protein) %>%
    dplyr::group_by(protein_function) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::top_n(1, protein_function) %>%
    dplyr::mutate(official_symbol = protein)
  
  table
}


process_table <- function(table = make_table()) {
  
  proteins <- unique(table$official_symbol)
  
  pb <- txtProgressBar(min = 1, max = length(proteins),
                       style = 3)
  
  protein_dfs <- lapply(seq_along(proteins), function(i) {
    setTxtProgressBar(pb, i)
    get_protein_function(table, proteins[i])
    })
  
  close(pb)
  
  Reduce(rbind, protein_dfs)
}

return_function_text <- function(proteins, protein) {
  
  protein_df <- get_protein_function(proteins, protein)
  prot_func <- protein_df[["protein_function"]]
  n <- protein_df[["n"]]
  
  protein_text <- paste0("Protein Function for ", protein,
                          ":")
  
  if(!is.na(prot_func) & n > 1) {
     function_text <- prot_func
  } else {
    function_text <- "not enough information to assign a function"
  }
  
  paste(protein_text, function_text)
}

read_protein <- function(proteins) {
  
  protein <- readline(prompt = "Enter protein: ") %>%
    tolower()
  
  return_function_text(proteins, protein)
}

start_shiny <- function() {
  
  start <- readline(prompt = "Start shiny app [y/n]?")
  
  if(tolower(start) %in% c("y", "yes")) {
    save(proteins, file = "shiny_data.rda")
    if(!require(shiny)) {
      install.packages("shiny")
    }
    shiny::runApp()
  }
}




