get_interactions <- function() {
  
  interactions <- read.table("BIOGRID-ORGANISM-Homo_sapiens-3.1.91.tab",
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
  
  ontologies <- readr::read_delim("gene_association.goa_human", "\t", quote = "",
                                  skip = 23) 
  
  # # get the headers from the README here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README
  # headers <- readLines("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README")
  # headers <-headers[(grep("^GAF2.1", headers) + 3):(grep("^GPAD1.1", headers) - 2)]
  # headers <- sub("\\t\\d{1,2}\\s{1,7}", "", headers)
  # headers <- gsub("\\(|\\)|:", "_", headers)
  # headers <- gsub(" ", "", headers)
  
  headers <- c("DB", "DB_Object_ID", "DB_Object_Symbol", "Qualifier", "GO_ID",
               "DB_Reference", "EvidenceCode", "With_or_From", "Aspect", "DB_Object_Name",
               "DB_Object_Synonym", "DB_Object_Type", "TaxonandInteractingtaxon",
               "Date", "Assigned_By", "Annotation_Extension", "Gene_Product_Form_ID")
  
  names(ontologies) <- headers
  
  ontologies[["DB_Object_Synonym"]] <- tolower(ontologies[["DB_Object_Synonym"]])

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
    dplyr::rename(aliases = DB_Object_Synonym, go_id = GO_ID) %>%
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


determine_protein_function <- function(interactions = get_interactions(), 
                                 lookup = make_lookup(), protein) {
  
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
    dplyr::group_by(go_id) %>%
    dplyr::summarize(number_of_interactions = n()) %>%
    dplyr::filter(number_of_interactions == max(number_of_interactions)) %>%
    dplyr::filter(number_of_interactions > 1)
    
  table
}


get_amigo_info <- function(go_id) {
  
  lines <- try(readLines(paste0("http://amigo.geneontology.org/amigo/term/", go_id),
                         n = 1000))
  
  if(class(lines) == "try-error") {
    return(data.frame(amigo_name = "[can't reach Amigo website]", 
                      amigo_definition = "[can't reach Amigo website]",
                      stringsAsFactors = FALSE))
  }
  
  name_start <- grep('<dt>Name</dt>', lines)[1]
  
  name <- trimws(gsub("</?dd>", "", lines[name_start + 1]), "both")
  
  if(name == "") {
    name <- "None"
  }
  
  definition_start <- grep('<dt>Definition</dt>', lines)[1]
  
  definition_stop <- grep('<dt>Comment</dt>', lines)
  definition_stop <- definition_stop[definition_stop > definition_start][1]
  
  definition <- lines[(definition_start + 1):(definition_stop - 1)]
  definition <- definition[!grepl("Source|cite|Comment", definition)]
  definition <- gsub("</?dd>|\\t", "", definition)
  definition <- trimws(grep("\\S", definition, value = TRUE), "both")
  
  data.frame(amigo_name = name, amigo_definition = definition,
             stringsAsFactors = FALSE)
}


attach_amigo_info <- function(protein_functions) {
  
  if(nrow(protein_functions) == 0) {
    return(protein_functions)
  }
  
  amigo_info <- lapply(protein_functions$go_id, get_amigo_info)
  amigo_info <- Reduce(rbind, amigo_info)
  
  cbind(protein_functions, amigo_info)
}


return_function_text <- function(interactions, lookup, protein) {
  
  cat(paste0("Protein Function(s) for ", protein,":\n"))
  
  if(tolower(protein) %in% lookup$official_symbol) {
    protein_df <- determine_protein_function(interactions, lookup, protein) %>%
      attach_amigo_info() %>%
      as.data.frame()
    
    if(nrow(protein_df) > 0) {
      return(protein_df)
    } else {
      cat("Not enough information to assign a function")
    }
  } else {
    cat("Can't match the protein as an Official Symbol")
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

