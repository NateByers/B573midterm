
get_interactions <- function() {
  interactions <- read.table("data/BIOGRID-ORGANISM-Homo_sapiens-3.1.91.tab.txt",
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
  #ontologies <- read.table("data/gene_association.goa_human.txt",
  #                         stringsAsFactors = FALSE, skip = 23, header = FALSE, sep = "\t",
  #                         comment.char = "", quote = ""
  #                         )
  
  ontologies <- readr::read_delim("data/gene_association.goa_human.txt", "\t", quote = "",
                                  skip = 23)
  
  for(i in names(ontologies)) {
    if(class(ontologies[[i]]) == "character") {
      ontologies[[i]] <- tolower(ontologies[[i]])
    }
  }
  
  # README is here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README
  
  headers <- readLines("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README")
  headers <-headers[(grep("^GAF2.1", headers) + 3):(grep("^GPAD1.1", headers) - 2)]
  headers <- sub("\\t\\d{1,2}\\s{1,7}", "", headers)
  headers <- gsub("\\(|\\)|:", "_", headers)
  headers <- gsub(" ", "", headers)
  
  names(ontologies) <- headers

  ontologies
}

