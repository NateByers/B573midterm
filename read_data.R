library(dplyr)

aliases <- read.table("data/BIOGRID-ORGANISM-Homo_sapiens-3.1.91.tab.txt",
                 stringsAsFactors = FALSE, skip = 35, header = TRUE, sep = "\t",
                 comment.char = "")


interactions <- read.table("data/gene_association.goa_human.txt",
                           stringsAsFactors = FALSE, skip = 23, header = FALSE, sep = "\t",
                           comment.char = "")

# README is here ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README

headers <- readLines("ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/README")
headers <-headers[(grep("^GAF2.1", headers) + 3):(grep("^GPAD1.1", headers) - 2)]
headers <- sub("\\t\\d{1,2}\\s{1,7}", "", headers)
headers <- gsub("\\(|\\)|:", "_", headers)
headers <- gsub(" ", "", headers)

names(interactions) <- headers
