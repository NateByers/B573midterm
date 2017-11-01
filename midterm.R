cat("[1] loading libraries (and installing if necessary)...")
source("libraries.R")
cat("\n[2] sourcing helper functions...")
source("helpers.R")

cat("\n[3] reading in interactions text file...")
interactions <- get_interactions() 

cat("\n[4] reading in ontologies text file...")
ontologies <- get_ontologies()

cat("\n[5] making a lookup table of proteins and functions...")
lookup <- make_lookup(interactions, ontologies)

cat("\n[6] making a useful table of interactions...")
proteins <- make_table(interactions, ontologies, lookup)

cat("\n[7] processing the table to get the associated function for each protein
     based on interactions...")
proteins <- process_table(proteins)
