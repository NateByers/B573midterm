library(dplyr)
library(parallel)
source("read_data.R")
source("make_lookup.R")
source("make_table.R")
source("get_protein_function.R")

process_table <- function(table = make_table()) {
  
  then <- Sys.time()
  protein_dfs <- lapply(unique(table$official_symbol)[1:10], get_protein_function,
                        table = table)
  Sys.time() - then
  
  protein_df <- Reduce(rbind, protein_dfs)
  
}
