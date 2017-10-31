library(dplyr)
library(parallel)
source("read_data.R")
source("make_lookup.R")
source("make_table.R")
source("get_protein_function.R")

process_table <- function(table = make_table(), parallel = FALSE, clusters) {

  table <- table %>%
    dplyr::group_by(interaction)
  
  if(parallel) {
    if(missing(clusters)) {
      clusters <- parallel::detectCores() - 1
    }
    cl <- parallel::makeCluster(clusters)
    clusterEvalQ(cl, {
      library(dplyr)
      NULL
    })
    protein_dfs <- parallel::parLapply(cl, unique(table$official_symbol), function(protein, table) {
      # protein <- unique(table$official_symbol)[1]
      table %>%
        dplyr::filter(protein %in% official_symbol) %>%
        dplyr::filter(official_symbol != protein) %>%
        dplyr::group_by(protein_function) %>%
        dplyr::summarize(n = n()) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::top_n(1, protein_function)
    }, table = table)
    stopCluster(cl)
  } else {
    then <- Sys.time()
    protein_dfs <- lapply(unique(table$official_symbol)[1:10], function(protein) {
      # protein <- unique(table$official_symbol)[1]
      table %>%
        dplyr::filter(protein %in% official_symbol) %>%
        dplyr::filter(official_symbol != protein) %>%
        dplyr::group_by(protein_function) %>%
        dplyr::summarize(n = n()) %>%
        dplyr::filter(n == max(n)) %>%
        dplyr::top_n(1, protein_function)
    })
    Sys.time() - then
  }
    
  Reduce(rbind, protein_dfs)
  
  
}


# process_table <- function(table = make_table(), lookup = make_lookup()) {
#   
#   table <- table %>%
#     dplyr::group_by(interaction)
# }