
get_protein_function <- function(table = make_table(), protein) {
  # protein <- "actn2"
  
  protein_df <- table %>%
    dplyr::filter(official_symbol == protein)
  
  table <- table %>%
    semi_join(protein_df, "interaction") %>%
    filter(official_symbol != protein) %>%
    dplyr::group_by(protein_function) %>%
    dplyr::summarize(n = n()) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::top_n(1, protein_function)
  
  table
}

