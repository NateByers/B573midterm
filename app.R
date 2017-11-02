library(shiny)
library(dplyr)

source("helpers.R")
load("shiny_data.rda")

symbols <- unique(c(interactions$OFFICIAL_SYMBOL_A, 
                    interactions$OFFICIAL_SYMBOL_B))
symbols <- symbols[order(symbols)]

server <- function(input, output) {
  
  output$interactions_table <- renderDataTable(interactions)
  
  get_user_selected_functions <- reactive({
    get_protein_function(interactions, lookup, input$protein) %>%
      dplyr::select(-official_symbol) %>%
      dplyr::mutate(protein_function = paste0('<a href="http://amigo.geneontology.org/amigo/term/',
                                              protein_function, '"  target="_blank">', 
                                              protein_function, '</a>'))
  })
  
  output$predicted_functions <- renderDataTable({
    get_user_selected_functions()
  }, escape = FALSE)
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("protein", "Select Protein:", symbols),
      textOutput("protein_function")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Predicted Functions", dataTableOutput("predicted_functions")),
        tabPanel("Interactions", dataTableOutput("interactions_table"))
        )
      )
    )
)

onStop(function() unlink("shiny_data.rda"))

shinyApp(ui = ui, server = server)