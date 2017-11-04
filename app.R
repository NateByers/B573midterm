library(shiny)
library(dplyr)

source("helpers.R")

if(!file.exists("shiny_data.rda")) {
  stop("must run midterm.R to open the shiny app")
}

load("shiny_data.rda")

symbols <- unique(c(interactions$OFFICIAL_SYMBOL_A, 
                    interactions$OFFICIAL_SYMBOL_B))
symbols <- symbols[order(symbols)]

server <- function(input, output) {
  
  get_user_selected_functions <- reactive({
    determine_protein_function(interactions, lookup, input$protein) %>%
      attach_amigo_info() %>%
      dplyr::mutate(go_id = paste0('<a href="http://amigo.geneontology.org/amigo/term/',
                                              go_id, '"  target="_blank">', 
                                              go_id, '</a>'))
  })
  
  output$predicted_functions <- renderDataTable({
    get_user_selected_functions()
  }, escape = FALSE)
  
  get_filtered_interactions <- reactive({
    interactions %>%
      dplyr::filter(OFFICIAL_SYMBOL_A == input$protein | OFFICIAL_SYMBOL_B == input$protein)
  })
  
  output$interactions_table <- renderDataTable({
    get_filtered_interactions()
  })
}

ui <- fluidPage(title = "Protein Function by Association",
  sidebarLayout(
    sidebarPanel(
      selectInput("protein", "Select Protein:", symbols),
      helpText("Choose from the list of proteins above.",
               "The table on the right will show the",
               "GO ids that are most associated with the",
               "protein from the interactions table.",
               "You can see the interactions table by",
               "clicking on the 'Interactions' tab on ",
               "the right."),
      width = 2
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Predicted Functions", dataTableOutput("predicted_functions")),
        tabPanel("Interactions", dataTableOutput("interactions_table"))
        )
      )
    )
)


shinyApp(ui = ui, server = server)