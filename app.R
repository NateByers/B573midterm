source("helpers.R")
load("shiny_data.rda")

symbols <- unique(c(interactions$OFFICIAL_SYMBOL_A, 
                    interactions$OFFICIAL_SYMBOL_B))
symbols <- symbols[order(symbols)]



server <- function(input, output) {
  
  output$interactions_table <- renderDataTable(interactions)
  
  output$predicted_functions <- renderDataTable({
    get_protein_function(interactions, lookup, input$protein)
  })
  
  
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

onStop(function() rm("shiny_data.rda"))

shinyApp(ui = ui, server = server)