source("helpers.R")
load("shiny_data.rda")

server <- function(input, output) {
  
  output$protein_table <- renderDataTable({
    proteins
  })
  
  output$protein_function <- renderText({
    return_function_text(proteins, input$protein)
  })
  
  
}

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      selectInput("protein", "Select Protein:",
                  unique(proteins$official_symbol)[order(unique(proteins$official_symbol))]),
      textOutput("protein_function")
    ),
    mainPanel(dataTableOutput("protein_table"))
  )
)

onStop(function() rm("shiny_data.rda"))

shinyApp(ui = ui, server = server)