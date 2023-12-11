library(shiny)

# Define UI for app that draws a histogram
ui <- fluidPage(

  # App title
  titlePanel("fastcpd - Fast Change Point Detection"),

  # Sidebar layout with input and output definitions
  sidebarLayout(

    # Sidebar panel for inputs
    sidebarPanel(fileInput(
      inputId = "upload",
      label = "Select a csv file",
      buttonLabel = "Select",
      multiple = FALSE,
      accept = c(".csv")
    )),

    # Main panel for displaying outputs
    mainPanel(textOutput(outputId = "summary"))
  )
)

server <- function(input, output) {
  output$summary <- renderText({
    uploaded_file <- input$upload
    if (is.null(uploaded_file)) {
      return(NULL)
    }
    data <- read.csv(file = uploaded_file$datapath, header = TRUE, sep = ",")
    fastcpd::fastcpd.lm(data)@cp_set
  })
}

shinyApp(ui = ui, server = server)
