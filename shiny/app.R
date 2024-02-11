source("ui.R")
source("server.R")

shinyApp(
  ui = ui, server = server, options = list(host = "0.0.0.0", port = 7469)
)
