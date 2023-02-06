source("tabs.R")

ui = fluidPage(
  navbarPage("Platform-of-1 Planner",
             tabPanel("Simulation", simulationTab),
             tabPanel("Power"),
             tabPanel("Help")
  )
)