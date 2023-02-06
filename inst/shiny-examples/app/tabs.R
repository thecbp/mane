simulationTab = sidebarLayout(
  sidebarPanel(
    # fluidRow(
    #   column(12, submitButton(text = h2("Simulate trials"),
    #                           width = "100%"))
    # ),
    wellPanel(
      h3("Trial Simulation Parameters"),
      numericInput(inputId = "n-trts",
                   label = "Number of Treatments",
                   value = 2),
      numericInput(inputId = "burn-in-lengths",
                   label = "Burn-in Period Lengths",
                   value = 1),
      numericInput(inputId = "adaptive-length",
                   label = "Adaptive Period Lengths",
                   value = 1),
      numericInput(inputId = "maximum-length",
                   label = "Maximum Number of Periods",
                   value = 50),
      numericInput(inputId = "n-sims",
                   label = "Number of Simulations",
                   value = 100),
      selectInput(inputId = "objective",
                  label = "Maximize or minimize reward?",
                  choices = c("Maximize", "Minimize"))
    ),
    wellPanel(
      h3("Model Parameters"),
      uiOutput("modelControls"),
      numericInput(inputId = "within-person-noise",
                   label = "Within-person noise",
                   value = 10)
    ),
    wellPanel(
      h3("Prior Parameters"),
      uiOutput("priorControls")
    ),
    wellPanel(
      h3("MCMC Parameters"),
      numericInput(inputId = "n-chains",
                   label = "Number of chains",
                   value = 4),
      numericInput(inputId = "samples-per-chain",
                   label = "Samples per chains",
                   value = 1000)
    )
  ),
  mainPanel(
    plotOutput("hist"),
    textOutput("test")
  )
)
