simulationTab = sidebarLayout(
  sidebarPanel(
    # fluidRow(
    #   column(12, actionButton(text = h2("Simulate trials"),
    #                           width = "100%"))
    # ),
    wellPanel(
      h3("Trial Simulation Parameters"),
      numericInput(inputId = "n-trts",
                   label = "Number of Treatments",
                   value = 2),
      numericInput(inputId = "n-burn-cycles",
                   label = "Number of Burn-in Cycles",
                   value = 1),
      numericInput(inputId = "burn-n-obvs-per-period",
                   label = "Number of Observations per Period (Burn-in)",
                   value = 1),
      numericInput(inputId = "adaptive-n-obvs-per-period",
                   label = "Number of Observations per Period (Adaptive)",
                   value = 1),
      numericInput(inputId = "maximum-duration",
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
      fluidRow(
        column(6, numericInput(inputId = "intercept-prior-mean",
                               label = paste0("Intercept Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "intercept-prior-variance",
                               label = paste0("Intercept Prior Variance"),
                               value = 100))

      ),
      fluidRow(
        column(6, numericInput(inputId = "treatment-prior-mean",
                               label = paste0("Treatment Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "treatment-prior-variance",
                               label = paste0("Treatment Prior Variance"),
                               value = 100))
      ),
      fluidRow(
        column(12, numericInput(inputId = "noise-prior-df",
                               label = paste0("Noise Prior Degrees of Freedom"),
                               value = 1))

      )
    ),
    wellPanel(
      h3("MCMC Parameters"),
      numericInput(inputId = "n-chains",
                   label = "Number of chains",
                   value = 4),
      numericInput(inputId = "samples-per-chain",
                   label = "Samples per chains",
                   value = 1000),
      numericInput(inputId = "adapt-delta",
                   label = "Adapt Delta",
                   value = 1000),
      numericInput(inputId = "max-treedepth",
                   label = "Maximum tree depth",
                   value = 15)
    )
  ),
  mainPanel(
    plotOutput("hist"),
    textOutput("test")
  )
)
