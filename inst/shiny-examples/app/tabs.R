simulationTab = sidebarLayout(
  sidebarPanel(
    fluidRow(
      column(12, actionButton(inputId = "simulate",
                              label = h2("Simulate trials"),
                              width = "100%"))
    ),
    wellPanel(
      h3("Trial Simulation Parameters"),
      numericInput(inputId = "n_trts",
                   label = "Number of Treatments",
                   value = 3),
      numericInput(inputId = "n_burn_cycles",
                   label = "Number of Burn-in Cycles",
                   value = 1),
      numericInput(inputId = "burn_n_obvs_per_period",
                   label = "Number of Observations per Period (Burn-in)",
                   value = 3),
      numericInput(inputId = "adaptive_n_obvs_per_period",
                   label = "Number of Observations per Period (Adaptive)",
                   value = 3),
      numericInput(inputId = "maximum_duration",
                   label = "Maximum Number of Periods",
                   value = 5),
      numericInput(inputId = "n_sims",
                   label = "Number of Simulations",
                   value = 1),
      selectInput(inputId = "objective",
                  label = "Maximize or minimize reward?",
                  choices = c("Maximize", "Minimize"))
    ),
    wellPanel(
      h3("Model Parameters"),
      uiOutput("modelControls"),
      numericInput(inputId = "within_person_noise",
                   label = "Within-person noise",
                   value = 10)
    ),
    wellPanel(
      h3("Prior Parameters"),
      fluidRow(
        column(6, numericInput(inputId = "intercept_prior_mean",
                               label = paste0("Intercept Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "intercept_prior_variance",
                               label = paste0("Intercept Prior Variance"),
                               value = 100))

      ),
      fluidRow(
        column(6, numericInput(inputId = "treatment_prior_mean",
                               label = paste0("Treatment Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = "treatment_prior_variance",
                               label = paste0("Treatment Prior Variance"),
                               value = 100))
      ),
      fluidRow(
        column(12, numericInput(inputId = "noise_prior_df",
                               label = paste0("Noise Prior Degrees of Freedom"),
                               value = 1))

      )
    ),
    wellPanel(
      h3("MCMC Parameters"),
      numericInput(inputId = "n_chains",
                   label = "Number of chains",
                   value = 4),
      numericInput(inputId = "samples_per_chain",
                   label = "Samples per chains",
                   value = 1000),
      numericInput(inputId = "adapt_delta",
                   label = "Adapt Delta",
                   value = 0.999),
      numericInput(inputId = "max_treedepth",
                   label = "Maximum tree depth",
                   value = 15)
    )
  ),
  mainPanel(
    plotOutput("hist"),
    verbatimTextOutput("test")
  )
)
