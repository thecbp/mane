simulationTab = sidebarLayout(
  sidebarPanel(
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
                   value = 5),
      numericInput(inputId = "adaptive_n_obvs_per_period",
                   label = "Number of Observations per Period (Adaptive)",
                   value = 5),
      numericInput(inputId = "maximum_duration",
                   label = "Maximum Number of Periods",
                   value = 10),
      numericInput(inputId = "n_sims",
                   label = "Number of Simulations",
                   value = 10),
      selectInput(inputId = "objective",
                  label = "Maximize or minimize reward?",
                  choices = c("Maximize", "Minimize"))
    ),
    wellPanel(
      h3("Model Parameters"),
      uiOutput("modelControls"),
      numericInput(inputId = "within_person_noise",
                   label = "Within-person noise",
                   value = 10),
      sliderInput(inputId = "serial_correlation",
                   label = "Degree of serial correlation",
                   min = 0, max = 0.99, value = 0, step = 0.01)
    ),
    wellPanel(
      h3("Decision Rule Parameters"),
      fluidRow(
        column(6, numericInput(inputId = "delta_eff",
                               label = paste0("Efficacy threshold"),
                               value = 0)),
        column(6, numericInput(inputId = "delta_fut",
                               label = paste0("Futility threshold"),
                               value = 0))
      )
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
                   value = 2000),
      numericInput(inputId = "adapt_delta",
                   label = "Adapt Delta",
                   value = 0.99),
      numericInput(inputId = "max_treedepth",
                   label = "Maximum tree depth",
                   value = 15)
    ),
    fluidRow(
      column(6, actionButton(inputId = "simulate",
                             label = "Simulate",
                             width = "100%")),
      column(6,  downloadButton("download", "Download Sims"))
    )
  ),
  mainPanel(
    conditionalPanel(condition = "input.simulate == 0",
                     wellPanel(
                       h2("No simulations created yet"),
                       p("No simulations have been created yet, so the plots have not been populated. Please run some of your own simulations.")),
                     ),
    textOutput("latex_text_output"),
    plotOutput("epp"),
    plotOutput("power"),
    plotOutput("fwer")
  )
)

howtoTab = fluidPage(
    fluidRow(
      column(1),
      column(10, h2("Step 1: Decide trial parameters")),
      column(1)
    ),
    fluidRow(
      column(2),
      column(8, uiOutput("step1")),
      column(2)
    ),
    fluidRow(
      column(1),
      column(10, h2("Step 2: Set treatment effects and reward model parameters")),
      column(1)
    ),
    fluidRow(
      column(2),
      column(8, uiOutput("step2")),
      column(2)
    ),
    fluidRow(
      column(1),
      column(10, h2("Step 3: Decide futility and decision rules")),
      column(1)
    ),
    fluidRow(
      column(2),
      column(8, uiOutput("step3")),
      column(2)
    ),
    fluidRow(
      column(1),
      column(10, h2("Optional: Set priors and MCMC parameters")),
      column(1)
    ),
    fluidRow(
      column(2),
      column(8, uiOutput("step4")),
      column(2)
    )
)

