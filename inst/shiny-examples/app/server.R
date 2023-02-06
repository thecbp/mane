source("helpers")

server <- function(input, output, session) {

  simulation_data = eventReactive(input[["simulate"]], {

    treatment_effect_input_names = paste0("trt-effect-", 1:input[["n-trts"]])
    betas = input[treatment_effect_input_names] %>% unlist()

    withProgress(message = "Simulating Platform-of-1 trials", {
      for (i in seq_len(input[["n-sims"]])) {

        sim = simulate(n_trts = input[["n-trts"]],
                       n_burn_cycles = input[["n-burn-cycles"]],
                       n_adaptive_periods = input[[]],
                       betas = betas,
                       y_sigma = input[["within-person-noise"]],
                       n_chains = input[["n-chains"]],
                       n_iter = input[["samples-per-chain"]],
                       lag = 1, # RECHECK LATER
                       stabilize = NULL,
                       objective = input[["objective"]],
                       adapt_delta = input[["adapt-delta"]],
                       max_treedepth = input[["max-treedepth"]])

        # Increment progress on the progress bar
        incProgress(1 / input$steps)
      }
      # ADD LATER: SOME STRUCTURE THAT PUTS TOGETHER ALL OF THE SIM DATA
    })



  })

  # Create UI to specify the true treatment effects
  output$modelControls = renderUI({

    model_parameter_panel = purrr::map(1:input[["n-trts"]], function(trt) {

      fluidRow(column(12, numericInput(inputId = paste0("trt-effect-", trt),
                               label = paste0("Treatment ", trt, " Effect"),
                               value = 0)))

    })

    model_parameter_panel

  })


  # Create UI to specify the priors for the treatment effects
  output$priorControls = renderUI({

    prior_parameter_panel = purrr::map(1:input[["n-trts"]], function(trt) {

      fluidRow(
        column(6, numericInput(inputId = paste0("trt-", trt, "-mean"),
                               label = paste0("Treatment ", trt, " Prior Mean"),
                               value = 0)),
        column(6, numericInput(inputId = paste0("trt-", trt, "-variance"),
                               label = paste0("Treatment ", trt, " Prior Variance"),
                               value = 100))

      )

    })

    prior_parameter_panel

  })

  output$hist = renderPlot({
    hist(rnorm(input[["burn-in-lengths"]]))
  })

  output$test = renderText({




  })





}
