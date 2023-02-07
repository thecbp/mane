source("helpers")

server <- function(input, output, session) {

  simulation_data = eventReactive(input[["simulate"]], {

    treatment_effect_input_names = paste0("trt-effect-", 1:input[["n-trts"]])
    betas = input[treatment_effect_input_names] %>% unlist()
    intercept_prior_string = paste0("normal(", input[["intercept-prior-mean"]],
                                    ",", input[["intercept-prior-variance"]], ")")
    treatment_prior_string = paste0("normal(", input[["treatment-prior-mean"]],
                                    ",", input[["treatment-prior-variance"]], ")")

    priors = list(
      "Intercept" = intercept_prior_string,
      "b" = treatment_prior_string
    )

    withProgress(message = "Simulating Platform-of-1 trials", {

      for (i in seq_len(input[["n-sims"]])) {

        sim = simulate(n_trts = input[["n-trts"]],
                       n_burn_cycles = input[["n-burn-cycles"]],
                       burn_obvs_per_period = input[["burn-n-obvs-per-period"]],
                       adaptive_obvs_per_period = input[["adaptive-n-obvs-per-period"]],
                       max_duration = input[["maximum-duration"]],
                       betas = betas,
                       y_sigma = input[["within-person-noise"]],
                       priors = priors,
                       n_chains = input[["n-chains"]],
                       n_iter = input[["samples-per-chain"]],
                       lag = 1, # RECHECK LATER
                       stabilize = NULL,
                       objective = input[["objective"]],
                       adapt_delta = input[["adapt-delta"]],
                       max_treedepth = input[["max-treedepth"]])

        # Increment progress on the progress bar
        incProgress(1 / input[["n-sims"]])
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

  output$hist = renderPlot({
    hist(rnorm(input[["burn-in-lengths"]]))
  })

  output$test = renderText({




  })





}
