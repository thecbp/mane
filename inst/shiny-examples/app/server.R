server <- function(input, output, session) {

  run_simulation = reactive({

    if (input$simulate > 0) {

      treatment_effect_input_names = paste0("trt-effect-", 1:input$n_trts)
      betas = c()
      for (effect in treatment_effect_input_names) {
        betas = c(betas, input[[effect]])
      }
      intercept_prior_string = paste0("normal(", input$intercept_prior_mean,
                                      ",", input$intercept_prior_variance, ")")
      treatment_prior_string = paste0("normal(", input$treatment_prior_mean,
                                      ",", input$treatment_prior_variance, ")")

      priors = list(
        "Intercept" = intercept_prior_string,
        "b" = treatment_prior_string
      )

      sims = list()

      withProgress(message = "Simulating Platform-of-1 trials...", {

        for (i in seq_len(input$n_sims)) {

          sim = simulate(n_trts = input$n_trts,
                         n_burn_cycles = input$n_burn_cycles,
                         burn_obvs_per_period = input$burn_n_obvs_per_period,
                         adaptive_obvs_per_period = input$adaptive_n_obvs_per_period,
                         max_duration = input$maximum_duration,
                         betas = betas,
                         y_sigma = input$within_person_noise,
                         priors = priors,
                         n_chains = input$n_chains,
                         n_iter = input$samples_per_chain,
                         lag = 1, # RECHECK LATER
                         stabilize = NULL,
                         objective = input$objective,
                         adapt_delta = input$adapt_delta,
                         max_treedepth = input$max_treedepth)

          sims[[i]] = sim

          # Increment progress on the progress bar
          incProgress(1 / input$n_sims)
        }

      })

      return(sims)
    }
  })

  # Create UI to specify the true treatment effects
  output$modelControls = renderUI({

    model_parameter_panel = purrr::map(1:input$n_trts, function(trt) {

      fluidRow(column(12, numericInput(inputId = paste0("trt-effect-", trt),
                               label = paste0("Treatment ", trt, " Effect"),
                               value = 0)))

    })

    model_parameter_panel

  })

  output$test = renderPrint({
    run_simulation()
  })

}
