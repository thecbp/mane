server <- function(input, output, session) {

  ##############

  # DYNAMIC UI #

  ##############

  # Create UI to specify the true treatment effects
  output$modelControls = renderUI({

    model_parameter_panel = purrr::map(1:input$n_trts, function(trt) {

      fluidRow(column(12, numericInput(inputId = paste0("trt-effect-", trt),
                                       label = paste0("Treatment ", trt, " Effect"),
                                       value = 0)))

    })

    model_parameter_panel

  })

  ############################

  # SIMULATION FUNCTIONALITY #

  ############################

  simulations = reactive({

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

  load_data = observeEvent("load", {
    print("Load button pressed")
  })

  save_data = observeEvent("save", {

    sims = simulations()
    path = paste0(lubridate::now(),"-platform-of-1-simulations.rds")

    saveRDS(sims, path)

    # TO-DO:
    # add some feedback to tell user that the data was saved
    # tell the user what folder the file is saved in
    # tell the user what the name of the field


  })

  output$allocation_probability_plot = renderPlot({

    sims = simulations()

    if (!is.null(sims)) {
      p = sims[[1]]$allocation_probs %>%
        ggplot(aes(x = period, y = X1)) +
        geom_point() +
        geom_line()
    } else {
      p = NULL
    }

    p

  })

  output$allocation_probability_plot = renderPlot({

    sims = simulations()

    if (!is.null(sims)) {
      p = sims[[1]]$allocation_probs %>%
        ggplot(aes(x = period, y = X1)) +
        geom_point() +
        geom_line()
    } else {
      p = NULL
    }

    p

  })

  output$epp_plot = renderPlot({

    sims = simulations()

    if (!is.null(sims)) {
      p = sims[[1]]$allocation_probs %>%
        ggplot(aes(x = period, y = X1)) +
        geom_point() +
        geom_line()
    } else {
      p = NULL
    }

    p

  })

  output$test = renderPrint({
    sims = simulations()
    sims
  })

}
