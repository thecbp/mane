server <- function(input, output, session) {

  ##############

  # DYNAMIC UI #

  ##############

  # Create UI to specify the true treatment effects
  output$modelControls = renderUI({

    model_parameter_panel = purrr::map(1:input$n_trts, function(trt) {

      if (trt == 1) {
        fluidRow(column(12, numericInput(inputId = paste0("trt-effect-", trt),
                                         label = paste0("Treatment ", trt, " Effect (Intercept)"),
                                         value = 100)))
      } else {
        fluidRow(column(12, numericInput(inputId = paste0("trt-effect-", trt),
                                         label = paste0("Treatment ", trt, " Effect"),
                                         value = 0)))
      }

    })

    model_parameter_panel

  })

  ############################

  # SIMULATION FUNCTIONALITY #

  ############################

  simulations = reactive({

    gammas_E = seq(0.50, 0.99, 0.01) # efficacy gammas
    gammas_F = seq(0.01, 0.50, 0.01) # futility gammas

    if (input$simulate > 0) {

      treatment_effect_input_names = paste0("trt-effect-", 1:input$n_trts)
      null_betas = c(input[["trt-effect-1"]], rep(0, input$n_trts - 1))
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

      sims = decision_storage = vector('list', length = input$n_sims)
      decision_storage = vector('list', length = input$n_sims)
      decision_storage = lapply(decision_storage, function(l) {
        vector('list', length = length(gammas_E))
      })

      # Simulate trials with the specified true treatment effect
      withProgress(message = "Simulating Platform-of-1 trials...", {

        for (i in seq_len(input$n_sims)) {

          sim = mane::simulate(n_trts = input$n_trts,
                               n_burn_cycles = input$n_burn_cycles,
                               burn_obvs_per_period = input$burn_n_obvs_per_period,
                               adaptive_obvs_per_period = input$adaptive_n_obvs_per_period,
                               max_duration = input$maximum_duration,
                               betas = betas,
                               y_sigma = input$within_person_noise,
                               priors = priors,
                               phi = input$serial_correlation,
                               n_chains = input$n_chains,
                               n_iter = input$samples_per_chain,
                               stabilize = NULL,
                               objective = input$objective,
                               adapt_delta = input$adapt_delta,
                               max_treedepth = input$max_treedepth)

          sims[[i]] = sim

          # Save the posterior samples throughout trial for all periods
          posteriors = generate_trial_posteriors(sim, i)

          # Generate decisions from posteriors based on sample size in gamma
          for (j in 1:length(gammas_E)) {

            # Percentage of posteriors that meet decision rule at given sample size + gamma
            # structure: list where values are period, values are decisions (T or F)
            decisions = posteriors %>%
              map(function(data) {
                if (is.null(data)) { NULL }
                else {
                  X2_sample = data %>% pull(b_X2)
                  X3_sample = data %>% pull(b_X3)

                  list(
                    futility_X2 = mean(X2_sample < input[["delta_fut"]]) < gammas_F[j],
                    futility_X3 = mean(X3_sample < input[["delta_fut"]]) < gammas_F[j],
                    efficacy_X2 = mean(X2_sample < input[["delta_eff"]]) > gammas_E[j],
                    efficacy_X3 = mean(X3_sample < input[["delta_eff"]]) > gammas_E[j]
                  )

                }

              })

            decision_storage[[i]][[j]] = decisions

          }

          # Increment progress on the progress bar
          incProgress(1 / input$n_sims)
        }

      })

      # Simulate trials with all null treatment effects
      null_sims = vector('list', length = input$n_sims)
      null_decision_storage = vector('list', length = input$n_sims)
      null_decision_storage = lapply(null_decision_storage, function(l) {
        vector('list', length = length(gammas_E))
      })

      withProgress(message = "Simulating null trials...", {

        for (i in seq_len(input$n_sims)) {

          sim = mane::simulate(n_trts = input$n_trts,
                               n_burn_cycles = input$n_burn_cycles,
                               burn_obvs_per_period = input$burn_n_obvs_per_period,
                               adaptive_obvs_per_period = input$adaptive_n_obvs_per_period,
                               max_duration = input$maximum_duration,
                               betas = null_betas,
                               y_sigma = input$within_person_noise,
                               priors = priors,
                               phi = input$serial_correlation,
                               n_chains = input$n_chains,
                               n_iter = input$samples_per_chain,
                               stabilize = NULL,
                               objective = input$objective,
                               adapt_delta = input$adapt_delta,
                               max_treedepth = input$max_treedepth)

          null_sims[[i]] = sim

          # Save the posterior samples throughout trial for all periods
          posteriors = generate_trial_posteriors(sim, i)

          # Generate decisions from posteriors based on sample size in gamma
          for (j in 1:length(gammas_E)) {

            # Percentage of posteriors that meet decision rule at given sample size + gamma
            # structure: list where values are period, values are decisions (T or F)
            decisions = posteriors %>%
              map(function(data) {
                if (is.null(data)) { NULL }
                else {
                  X2_sample = data %>% pull(b_X2)
                  X3_sample = data %>% pull(b_X3)

                  list(
                    futility_X2 = mean(X2_sample < input[["delta_fut"]]) < gammas_F[j],
                    futility_X3 = mean(X3_sample < input[["delta_fut"]]) < gammas_F[j],
                    efficacy_X2 = mean(X2_sample < input[["delta_eff"]]) > gammas_E[j],
                    efficacy_X3 = mean(X3_sample < input[["delta_eff"]]) > gammas_E[j]
                  )

                }

              })

            null_decision_storage[[i]][[j]] = decisions

          }

          # Increment progress on the progress bar
          incProgress(1 / input$n_sims)
        }

      })

      out = list()
      out[["sims"]] = sims
      out[["decision_storage"]] = decision_storage
      out[["null_sims"]] = null_sims
      out[["null_decision_storage"]] = null_decision_storage

      saveRDS(out, paste0(lubridate::now(), "-sims.rds"))

      return(out)

    }
  })

  optimal_arm = reactive({

    treatment_effect_input_names = paste0("trt-effect-", 1:input$n_trts)
    effects = c()

    # Add mean changes to intercept to get total effect for non-reference arms
    for (i in 1:input$n_trts) {
      if (i == 1)  {
        effects = c(effects, input[[treatment_effect_input_names[i]]])
      } else {
        effects = c(effects, input[[treatment_effect_input_names[1]]] + input[[treatment_effect_input_names[i]]])
      }
    }

    if (input$objective == "Maximize") {
      optimal_index = which(effects == max(effects))
    } else {
      optimal_index = which(effects == min(effects))
    }

    # Give back the name of the optimal arm
    if (length(optimal_index) > 1) {
      paste0("X1")
    } else {
      paste0("X", optimal_index)
    }

  })

  output$epp = renderPlot({

    sims = simulations()
    effect_sims = sims[["sims"]]

    # Placeholder if the simulations have not been made yet
    if (is.null(sims)) {

      p = ggplot() + theme_bw() +
        labs(title = "Expected proportion of periods on each treatment",
             x = "Period",
             y = "EPP") +
        theme(plot.title = element_text(hjust = 0.5))

    } else {

      allocs = lapply(effect_sims, function(s) {s$allocation_probs}) %>%
        bind_rows() %>%
        pivot_longer(cols = starts_with("X"),
                     names_to = "treatment",
                     values_to = "probability"
        ) %>%
        group_by(treatment, period) %>%
        summarize(
          prob = mean(probability)
        )

      p = allocs %>%
        ggplot(aes(x = period, y = prob, color = treatment)) +
        geom_point() +
        geom_line() +
        theme_bw() +
        labs(title = "Expected proportion of periods on each treatment",
             x = "Period",
             y = "EPP") +
        theme(plot.title = element_text(hjust = 0.5))

    }

    p

  })

  output$fwer = renderPlot({

    sims = simulations()
    optimal_arm = optimal_arm()

    # Placeholder if the simulations have not been made yet
    if (is.null(sims)) {

      p = ggplot() +
        theme_bw() +
        labs(title = "FWER Heatmap",
             x = "Period",
             y = "Efficacy cutoff") +
        theme(plot.title = element_text(hjust = 0.5))

    } else {

      p = create_power_plot(sims[["null_decision_storage"]],
                            nsims = input$n_sims,
                            duration = input$maximum_duration,
                            ntrts = input$n_trts,
                            optimal = optimal_arm,
                            fwer = T) +
        labs(title = "FWER Heatmap",
             x = "Period",
             y = "Efficacy cutoff")

    }

    p


  })

  output$power = renderPlot({

    sims = simulations()
    optimal_arm = optimal_arm()

    # Placeholder if the simulations have not been made yet
    if (is.null(sims)) {

      p = ggplot() + theme_bw() +
        labs(title = "Power Heatmap",
             x = "Period",
             y = "Efficacy cutoff") +
        theme(plot.title = element_text(hjust = 0.5))

    } else {

      p = create_power_plot(sims[["decision_storage"]],
                            nsims = input$n_sims,
                            duration = input$maximum_duration,
                            ntrts = input$n_trts,
                            optimal = optimal_arm) +
        labs(title = "Power Heatmap",
             x = "Period",
             y = "Efficacy cutoff")

    }

    p

  })

  output$download= downloadHandler(

    filename = function() {
      paste0(lubridate::now(), "-sims", ".rds")
    },
    content = function(file) {
      saveRDS(simulations(), file)
    }
  )

  howto = reactive({
    readLines("howto.txt")
  })

  # Text for helping guide the user to use the app
  output$step1 = renderUI({ withMathJax(howto()[1]) })
  output$step2 = renderUI({ withMathJax(howto()[2]) })
  output$step3 = renderUI({ withMathJax(howto()[3]) })
  output$step4 = renderUI({ withMathJax(howto()[4]) })

}
