server <- function(input, output, session) {

  # Create UI to specify the true treatment effects
  output$modelControls = renderUI({

    model_parameter_panel = purrr::map(1:input[["n-trts"]], function(trt) {

      fluidRow(
        column(12, numericInput(inputId = paste0("trt-", trt, "-mean"),
                               label = paste0("Treatment ", trt, " Effect"),
                               value = 0))

      )

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
    input[["burn-in-lengths"]]
  })



}
