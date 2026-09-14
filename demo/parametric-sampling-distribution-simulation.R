library(shiny)

SIMULATION_INTERVALS <- c(
  Slow = 10000,
  Normal = 3000,
  Fast = 200,
  "Super Fast" = 50
)

ui <- fluidPage(
  titlePanel("STAT 131A Parallel Universe Simulator"),
  hr(style = "border-color: grey;"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "dist",
        label = "Data distribution",
        choices = c("Bernoulli", "Uniform", "Normal")
      ),
      conditionalPanel(
        condition = "input.dist == 'Bernoulli'",
        numericInput(
          inputId = "p",
          label = "Probability of success ( p )",
          value = 0.5,
          min = 0,
          max = 1
        )
      ),
      conditionalPanel(
        condition = "input.dist == 'Uniform'",
        numericInput(
          inputId = "a",
          label = "Minimum ( a )",
          value = 0
        )
      ),
      conditionalPanel(
        condition = "input.dist == 'Uniform'",
        numericInput(
          inputId = "b",
          label = "Maximum ( b )",
          value = 1
        )
      ),
      conditionalPanel(
        condition = "input.dist == 'Normal'",
        numericInput(
          inputId = "m",
          label = "Mean ( \u03bc )",
          value = 0
        )
      ),
      conditionalPanel(
        condition = "input.dist == 'Normal'",
        numericInput(
          inputId = "sd",
          label = "Standard Deviation ( \u03c3 )",
          value = 1,
          min = 0
        )
      ),
      numericInput(
        inputId = "n",
        label = "Sample size",
        value = 30,
        min = 1
      ),
      selectInput(
        inputId = "speed",
        label = "Simulation Speed",
        choices = names(SIMULATION_INTERVALS),
        selected = "Normal"
      ),
      hr(style = "border-color: grey;"),
      fluidRow(
        column(width = 4, actionButton("reset", "Reset")),
        column(width = 4, actionButton("stop", "Stop")),
        column(width = 4, actionButton("play", "Play"))
      )
    ),
    mainPanel(
      plotOutput(outputId = "data_dist"),
      plotOutput(outputId = "sampling_dist")
    )
  )
)

server <- function(input, output, session) {
  state <- reactiveValues(
    playing = FALSE,
    curr_sample = NULL,
    curr_mean = NA_real_,
    mean_vec = numeric()
  )

  simulation_number <- reactive(length(state$mean_vec))

  plot_limits <- reactive({
    switch(
      input$dist,
      Bernoulli = c(0, 1),
      Uniform = {
        req(input$a < input$b)
        c(input$a, input$b)
      },
      Normal = {
        req(input$sd > 0)
        input$m + c(-3, 3) * input$sd
      }
    )
  })

  draw_sample <- function() {
    req(input$dist, input$n, input$n >= 1)

    sample_values <- switch(
      input$dist,
      Bernoulli = {
        req(input$p >= 0, input$p <= 1)
        rbinom(input$n, size = 1, prob = input$p)
      },
      Uniform = {
        req(input$a < input$b)
        runif(input$n, min = input$a, max = input$b)
      },
      Normal = {
        req(input$sd > 0)
        rnorm(input$n, mean = input$m, sd = input$sd)
      }
    )

    sample_mean <- mean(sample_values)
    state$curr_sample <- sample_values
    state$curr_mean <- sample_mean
    state$mean_vec <- c(state$mean_vec, sample_mean)
  }

  observeEvent(input$play, {
    state$playing <- TRUE
  })

  observeEvent(input$stop, {
    state$playing <- FALSE
  })

  observeEvent(input$reset, {
    state$playing <- FALSE
    state$curr_sample <- NULL
    state$curr_mean <- NA_real_
    state$mean_vec <- numeric()
  })

  observe({
    req(state$playing)
    req(input$speed %in% names(SIMULATION_INTERVALS))

    invalidateLater(SIMULATION_INTERVALS[[input$speed]], session)
    isolate(draw_sample())
  })

  output$data_dist <- renderPlot({
    req(state$curr_sample)

    sample_values <- state$curr_sample
    sample_mean <- state$curr_mean
    sample_sd <- sd(sample_values)
    limits <- plot_limits()

    hist(
      sample_values,
      xlim = limits,
      main = paste0("Distribution of random sample #", simulation_number()),
      xlab = "Values from a single random sample"
    )

    abline(v = sample_mean, col = "red")

    mtext(
      side = 3,
      text = paste0(
        "Sample mean: ", round(sample_mean, 3),
        "    Sample standard deviation (SD): ", round(sample_sd, 3),
        "    Sample size (n): ", length(sample_values)
      )
    )
  })

  output$sampling_dist <- renderPlot({
    req(length(state$mean_vec) > 0)

    means <- state$mean_vec
    sampling_mean <- mean(means)
    sampling_sd <- sd(means)
    limits <- plot_limits()

    hist(
      means,
      xlim = limits,
      main = "A growing sampling distribution",
      xlab = paste("Sample means from", length(means), "random samples")
    )

    abline(v = sampling_mean, col = "blue", lwd = 3)

    mtext(
      side = 3,
      text = paste0(
        "Mean of sampling distribution: ", round(sampling_mean, 3),
        "    Std. dev. of sampling distribution (Std. error): ",
        round(sampling_sd, 3),
        "    Number of sample means plotted: ", length(means)
      )
    )
  })
}

runApp(shinyApp(ui, server), launch.browser = TRUE)
