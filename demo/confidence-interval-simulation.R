library(shiny)
library(tidyverse)

INF_EQUIVALENT = 10000
MAX_SIM = 10000

ui = fluidPage(
  titlePanel("STAT 131A Confidence interval simulation"),
  hr(style="border-color: grey;"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "dist",
        label = "Data distribution",
        choices = c("Bernoulli", "Uniform", "Normal")
        # choices = c("Bernoulli", "Uniform", "Normal", "Binomial", "Geometric")
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
          value = 1
        )
      ),
      numericInput(
        inputId = "n",
        label = "Sample size",
        value = 30,
        min = 1
      ),
      numericInput(
        inputId = "cl",
        label = "Confidence level (1 - \u03b1)",
        value = 0.95,
        min = 0.001,
        max = 0.999
      ),
      selectInput(
        inputId = "speed",
        label = "Simulation Speed",
        choices = c("Standard", "Fast", "Super fast"),
        selected = "Standard"
      ),
     hr(style="border-color: grey;"),
     fluidRow(
       column(width=7,actionButton("reset","Reset")),
       column(width=7,actionButton("stop","Stop")),
       column(width=7,actionButton("play","Play"))
     )
    ),
    
    # plot panel
    mainPanel(
      
      textOutput(outputId='running_coverage'),
      plotOutput(outputId='data_dist', height='200px'),
      plotOutput(outputId='normal_dist', height='200px'),
      plotOutput(outputId='sampling_dist', height='200px'),
      
    )
  )
  
)

reset_vars = function(var_list) {
  # current simulation (i.e., sample) number
  var_list$curr_sim = 0
  
  var_list$true_mean = NA
  
  # stores the se from the simulated sampling distribution
  var_list$se = NA
  
  # stores the estimates of the simulated sampling distribution
  var_list$real_estimates = NA
  
  # cached density of real_estimates
  var_list$real_density = NA
  
  # stored the current sample
  var_list$curr_sample = NA
  
  var_list$means = NA
  
  # stores the current mean of the sample
  var_list$curr_mean = NA
  
  var_list$ci_lower = NA
  var_list$ci_upper = NA
  
  # stores the current confidence interval
  var_list$curr_ci = NA
  
  var_list$samples_to_iter = NA
  
  var_list$contains_true_mean = NA
  
  var_list$contains_true_mean_vec = rep(NA, 10000)
}

forward = function(var_list, input) {
  
  req(input$dist)
  req(input$p)
  req(input$a)
  req(input$b)
  req(input$m)
  req(input$sd)
  req(input$n)
  req(input$cl)
  req(input$speed)
  
  # initialize true sampling distribution
  if (is.na(var_list$real_estimates[1])) {
    
    if (input$dist=="Bernoulli") {
      
      var_list$real_estimates = colMeans(matrix(
        rbinom(INF_EQUIVALENT * input$n, size=1, prob=input$p),
        nrow = input$n
      ))
      
      var_list$true_mean = input$p
      
      var_list$samples_to_iter = matrix(
        rbinom(MAX_SIM * input$n, size=1, prob=input$p),
        nrow = input$n
      )
      
    } else if (input$dist=="Uniform") { 
      
      var_list$real_estimates = colMeans(matrix(
        runif(INF_EQUIVALENT * input$n, min=input$a, max=input$b),
        nrow = input$n
      ))
      
      var_list$true_mean = (input$a + input$b) / 2
      
      var_list$samples_to_iter = matrix(
        runif(MAX_SIM * input$n, min=input$a, max=input$b),
        nrow = input$n
      )
      
    } else if (input$dist=="Normal") {
      var_list$real_estimates = colMeans(matrix(
        rnorm(INF_EQUIVALENT * input$n, mean=input$m, sd=input$sd),
        nrow = input$n
      ))
      
      var_list$true_mean = input$m
      
      var_list$samples_to_iter = matrix(
        rnorm(MAX_SIM * input$n, mean=input$m, sd=input$sd),
        nrow = input$n
      )
    }
    
    var_list$means = colMeans(var_list$samples_to_iter)
    
    var_list$se = sd(var_list$real_estimates)
    
    # cache density of sampling distribution for rendering
    var_list$real_density = density(var_list$real_estimates, adjust = 5)
    
    alpha = 1 - input$cl
    var_list$ci_lower = var_list$means - qnorm(1 - alpha/2) * var_list$se
    var_list$ci_upper = var_list$means + qnorm(1 - alpha/2) * var_list$se
    var_list$contains_true_mean_vec = 
      (var_list$true_mean > var_list$ci_lower) &
      (var_list$true_mean < var_list$ci_upper)
  }
  
  if (input$speed=="Standard") {
    samples_per_iter = 1
  } else if (input$speed=="Fast") {
    samples_per_iter = 1
  } else if (input$speed=="Super fast") {
    samples_per_iter = 10
  }
  
  var_list$curr_sample = var_list$samples_to_iter[, var_list$curr_sim]
  var_list$curr_mean = var_list$means[var_list$curr_sim]
  var_list$curr_ci = c(
    var_list$ci_lower[var_list$curr_sim], 
    var_list$ci_upper[var_list$curr_sim]
  )
  var_list$contains_true_mean = var_list$contains_true_mean_vec[var_list$curr_sim]
  
  var_list$curr_sim = var_list$curr_sim + samples_per_iter
}
  
  

server = function(input,output){
  
  # reactive to store all reactive variables
  var_list = reactiveValues() 
  
  reset_vars(var_list)

  session = reactiveValues()
  session$timer = reactiveTimer(Inf)
  
  # handles the time steps of the animation
  observeEvent(
    eventExpr=input$play,
    handlerExpr={
      session$timer=reactiveTimer(
        intervalMs = 
          if (input$speed == 'Standard') {
            3000
          } else if (input$speed == 'Fast') {
            300
          } else if (input$speed == 'Super fast') {
            300
          }
      )
      observeEvent(
        eventExpr=session$timer(),
        handlerExpr={
          forward(var_list, input)
        }
      )
    }
  )

  # handles the stop button
  observeEvent(
    eventExpr=input$stop,
    handlerExpr={
      # resets the timer
      session$timer = reactiveTimer(Inf)
    }
  )

  ## handles the reset button (sets everything to original values)
  observeEvent(
    eventExpr=input$reset,
    handlerExpr={
      
      session$timer = reactiveTimer(Inf)
      
      reset_vars(var_list)
    }
  )
  
  # data distribution
  output$data_dist = renderPlot({
    
    if (is.na(var_list$curr_sample[1])) {
      
      # don't render if no data
      return()
      
    } else {
      
      xmin = if (input$dist == "Bernoulli") {
        0
      } else if (input$dist == "Uniform") {
        input$a
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli") {
        1
      } else if (input$dist == "Uniform") {
        input$b
      } else if (input$dist == "Normal") {
        input$m + 3*input$sd
      }
      
      hist(
        var_list$curr_sample,
        xlim = c(xmin, xmax),
        main = paste0('Distribution of random sample #', var_list$curr_sim),
        xlab = 'Possible values of random variable'
      )
      
      # add a vertical line at the true mean
      abline(v=var_list$true_mean, col="red", lwd=3)
      
      # add vertical lines at confidence interval
      abline(v=var_list$curr_ci[1], col="blue", lwd=3, lty=2)
      abline(v=var_list$curr_ci[2], col="blue", lwd=3, lty=2)
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "Sample mean: ", 
        round(var_list$curr_mean, 3),
        "    ",
        "Sample size (n): ",
        input$n,
        "    ",
        "Confidence interval (CI): [",
        paste0(round(var_list$curr_ci, 3), collapse=","),
        "]",
        "    ",
        "CI captures true mean: ",
        var_list$contains_true_mean
      ))
      
    }
    
  })
  
  output$sampling_dist = renderPlot({
    
    if (is.na(var_list$real_estimates[1])) {
      # don't render if no data
      return()
    } else {
      
      xmin = if (input$dist == "Bernoulli") {
        0
      } else if (input$dist == "Uniform") {
        input$a
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli") {
        1
      } else if (input$dist == "Uniform") {
        input$b
      } else if (input$dist == "Normal") {
        input$m + 3*input$sd
      }
      
      plot(
        var_list$real_density,
        xlim = c(xmin, xmax),
        main = paste0('Sampling distribution of estimator given n=', input$n),
        xlab = 'Possible values of estimator',
        lwd = 2
      )
      
      # add a vertical line at the true mean
      abline(v=var_list$true_mean, col="red", lwd=3)
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "Mean of sampling distribution: ",
        round(var_list$true_mean, 3),
        "    ",
        "SD of sampling distribution (SE): ", 
        round(var_list$se, 3)
        # "    ",
        # "Number of true estimates plotted: ",
        # INF_EQUIVALENT
      ))
    }
    
  })
  
  output$normal_dist = renderPlot({
    
    curve(
      dnorm(x, 0, 1), 
      from=-3.5, 
      to=3.5,
      xname='x',
      ylab='Density',
      main='Standard normal distribution: N(μ=0, σ=1)'
    )
    
    alpha = 1 - input$cl
    z_alpha2 = qnorm(1 - alpha/2)
    
    abline(v=z_alpha2, col='black')
    abline(v=-1 * z_alpha2, col='black')
    # quote('x = \u03a6'^{"-1"}~'(Area left of x)'),
    mtext(side=3, text=paste0(
      "Z(1 - \u03b1/2) = ",
      round(z_alpha2, 3),
      "    ",
      "CL = ",
      input$cl,
      " = ",
      "Area between Z(\u03b1/2) and Z(1 - \u03b1/2)"
    ))
    
    
    
  })
  
  output$running_coverage = renderText({
    
    if (is.na(var_list$contains_true_mean_vec[1])) {
      # don't render if no data
      return()
    } else {
      running_coverage = mean(var_list$contains_true_mean_vec[1:var_list$curr_sim])
    
      paste0(
        "Proportion of ",
        var_list$curr_sim,
        " random confidence intervals that capture the true fixed mean: ",
        round(running_coverage, 3)
      )
    }
    
  })
  
}

runApp(shinyApp(ui,server),launch.browser = TRUE)