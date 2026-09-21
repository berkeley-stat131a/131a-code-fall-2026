library(shiny)
library(tidyverse)

INF_EQUIVALENT = 20000

ui = fluidPage(
  titlePanel("STAT 131A Parametric bootstrap simulation"),
  hr(style="border-color: grey;"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        inputId = "dist",
        label = "Data distribution",
        choices = c("Bernoulli", "Uniform", "Normal", "Fair coin with fixed sample")
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
      conditionalPanel(
        condition = "input.dist == 'Normal'", 
        numericInput(
          inputId = "sd",
          label = "Standard Deviation ( \u03c3 )",
          value = 1
        )
      ),
      conditionalPanel(
        condition = "input.dist == 'Fair coin with fixed sample'", 
        numericInput(
          inputId = "fixed_n_heads",
          label = "Number of heads in sample",
          value = 5,
          min = 0,
        )
      ),
      # conditionalPanel(
      #   condition = "input.dist == 'Binomial'", 
      #   numericInput(
      #     inputId = "size",
      #     label = "Number of trials",
      #     value = 10,
      #     min = 1
      #   )
      # ),
      # conditionalPanel(
      #   condition = "input.dist == 'Binomial'", 
      #   numericInput(
      #     inputId = "binom_prob",
      #     label = "Probability of success",
      #     value = 0.5,
      #     min = 0,
      #     max = 1
      #   )
      # ),
      # conditionalPanel(
      #   condition = "input.dist == 'Geometric'", 
      #   numericInput(
      #     inputId = "geom_prob",
      #     label = "Probability of success",
      #     value = 0.5,
      #     min = 0,
      #     max = 1
      #   )
      # ),
      numericInput(
        inputId = "n",
        label = "Sample size (n)",
        value = 30,
        min = 2
      ),
      selectInput(
        inputId = "theta_hat",
        label = "Estimator",
        choices = c("Mean", "Median", "Minimum", "Maximum", "Quantile")
        # choices = c("Bernoulli", "Uniform", "Normal", "Binomial", "Geometric")
      ),
      conditionalPanel(
        condition = "input.theta_hat == 'Quantile'", 
        numericInput(
          inputId = "q",
          label = "Quantile",
          value = 0.5,
          min = 0,
          max = 1
        )
      ),
      numericInput(
        inputId = "B",
        label = "Number of bootstrap resamples (B)",
        value = 1000,
        min = 1
      ),
      # selectInput(
      #   inputId = "speed",
      #   label = "Simulation Speed",
      #   choices = c("Slow", "Normal", "Fast"),
      #   selected = "Normal"
      # ),
      
     hr(style="border-color: grey;"),
     fluidRow(
       # column(2, actionButton("reset","Reset")),
       # column(width=7,actionButton("stop","Stop")),
       column(2, actionButton("play","Play"))
     )
    ),
    
    # plot panel
    mainPanel(
      
      # fluidRow(
      #  column(6,plotOutput('data_dist')),
      #  column(6,plotOutput('resampled_data_dist'))
      # ),
      
      plotOutput('data_dist', height='200px'),
      plotOutput('sampling_dist', height='200px'),
      plotOutput('synth_sampling_dist', height='200px'),
      
    )
  )
  
)

server = function(input,output){
  
  # reactive to store all reactive variables
  var_list = reactiveValues() 
  
  # current simulation (i.e., sample) number
  # var_list$curr_sim = 1
  
  # var_list$resetindicator = 0  
  
  var_list$orig_sample = NA
  
  var_list$estimator_func = NA
  
  # stored the current sample
  # var_list$curr_synth_sample = NA
  
  # stores the current mean of the resample
  # var_list$curr_synth_estimate = NA
  
  # initialize a big vector to store real estimates
  var_list$real_estimates = rep(NA, INF_EQUIVALENT)
  
  # initialize a big vector to store synth_estimates
  var_list$synth_estimates = rep(NA, INF_EQUIVALENT)

  
  forward = function() {
    
    # distribution choice
    req(input$dist)
    
    # bernoulli
    req(input$p)
    
    # uniform
    req(input$a)
    req(input$b)
    
    # normal
    req(input$m)
    req(input$sd)
    
    # sample size
    req(input$n)
    
    req(input$theta_hat)
    
    # If quantile selected as estimator
    req(input$q)
    
    # req(input$speed)
    
    if (input$theta_hat == "Mean") {
      
      var_list$estimator_func = mean
      
    } else if (input$theta_hat == "Median") {
      
      var_list$estimator_func = median
      
    } else if (input$theta_hat == "Quantile") {
      
      var_list$estimator_func = function(x) quantile(x, probs=input$q)
      
    } else if (input$theta_hat == "Minimum") {
      
      var_list$estimator_func = min
      
    } else if (input$theta_hat == "Maximum") {
      
      var_list$estimator_func = max
      
    }
    
    # initialize the real sample and real sampling distribution
    if (is.na(var_list$real_estimates[1])) {
      if (input$dist=="Bernoulli") {
        
        var_list$orig_sample = rbinom(n=input$n, size=1, prob=input$p)
        
        var_list$real_estimates = replicate(
          INF_EQUIVALENT, rbinom(n=input$n, size=1, prob=input$p)
          ) %>% 
          apply(2, var_list$estimator_func)
        
      } else if (input$dist=="Uniform"){ 
        
        var_list$orig_sample = runif(n=input$n, min=input$a, max=input$b)
        
        var_list$real_estimates = replicate(
          INF_EQUIVALENT, runif(n=input$n, min=input$a, max=input$b)
          ) %>% 
          apply(2, var_list$estimator_func)
        
      } else if (input$dist=="Normal") {
        
        var_list$orig_sample = rnorm(n=input$n, mean=input$m, sd=input$sd)
        
        var_list$real_estimates = replicate(
          INF_EQUIVALENT, rnorm(n=input$n, mean=input$m, sd=input$sd)
          ) %>% 
          apply(2, var_list$estimator_func)
      } else if (input$dist=="Fair coin with fixed sample") {
      
        var_list$orig_sample = rep(c(1, 0), c(input$fixed_n_heads, input$n - input$fixed_n_heads))
        
        var_list$real_estimates = replicate(
          INF_EQUIVALENT, rbinom(n=input$n, size=1, prob=0.5)
          ) %>% 
          apply(2, var_list$estimator_func)
      }
      
      var_list$estimator_name = if (input$theta_hat == "Mean") {
        "Sample mean"
      } else if (input$theta_hat == "Median") {
        "Sample median"
      } else if (input$theta_hat == "Mode") {
        "Sample mode"
      } else if (input$theta_hat == "Quantile") {
        paste0(input$q, " quantile")
      } else if (input$theta_hat == "Maximum") {
        "Sample maximum"
      } else if (input$theta_hat == "Minimum") {
        "Sample minimum"
      }
      
    }
    
    var_list$synth_estimates = replicate(
        input$B, 
        sample(var_list$orig_sample, size=length(var_list$orig_sample), replace=TRUE)
      ) %>% 
      apply(2, var_list$estimator_func)
    
    # var_list$curr_synth_sample = 
    #   sample(var_list$orig_sample, size=input$n, replace=TRUE)
    # 
    # var_list$curr_synth_estimate = var_list$estimator_func(var_list$curr_synth_sample)
    # 
    # var_list$synth_estimates[var_list$curr_sim] = var_list$curr_synth_estimate
    # 
    # var_list$curr_sim = var_list$curr_sim + 1
  }

  # session = reactiveValues()
  # session$timer = reactiveTimer(Inf)
  
  # handles the time steps of the animation
  observeEvent(
    eventExpr=input$play,
    handlerExpr={ 
      
      var_list$orig_sample = NA
      
      # initialize a big vector to store real estimates
      var_list$real_estimates = rep(NA, INF_EQUIVALENT)
      
      # initialize a big vector to store synth_estimates
      var_list$synth_estimates = rep(NA, INF_EQUIVALENT)
      
      forward() 
    }
    # handlerExpr={
    #   session$timer=reactiveTimer(
    #     intervalMs = 
    #       if (input$speed == 'Slow') {
    #         20000
    #       } else if (input$speed == 'Normal') {
    #         3000
    #       } else if (input$speed == 'Fast') {
    #         200
    #       }
    #   )
    #   observeEvent(
    #     eventExpr=session$timer(),
    #     handlerExpr={
    #       forward()
    #     }
    #   )
    # }
  )

  # # handles the stop button
  # observeEvent(
  #   eventExpr=input$stop,
  #   handlerExpr={
  #     # resets the timer
  #     session$timer = reactiveTimer(Inf)
  #   }
  # )

  # ## handles the reset button (sets everything to original values)
  # observeEvent(
  #   eventExpr=input$reset,
  #   handlerExpr={
  #     
  #     # session$timer = reactiveTimer(Inf)
  #     
  #     # current simulation (i.e., sample) number
  #     # var_list$curr_sim = 1
  #     
  #     # var_list$resetindicator = 0  
  #     
  #     var_list$orig_sample = NA
  #     
  #     # stored the current sample
  #     # var_list$curr_synth_sample = NA
  #     
  #     # stores the current mean of the resample
  #     # var_list$curr_synth_estimate = NA
  #     
  #     # initialize a big vector to store real estimates
  #     var_list$real_estimates = rep(NA, INF_EQUIVALENT)
  #     
  #     # initialize a big vector to store synth_estimates
  #     var_list$synth_estimates = rep(NA, INF_EQUIVALENT)
  #   }
  # )
  
  # data distribution
  output$data_dist = renderPlot({
    
    if (is.na(var_list$orig_sample[1])) {
      
      # don't render if no data
      return()
      
    } else {
      
      xmin = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        0
      } else if (input$dist == "Uniform") {
        input$a
      } else if (input$dist == "Geometric") {
        0
      } else if (input$dist == "Binomial") {
        0
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        1
      } else if (input$dist == "Uniform") {
        input$b
      } else if (input$dist == "Geometric") {
        10
      } else if (input$dist == "Binomial") {
        input$size
      } else if (input$dist == "Normal") {
        input$m + 3*input$sd
      }
      
      hist(
        var_list$orig_sample,
        xlim = c(xmin, xmax),
        main = paste0('Distribution of data in original sample (observed)'),
        xlab = 'Possible values of random variable'
      )
      
      estimate = var_list$estimator_func(var_list$orig_sample)
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        var_list$estimator_name, " (point estimate; theta hat): ", round(estimate, 3)
      ))
      
      abline(v=estimate, col='red')
    
    }
    
  })
  
  # # resampled data distribution
  # output$resampled_data_dist = renderPlot({
  #   
  #   if (is.na(var_list$curr_synth_sample[1])) {
  #     
  #     # don't render if no resampled data yet
  #     return()
  #     
  #   } else {
  #     
  #     xmin = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
  #       0
  #     } else if (input$dist == "Uniform") {
  #       input$a
  #     } else if (input$dist == "Geometric") {
  #       0
  #     } else if (input$dist == "Binomial") {
  #       0
  #     } else if (input$dist == "Normal") {
  #       input$m - 3*input$sd
  #     }
  #     
  #     xmax = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
  #       1
  #     } else if (input$dist == "Uniform") {
  #       input$b
  #     } else if (input$dist == "Geometric") {
  #       10
  #     } else if (input$dist == "Binomial") {
  #       input$size
  #     } else if (input$dist == "Normal") {
  #       input$m + 3*input$sd
  #     }
  #     
  #     hist(
  #       var_list$curr_synth_sample,
  #       xlim = c(xmin, xmax),
  #       main = paste0('Distribution of resample ', var_list$curr_sim),
  #       xlab = paste0('Values from resample ', var_list$curr_sim)
  #     )
  #     
  #   }
  #   
  # })
  
  output$sampling_dist = renderPlot({
    
    if (is.na(var_list$real_estimates[1])) {
      # don't render if no data
      return()
    } else {
      
      xmin = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        0
      } else if (input$dist == "Uniform") {
        input$a
      } else if (input$dist == "Geometric") {
        0
      } else if (input$dist == "Binomial") {
        0
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        1
      } else if (input$dist == "Uniform") {
        input$b
      } else if (input$dist == "Geometric") {
        10
      } else if (input$dist == "Binomial") {
        input$size
      } else if (input$dist == "Normal") {
        input$m + 3*input$sd
      }
      
      true_sampling_kde = density(
        var_list$real_estimates,
        na.rm = TRUE,
        adjust = 3,
        from = xmin,
        to = xmax
      )

      plot(
        true_sampling_kde,
        xlim = c(xmin, xmax),
        main = paste0('True sampling distribution of ', var_list$estimator_name, ' (unobserved)'),
        xlab = paste('Possible values of estimator'),
        ylab = 'Density',
        type = 'l',
        lwd = 2
      )
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "SD of true sampling distribution (true SE): ", 
        round(sd(var_list$real_estimates, na.rm=TRUE), 3)
        # "    ",
        # "Number of true estimates plotted: ",
        # INF_EQUIVALENT
      ))
    }
    
  })
  
  output$synth_sampling_dist = renderPlot({
    
    if (is.na(var_list$synth_estimates[1])) {
      # don't render if no data
      return()
    } else {
      
      xmin = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        0
      } else if (input$dist == "Uniform") {
        input$a
      } else if (input$dist == "Geometric") {
        0
      } else if (input$dist == "Binomial") {
        0
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli" | input$dist == "Fair coin with fixed sample") {
        1
      } else if (input$dist == "Uniform") {
        input$b
      } else if (input$dist == "Geometric") {
        10
      } else if (input$dist == "Binomial") {
        input$size
      } else if (input$dist == "Normal") {
        input$m + 3*input$sd
      }
      
      hist(
        var_list$synth_estimates,
        xlim = c(xmin, xmax),
        main = paste0('Synthetic sampling distribution of ', var_list$estimator_name, ' (observed)'),
        xlab = paste('Synthetic point estimates (synth. theta hats)')
      )
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "SD of synthetic sampling distribution (synth. SE): ", 
        round(sd(var_list$synth_estimates, na.rm=TRUE), 3),
        "    ",
        "Number of synthetic estimates plotted (B): ",
        input$B
      ))
    }
    
  })
  
  # ## visual data outputs
  # output$curr_sim = renderText({
  #   paste("Simulation number: ", var_list$curr_sim)
  # })
  
  
  
  # ## pie plot output
  # output$pie<-renderPlot({
  #   if(sum(waits$slices)==0){
  #     return() # dont error if no data
  #   }
  #   par(xpd=TRUE) # allow legend outside plot margins
  #   pie(waits$slices,labels = waits$slices,col=c("green3","yellow2","red2"))
  #   legend(-1.3,-0.5,legend=c("<26","26 to 25","36 +"),fill=c("green3","yellow2","red2"),bty='n',cex=0.8)
  # })
  # 
  # output$voltrack<-renderPlot({
  #   ytop<-max(waits$trackvol)+100 # upper limit for y axis
  #   plot(waits$trackvol,type='l',ylab="Waitlist Volume",xaxt='n',xlab="Week Number",bty='n',lwd=1.5,ylim=c(0,ytop),main="20 Week Volume")
  #   points(20,tail(waits$trackvol,1),pch=21,bg="blue")
  #   axis(1,at=1:20,labels = waits$xdates)
  # })
  
}

runApp(shinyApp(ui,server),launch.browser = TRUE)
