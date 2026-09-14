library(shiny)
library(tidyverse)

ui = fluidPage(
  titlePanel("STAT 131A Parametric sampling distribution simulation"),
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
        label = "Sample size",
        value = 30,
        min = 1
      ),
      selectInput(
        inputId = "speed",
        label = "Simulation Speed",
        choices = c("Slow", "Normal", "Fast"),
        selected = "Normal"
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
      
      plotOutput(outputId='data_dist'),
      plotOutput(outputId='sampling_dist'),
      
      # # tab layout
      # tabsetPanel(
      #   
      #   tabPanel(
      #     title="Distribution",
      #     plotOutput(outputId='mygraph')
      #   ),
      #   # tabPanel("Summary",
      #   #          
      #   #          # plots on same row
      #   #          fluidRow(
      #   #            
      #   #            column(6,
      #   #                   
      #   #                   plotOutput('voltrack')
      #   #            ),
      #   #            column(6,
      #   #                   plotOutput('pie')  
      #   #            )
      #   #            
      #   #          )
      #   #          
      #   # )
      # ),
      
      # # visual data on same row
      # fluidRow(
      #   column(width=1,textOutput("curr_sim")),
      #   column(width=1,textOutput("curr_mean")),
      #   column(width=9,textOutput("curr_sample"))
      # )
      
    )
  )
  
)

server = function(input,output){
  
  # reactive to store all reactive variables
  var_list = reactiveValues() 
  
  # current simulation (i.e., sample) number
  var_list$curr_sim = 1
  
  var_list$resetindicator = 0  
  
  # stored the current sample
  var_list$curr_sample = NA
  
  # stores the current mean of the sample
  var_list$curr_mean = NA
  
  # initialize a big vector to store means
  var_list$mean_vec = rep(NA, 10000)

  
  forward = function() {
    
    req(input$dist)
    req(input$p)
    req(input$a)
    req(input$b)
    req(input$n)
    req(input$speed)

    if (input$dist=="Bernoulli") {
      
      var_list$curr_sample = rbinom(n=input$n, size=1, prob=input$p)
      
    } else if (input$dist=="Uniform"){ 
      
      var_list$curr_sample = runif(n=input$n, min=input$a, max=input$b)
      
    } else if (input$dist=="Normal") {
      
      var_list$curr_sample = rnorm(n=input$n, mean=input$m, sd=input$sd)
      
    } else if (input$dist=="Binomial") {
      
      var_list$curr_sample = rbinom(n=input$n, size=input$size, prob=input$binom_prob)
      
    } else if (input$dist=="Geometric") {
      
      var_list$curr_sample = rgeom(n=input$n, prob=input$geom_prob)
      
    }
    
    var_list$curr_mean = mean(var_list$curr_sample)
    
    var_list$mean_vec[var_list$curr_sim] = var_list$curr_mean
    
    var_list$curr_sim = var_list$curr_sim + 1
    
  }

  session = reactiveValues()
  session$timer = reactiveTimer(Inf)
  
  # handles the time steps of the animation
  observeEvent(
    eventExpr=input$play,
    handlerExpr={
      session$timer=reactiveTimer(
        intervalMs = 
          if (input$speed == 'Slow') {
            10000
          } else if (input$speed == 'Normal') {
            3000
          } else if (input$speed == 'Fast') {
            200
          }
      )
      observeEvent(
        eventExpr=session$timer(),
        handlerExpr={
          forward()
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
    
      var_list$curr_sim = 1
      
      var_list$resetindicator = 0   # used to change button labels
      
      var_list$curr_sample = NA
      var_list$curr_mean = NA
      
      var_list$mean_vec = rep(NA, 10000)
      
      session$timer = reactiveTimer(Inf)
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
      } else if (input$dist == "Geometric") {
        0
      } else if (input$dist == "Binomial") {
        0
      } else if (input$dist == "Normal") {
        input$m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli") {
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
        var_list$curr_sample,
        xlim = c(xmin, xmax),
        main = paste0('Distribution of random sample #', var_list$curr_sim),
        xlab = 'Values from a single random sample'
      )
      
      # add a vertical line at the sample mean
      abline(v=var_list$curr_mean, col="red")
      
      sample_sd = sd(var_list$curr_sample)
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "Sample mean: ", 
        round(var_list$curr_mean, 3),
        "    ",
        "Sample standard deviation (SD): ", 
        round(sample_sd, 3),
        "    ",
        "Sample size (n): ",
        input$n
        # "    ",
        # "SD / √n: ",
        # round(sample_sd / sqrt(input$n), 3)
      ))
    
    }
    
  })
  
  output$sampling_dist<-renderPlot({
    
    if (is.na(var_list$mean_vec[1])) {
      # don't render if no data
      return()
    } else {
      
      xmin = if (input$dist == "Bernoulli") {
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
      
      xmax = if (input$dist == "Bernoulli") {
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
        var_list$mean_vec,
        xlim = c(xmin, xmax),
        main = 'A growing sampling distribution',
        xlab = paste('Sample means from', var_list$curr_sim, 'random samples')
      )
      
      # add a vertical line at the sampling distrbution mean
      abline(v=mean(var_list$mean_vec, na.rm=TRUE), col="blue", lwd=3)
      
      # write the sample mean and SD above the plot
      mtext(side=3, text=paste0(
        "Mean of sampling distribution: ", 
        round(mean(var_list$mean_vec, na.rm=TRUE), 3),
        "    ",
        "Std. dev. of sampling distribution (Std. error): ", 
        round(sd(var_list$mean_vec, na.rm=TRUE), 3),
        "    ",
        "Number of sample means plotted: ",
        var_list$curr_sim
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