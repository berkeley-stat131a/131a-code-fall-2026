library(shiny)

INF_EQUIVALENT = 10000
MAX_SIM = 10000

ui = fluidPage(
  withMathJax(),
  tags$head(tags$style(HTML("
    .well { padding: 8px 12px; margin-bottom: 8px; }
    .well h4 { margin-top: 0; margin-bottom: 6px; }
    .well .form-group { margin-bottom: 6px; }
    .well label { margin-bottom: 2px; font-weight: normal; }
    .well .form-control { height: 28px; padding: 2px 8px; }
    .well select.form-control { height: 30px; }
    .well .selectize-input { min-height: 28px; padding: 3px 8px; }
    .well .row .btn { padding: 3px 10px; }
  "))),
  titlePanel("STAT 131A Hypothesis test simulation"),
  hr(style="border-color: grey;"),
  sidebarLayout(
    sidebarPanel(
      wellPanel(
        h4("Data-generation"),
        selectInput(
          inputId = "dist",
          label = "Data distribution",
          choices = c("Bernoulli", "Uniform", "Normal")
          # choices = c("Bernoulli", "Uniform", "Normal", "Binomial", "Geometric")
        ),
        conditionalPanel(
          condition = "input.dist == 'Bernoulli'",
          numericInput(
            inputId = "true_p",
            label = "True probability of success \\(p\\)",
            value = 0.5,
            min = 0,
            max = 1
          )
        ),
        conditionalPanel(
          condition = "input.dist == 'Uniform'",
          # minimum is fixed at 0: shown for reference but disabled
          div(
            class = "form-group shiny-input-container",
            tags$label("Minimum \\(a\\)"),
            tags$input(type = "number", class = "form-control", value = 0, disabled = NA)
          ),
          numericInput(
            inputId = "true_u_mean",
            label = "True mean \\(\\mu\\)",
            value = 1,
            min = 0.001
          ),
          # maximum is determined by the mean (b = 2 * mean): autofilled and disabled
          htmltools::tagQuery(numericInput(
            inputId = "true_b",
            label = "True maximum \\(b\\)",
            value = 2,
            min = 0
          ))$find("input")$addAttrs(disabled = NA)$allTags()
        ),
        conditionalPanel(
          condition = "input.dist == 'Normal'",
          numericInput(
            inputId = "true_m",
            label = "True mean \\(\\mu\\)",
            value = 0
          ),
          numericInput(
            inputId = "sd",
            label = "Standard deviation \\(\\sigma\\)",
            value = 1
          )
        ),
        numericInput(
          inputId = "n",
          label = "Sample size",
          value = 30,
          min = 1
        )
      ),
      wellPanel(
        h4("Testing choices"),
        conditionalPanel(
          condition = "input.dist == 'Bernoulli'",
          numericInput(
            inputId = "null_p",
            label = "Null probability of success \\(p_{H_0}\\)",
            value = 0.5,
            min = 0,
            max = 1
          )
        ),
        conditionalPanel(
          condition = "input.dist == 'Uniform'",
          # minimum is fixed at 0: shown for reference but disabled
          div(
            class = "form-group shiny-input-container",
            tags$label("Null minimum \\(a_{H_0}\\)"),
            tags$input(type = "number", class = "form-control", value = 0, disabled = NA)
          ),
          numericInput(
            inputId = "null_u_mean",
            label = "Null mean \\(\\mu_{H_0}\\)",
            value = 1,
            min = 0.001
          ),
          # maximum is determined by the mean (b = 2 * mean): autofilled and disabled
          htmltools::tagQuery(numericInput(
            inputId = "null_b",
            label = "Null maximum \\(b_{H_0}\\)",
            value = 2,
            min = 0
          ))$find("input")$addAttrs(disabled = NA)$allTags()
        ),
        conditionalPanel(
          condition = "input.dist == 'Normal'",
          numericInput(
            inputId = "null_m",
            label = "Null mean \\(\\mu_{H_0}\\)",
            value = 0
          )
        ),
        numericInput(
          inputId = "alpha",
          label = "Significance level \\(\\alpha\\)",
          value = 0.05,
          min = 0.001,
          max = 0.999
        ),
        # confidence level is derived from the significance level
        div(
          strong("Confidence level \\(1 - \\alpha\\): "),
          textOutput("conf_level", inline = TRUE)
        )
      ),
      wellPanel(
        selectInput(
          inputId = "speed",
          label = "Simulation Speed",
          choices = c("Standard", "Fast", "Super fast"),
          selected = "Standard"
        ),
        fluidRow(
          column(width=4,actionButton("reset","Reset")),
          column(width=4,actionButton("stop","Stop")),
          column(width=4,actionButton("play","Play"))
        )
      )
    ),
    
    # plot panel
    mainPanel(
      
      uiOutput(outputId='running_summary'),
      plotOutput(outputId='data_dist', height='200px'),
      # for Bernoulli the estimated null sampling distribution is exact, so only show that one
      conditionalPanel(
        condition = "input.dist != 'Bernoulli'",
        plotOutput(outputId='sampling_dist', height='200px')
      ),
      plotOutput(outputId='est_sampling_dist', height='200px'),
      
    )
  )
  
)

# Bound temporary allocations while simulating the reference sample means.
# Consecutive batches retain the same random-draw order as one large matrix.
simulate_means = function(draw, n, count, max_values = 1000000) {
  batch_size = max(1, floor(max_values / n))
  means = numeric(count)
  for (first in seq.int(1, count, by = batch_size)) {
    last = min(first + batch_size - 1, count)
    means[first:last] = colMeans(matrix(draw(n * (last - first + 1)), nrow = n))
  }
  means
}

# x-axis range covering both the true and null sampling distributions
sampling_xlim = function(input) {
  if (input$dist == "Bernoulli") {
    c(0, 1)
  } else if (input$dist == "Uniform") {
    c(0, 2 * max(input$true_u_mean, input$null_u_mean))
  } else if (input$dist == "Normal") {
    c(
      min(input$true_m, input$null_m) - 3*input$sd,
      max(input$true_m, input$null_m) + 3*input$sd
    )
  }
}

# small arrow below the x axis of the current plot, pointing up at x
arrow_below_axis = function(x, col) {
  usr = par("usr")
  y_per_inch = (usr[4] - usr[3]) / par("pin")[2]
  arrows(
    x0 = x, y0 = usr[3] - 0.14 * y_per_inch,
    x1 = x, y1 = usr[3] - 0.03 * y_per_inch,
    col = col, lwd = 3, length = 0.08, xpd = TRUE
  )
}

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
  
  # null mean and cached density of the null sampling distribution
  var_list$null_mean = NA
  var_list$null_density = NA
  var_list$null_estimates = NA
  
  # running proportion of rejected tests, using the true null distribution
  # and the estimated (s / sqrt(n)) null distribution
  var_list$reject_true_prop = NA
  var_list$reject_est_prop = NA
  
  # stores the se from the simulated null sampling distribution
  var_list$null_se = NA
  
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

  var_list$contains_true_mean_vec = NA
  var_list$coverage = NA
  
  # confidence intervals built with the SE estimated from each sample (s / sqrt(n))
  var_list$ci_est_lower = NA
  var_list$ci_est_upper = NA
  var_list$contains_true_mean_est_vec = NA
  var_list$coverage_est = NA
  
  # estimated SEs (s / sqrt(n)) of every simulated sample
  var_list$est_ses = NA
}

forward = function(var_list, input) {
  
  req(input$dist)
  req(input$sd)
  req(input$n)
  req(input$alpha)
  req(input$speed)

  if (input$dist == "Bernoulli") {
    req(input$true_p, input$null_p)
  } else if (input$dist == "Uniform") {
    req(input$true_u_mean, input$null_u_mean)
    # the minimum is fixed at 0, so the mean (= maximum / 2) must be positive
    req(input$true_u_mean > 0, input$null_u_mean > 0)
  } else if (input$dist == "Normal") {
    req(input$true_m, input$null_m)
  }
  
  if (var_list$curr_sim >= MAX_SIM) return(invisible(NULL))

  # initialize true sampling distribution
  if (is.na(var_list$real_estimates[1])) {
    
    if (input$dist=="Bernoulli") {
      
      var_list$real_estimates = simulate_means(
        function(count) rbinom(count, size=1, prob=input$true_p),
        input$n, INF_EQUIVALENT
      )
      
      var_list$true_mean = input$true_p
      
      var_list$samples_to_iter = matrix(
        rbinom(MAX_SIM * input$n, size=1, prob=input$true_p),
        nrow = input$n
      )
      
    } else if (input$dist=="Uniform") { 
      
      var_list$real_estimates = simulate_means(
        function(count) runif(count, min=0, max=2 * input$true_u_mean),
        input$n, INF_EQUIVALENT
      )
      
      var_list$true_mean = input$true_u_mean
      
      var_list$samples_to_iter = matrix(
        runif(MAX_SIM * input$n, min=0, max=2 * input$true_u_mean),
        nrow = input$n
      )
      
    } else if (input$dist=="Normal") {
      var_list$real_estimates = simulate_means(
        function(count) rnorm(count, mean=input$true_m, sd=input$sd),
        input$n, INF_EQUIVALENT
      )
      
      var_list$true_mean = input$true_m
      
      var_list$samples_to_iter = matrix(
        rnorm(MAX_SIM * input$n, mean=input$true_m, sd=input$sd),
        nrow = input$n
      )
    }
    
    var_list$means = colMeans(var_list$samples_to_iter)
    
    var_list$se = sd(var_list$real_estimates)
    
    # cache density of sampling distribution for rendering
    # sample proportions must stay within [0, 1], so don't let the density extend past it
    density_range = if (input$dist == "Bernoulli") list(from = 0, to = 1) else list()
    var_list$real_density = do.call(
      density, c(list(var_list$real_estimates, adjust = 5), density_range)
    )
    
    # simulate the sampling distribution under the null hypothesis
    if (input$dist=="Bernoulli") {
      null_estimates = simulate_means(
        function(count) rbinom(count, size=1, prob=input$null_p),
        input$n, INF_EQUIVALENT
      )
      var_list$null_mean = input$null_p
    } else if (input$dist=="Uniform") {
      null_estimates = simulate_means(
        function(count) runif(count, min=0, max=2 * input$null_u_mean),
        input$n, INF_EQUIVALENT
      )
      var_list$null_mean = input$null_u_mean
    } else if (input$dist=="Normal") {
      null_estimates = simulate_means(
        function(count) rnorm(count, mean=input$null_m, sd=input$sd),
        input$n, INF_EQUIVALENT
      )
      var_list$null_mean = input$null_m
    }
    var_list$null_estimates = null_estimates
    var_list$null_se = sd(null_estimates)
    var_list$null_density = do.call(
      density, c(list(null_estimates, adjust = 5), density_range)
    )
    
    # run the hypothesis test on every simulated sample
    dev_from_null = abs(var_list$means - var_list$null_mean)
    null_devs = sort(abs(null_estimates - var_list$null_mean))
    
    # p-value from the simulated null distribution: P(|null estimate - null mean| >= deviation)
    p_true = 1 - findInterval(dev_from_null, null_devs, left.open = TRUE) / length(null_devs)
    reject_true = p_true < input$alpha
    
    # p-value from N(null mean, s / sqrt(n)), with s from each sample
    sample_sds = sqrt(pmax(
      colSums(var_list$samples_to_iter^2) - input$n * var_list$means^2, 0
    ) / (input$n - 1))
    p_est = 2 * pnorm(-dev_from_null / (sample_sds / sqrt(input$n)))
    reject_est = p_est < input$alpha
    reject_est[is.na(reject_est)] = FALSE
    # for Bernoulli the null sampling distribution is known exactly, so the
    # estimated null distribution is identical to the true one
    if (input$dist == "Bernoulli") reject_est = reject_true
    
    var_list$reject_true_prop = cumsum(reject_true) / seq_along(reject_true)
    var_list$reject_est_prop = cumsum(reject_est) / seq_along(reject_est)
    
    alpha = input$alpha
    margin = qnorm(1 - alpha/2) * var_list$se
    
    # for Uniform and Normal, also build the CI with the SE estimated from each sample
    var_list$est_ses = sample_sds / sqrt(input$n)
    margin_est = qnorm(1 - alpha/2) * var_list$est_ses
    var_list$ci_est_lower = var_list$means - margin_est
    var_list$ci_est_upper = var_list$means + margin_est
    var_list$contains_true_mean_est_vec =
      (var_list$true_mean > var_list$ci_est_lower) &
      (var_list$true_mean < var_list$ci_est_upper)
    var_list$contains_true_mean_est_vec[is.na(var_list$contains_true_mean_est_vec)] = FALSE
    var_list$coverage_est = cumsum(var_list$contains_true_mean_est_vec) /
      seq_along(var_list$contains_true_mean_est_vec)
    
    var_list$ci_lower = var_list$means - margin
    var_list$ci_upper = var_list$means + margin
    var_list$contains_true_mean_vec =
      (var_list$true_mean > var_list$ci_lower) &
      (var_list$true_mean < var_list$ci_upper)
    var_list$coverage = cumsum(var_list$contains_true_mean_vec) /
      seq_along(var_list$contains_true_mean_vec)
  }
  
  if (input$speed=="Standard") {
    samples_per_iter = 1
  } else if (input$speed=="Fast") {
    samples_per_iter = 1
  } else if (input$speed=="Super fast") {
    samples_per_iter = 50
  }
  
  var_list$curr_sim = min(var_list$curr_sim + samples_per_iter, MAX_SIM)

  var_list$curr_sample = var_list$samples_to_iter[, var_list$curr_sim]
  var_list$curr_mean = var_list$means[var_list$curr_sim]
  if (input$dist == "Bernoulli") {
    # the SE is known exactly, so there is only one confidence interval
    var_list$curr_ci = c(
      var_list$ci_lower[var_list$curr_sim],
      var_list$ci_upper[var_list$curr_sim]
    )
    var_list$contains_true_mean = var_list$contains_true_mean_vec[var_list$curr_sim]
  } else {
    # only show the CI constructed with the (incorrectly) estimated SE
    var_list$curr_ci = c(
      var_list$ci_est_lower[var_list$curr_sim],
      var_list$ci_est_upper[var_list$curr_sim]
    )
    var_list$contains_true_mean = var_list$contains_true_mean_est_vec[var_list$curr_sim]
  }
  
}
  
  

server = function(input, output, session){
  
  # reactive to store all reactive variables
  var_list = reactiveValues() 
  
  reset_vars(var_list)

  running = reactiveVal(FALSE)

  # the uniform maximum is determined by the mean since the minimum is fixed at 0
  observeEvent(input$true_u_mean, {
    req(input$true_u_mean)
    updateNumericInput(session, "true_b", value = 2 * input$true_u_mean)
  })

  observeEvent(input$null_u_mean, {
    req(input$null_u_mean)
    updateNumericInput(session, "null_b", value = 2 * input$null_u_mean)
  })

  observeEvent(input$play, {
    if (input$dist == "Uniform" && !isTRUE(input$true_u_mean > 0 && input$null_u_mean > 0)) {
      showNotification(
        "Uniform means must be positive: the minimum is 0, so the maximum is 2 \u00d7 mean.",
        type = "error"
      )
      return()
    }
    if (var_list$curr_sim < MAX_SIM) running(TRUE)
  })

  observeEvent(input$stop, {
    running(FALSE)
  })

  observeEvent(input$reset, {
    running(FALSE)
    reset_vars(var_list)
  })

  # One observer for the session; Play never creates another timer observer.
  observe({
    req(running())
    interval = isolate({
      forward(var_list, input)
      if (var_list$curr_sim >= MAX_SIM) {
        running(FALSE)
        NULL
      } else {
        if (input$speed == "Standard") 3000 else 300
      }
    })
    if (!is.null(interval)) invalidateLater(interval, session)
  })

  # data distribution
  output$data_dist = renderPlot({
    
    if (is.na(var_list$curr_sample[1])) {
      
      # don't render if no data
      return()
      
    } else {
      
      xmin = if (input$dist == "Bernoulli") {
        0
      } else if (input$dist == "Uniform") {
        0
      } else if (input$dist == "Normal") {
        input$true_m - 3*input$sd
      }
      
      xmax = if (input$dist == "Bernoulli") {
        1
      } else if (input$dist == "Uniform") {
        2 * input$true_u_mean
      } else if (input$dist == "Normal") {
        input$true_m + 3*input$sd
      }
      
      hist(
        var_list$curr_sample,
        xlim = c(xmin, xmax),
        main = paste0('Distribution of random sample #', var_list$curr_sim),
        xlab = 'Possible values of random variable'
      )
      
      # add a small arrow below the x axis pointing at the true mean (estimand)
      arrow_below_axis(var_list$true_mean, "red")
      
      # add a vertical line at the point estimate (sample mean)
      abline(v=var_list$curr_mean, col="darkgreen", lwd=3)
      
      # add vertical lines at confidence interval
      abline(v=var_list$curr_ci[1], col="black", lwd=1, lty=2)
      abline(v=var_list$curr_ci[2], col="black", lwd=1, lty=2)

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
      
      xlim = sampling_xlim(input)
      
      # true null sampling distribution
      plot(
        var_list$null_density,
        xlim = xlim,
        main = paste0('True null sampling distribution of estimator given n=', input$n),
        xlab = 'Possible values of estimator',
        col = "blue",
        lwd = 2
      )
      
      # add a small arrow below the x axis pointing at the null mean
      arrow_below_axis(var_list$null_mean, "blue")
      
      # very light dashed lines at the rejection region boundaries
      reject_dist = quantile(abs(var_list$null_estimates - var_list$null_mean), 1 - input$alpha)
      abline(
        v = var_list$null_mean + c(-1, 1) * reject_dist,
        col = adjustcolor("black", alpha.f = 0.25), lty = 2
      )
      
      p_value = NA
      if (!is.na(var_list$curr_mean)) {
        # distance of the point estimate from the null mean
        dist_from_null = abs(var_list$curr_mean - var_list$null_mean)
        
        # lightly shade the two-sided area under the null sampling distribution
        dens = var_list$null_density
        shade_col = adjustcolor("blue", alpha.f = 0.2)
        for (tail_idx in list(
          which(dens$x <= var_list$null_mean - dist_from_null),
          which(dens$x >= var_list$null_mean + dist_from_null)
        )) {
          if (length(tail_idx) > 1) {
            polygon(
              c(dens$x[tail_idx], rev(dens$x[tail_idx])),
              c(dens$y[tail_idx], rep(0, length(tail_idx))),
              col = shade_col, border = NA
            )
          }
        }
        
        # add a vertical line at the point estimate
        abline(v=var_list$curr_mean, col="darkgreen", lwd=3)
        
        # two-sided p-value from the simulated null sampling distribution
        p_value = mean(
          abs(var_list$null_estimates - var_list$null_mean) >= dist_from_null
        )
      }

      # write the sample mean and SD above the plot
      # (base graphics use plotmath, R's equivalent of LaTeX, for H_0)
      mtext(side=3, text=bquote(
        .(paste0(
          "SE of null sampling distribution: ",
          round(var_list$null_se, 3),
          "    two-sided p-value: ",
          round(p_value, 3),
          "    "
        )) * H[0] * .(paste0(" rejected? ", p_value < input$alpha))
      ))
    }
    
  })
  
  # null sampling distribution estimated from the current sample: N(null mean, s/sqrt(n)),
  # except for Bernoulli where it is identical to the true null sampling distribution
  output$est_sampling_dist = renderPlot({
    
    if (is.na(var_list$curr_sample[1])) {
      # don't render if no data
      return()
    }
    
    is_bernoulli = input$dist == "Bernoulli"
    est_se = if (is_bernoulli) {
      var_list$null_se
    } else {
      sd(var_list$curr_sample) / sqrt(input$n)
    }
    if (is.na(est_se) || est_se == 0) {
      plot.new()
      title("Estimated null sampling distribution: sample SD is 0 or undefined")
      return()
    }
    
    xlim = sampling_xlim(input)
    null_mean = var_list$null_mean
    if (is_bernoulli) {
      x = var_list$null_density$x
      y = var_list$null_density$y
    } else {
      x = seq(xlim[1], xlim[2], length.out = 1000)
      y = dnorm(x, mean = null_mean, sd = est_se)
    }
    
    plot(
      x, y,
      type = "l",
      xlim = xlim,
      main = paste0(
        if (is_bernoulli) 'Null sampling distribution' else 'Estimated null sampling distribution',
        ' of estimator given n=', input$n
      ),
      xlab = 'Possible values of estimator',
      ylab = 'Density',
      col = "blue",
      lwd = 2
    )
    
    # add a small arrow below the x axis pointing at the null mean
    arrow_below_axis(null_mean, "blue")
    
    # very light dashed lines at the rejection region boundaries
    reject_dist = if (is_bernoulli) {
      quantile(abs(var_list$null_estimates - null_mean), 1 - input$alpha)
    } else {
      qnorm(1 - input$alpha/2) * est_se
    }
    abline(
      v = null_mean + c(-1, 1) * reject_dist,
      col = adjustcolor("black", alpha.f = 0.25), lty = 2
    )
    
    # lightly shade the two-sided area under the estimated null distribution
    dist_from_null = abs(var_list$curr_mean - null_mean)
    shade_col = adjustcolor("blue", alpha.f = 0.2)
    for (tail_idx in list(
      which(x <= null_mean - dist_from_null),
      which(x >= null_mean + dist_from_null)
    )) {
      if (length(tail_idx) > 1) {
        polygon(
          c(x[tail_idx], rev(x[tail_idx])),
          c(y[tail_idx], rep(0, length(tail_idx))),
          col = shade_col, border = NA
        )
      }
    }
    
    # add a vertical line at the point estimate
    abline(v=var_list$curr_mean, col="darkgreen", lwd=3)
    
    est_p_value = if (is_bernoulli) {
      mean(abs(var_list$null_estimates - null_mean) >= dist_from_null)
    } else {
      2 * pnorm(-dist_from_null / est_se)
    }
    
    # plotmath is used so the square root is drawn with a radical sign
    se_value = round(
      if (is_bernoulli) sqrt(input$null_p * (1 - input$null_p) / input$n) else est_se,
      3
    )
    se_label = if (is_bernoulli) {
      bquote("SE: "*sqrt(.(input$null_p)*(1 - .(input$null_p))/.(input$n))*" = "*.(se_value))
    } else {
      bquote("Estimated SE ("*s/sqrt(n)*"): "*.(se_value))
    }
    
    mtext(side=3, text=bquote(
      .(se_label) * .(paste0(
        "    two-sided p-value: ",
        round(est_p_value, 3),
        "    "
      )) * H[0] * .(paste0(" rejected? ", est_p_value < input$alpha))
    ))
    
  })
  
  output$conf_level = renderText({
    req(input$alpha)
    1 - input$alpha
  })
  
  output$running_summary = renderUI({
    req(var_list$curr_sim >= 1, !is.na(var_list$coverage[1]))
    i = var_list$curr_sim
    tagList(
      paste0("Proportion of ", i, " random samples with:"),
      tags$ul(
        if (input$dist == "Bernoulli") {
          tags$li(paste0(
            "A confidence interval that captures the true fixed mean: ",
            round(var_list$coverage[i], 3)
          ))
        } else {
          tagList(
            tags$li(paste0(
              "A confidence interval that SHOULD capture the true fixed mean: ",
              round(var_list$coverage[i], 3)
            )),
            tags$li(paste0(
              "A confidence interval that ACTUALLY captures the true fixed mean: ",
              round(var_list$coverage_est[i], 3)
            ))
          )
        },
        if (input$dist == "Bernoulli") {
          tags$li(paste0(
            "A hypothesis test that rejects the null hypothesis: ",
            round(var_list$reject_est_prop[i], 3)
          ))
        } else {
          tagList(
            tags$li(paste0(
              "A hypothesis test that SHOULD reject the null hypothesis: ",
              round(var_list$reject_true_prop[i], 3)
            )),
            tags$li(paste0(
              "A hypothesis test that ACTUALLY rejects the null hypothesis: ",
              round(var_list$reject_est_prop[i], 3)
            ))
          )
        }
      )
    )
  })
  
}

runApp(shinyApp(ui,server),launch.browser = TRUE)