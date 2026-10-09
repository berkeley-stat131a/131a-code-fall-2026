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
  titlePanel("STAT 131A Multiple testing simulation"),
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
      
      uiOutput(outputId='running_ci_summary'),
      plotOutput(outputId='data_dist', height='200px'),
      uiOutput(outputId='running_test_summary'),
      plotOutput(outputId='sampling_dist', height='200px'),
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

# exact null distribution of a Bernoulli sample proportion: each possible value k/n,
# its probability under H0, and its exact two-sided p-value
bernoulli_null_dist = function(p0, n) {
  support = (0:n) / n
  prob = dbinom(0:n, size = n, prob = p0)
  dev = abs(support - p0)
  # tolerance guards against floating point error when comparing distances
  p_value = vapply(dev, function(d) sum(prob[dev >= d - 1e-9]), numeric(1))
  list(support = support, prob = prob, p_value = p_value)
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
  
  # running proportion of rejected tests, using the true null distribution
  # and the estimated (s / sqrt(n)) null distribution
  var_list$reject_true_prop = NA
  var_list$reject_est_prop = NA
  
  # theoretical SE of the null sampling distribution
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
    # (Bernoulli means sit on a 1/n grid, so they need extra smoothing; the
    # continuous cases use the default bandwidth so the curve's width matches the SE)
    density_range = if (input$dist == "Bernoulli") {
      list(from = 0, to = 1, adjust = 5)
    } else {
      list(adjust = 1)
    }
    var_list$real_density = do.call(
      density, c(list(var_list$real_estimates), density_range)
    )
    
    # theoretical mean and SE of the sample mean under the null hypothesis
    if (input$dist=="Bernoulli") {
      var_list$null_mean = input$null_p
      var_list$null_se = sqrt(input$null_p * (1 - input$null_p) / input$n)
    } else if (input$dist=="Uniform") {
      # Uniform(0, b) has SD b / sqrt(12), with b = 2 * mean
      var_list$null_mean = input$null_u_mean
      var_list$null_se = 2 * input$null_u_mean / sqrt(12 * input$n)
    } else if (input$dist=="Normal") {
      var_list$null_mean = input$null_m
      var_list$null_se = input$sd / sqrt(input$n)
    }
    
    # run the hypothesis test on every simulated sample
    dev_from_null = abs(var_list$means - var_list$null_mean)
    
    # p-value from the theoretical null distribution N(null mean, true null SE)
    p_true = 2 * pnorm(-dev_from_null / var_list$null_se)
    reject_true = p_true < input$alpha
    reject_true[is.na(reject_true)] = FALSE
    
    # p-value from N(null mean, s / sqrt(n)), with s from each sample
    sample_sds = sqrt(pmax(
      colSums(var_list$samples_to_iter^2) - input$n * var_list$means^2, 0
    ) / (input$n - 1))
    p_est = 2 * pnorm(-dev_from_null / (sample_sds / sqrt(input$n)))
    reject_est = p_est < input$alpha
    reject_est[is.na(reject_est)] = FALSE
    if (input$dist == "Bernoulli") {
      # the estimated test uses the normal approximation N(p_H0, sqrt(p_H0 (1 - p_H0) / n)),
      # which is what reject_true holds so far; the true test uses the exact Binomial
      reject_est = reject_true
      exact = bernoulli_null_dist(input$null_p, input$n)
      reject_true = exact$p_value[round(var_list$means * input$n) + 1] < input$alpha
    }
    
    var_list$reject_true_prop = cumsum(reject_true) / seq_along(reject_true)
    var_list$reject_est_prop = cumsum(reject_est) / seq_along(reject_est)
    
    alpha = input$alpha
    margin = qnorm(1 - alpha/2) * var_list$se
    
    # also build the CI with the SE estimated from each sample:
    # sqrt(p_hat (1 - p_hat) / n) for Bernoulli, s / sqrt(n) otherwise
    var_list$est_ses = if (input$dist == "Bernoulli") {
      sqrt(var_list$means * (1 - var_list$means) / input$n)
    } else {
      sample_sds / sqrt(input$n)
    }
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
  # only show the CI constructed with the (incorrectly) estimated SE
  var_list$curr_ci = c(
    var_list$ci_est_lower[var_list$curr_sim],
    var_list$ci_est_upper[var_list$curr_sim]
  )
  var_list$contains_true_mean = var_list$contains_true_mean_est_vec[var_list$curr_sim]
  
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

  # samples are generated once per run, so changing any data-generation or testing
  # input stops and resets the run rather than mixing old samples with new settings
  observeEvent(
    list(
      input$dist, input$true_p, input$true_u_mean, input$true_m, input$sd, input$n,
      input$null_p, input$null_u_mean, input$null_m, input$alpha
    ),
    {
      running(FALSE)
      reset_vars(var_list)
    },
    ignoreInit = TRUE
  )

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
      null_mean = var_list$null_mean
      main = paste0('True null sampling distribution of estimator given n=', input$n)
      
      # distance of the point estimate from the null mean
      dist_from_null = abs(var_list$curr_mean - null_mean)
      shade_col = adjustcolor("blue", alpha.f = 0.2)
      boundary_col = adjustcolor("black", alpha.f = 0.25)
      
      if (input$dist == "Bernoulli") {
        # exact Binomial distribution: one bar for each possible sample proportion
        exact = bernoulli_null_dist(input$null_p, input$n)
        in_tail = abs(exact$support - null_mean) >= dist_from_null - 1e-9
        bar_half_width = 0.4 / input$n
        
        plot(
          NA,
          xlim = xlim,
          ylim = c(0, max(exact$prob)),
          main = main,
          xlab = 'Possible values of estimator',
          ylab = 'Probability'
        )
        # lightly shade the bars at least as far from the null mean as the estimate
        rect(
          exact$support - bar_half_width, 0,
          exact$support + bar_half_width, exact$prob,
          col = ifelse(in_tail, shade_col, NA),
          border = "blue"
        )
        
        # very light dashed lines between the last kept and first rejected bar on each side
        rejected = exact$p_value < input$alpha
        for (side in c(-1, 1)) {
          on_side = which(sign(exact$support - null_mean) == side & rejected)
          if (length(on_side) > 0) {
            innermost = on_side[which.min(abs(exact$support[on_side] - null_mean))]
            abline(v = exact$support[innermost] - side * 0.5 / input$n, col = boundary_col, lty = 2)
          }
        }
        
        # exact two-sided p-value
        p_value = sum(exact$prob[in_tail])
      } else {
        # theoretical N(null mean, true null SE)
        # (exact for Normal data; the CLT approximation for Uniform data)
        x = seq(xlim[1], xlim[2], length.out = 1000)
        y = dnorm(x, mean = null_mean, sd = var_list$null_se)
        plot(
          x, y,
          type = "l",
          xlim = xlim,
          main = main,
          xlab = 'Possible values of estimator',
          ylab = 'Density',
          col = "blue",
          lwd = 2
        )
        
        # very light dashed lines at the rejection region boundaries
        reject_dist = qnorm(1 - input$alpha/2) * var_list$null_se
        abline(v = null_mean + c(-1, 1) * reject_dist, col = boundary_col, lty = 2)
        
        # lightly shade the two-sided area under the null sampling distribution
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
        
        # two-sided p-value from the theoretical null sampling distribution
        p_value = 2 * pnorm(-dist_from_null / var_list$null_se)
      }
      
      # add a small arrow below the x axis pointing at the null mean
      arrow_below_axis(null_mean, "blue")
      
      # add a vertical line at the point estimate
      abline(v=var_list$curr_mean, col="darkgreen", lwd=3)

      # write the SE and test result above the plot
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
  
  # null sampling distribution estimated with the normal approximation:
  # N(p_H0, sqrt(p_H0 (1 - p_H0) / n)) for Bernoulli, N(null mean, s / sqrt(n)) otherwise
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
    main = paste0('Normally-approximated null sampling distribution of estimator given n=', input$n)
    dist_from_null = abs(var_list$curr_mean - null_mean)
    shade_col = adjustcolor("blue", alpha.f = 0.2)
    boundary_col = adjustcolor("black", alpha.f = 0.25)
    
    x = seq(xlim[1], xlim[2], length.out = 1000)
    y = dnorm(x, mean = null_mean, sd = est_se)
    
    plot(
      x, y,
      type = "l",
      xlim = xlim,
      main = main,
      xlab = 'Possible values of estimator',
      ylab = 'Density',
      col = "blue",
      lwd = 2
    )
    
    # very light dashed lines at the rejection region boundaries
    reject_dist = qnorm(1 - input$alpha/2) * est_se
    abline(v = null_mean + c(-1, 1) * reject_dist, col = boundary_col, lty = 2)
    
    # lightly shade the two-sided area under the estimated null distribution
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
    
    est_p_value = 2 * pnorm(-dist_from_null / est_se)
    
    # add a small arrow below the x axis pointing at the null mean
    arrow_below_axis(null_mean, "blue")
    
    # add a vertical line at the point estimate
    abline(v=var_list$curr_mean, col="darkgreen", lwd=3)
    
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
  
  # running confidence interval results, above the plot of the current sample
  output$running_ci_summary = renderUI({
    req(var_list$curr_sim >= 1, !is.na(var_list$coverage[1]))
    i = var_list$curr_sim
    tagList(
      paste0("Proportion of ", i, " random samples with:"),
      tags$ul(
        tags$li(paste0(
          "A confidence interval that SHOULD capture the true fixed mean: ",
          round(var_list$coverage[i], 3)
        )),
        tags$li(paste0(
          "A confidence interval that ACTUALLY captures the true fixed mean: ",
          round(var_list$coverage_est[i], 3)
        ))
      )
    )
  })
  
  # running hypothesis test results, above the null sampling distribution plots
  output$running_test_summary = renderUI({
    req(var_list$curr_sim >= 1, !is.na(var_list$reject_true_prop[1]))
    i = var_list$curr_sim
    tagList(
      paste0("Proportion of ", i, " random samples with:"),
      tags$ul(
        tags$li(paste0(
          "A hypothesis test that SHOULD reject the null hypothesis: ",
          round(var_list$reject_true_prop[i], 3)
        )),
        tags$li(paste0(
          "A hypothesis test that ACTUALLY rejects the null hypothesis: ",
          round(var_list$reject_est_prop[i], 3)
        ))
      )
    )
  })
  
}

runApp(shinyApp(ui,server),launch.browser = TRUE)