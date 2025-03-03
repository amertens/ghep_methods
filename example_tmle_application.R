# Load required packages
library(tidyverse)
library(survival)
library(tmle3)
library(SuperLearner)
library(survminer)
library(gridExtra)

#' 1. Data Simulation ----
simulate_hcv_data <- function(n = 1000, seed = 123) {
  set.seed(seed)
  
  # Baseline covariates
  data <- tibble(
    # Demographics
    age = rnorm(n, mean = 50, sd = 10),
    sex = rbinom(n, 1, 0.6),
    
    # Comorbidities
    diabetes = rbinom(n, 1, 0.2),
    hypertension = rbinom(n, 1, 0.3),
    liver_disease = rbinom(n, 1, 0.4),
    
    # Baseline kidney function
    baseline_egfr = rnorm(n, mean = 90, sd = 15),
    
    # Create treatment assignment based on covariates
    ps_true = plogis(-1 + 0.02 * age + 0.5 * diabetes + 
                       0.3 * hypertension - 0.01 * baseline_egfr)
  )
  
  # Assign treatment (SOF vs non-SOF)
  data$treatment <- rbinom(n, 1, data$ps_true)
  
  # Simulate survival times
  lambda <- 0.001 * exp(
    0.3 * data$treatment +
      0.02 * scale(data$age) +
      0.5 * data$diabetes +
      -0.01 * scale(data$baseline_egfr)
  )
  
  # Generate survival times and censoring
  data <- data %>%
    mutate(
      surv_time = rexp(n, lambda),
      cens_time = rexp(n, 0.0005),
      time = pmin(surv_time, cens_time),
      event = as.numeric(surv_time <= cens_time)
    )
  
  return(data)
}

#' 2. Enhanced Diagnostics ----

#' 2.1 Cox Proportional Hazards Assumptions
check_ph_assumptions <- function(fit, data) {
  # Test proportional hazards
  ph_test <- cox.zph(fit)
  
  # Create diagnostic plots
  ph_plots <- list()
  
  # Schoenfeld residuals plot
  ph_plots$schoenfeld <- ggcoxzph(ph_test)
  
  # Log-log plots
  ph_plots$loglog <- ggsurvplot(
    survfit(Surv(time, event) ~ treatment, data = data),
    fun = "cloglog",
    conf.int = TRUE,
    title = "Log-log plot"
  )
  
  return(list(test = ph_test, plots = ph_plots))
}

#' 2.2 Extended Residual Plots
create_residual_plots <- function(fit, data) {
  # Calculate different types of residuals
  res_data <- data %>%
    mutate(
      mart_res = residuals(fit, type = "martingale"),
      dev_res = residuals(fit, type = "deviance"),
      score_res = residuals(fit, type = "score")[,1], # First column for treatment
      schoen_res = residuals(fit, type = "schoenfeld")[,1] # First column for treatment
    )
  
  # Create plots
  plots <- list()
  
  # Martingale residuals vs. linear predictor
  plots$martingale <- ggplot(res_data, aes(x = predict(fit), y = mart_res)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess") +
    theme_minimal() +
    labs(title = "Martingale Residuals vs. Linear Predictor",
         x = "Linear Predictor",
         y = "Martingale Residuals")
  
  # Deviance residuals vs. linear predictor
  plots$deviance <- ggplot(res_data, aes(x = predict(fit), y = dev_res)) +
    geom_point(alpha = 0.5) +
    geom_smooth(method = "loess") +
    theme_minimal() +
    labs(title = "Deviance Residuals vs. Linear Predictor",
         x = "Linear Predictor",
         y = "Deviance Residuals")
  
  # Score residuals vs time
  plots$score <- ggplot(res_data, aes(x = time, y = score_res)) +
    geom_point(alpha = 0.5) +
    theme_minimal() +
    labs(title = "Score Residuals vs. Time",
         x = "Time",
         y = "Score Residuals")
  
  return(plots)
}
check_tmle_convergence <- function(tmle_fit) {
  # Extract key convergence metrics
  metrics <- list(
    # Initial vs updated estimates
    initial_vs_updated = data.frame(
      stage = c("Initial", "Updated"),
      estimate = c(tmle_fit$initial_estimate, tmle_fit$estimates$ATE),
      se = c(tmle_fit$initial_se, tmle_fit$estimates$se)
    ),
    
    # EIC mean (should be approximately zero)
    eic_mean = mean(tmle_fit$eic),
    
    # EIC standard error
    eic_se = sd(tmle_fit$eic)/sqrt(length(tmle_fit$eic)),
    
    # Convergence status
    converged = tmle_fit$converged,
    
    # Number of iterations
    n_iterations = tmle_fit$n_iterations,
    
    # Variance of influence curve
    ic_variance = var(tmle_fit$eic),
    
    # Relative efficiency vs initial
    rel_efficiency = tmle_fit$estimates$se / tmle_fit$initial_se
  )
  
  # Create diagnostic plots
  plots <- list()
  
  # EIC distribution
  plots$eic_dist <- ggplot(data.frame(eic = tmle_fit$eic), aes(x = eic)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Distribution of Efficient Influence Curve",
         x = "EIC Values",
         y = "Count")
  
  # Convergence plot (if iteration history available)
  if (!is.null(tmle_fit$iteration_history)) {
    plots$convergence <- ggplot(
      data.frame(
        iteration = seq_along(tmle_fit$iteration_history),
        estimate = tmle_fit$iteration_history
      ),
      aes(x = iteration, y = estimate)
    ) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      labs(title = "TMLE Convergence Path",
           x = "Iteration",
           y = "Estimate")
  }
  
  # Q-function comparison plot
  plots$q_functions <- ggplot(
    data.frame(
      initial = tmle_fit$initial_preds,
      updated = tmle_fit$targeted_preds
    ),
    aes(x = initial, y = updated)
  ) +
    geom_point(alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Initial vs Targeted Q-function Predictions",
         x = "Initial Predictions",
         y = "Updated Predictions")
  
  # Clever covariate distribution
  plots$clever_cov <- ggplot(
    data.frame(clever_cov = tmle_fit$clever_covariates),
    aes(x = clever_cov)
  ) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Distribution of Clever Covariates",
         x = "Clever Covariate Values",
         y = "Count")
  
  # Compute additional diagnostics
  diagnostics <- list(
    # Check positivity
    positivity_violations = check_positivity(tmle_fit),
    
    # Assess balance improvement
    balance_metrics = compute_balance_metrics(tmle_fit),
    
    # Numerical stability checks
    numerical_stability = check_numerical_stability(tmle_fit)
  )
  
  return(list(
    metrics = metrics,
    plots = plots,
    diagnostics = diagnostics
  ))
}

# Helper functions for additional diagnostics
check_positivity <- function(tmle_fit) {
  ps_bounds <- range(tmle_fit$ps_scores)
  near_bounds <- sum(tmle_fit$ps_scores < 0.05 | tmle_fit$ps_scores > 0.95)
  
  list(
    ps_range = ps_bounds,
    n_near_bounds = near_bounds,
    prop_near_bounds = near_bounds / length(tmle_fit$ps_scores)
  )
}

compute_balance_metrics <- function(tmle_fit) {
  # Compare standardized differences before/after targeting
  initial_balance <- compute_std_diff(tmle_fit$initial_preds, tmle_fit$treatment)
  final_balance <- compute_std_diff(tmle_fit$targeted_preds, tmle_fit$treatment)
  
  list(
    initial_max_std_diff = max(abs(initial_balance)),
    final_max_std_diff = max(abs(final_balance)),
    balance_improvement = mean(abs(initial_balance)) - mean(abs(final_balance))
  )
}

check_numerical_stability <- function(tmle_fit) {
  list(
    # Check for extreme weights
    max_weight = max(abs(tmle_fit$weights)),
    prop_extreme_weights = mean(abs(tmle_fit$weights) > 10),
    
    # Check for numerical instability in targeting step
    targeting_stability = var(diff(tmle_fit$iteration_history)),
    
    # Check condition number of design matrix if available
    condition_number = if(!is.null(tmle_fit$design_matrix)) {
      kappa(tmle_fit$design_matrix)
    } else {
      NA
    }
  )
}

# Helper function to compute standardized differences
compute_std_diff <- function(x, treatment) {
  treated_mean <- mean(x[treatment == 1])
  control_mean <- mean(x[treatment == 0])
  pooled_sd <- sqrt((var(x[treatment == 1]) + var(x[treatment == 0])) / 2)
  
  (treated_mean - control_mean) / pooled_sd
}


# Load required packages
library(tidyverse)
library(survival)
library(tmle)
library(SuperLearner)
library(gridExtra)

#' Run diagnostics and interpret results
interpret_tmle_analysis <- function() {
  # 1. Generate example data
  set.seed(123)
  data <- simulate_hcv_data(n = 1000)
  
  # 2. Fit TMLE for survival outcome
  tmle_fit <- tmle(
    Y = data$time,
    A = data$treatment,
    W = data[, c("age", "sex", "diabetes", "hypertension", 
                 "liver_disease", "baseline_egfr")],
    Delta = data$event,  # Censoring indicator
    g.SL.library = c("SL.glm", "SL.mean"),
    Q.SL.library = c("SL.glm", "SL.mean"),
    family = "survival"
  )
  
  # 3. Run diagnostics
  diagnostics <- check_tmle_convergence(tmle_fit)
  
  # 4. Interpret results
  
  # 4.1 Basic convergence checks
  cat("\n=== Basic Convergence Checks ===\n")
  cat("Number of updates:", length(tmle_fit$epsilon), "\n")
  
  # Interpretation function for convergence
  interpret_convergence <- function(tmle_fit) {
    n_updates <- length(tmle_fit$epsilon)
    if(n_updates < 5) {
      return("✓ Quick convergence achieved")
    } else if(n_updates < 10) {
      return("✓ Normal convergence achieved")
    } else {
      return("⚠ Many updates needed - check for stability issues")
    }
  }
  
  cat("Interpretation:", interpret_convergence(tmle_fit), "\n")
  
  # 4.2 EIC Analysis
  cat("\n=== EIC Analysis ===\n")
  eic_mean <- mean(tmle_fit$IC)
  eic_se <- sd(tmle_fit$IC)/sqrt(length(tmle_fit$IC))
  cat("EIC mean:", round(eic_mean, 6), "\n")
  cat("EIC SE:", round(eic_se, 6), "\n")
  
  # Interpretation function for EIC
  interpret_eic <- function(eic_mean, eic_se) {
    eic_standardized <- abs(eic_mean / eic_se)
    if(eic_standardized < 0.1) {
      return("✓ EIC centered well around zero")
    } else if(eic_standardized < 0.5) {
      return("⚠ EIC slightly off-center - check for systematic bias")
    } else {
      return("❌ EIC significantly off-center - results may be biased")
    }
  }
  
  cat("Interpretation:", interpret_eic(eic_mean, eic_se), "\n")
  
  # 4.3 Treatment Mechanism Assessment
  cat("\n=== Treatment Mechanism Assessment ===\n")
  g1W <- tmle_fit$g$g1W
  g0W <- 1 - g1W
  near_bounds <- sum(g1W < 0.05 | g1W > 0.95)
  prop_near_bounds <- near_bounds / length(g1W)
  
  cat("Proportion near bounds:", round(prop_near_bounds, 3), "\n")
  
  # Interpretation function for treatment mechanism
  interpret_tx_mechanism <- function(prop_near_bounds) {
    if(prop_near_bounds < 0.01) {
      return("✓ Good treatment mechanism estimation")
    } else if(prop_near_bounds < 0.05) {
      return("⚠ Some extreme values - results may be less stable")
    } else {
      return("❌ Many extreme values - results may be unreliable")
    }
  }
  
  cat("Interpretation:", interpret_tx_mechanism(prop_near_bounds), "\n")
  
  # 4.4 Estimate Assessment
  cat("\n=== Estimate Assessment ===\n")
  cat("Initial estimate:", round(tmle_fit$Qinit$Q, 3), "\n")
  cat("Updated estimate:", round(tmle_fit$estimates$ATE$psi, 3), "\n")
  cat("Standard error:", round(tmle_fit$estimates$ATE$var.psi^0.5, 3), "\n")
  
  # Create diagnostic plots
  plots <- list()
  
  # EIC Distribution
  plots$eic <- ggplot(data.frame(ic = tmle_fit$IC), aes(x = ic)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = "Influence Curve Distribution",
         subtitle = paste("Mean:", round(mean(tmle_fit$IC), 6)))
  
  # Treatment mechanism
  plots$tx_mechanism <- ggplot(data.frame(g1W = g1W), aes(x = g1W)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Treatment Mechanism Distribution",
         subtitle = "Check for values near 0 or 1")
  
  # Clever covariate plot
  plots$clever_cov <- ggplot(data.frame(H = tmle_fit$H), aes(x = H)) +
    geom_histogram(bins = 30, fill = "blue", alpha = 0.5) +
    theme_minimal() +
    labs(title = "Clever Covariate Distribution")
  
  # Display plots in grid
  do.call(grid.arrange, c(plots, ncol = 2))
  
  # Overall Assessment
  cat("\n=== Overall Assessment ===\n")
  quality_score <- mean(c(
    length(tmle_fit$epsilon) < 10,  # Convergence
    abs(eic_mean/eic_se) < 0.5,     # EIC
    prop_near_bounds < 0.05         # Treatment mechanism
  ))
  
  cat("Overall quality score:", round(quality_score * 100), "%\n")
  
  quality_interpretation <- case_when(
    quality_score > 0.8 ~ "✓ High quality results - proceed with confidence",
    quality_score > 0.6 ~ "⚠ Moderate quality - interpret with some caution",
    TRUE ~ "❌ Low quality - results may not be reliable"
  )
  
  cat("Final interpretation:", quality_interpretation, "\n")
  
  return(list(
    tmle_fit = tmle_fit,
    diagnostics = list(
      eic_mean = eic_mean,
      eic_se = eic_se,
      prop_near_bounds = prop_near_bounds,
      quality_score = quality_score
    ),
    plots = plots
  ))
}

# Example usage
if(TRUE) {  # Set to TRUE to run
  results <- interpret_tmle_analysis()
  
  # Access components
  print(results$tmle_fit$estimates$ATE)  # Final estimates
  print(results$diagnostics)             # Diagnostic metrics
  print(results$plots$eic)               # EIC distribution plot
}