# fit_indgns.r
# Diagonal-constrained multivariate MGNS estimation (dddd model)
# Estimate y1, y2 separately with sacsarlm()
# - R matrix: diagonal (ρ₁₁ only, ρ₂₂ only)
# - Λ matrix: diagonal (λ₁₁ only, λ₂₂ only)  
# - A matrix: diagonal (time lag also independent, no cross terms)
# - Σ matrix: diagonal (errors are independent)

#' Diagonal-constrained multivariate MGNS estimation (dddd model)
#' 
#' Estimate each variable separately with sacsarlm() and merge the results
#' 
#' @param data_file path to the data file
#' @param weight_file path to the spatial-weight-matrix file
#' @param y_vars vector of dependent variable names
#' @param x_vars list of regressor names (per dependent variable)
#' @param time_var time variable name
#' @param time_point time point to use (latest if NULL)
#' @param region_var region variable name
#' @param include_intercept whether to include an intercept
#' @param include_time_lag whether to include a time lag
#' @param verbose verbose output
#' @return multivar_indgns object
#' 
fit_indgns <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  data = NULL,
  W_matrix = NULL,
  verbose = FALSE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained MGNS estimation (dddd model) ===\n")
    cat(sprintf("    Time lag: %s (no cross terms)\n", ifelse(include_time_lag, "included", "not included")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("the spatialreg package is required: install.packages('spatialreg')")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("the spdep package is required: install.packages('spdep')")
  }
  
  k <- length(y_vars)
  
  # S0: data preparation
  
  if (verbose) cat("\nS0: data preparation...\n")
  
  if (is.null(data)) data <- read.csv(data_file)
  if (is.null(W_matrix)) W_matrix <- as.matrix(read.csv(weight_file, header = TRUE))
  
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  if (include_time_lag) {
    data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
    if (nrow(data_t_lag) != n) {
      stop("the number of regions in the time-lag data does not match")
    }
  }
  
  W_listw <- spdep::mat2listw(W_matrix, style = "W")
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Time point used: t = %d\n", time_point))
    cat(sprintf("  Time lag: %s\n", ifelse(include_time_lag, "included (no cross terms)", "not included")))
  }
  
  # S1-S2: estimate each variable separately with sacsarlm()
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Per-variable MGNS estimation (sacsarlm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  R_diag <- matrix(0, k, k)
  Lambda_diag <- matrix(0, k, k)
  Sigma_diag <- matrix(0, k, k)
  
  individual_models <- list()
  beta_list <- list()
  residuals_list <- list()
  fitted_list <- list()
  
  total_loglik <- 0
  total_params <- 0
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    if (verbose) {
      cat(sprintf("\n--- Estimating variable %s ---\n", var_name))
    }
    
    df_i <- data.frame(
      y = data_t[[var_name]]
    )
    
    for (xvar in xi_vars) {
      df_i[[xvar]] <- data_t[[xvar]]
    }
    
    # Add time-lag variables (own lag only, no cross terms)
    if (include_time_lag) {
      lag_var_name <- paste0(var_name, "_lag")
      df_i[[lag_var_name]] <- data_t_lag[[var_name]]
    }
    
    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)
    
    tryCatch({
      model_i <- spatialreg::sacsarlm(
        formula = formula_i,
        data = df_i,
        listw = W_listw,
        method = "eigen",
        quiet = !verbose
      )
      
      R_diag[i, i] <- model_i$rho
      Lambda_diag[i, i] <- model_i$lambda
      Sigma_diag[i, i] <- model_i$s2
      
      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef
      
      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))
      
      s <- summary(model_i)
      
      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  #  add σ² only (the spatial parameter is included in coef)
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i
      
      rho <- model_i$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      
      lambda <- model_i$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      y_i <- df_i$y
      fitted_i <- fitted(model_i)
      residuals_i <- residuals(model_i)
      
      r2_results <- compute_r_squared(
        y = y_i,
        fitted = as.numeric(fitted_i),
        residuals = as.numeric(residuals_i),
        loglik = loglik_i,
        num_params = num_params_i
      )
      
      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = get_signif_code(rho_p),
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = get_signif_code(lambda_p)
        ),
        
        coefficients = coef_table,
        sigma2 = model_i$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params_i,
          num_obs = length(y_i)
        )
      )
      
      if (verbose) {
        cat("  [Spatial parameters]\n")
        cat(sprintf("    ρ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    rho, rho_se, rho_z, rho_p, get_signif_code(rho_p)))
        cat(sprintf("    λ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    lambda, lambda_se, lambda_z, lambda_p, get_signif_code(lambda_p)))
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f\n", model_i$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²     = %.4f\n", r2_results$R2_adj))
        cat(sprintf("    Log-likelihood = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: sacsarlm() failed when estimating variable %s\n", var_name))
      cat("Error message:", e$message, "\n")
      stop("Estimation failed")
    })
  }
  
  # S3: merge the results
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Merging the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  ic <- compute_information_criteria(total_loglik, total_params, k * n)
  
  # Vectorize residuals and fitted values (in y1, y2 order)
  residuals_raw <- unlist(residuals_list)
  fitted_vals <- unlist(fitted_list)
  
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_list[[y_vars[i]]] / sqrt(Sigma_diag[i, i])
  }
  
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- data_t[[y_vars[i]]]
  }
  
  r2_list <- sapply(individual_models, function(m) m$fit$R2)
  r2_adj_list <- sapply(individual_models, function(m) m$fit$R2_adj)
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  coefficients <- list(
    R = R_diag,
    Lambda = Lambda_diag,
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  R_se <- matrix(0, k, k)
  Lambda_se <- matrix(0, k, k)
  for (i in 1:k) {
    var_name <- y_vars[i]
    R_se[i, i] <- individual_models[[var_name]]$spatial_params$rho_se
    Lambda_se[i, i] <- individual_models[[var_name]]$spatial_params$lambda_se
  }
  
  std_errors <- list(
    R = R_se,
    Lambda = Lambda_se
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  # Output results
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[Diagonal R matrix (spatial lag)]\n")
    print(round(R_diag, 4))
    
    cat("\n[Diagonal Λ matrix (spatial error)]\n")
    print(round(Lambda_diag, 4))
    
    cat("\n[Diagonal Σ matrix (error variances)]\n")
    print(round(Sigma_diag, 4))
    
    cat("\n[Goodness of fit (merged)]\n")
    cat(sprintf("  Total log-likelihood: %.4f\n", total_loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  Mean R²: %.4f\n", r2_mean))
    cat(sprintf("  Mean Adj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Build the result object (unified output format)
  
  dl <- build_data_list_from_parts(
    y = y_vec, W = W_matrix, W_listw = W_listw,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = include_time_lag,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "IndGNS",
    R                 = R_diag,
    Lambda_mat             = Lambda_diag,
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "sacsarlm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    std_errors_R      = R_se,
    std_errors_Lambda      = Lambda_se,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}

# Standard methods

#' print method
print.multivar_indgns <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained MGNS (dddd model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time lag: %s\n", ifelse(x$data_info$include_time_lag, "included (no cross terms)", "not included")))
  
  cat("\n[Diagonal R matrix (spatial lag ρ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$R)), collapse=", ")))
  
  cat("\n[Diagonal Λ matrix (spatial error λ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$Lambda)), collapse=", ")))
  
  cat("\n[Goodness of fit]\n")
  cat(sprintf("  Log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  Mean R²: %.4f\n", x$fit$R2))
  cat(sprintf("  Mean Adj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for details\n\n")
  
  invisible(x)
}

#' summary method
summary.multivar_indgns <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained MGNS - detailed results (dddd model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Dependent variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  cat(sprintf("  Time lag: %s\n", ifelse(x$data_info$include_time_lag, "included (no cross terms)", "not included")))
  
  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. Estimation results for %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    cat("\nSpatial parameters:\n")
    sp <- model_i$spatial_params
    cat(sprintf("  ρ (spatial lag) = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$rho, sp$rho_se, sp$rho_z, sp$rho_p, sp$rho_signif))
    cat(sprintf("  λ (spatial error) = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$lambda, sp$lambda_se, sp$lambda_z, sp$lambda_p, sp$lambda_signif))
    
    cat(sprintf("\nError variance: σ² = %.6f\n", model_i$sigma2))
    
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  Log-likelihood = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))
  }
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. Merged goodness of fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\nDiagonal R matrix:\n")
  print(round(x$coefficients$R, digits))
  
  cat("\nDiagonal Λ matrix:\n")
  print(round(x$coefficients$Lambda, digits))
  
  cat("\nDiagonal Σ matrix:\n")
  print(round(x$coefficients$Sigma, digits))
  
  cat("\nMerged goodness of fit:\n")
  cat(sprintf("  Total log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  Mean R²       = %.4f\n", x$fit$R2))
  cat(sprintf("  Mean Adj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}

#' coef method
coef.multivar_indgns <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(object$coefficients$R)
  } else if (type == "Lambda") {
    return(object$coefficients$Lambda)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}

#' logLik method
logLik.multivar_indgns <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}

#' AIC method
AIC.multivar_indgns <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}

#' BIC method
BIC.multivar_indgns <- function(object, ...) {
  return(object$fit$BIC)
}

#' residuals method
residuals.multivar_indgns <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}

#' fitted method
fitted.multivar_indgns <- function(object, ...) {
  
  fitted_vals <- numeric(object$data_info$k * object$data_info$n)
  
  for (i in 1:object$data_info$k) {
    var_name <- object$data_info$y_vars[i]
    idx <- ((i - 1) * object$data_info$n + 1):(i * object$data_info$n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}
