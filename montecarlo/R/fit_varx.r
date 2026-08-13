# fit_varx.r
# VARX model estimation (0011 model)
# No spatial correlation, with time-lag cross terms, with error correlation
# Uses SUR estimation from the systemfit package
# Model equations:
#   y1_t = α11 * y1_{t-1} + α12 * y2_{t-1} + X1 * β1 + ε1
#   y2_t = α21 * y1_{t-1} + α22 * y2_{t-1} + X2 * β2 + ε2
#   
#   ε ~ N(0, Σ)  where Σ is non-diagonal (correlated errors)
# Parameters:
#   - R matrix: 0 (no spatial lag)
#   - Λ matrix: 0 (no spatial error)
#   - A matrix: non-diagonal (α11, α12, α21, α22)
#   - Σ matrix: non-diagonal (correlated errors)

#' VARX model estimation (0011 model)
#' 
#' Uses SUR estimation from the systemfit package
#' 
#' @param data_file path to the data file
#' @param y_vars vector of dependent variable names
#' @param x_vars list of regressor names (per dependent variable)
#' @param time_var time variable name
#' @param time_point time point to use (latest if NULL)
#' @param region_var region variable name
#' @param include_intercept whether to include an intercept
#' @param method estimation method ("SUR", "OLS", "WLS", "3SLS")
#' @param verbose verbose output
#' @return multivar_varx object
#' 
fit_varx <- function(
  data_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  method = "SUR",
  data = NULL,
  verbose = FALSE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== VARX model estimation (0011 model) ===\n")
    cat(sprintf("    Estimation method: %s\n", method))
    cat(sprintf("    Time lag: with cross terms\n"))
    cat(sprintf("    Error correlation: present\n"))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("systemfit", quietly = TRUE)) {
    stop("the systemfit package is required: install.packages('systemfit')")
  }
  
  k <- length(y_vars)
  
  # S0: data preparation
  
  if (verbose) cat("\nS0: data preparation...\n")
  
  if (is.null(data)) data <- read.csv(data_file)
  
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
  if (nrow(data_t_lag) != n) {
    stop("the number of regions in the time-lag data does not match")
  }
  
  if (region_var %in% names(data_t)) {
    data_t <- data_t[order(data_t[[region_var]]), ]
    data_t_lag <- data_t_lag[order(data_t_lag[[region_var]]), ]
  }
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Time point used: t = %d\n", time_point))
  }
  
  # S1: build the data frame for SUR
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Building the data frame (for SUR)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  df_sur <- data.frame(row.names = 1:n)
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    df_sur[[var_name]] <- data_t[[var_name]]
  }
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    lag_name <- paste0(var_name, "_lag")
    df_sur[[lag_name]] <- data_t_lag[[var_name]]
  }
  
  all_x_vars <- unique(unlist(x_vars))
  for (xvar in all_x_vars) {
    df_sur[[xvar]] <- data_t[[xvar]]
  }
  
  if (verbose) {
    cat(sprintf("  Data frame: %d rows × %d cols\n", nrow(df_sur), ncol(df_sur)))
    cat(sprintf("  Variables: %s\n", paste(names(df_sur), collapse=", ")))
  }
  
  # S2: build the formulas
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Building the formulas\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  formula_list <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    lag_vars <- paste0(y_vars, "_lag")
    
    rhs_vars <- c(lag_vars, xi_vars)
    
    if (include_intercept) {
      formula_str <- paste(var_name, "~", paste(rhs_vars, collapse = " + "))
    } else {
      formula_str <- paste(var_name, "~ -1 +", paste(rhs_vars, collapse = " + "))
    }
    
    formula_list[[var_name]] <- as.formula(formula_str)
    
    if (verbose) {
      cat(sprintf("  %s: %s\n", var_name, formula_str))
    }
  }
  
  # S3: SUR estimation
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("SUR estimation (%s)\n", method))
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  tryCatch({
    sur_result <- systemfit::systemfit(
      formula = formula_list,
      method = method,
      data = df_sur
    )
    
    if (verbose) {
      cat("  Estimation complete\n")
      cat("  Class of the result object:", class(sur_result), "\n")
      cat("  Elements of the result object:", paste(names(sur_result), collapse=", "), "\n")
    }
    
  }, error = function(e) {
    cat(sprintf("\nError: systemfit() failed\n"))
    cat("Error message:", e$message, "\n")
    stop("Estimation failed")
  })
  
  # S4: extract the results
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Extracting the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # residuals() and fitted() return data.frames
  res_df <- residuals(sur_result)
  fit_df <- fitted(sur_result)
  
  residuals_mat <- as.matrix(res_df)
  fitted_mat <- as.matrix(fit_df)
  
  res_colnames <- colnames(residuals_mat)
  
  if (verbose) {
    cat(sprintf("  Residual matrix: %d x %d\n", nrow(residuals_mat), ncol(residuals_mat)))
    cat(sprintf("  Column names: %s\n", paste(res_colnames, collapse=", ")))
  }
  
  # Reorder if the column names do not match y_vars
  if (!is.null(res_colnames) && all(y_vars %in% res_colnames)) {
    residuals_mat <- residuals_mat[, y_vars, drop = FALSE]
    fitted_mat <- fitted_mat[, y_vars, drop = FALSE]
  }
  
  Sigma <- (t(residuals_mat) %*% residuals_mat) / n
  rownames(Sigma) <- y_vars
  colnames(Sigma) <- y_vars
  
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  loglik <- -0.5 * k * n * log(2 * pi) - 0.5 * n * log_det_Sigma - 0.5 * n * k
  
  num_params_per_eq <- sapply(1:k, function(i) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    # intercept + k time lags + regressors
    n_params <- ifelse(include_intercept, 1, 0) + k + length(xi_vars)
    return(n_params)
  })
  num_beta <- sum(num_params_per_eq)
  num_sigma <- k * (k + 1) / 2  # number of independent parameters of Σ
  total_params <- num_beta + num_sigma
  
  ic <- compute_information_criteria(loglik, total_params, k * n)
  
  if (verbose) cat("  Information criteria computed\n")
  
  sur_summary <- summary(sur_result)
  
  if (verbose) cat("  summary obtained\n")
  
  A_mat <- matrix(0, k, k)
  A_se <- matrix(0, k, k)
  A_z <- matrix(0, k, k)
  A_p <- matrix(0, k, k)
  rownames(A_mat) <- y_vars
  colnames(A_mat) <- y_vars
  rownames(A_se) <- y_vars
  colnames(A_se) <- y_vars
  
  # coef(sur_summary) may return a list of summary.systemfit.equation objects
  # Use sur_result$coefficients directly
  all_coef_vec <- coef(sur_result)
  all_se_vec <- sqrt(diag(vcov(sur_result)))
  
  if (verbose) {
    cat("  Length of the coefficient vector:", length(all_coef_vec), "\n")
    cat("  Coefficient names:\n")
    print(names(all_coef_vec))
  }
  
  coef_names <- names(all_coef_vec)
  
  for (i in 1:k) {
    eq_prefix <- y_vars[i]
    
    for (j in 1:k) {
      lag_name <- paste0(y_vars[j], "_lag")
      # systemfit coefficient names have the form "eqname_varname"
      coef_name <- paste0(eq_prefix, "_", lag_name)
      
      if (coef_name %in% coef_names) {
        A_mat[i, j] <- all_coef_vec[coef_name]
        A_se[i, j] <- all_se_vec[coef_name]
        A_z[i, j] <- all_coef_vec[coef_name] / all_se_vec[coef_name]
        A_p[i, j] <- 2 * pnorm(-abs(A_z[i, j]))
      }
    }
  }
  
  if (verbose) cat("  A matrix extracted\n")
  
  beta_list <- list()
  coef_tables <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    eq_prefix <- var_name
    
    eq_coef_idx <- grep(paste0("^", eq_prefix, "_"), coef_names)
    eq_coef_names <- coef_names[eq_coef_idx]
    eq_coef_vals <- all_coef_vec[eq_coef_idx]
    eq_se_vals <- all_se_vec[eq_coef_idx]
    
    short_names <- sub(paste0("^", eq_prefix, "_"), "", eq_coef_names)
    
    lag_names <- paste0(y_vars, "_lag")
    non_lag_mask <- !short_names %in% lag_names
    
    beta_coef <- eq_coef_vals[non_lag_mask]
    names(beta_coef) <- short_names[non_lag_mask]
    beta_list[[var_name]] <- beta_coef
    
    z_vals <- eq_coef_vals / eq_se_vals
    p_vals <- 2 * pnorm(-abs(z_vals))
    
    coef_table <- data.frame(
      parameter = short_names,
      estimate = eq_coef_vals,
      std_error = eq_se_vals,
      z_value = z_vals,
      p_value = p_vals,
      stringsAsFactors = FALSE
    )
    coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
    rownames(coef_table) <- NULL
    coef_tables[[var_name]] <- coef_table
  }
  
  if (verbose) cat("  Regression coefficients extracted\n")
  
  residuals_raw <- as.vector(residuals_mat)
  fitted_vals <- as.vector(fitted_mat)
  
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_mat[, i] / sqrt(Sigma[i, i])
  }
  
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- df_sur[[y_vars[i]]]
  }
  
  r2_list <- numeric(k)
  r2_adj_list <- numeric(k)
  names(r2_list) <- y_vars
  names(r2_adj_list) <- y_vars
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    y_i <- df_sur[[var_name]]
    fitted_i <- fitted_mat[, i]
    residuals_i <- residuals_mat[, i]
    
    r2_results <- compute_r_squared(
      y = y_i,
      fitted = fitted_i,
      residuals = residuals_i,
      loglik = NULL,
      num_params = num_params_per_eq[i]
    )
    
    r2_list[i] <- r2_results$R2
    r2_adj_list[i] <- r2_results$R2_adj
  }
  
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  # S5: build the per-variable model results
  
  individual_models <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    
    individual_models[[var_name]] <- list(
      time_lag_params = list(
        alpha = A_mat[i, ],
        alpha_se = A_se[i, ],
        alpha_z = A_z[i, ],
        alpha_p = A_p[i, ],
        alpha_signif = sapply(A_p[i, ], get_signif_code)
      ),
      
      coefficients = coef_tables[[var_name]],
      
      sigma2 = Sigma[i, i],
      
      residuals = residuals_mat[, i],
      fitted = fitted_mat[, i],
      
      fit = list(
        R2 = r2_list[i],
        R2_adj = r2_adj_list[i],
        num_params = num_params_per_eq[i],
        num_obs = n
      )
    )
  }
  
  # Output results
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[A matrix (time-lag coefficients)]\n")
    print(round(A_mat, 4))
    
    cat("\n[Standard errors of the A matrix]\n")
    print(round(A_se, 4))
    
    cat("\n[p-values of the A matrix]\n")
    print(round(A_p, 4))
    
    cat("\n[Σ matrix (error covariance)]\n")
    print(round(Sigma, 6))
    
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    rownames(corr_mat) <- y_vars
    colnames(corr_mat) <- y_vars
    cat("\n[Error correlation matrix]\n")
    print(round(corr_mat, 4))
    
    cat("\n[Goodness of fit (merged)]\n")
    cat(sprintf("  Log-likelihood: %.4f\n", loglik))
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
    y = y_vec, W = NULL, W_listw = NULL,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = TRUE,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "VARX",
    Sigma             = Sigma,
    loglik            = loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = method,
    data_list         = dl,
    beta0             = beta_list,
    alpha             = A_mat,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  result$coefficients$A <- A_mat
  result$std_errors$A <- A_se
  result$inference <- list(A_z = A_z, A_p = A_p)
  result$sur_result <- sur_result
  result$data_info$time_lag_cross <- TRUE
  
  return(result)
}

# Standard methods

#' print method
print.multivar_varx <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("VARX model (0011 model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time lag: with cross terms\n"))
  cat(sprintf("  Error correlation: present\n"))
  
  cat("\n[A matrix (time-lag coefficients)]\n")
  print(round(x$coefficients$A, digits))
  
  cat("\n[Σ matrix (error covariance) diagonal elements]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$Sigma)), collapse=", ")))
  
  Sigma <- x$coefficients$Sigma
  if (nrow(Sigma) > 1) {
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    cat(sprintf("  Error correlation (1,2): %.4f\n", corr_mat[1,2]))
  }
  
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
summary.multivar_varx <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("VARX model - detailed results (0011 model)\n")
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
  cat(sprintf("  Estimation method: %s\n", x$convergence$method))
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[2. Time-lag coefficients (A matrix)]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\nCoefficient matrix:\n")
  print(round(x$coefficients$A, digits))
  
  cat("\nStandard errors:\n")
  print(round(x$std_errors$A, digits))
  
  cat("\np-values:\n")
  print(round(x$inference$A_p, digits))
  
  cat("\nDetails (rows: dependent variables, columns: lag variables):\n")
  k <- x$data_info$k
  y_vars <- x$data_info$y_vars
  
  for (i in 1:k) {
    for (j in 1:k) {
      cat(sprintf("  α[%s,%s_lag] = %7.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  y_vars[i], y_vars[j],
                  x$coefficients$A[i, j],
                  x$std_errors$A[i, j],
                  x$inference$A_z[i, j],
                  x$inference$A_p[i, j],
                  get_signif_code(x$inference$A_p[i, j])))
    }
  }
  
  for (i in 1:k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. Details of the %s equation]\n", i + 2, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    cat(sprintf("\nError variance: σ² = %.6f\n", model_i$sigma2))
    
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²       = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²   = %.4f\n", model_i$fit$R2_adj))
  }
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. Error covariance matrix Σ]\n", k + 3))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\nCovariance matrix:\n")
  print(round(x$coefficients$Sigma, digits + 2))
  
  Sigma <- x$coefficients$Sigma
  D <- diag(1/sqrt(diag(Sigma)))
  corr_mat <- D %*% Sigma %*% D
  rownames(corr_mat) <- y_vars
  colnames(corr_mat) <- y_vars
  cat("\nCorrelation matrix:\n")
  print(round(corr_mat, digits))
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. Merged goodness of fit]\n", k + 4))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  Log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  Mean R²       = %.4f\n", x$fit$R2))
  cat(sprintf("  Mean Adj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}

#' coef method
coef.multivar_varx <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(NULL)
  } else if (type == "Lambda") {
    return(NULL)
  } else if (type == "A") {
    return(object$coefficients$A)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}

#' logLik method
logLik.multivar_varx <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}

#' AIC method
AIC.multivar_varx <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}

#' BIC method
BIC.multivar_varx <- function(object, ...) {
  return(object$fit$BIC)
}

#' residuals method
residuals.multivar_varx <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}

#' fitted method
fitted.multivar_varx <- function(object, ...) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars
  
  fitted_vals <- numeric(k * n)
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    idx <- ((i - 1) * n + 1):(i * n)
    fitted_vals[idx] <- object$individual_models[[var_name]]$fitted
  }
  
  return(fitted_vals)
}
