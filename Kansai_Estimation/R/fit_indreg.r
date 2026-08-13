# fit_indreg.r
# Diagonal-constrained OLS estimation (000d model): each variable estimated
# separately with lm(); no spatial lag/error, no time lag, diagonal Σ.

#' Diagonal-constrained OLS estimation (000d model)
#' 
#' Estimate each variable separately with lm() and merge the results
#' 
#' @param data_file path to the data file
#' @param y_vars vector of dependent variable names
#' @param x_vars list of regressor names (per dependent variable)
#' @param time_var time variable name
#' @param time_point time point to use (latest if NULL)
#' @param region_var region variable name
#' @param include_intercept whether to include an intercept
#' @param verbose verbose output
#' @return multivar_indreg object
#' 
fit_indreg <- function(
  data_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  verbose = FALSE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained OLS estimation (000d model) ===\n")
    cat(sprintf("    Time lag: none\n"))
    cat(sprintf("    Spatial correlation: none\n"))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  k <- length(y_vars)

  # S0: data preparation
  if (verbose) cat("\nS0: data preparation...\n")

  data <- read.csv(data_file)

  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }

  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)

  if (region_var %in% names(data_t)) {
    data_t <- data_t[order(data_t[[region_var]]), ]
  }
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Time point used: t = %d\n", time_point))
  }
  
  # S1-S2: estimate each variable separately with lm()
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Per-variable OLS estimation (lm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
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

    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)

    tryCatch({
      model_i <- lm(
        formula = formula_i,
        data = df_i
      )

      s <- summary(model_i)
      sigma2 <- s$sigma^2
      Sigma_diag[i, i] <- sigma2

      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef

      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))

      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  # β + σ²
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i

      coef_mat <- s$coefficients
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, "Estimate"],
        std_error = coef_mat[, "Std. Error"],
        t_value = coef_mat[, "t value"],
        p_value = coef_mat[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL

      R2 <- s$r.squared
      R2_adj <- s$adj.r.squared

      f_stat <- s$fstatistic
      if (!is.null(f_stat)) {
        f_value <- f_stat[1]
        f_df1 <- f_stat[2]
        f_df2 <- f_stat[3]
        f_p <- pf(f_value, f_df1, f_df2, lower.tail = FALSE)
      } else {
        f_value <- NA
        f_df1 <- NA
        f_df2 <- NA
        f_p <- NA
      }

      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        coefficients = coef_table,
        sigma2 = sigma2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = R2,
          R2_adj = R2_adj,
          F_statistic = f_value,
          F_df1 = f_df1,
          F_df2 = f_df2,
          F_p_value = f_p,
          num_params = num_params_i,
          num_obs = length(df_i$y),
          df_residual = s$df[2]
        )
      )

      if (verbose) {
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f (σ = %.4f)\n", sigma2, sqrt(sigma2)))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, t: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$t_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", R2))
        cat(sprintf("    Adj.R²     = %.4f\n", R2_adj))
        cat(sprintf("    Log-likelihood = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
        
        if (!is.na(f_value)) {
          cat(sprintf("    F statistic = %.2f (df1=%d, df2=%d, p=%.4f)\n",
                      f_value, f_df1, f_df2, f_p))
        }
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: lm() failed when estimating variable %s\n", var_name))
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

  # For OLS, R, Λ, A are all absent
  coefficients <- list(
    R = NULL,
    Lambda = NULL,
    A = NULL,
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")

  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[Diagonal Σ matrix (error variances)]\n")
    print(round(Sigma_diag, 6))
    
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
    y = y_vec,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = FALSE,
    include_intercept = include_intercept
  )
  
  result <- build_result_object(
    model_type        = "IndReg",
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "lm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}


# Standard methods

#' print method
print.multivar_indreg <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained OLS (000d model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time lag: none\n"))
  cat(sprintf("  Spatial correlation: none\n"))
  
  cat("\n[Diagonal Σ matrix (error variances)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$Sigma)), collapse=", ")))
  
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
summary.multivar_indreg <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained OLS - detailed results (000d model)\n")
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

  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. Estimation results for %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")

    cat(sprintf("\nError variance: σ² = %.6f (σ = %.4f)\n",
                model_i$sigma2, sqrt(model_i$sigma2)))

    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, t: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$t_value, row$p_value, row$signif))
    }

    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  Log-likelihood = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))

    if (!is.na(model_i$fit$F_statistic)) {
      cat(sprintf("  F statistic = %.2f (df1=%d, df2=%d, p=%.4f)\n",
                  model_i$fit$F_statistic, 
                  model_i$fit$F_df1, 
                  model_i$fit$F_df2, 
                  model_i$fit$F_p_value))
    }
  }

  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. Merged goodness of fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\nDiagonal Σ matrix:\n")
  print(round(x$coefficients$Sigma, digits + 2))
  
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
coef.multivar_indreg <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(NULL)
  } else if (type == "Lambda") {
    return(NULL)
  } else if (type == "A") {
    return(NULL)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' logLik method
logLik.multivar_indreg <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' AIC method
AIC.multivar_indreg <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' BIC method
BIC.multivar_indreg <- function(object, ...) {
  return(object$fit$BIC)
}


#' residuals method
residuals.multivar_indreg <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' fitted method
fitted.multivar_indreg <- function(object, ...) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars
  
  fitted_vals <- numeric(k * n)
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    idx <- ((i - 1) * n + 1):(i * n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}

