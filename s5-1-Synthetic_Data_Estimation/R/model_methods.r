# model_methods.r
# Unified output functions for the multivariate spatial regression models (MSAR/MSEM/MGNS)
# Provides the S3 methods of the multivar_spatial class (print/summary/coef/vcov/fitted/
# residuals/logLik/AIC/BIC/confint) and significance tests (add_inference).

# 1. print function

#' Basic display of a multivar_spatial object
#' 
#' @param x multivar_spatial object
#' @param digits number of decimal places
#' @param ... other arguments
#' @export
print.multivar_spatial <- function(x, digits = 4, ...) {
  
  model_name <- switch(x$model_type,
    "MSAR" = "Spatial Lag of Y Model (full)",
    "MSEM" = "Spatial Error Model (full)",
    "MGNS" = "General Nesting Spatial Model (full)",
    "IndSAR" = "Spatial Lag of Y Model (diagonal)",
    "IndSEM" = "Spatial Error Model (diagonal)",
    "IndGNS" = "General Nesting Spatial Model (diagonal)",
    "VARX" = "VARX Model",
    "IndReg" = "Independent Regression Model (IndReg)",
    x$model_type  # fallback
  )
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat(sprintf("Multivariate spatial regression model (%s)\n", model_name))
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  
  if (!is.null(x$data_info$include_time_lag)) {
    cat(sprintf("  Time lag: %s\n", 
                ifelse(x$data_info$include_time_lag, "included", "not included")))
  }
  
  cat("\n[Convergence status]\n")
  cat(sprintf("  Converged: %s\n", ifelse(x$convergence$converged, "yes", "no")))
  cat(sprintf("  Number of iterations: %d\n", x$convergence$iterations))
  cat(sprintf("  Optimization method: %s\n", x$convergence$method))
  
  if (x$model_type %in% c("MSAR", "MGNS", "IndSAR", "IndGNS") && !is.null(x$coefficients$R)) {
    cat("\n[Spatial lag matrix R]\n")
    cat(sprintf("  Diagonal elements: %s\n", 
                paste(sprintf("%.4f", diag(x$coefficients$R)), collapse=", ")))
    
    R <- x$coefficients$R
    offdiag <- R[row(R) != col(R)]
    if (any(abs(offdiag) > 0.001)) {
      cat(sprintf("  Max absolute off-diagonal value: %.4f\n", max(abs(offdiag))))
    }
  }
  
  if (x$model_type %in% c("MSEM", "MGNS", "IndSEM", "IndGNS") && !is.null(x$coefficients$Lambda)) {
    cat("\n[Spatial error matrix Λ]\n")
    cat(sprintf("  Diagonal elements: %s\n", 
                paste(sprintf("%.4f", diag(x$coefficients$Lambda)), collapse=", ")))
    
    Lambda_mat <- x$coefficients$Lambda
    offdiag <- Lambda_mat[row(Lambda_mat) != col(Lambda_mat)]
    if (any(abs(offdiag) > 0.001)) {
      cat(sprintf("  Max absolute off-diagonal value: %.4f\n", max(abs(offdiag))))
    }
  }
  
  cat("\n[Goodness of fit]\n")
  cat(sprintf("  Log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  
  cat("\n")
  cat("Use summary() for details\n")
  cat("\n")
  
  invisible(x)
}

# Aliases for MSAR, MSEM, MGNS
print.multivar_msar <- function(x, ...) print.multivar_spatial(x, ...)
print.multivar_msem <- function(x, ...) print.multivar_spatial(x, ...)
print.multivar_mgns <- function(x, ...) print.multivar_spatial(x, ...)

# 2. summary function

#' Detailed display of a multivar_spatial object
#' 
#' @param object multivar_spatial object
#' @param show_s1 whether to display the S1 initial-estimation results
#' @param show_inference whether to display significance-test results
#' @param digits number of decimal places
#' @param ... other arguments
#' @export
summary.multivar_spatial <- function(object, show_s1 = FALSE, show_inference = TRUE, 
                                      digits = 4, ...) {
  
  x <- object
  
  model_name <- switch(x$model_type,
    "MSAR" = "Spatial Lag of Y Model (MSAR, full)",
    "MSEM" = "Spatial Error Model (MSEM, full)",
    "MGNS" = "General Nesting Spatial Model (MGNS, full)",
    "IndSAR" = "Spatial Lag of Y Model (MSAR, diagonal)",
    "IndSEM" = "Spatial Error Model (MSEM, diagonal)",
    "IndGNS" = "General Nesting Spatial Model (MGNS, diagonal)",
    "VARX" = "VARX Model",
    "IndReg" = "Independent Regression Model (IndReg)",
    x$model_type  # fallback
  )
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat(sprintf("Multivariate spatial regression model - detailed results\n"))
  cat(sprintf("Model type: %s\n", model_name))
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Basic information
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Dependent variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  Number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  Number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  
  if (!is.null(x$data_info$time_point_used)) {
    cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  }
  
  if (!is.null(x$data_info$include_time_lag)) {
    cat(sprintf("  Time lag: %s\n", 
                ifelse(x$data_info$include_time_lag, "included (y_{t-1})", "not included")))
  }
  
  # Spatial parameters
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[2. Spatial parameters]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (x$model_type %in% c("MSAR", "MGNS", "IndSAR", "IndGNS") && !is.null(x$coefficients$R)) {
    cat("\nSpatial lag matrix R:\n")
    print(round(x$coefficients$R, digits))
    
    if (show_inference && !is.null(x$std_errors$R)) {
      cat("\nStandard errors of the R matrix:\n")
      print(round(x$std_errors$R, digits))
    }
  }
  
  if (x$model_type %in% c("MSEM", "MGNS", "IndSEM", "IndGNS") && !is.null(x$coefficients$Lambda)) {
    cat("\nSpatial error matrix Λ:\n")
    print(round(x$coefficients$Lambda, digits))
    
    if (show_inference && !is.null(x$std_errors$Lambda)) {
      cat("\nStandard errors of the Λ matrix:\n")
      print(round(x$std_errors$Lambda, digits))
    }
  }
  
  # Regression coefficients
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[3. Regression coefficients]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (show_inference && !is.null(x$inference$coefficients_table)) {
    cat("\nCoefficient table (with significance tests):\n")
    print_coef_table(x$inference$coefficients_table, digits)
  } else {
    for (i in 1:x$data_info$k) {
      var_name <- x$data_info$y_vars[i]
      cat(sprintf("\n%s coefficients:\n", var_name))
      
      if (!is.null(x$coefficients$beta0)) {
        coef_i <- x$coefficients$beta0[[var_name]]
        if (is.null(coef_i)) {
          coef_i <- x$coefficients$beta0[[paste0("y", i)]]
        }
        
        if (!is.null(coef_i)) {
          for (j in seq_along(coef_i)) {
            cat(sprintf("  %-20s: %10.4f\n", names(coef_i)[j], coef_i[j]))
          }
        }
      }
    }
  }
  
  if (!is.null(x$coefficients$alpha) && 
      !is.null(x$data_info$include_time_lag) && 
      x$data_info$include_time_lag) {
    cat("\nTime-lag coefficient matrix α:\n")
    print(round(x$coefficients$alpha, digits))
  }
  
  if (!is.null(x$coefficients$A)) {
    cat("\nTime-lag coefficient matrix A:\n")
    print(round(x$coefficients$A, digits))
  }
  
  # Error covariance matrix
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[4. Error covariance matrix Σ]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\n")
  print(round(x$coefficients$Sigma, digits))
  
  Sigma <- x$coefficients$Sigma
  if (nrow(Sigma) > 1) {
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    cat("\nError correlation matrix:\n")
    print(round(corr_mat, digits))
  }
  
  # Goodness of fit
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[5. Goodness-of-fit statistics]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  Log-likelihood:     %12.4f\n", x$fit$loglik))
  
  if (!is.null(x$fit$profile_loglik)) {
    cat(sprintf("  Profile likelihood: %12.4f\n", x$fit$profile_loglik))
  }
  
  cat(sprintf("  AIC:                %12.4f\n", x$fit$AIC))
  cat(sprintf("  BIC:                %12.4f\n", x$fit$BIC))
  
  if (!is.null(x$fit$R2)) {
    cat(sprintf("  R²:                 %12.4f\n", x$fit$R2))
  }
  if (!is.null(x$fit$R2_adj)) {
    cat(sprintf("  Adjusted R²:        %12.4f\n", x$fit$R2_adj))
  }
  
  # Convergence diagnostics
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[6. Convergence diagnostics]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  Convergence status: %s\n", 
              ifelse(x$convergence$converged, "converged", "not converged")))
  cat(sprintf("  Number of iterations: %d\n", x$convergence$iterations))
  cat(sprintf("  Optimization method: %s\n", x$convergence$method))
  
  if (!is.null(x$convergence$final_gradient_norm)) {
    cat(sprintf("  Final gradient norm: %.2e\n", x$convergence$final_gradient_norm))
  }
  if (!is.null(x$convergence$final_param_change)) {
    cat(sprintf("  Final parameter change: %.2e\n", x$convergence$final_param_change))
  }
  if (!is.null(x$convergence$final_lik_change)) {
    cat(sprintf("  Final likelihood change: %.2e\n", x$convergence$final_lik_change))
  }
  
  cat(sprintf("  Message: %s\n", x$convergence$message))
  
  # Residual summary
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[7. Residual summary]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (!is.null(x$residuals$standardized)) {
    cat("\nDistribution of standardized residuals:\n")
    res_summary <- summary(x$residuals$standardized)
    print(res_summary)
  } else if (!is.null(x$residuals$raw)) {
    cat("\nDistribution of residuals:\n")
    res_summary <- summary(x$residuals$raw)
    print(res_summary)
  }
  
  # S1 initial-estimation results (optional)
  
  if (show_s1 && !is.null(x$initial_values$individual_estimates)) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("[8. S1 initial-estimation results]\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    for (var_name in names(x$initial_values$individual_estimates)) {
      s1 <- x$initial_values$individual_estimates[[var_name]]
      cat(sprintf("\n%s:\n", var_name))
      
      if (!is.null(s1$spatial_params$rho)) {
        cat(sprintf("  ρ = %.4f (SE: %.4f, p: %.4f %s)\n",
                    s1$spatial_params$rho,
                    s1$spatial_params$rho_se,
                    s1$spatial_params$rho_p,
                    s1$spatial_params$rho_signif))
      }
      if (!is.null(s1$spatial_params$lambda)) {
        cat(sprintf("  λ = %.4f (SE: %.4f, p: %.4f %s)\n",
                    s1$spatial_params$lambda,
                    s1$spatial_params$lambda_se,
                    s1$spatial_params$lambda_p,
                    s1$spatial_params$lambda_signif))
      }
      
      if (!is.null(s1$fit)) {
        cat(sprintf("  R² = %.4f, Adj.R² = %.4f\n", 
                    s1$fit$R2, s1$fit$R2_adj))
        cat(sprintf("  Log-likelihood = %.4f, AIC = %.4f\n",
                    s1$fit$loglik, s1$fit$AIC))
      }
    }
  }
  
  # Execution information
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[Execution information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (!is.null(x$execution$time)) {
    cat(sprintf("\n  Execution time: %.2f s\n", as.numeric(x$execution$time)))
  }
  if (!is.null(x$execution$R_version)) {
    cat(sprintf("  R version: %s\n", x$execution$R_version))
  }
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Significance levels: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("\n")
  
  invisible(x)
}

# Aliases for MSAR, MSEM, MGNS
summary.multivar_msar <- function(object, ...) summary.multivar_spatial(object, ...)
summary.multivar_msem <- function(object, ...) summary.multivar_spatial(object, ...)
summary.multivar_mgns <- function(object, ...) summary.multivar_spatial(object, ...)

# 3. Formatted output of the coefficient table

#' Formatted output of the coefficient table
#' 
#' @param coef_table data.frame (parameter, estimate, std_error, z_value, p_value, signif)
#' @param digits number of decimal places
print_coef_table <- function(coef_table, digits = 4) {
  
  cat(paste(rep("-", 75), collapse=""), "\n")
  cat(sprintf("%-25s %10s %10s %10s %10s %5s\n",
              "Parameter", "Estimate", "Std.Error", "z-value", "p-value", ""))
  cat(paste(rep("-", 75), collapse=""), "\n")
  
  for (i in 1:nrow(coef_table)) {
    row <- coef_table[i, ]
    
    est <- if (is.na(row$estimate)) "NA" else sprintf("%10.4f", row$estimate)
    se <- if (is.na(row$std_error)) "NA" else sprintf("%10.4f", row$std_error)
    z <- if (is.na(row$z_value)) "NA" else sprintf("%10.4f", row$z_value)
    p <- if (is.na(row$p_value)) "NA" else sprintf("%10.4f", row$p_value)
    sig <- if (is.na(row$signif)) "" else row$signif
    
    cat(sprintf("%-25s %10s %10s %10s %10s %5s\n",
                row$parameter, est, se, z, p, sig))
  }
  
  cat(paste(rep("-", 75), collapse=""), "\n")
}

# 4. Standard method functions

#' extract coefficients
#' 
#' @param object multivar_spatial object
#' @param type extraction type ("all", "R", "Lambda", "beta", "alpha", "Sigma")
#' @param ... other arguments
#' @export
coef.multivar_spatial <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(object$coefficients$R)
  } else if (type == "Lambda") {
    return(object$coefficients$Lambda)
  } else if (type == "beta") {
    return(object$coefficients$beta)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else if (type == "alpha") {
    return(object$coefficients$alpha)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

coef.multivar_msar <- function(object, ...) coef.multivar_spatial(object, ...)
coef.multivar_msem <- function(object, ...) coef.multivar_spatial(object, ...)
coef.multivar_mgns <- function(object, ...) coef.multivar_spatial(object, ...)

#' extract the variance-covariance matrix
#' 
#' @param object multivar_spatial object
#' @param type extraction type ("beta", "spatial", "full")
#' @param ... other arguments
#' @export
vcov.multivar_spatial <- function(object, type = "beta", ...) {
  
  if (is.null(object$vcov)) {
    warning("The variance-covariance matrix has not been computed")
    return(NULL)
  }
  
  if (type == "beta") {
    return(object$vcov$beta)
  } else if (type == "spatial") {
    return(object$vcov$spatial)
  } else if (type == "full") {
    return(object$vcov$full)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

vcov.multivar_msar <- function(object, ...) vcov.multivar_spatial(object, ...)
vcov.multivar_msem <- function(object, ...) vcov.multivar_spatial(object, ...)
vcov.multivar_mgns <- function(object, ...) vcov.multivar_spatial(object, ...)

#' extract fitted values
#' 
#' @param object multivar_spatial object
#' @param ... other arguments
#' @export
fitted.multivar_spatial <- function(object, ...) {
  
  if (!is.null(object$residuals$fitted)) {
    return(object$residuals$fitted)
  }
  
  if (!is.null(object$model_data$y) && !is.null(object$residuals$raw)) {
    return(object$model_data$y - object$residuals$raw)
  }
  
  warning("Cannot compute fitted values")
  return(NULL)
}

fitted.multivar_msar <- function(object, ...) fitted.multivar_spatial(object, ...)
fitted.multivar_msem <- function(object, ...) fitted.multivar_spatial(object, ...)
fitted.multivar_mgns <- function(object, ...) fitted.multivar_spatial(object, ...)

#' extract residuals
#' 
#' @param object multivar_spatial object
#' @param type residual type ("raw", "standardized")
#' @param ... other arguments
#' @export
residuals.multivar_spatial <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

residuals.multivar_msar <- function(object, ...) residuals.multivar_spatial(object, ...)
residuals.multivar_msem <- function(object, ...) residuals.multivar_spatial(object, ...)
residuals.multivar_mgns <- function(object, ...) residuals.multivar_spatial(object, ...)

#' extract the log-likelihood
#' 
#' @param object multivar_spatial object
#' @param ... other arguments
#' @export
logLik.multivar_spatial <- function(object, ...) {
  
  ll <- object$fit$loglik
  attr(ll, "df") <- object$fit$num_params
  attr(ll, "nobs") <- object$fit$num_obs
  class(ll) <- "logLik"
  
  return(ll)
}

logLik.multivar_msar <- function(object, ...) logLik.multivar_spatial(object, ...)
logLik.multivar_msem <- function(object, ...) logLik.multivar_spatial(object, ...)
logLik.multivar_mgns <- function(object, ...) logLik.multivar_spatial(object, ...)

#' extract the AIC
#' 
#' @param object multivar_spatial object
#' @param ... other arguments
#' @param k penalty (default: 2)
#' @export
AIC.multivar_spatial <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}

AIC.multivar_msar <- function(object, ...) AIC.multivar_spatial(object, ...)
AIC.multivar_msem <- function(object, ...) AIC.multivar_spatial(object, ...)
AIC.multivar_mgns <- function(object, ...) AIC.multivar_spatial(object, ...)

#' extract the BIC
#' 
#' @param object multivar_spatial object
#' @param ... other arguments
#' @export
BIC.multivar_spatial <- function(object, ...) {
  return(object$fit$BIC)
}

BIC.multivar_msar <- function(object, ...) BIC.multivar_spatial(object, ...)
BIC.multivar_msem <- function(object, ...) BIC.multivar_spatial(object, ...)
BIC.multivar_mgns <- function(object, ...) BIC.multivar_spatial(object, ...)

# 5. Confidence intervals

#' Compute confidence intervals for the parameters
#' 
#' @param object multivar_spatial object
#' @param parm parameter name (optional)
#' @param level confidence level (default: 0.95)
#' @param type target parameters ("beta", "R", "Lambda", "all")
#' @param ... other arguments
#' @export
confint.multivar_spatial <- function(object, parm = NULL, level = 0.95, 
                                      type = "beta", ...) {
  
  alpha <- 1 - level
  z_crit <- qnorm(1 - alpha/2)
  
  if (type == "beta" || type == "all") {
    if (is.null(object$std_errors$beta)) {
      warning("The standard errors of β have not been computed")
      return(NULL)
    }
    
    beta <- object$coefficients$beta
    se_beta <- object$std_errors$beta
    
    ci <- data.frame(
      estimate = beta,
      lower = beta - z_crit * se_beta,
      upper = beta + z_crit * se_beta
    )
    
    rownames(ci) <- names(beta)
    colnames(ci) <- c("Estimate", sprintf("%.1f%%", 100*alpha/2), 
                       sprintf("%.1f%%", 100*(1-alpha/2)))
    
    return(ci)
  }
  
  warning("Currently only confidence intervals for β are supported")
  return(NULL)
}

confint.multivar_msar <- function(object, ...) confint.multivar_spatial(object, ...)
confint.multivar_msem <- function(object, ...) confint.multivar_spatial(object, ...)
confint.multivar_mgns <- function(object, ...) confint.multivar_spatial(object, ...)

# 6. Model comparison

#' Build a comparison table of multiple models
#' 
#' @param ... multivar_spatial objects
#' @param names model names (optional)
#' @return comparison table (data.frame)
compare_models <- function(..., names = NULL) {
  
  models <- list(...)
  n_models <- length(models)
  
  if (is.null(names)) {
    names <- paste0("Model_", 1:n_models)
  }
  
  comparison <- data.frame(
    Model = names,
    Type = character(n_models),
    k = integer(n_models),
    n = integer(n_models),
    Params = integer(n_models),
    LogLik = numeric(n_models),
    AIC = numeric(n_models),
    BIC = numeric(n_models),
    Converged = logical(n_models),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_models) {
    m <- models[[i]]
    comparison$Type[i] <- m$model_type
    comparison$k[i] <- m$data_info$k
    comparison$n[i] <- m$data_info$n
    comparison$Params[i] <- m$fit$num_params
    comparison$LogLik[i] <- m$fit$loglik
    comparison$AIC[i] <- m$fit$AIC
    comparison$BIC[i] <- m$fit$BIC
    comparison$Converged[i] <- m$convergence$converged
  }
  
  comparison <- comparison[order(comparison$AIC), ]
  
  comparison$Delta_AIC <- comparison$AIC - min(comparison$AIC)
  comparison$Delta_BIC <- comparison$BIC - min(comparison$BIC)
  
  class(comparison) <- c("model_comparison", "data.frame")
  return(comparison)
}

#' Display the model-comparison table
#' 
#' @param x model_comparison object
#' @param ... other arguments
print.model_comparison <- function(x, ...) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  cat("Model comparison\n")
  cat(paste(rep("=", 80), collapse=""), "\n\n")
  
  x$LogLik <- round(x$LogLik, 2)
  x$AIC <- round(x$AIC, 2)
  x$BIC <- round(x$BIC, 2)
  x$Delta_AIC <- round(x$Delta_AIC, 2)
  x$Delta_BIC <- round(x$Delta_BIC, 2)
  
  print(as.data.frame(x), row.names = FALSE)
  
  cat("\nBest model (AIC):", x$Model[1], "\n")
  cat("Best model (BIC):", x$Model[which.min(x$BIC)], "\n")
  cat("\n")
  
  invisible(x)
}

# 7. Significance-test computation (run after estimation)

#' Add significance tests to the estimation result
#' 
#' Call once after estimation completes to compute vcov, std_errors, inference
#' 
#' @param object multivar_spatial object
#' @param compute_spatial_se whether to compute spatial-parameter standard errors (uses the numerical Hessian)
#' @param verbose verbose output
#' @return the object with significance-test info added
add_inference <- function(object, compute_spatial_se = TRUE, gamma = 0, verbose = FALSE) {
  
  # diagonal/VARX/OLS models already have std_errors from per-variable estimation
  if (object$model_type %in% c("IndSAR", "IndSEM", "IndGNS",
                                "VARX", "IndReg")) {
    if (verbose) {
      cat(sprintf("\n%s: using std_errors from per-variable estimation (skipping add_inference)\n", object$model_type))
    }
    return(object)
  }
  
  if (verbose) {
    cat("\nComputing significance tests...\n")
  }
  
  k <- object$data_info$k
  n <- object$data_info$n
  X <- object$model_data$X
  Sigma <- object$coefficients$Sigma
  
  # Variance-covariance matrix of β (analytical)
  
  if (verbose) cat("  Computing the variance-covariance matrix of β...\n")
  
  # For MSEM/MGNS, pass the Λ matrix and W
  Lambda_for_vcov <- if (object$model_type %in% c("MSEM", "MGNS")) object$coefficients$Lambda else NULL
  W_for_vcov <- if (!is.null(Lambda_for_vcov)) object$model_data$W else NULL
  
  Psi <- tryCatch({
    compute_vcov_beta(X, Sigma, k, n, Lambda_mat = Lambda_for_vcov, W = W_for_vcov)
  }, error = function(e) {
    warning("Error computing the variance-covariance matrix of β: ", e$message)
    NULL
  })
  
  se_beta <- NULL
  if (!is.null(Psi)) {
    se_beta <- sqrt(diag(Psi))
  }
  
  # Spatial-parameter standard errors (numerical Hessian)
  
  se_R <- NULL
  se_Lambda <- NULL
  hessian_result <- NULL
  vcov_spatial <- NULL
  
  if (compute_spatial_se) {
    if (verbose) cat("  Computing spatial-parameter standard errors (numerical Hessian)...\n")
    
    if (object$model_type == "MSAR" && !is.null(object$coefficients$R)) {
      hessian_result <- compute_hessian_for_R(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          se_R <- matrix(sqrt(diag(vcov_spatial)), nrow = k, ncol = k, byrow = TRUE)
        }
      }
      
    } else if (object$model_type == "MSEM" && !is.null(object$coefficients$Lambda)) {
      hessian_result <- compute_hessian_for_Lambda(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          se_Lambda <- matrix(sqrt(diag(vcov_spatial)), nrow = k, ncol = k, byrow = TRUE)
        }
      }
      
    } else if (object$model_type == "MGNS") {
      hessian_result <- compute_hessian_for_RLambda(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          n_R <- k^2
         se_vec <- sqrt(diag(vcov_spatial))
         se_R <- matrix(se_vec[1:n_R], nrow = k, ncol = k, byrow = TRUE)
         se_Lambda <- matrix(se_vec[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
        }
      }
    }
  }
  
  if (verbose) cat("  Building the coefficient table...\n")
  
  coef_table <- create_coefficient_table(object, se_beta, se_R, se_Lambda)
  
  object$vcov <- list(
    beta = Psi,
    spatial = vcov_spatial,
    full = NULL
  )
  
  object$hessian <- list(
    matrix = hessian_result,
    method = "numerical"
  )
  
  object$std_errors <- list(
    beta = se_beta,
    R = se_R,
    Lambda = se_Lambda
  )
  
  object$inference <- list(
    coefficients_table = coef_table
  )
  
  if (verbose) {
    cat("  Done\n\n")
  }
  
  return(object)
}

#' Hessian computation for the R matrix (for MSAR)
#' 
#' @param object multivar_spatial object
#' @param eps step size for numerical differentiation
#' @param verbose verbose output
#' @return Hessian matrix
compute_hessian_for_R <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  R <- object$coefficients$R
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  R_vec <- as.vector(t(R))
  n_params <- length(R_vec)
  
  neg_loglik <- function(R_vec_temp) {
    R_temp <- matrix(R_vec_temp, nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      compute_log_likelihood(
        R = R_temp, beta = beta, Sigma = Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W,
        k = k, n = n, verbose = 0
      )
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  H <- compute_hessian_numerical(R_vec, neg_loglik, eps = eps)
  
  return(H)
}

#' Hessian computation for the Λ matrix (for MSEM)
#' 
#' @param object multivar_spatial object
#' @param eps step size for numerical differentiation
#' @param verbose verbose output
#' @return Hessian matrix
compute_hessian_for_Lambda <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  Lambda_mat <- object$coefficients$Lambda
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  Lambda_vec <- as.vector(t(Lambda_mat))
  n_params <- length(Lambda_vec)
  
  neg_loglik <- function(Lambda_vec_temp) {
    Lambda_temp <- matrix(Lambda_vec_temp, nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      if (exists("compute_log_likelihood_msem")) {
        compute_log_likelihood_msem(
          Lambda_mat = Lambda_temp, beta = beta, Sigma = Sigma,
          y = y, X = X, W = W, eigen_W = eigen_W,
          k = k, n = n, verbose = 0
        )
      } else {
        -Inf
      }
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  H <- compute_hessian_numerical(Lambda_vec, neg_loglik, eps = eps)
  
  return(H)
}

#' Hessian computation for the R, Λ matrices (for MGNS)
#' 
#' @param object multivar_spatial object
#' @param eps step size for numerical differentiation
#' @param verbose verbose output
#' @return Hessian matrix
compute_hessian_for_RLambda <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  R <- object$coefficients$R
  Lambda_mat <- object$coefficients$Lambda
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  R_vec <- as.vector(t(R))
  Lambda_vec <- as.vector(t(Lambda_mat))
  param_vec <- c(R_vec, Lambda_vec)
  
  neg_loglik <- function(param_temp) {
    n_R <- k^2
    R_temp <- matrix(param_temp[1:n_R], nrow = k, ncol = k, byrow = TRUE)
    Lambda_temp <- matrix(param_temp[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      if (exists("compute_log_likelihood_mgns")) {
        compute_log_likelihood_mgns(
          R = R_temp, Lambda_mat = Lambda_temp, beta = beta, Sigma = Sigma,
          y = y, X = X, W = W, eigen_W = eigen_W,
          k = k, n = n, verbose = 0
        )
      } else {
        -Inf
      }
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  H <- compute_hessian_numerical(param_vec, neg_loglik, eps = eps)
  
  return(H)
}

#' Build the coefficient table
#' 
#' @param object multivar_spatial object
#' @param se_beta standard errors of β
#' @param se_R standard-error matrix of R
#' @param se_Lambda standard-error matrix of Λ
#' @return data.frame
create_coefficient_table <- function(object, se_beta, se_R, se_Lambda) {
  
  k <- object$data_info$k
  rows <- list()
  
  if (!is.null(object$coefficients$R)) {
    R <- object$coefficients$R
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("R[%d,%d]", i, j)
        est <- R[i, j]
        se <- if (!is.null(se_R)) se_R[i, j] else NA
        z <- if (!is.na(se) && se > 0) est / se else NA
        p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
        sig <- get_signif_code(p)
        
        rows[[length(rows) + 1]] <- data.frame(
          parameter = param_name,
          estimate = est,
          std_error = se,
          z_value = z,
          p_value = p,
          signif = sig,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (!is.null(object$coefficients$Lambda)) {
    Lambda_mat <- object$coefficients$Lambda
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("Lambda[%d,%d]", i, j)
        est <- Lambda_mat[i, j]
        se <- if (!is.null(se_Lambda)) se_Lambda[i, j] else NA
        z <- if (!is.na(se) && se > 0) est / se else NA
        p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
        sig <- get_signif_code(p)
        
        rows[[length(rows) + 1]] <- data.frame(
          parameter = param_name,
          estimate = est,
          std_error = se,
          z_value = z,
          p_value = p,
          signif = sig,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  beta <- object$coefficients$beta
  if (!is.null(beta) && !is.null(se_beta)) {
    beta_names <- names(beta)
    if (is.null(beta_names)) {
      beta_names <- paste0("beta[", 1:length(beta), "]")
    }
    
    for (i in seq_along(beta)) {
      est <- beta[i]
      se <- se_beta[i]
      z <- if (!is.na(se) && se > 0) est / se else NA
      p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
      sig <- get_signif_code(p)
      
      rows[[length(rows) + 1]] <- data.frame(
        parameter = beta_names[i],
        estimate = est,
        std_error = se,
        z_value = z,
        p_value = p,
        signif = sig,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0) {
    return(NULL)
  }
  
  coef_table <- do.call(rbind, rows)
  rownames(coef_table) <- NULL
  
  return(coef_table)
}
