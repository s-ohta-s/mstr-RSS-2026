# initial_estimation_extended.r
# Extended version of the S1 initial estimation

# 1. Extended initial estimation for MSAR

#' S1: initial estimation via per-variable MSAR models (extended version)
#' 
#' Returns detailed estimation results including significance info and R²
#' 
#' @param data_list output of prepare_data_extended()
#' @param verbose verbosity level (0: none, 1: basic, 2: detailed)
#' @return list(R_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sly_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation via per-variable MSAR models (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("the spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  R_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d:\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    } else {
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[2:(1+length(xi_names))] <- xi_names
      }
    }
    
    tryCatch({
      if (verbose >= 2) {
        cat("  Calling lagsarlm()...\n")
      }
      
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::lagsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::lagsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      R_init[i, i] <- model$rho
      Sigma_init[i, i] <- model$s2
      
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      s <- summary(model)
      
      # === Significance of the spatial parameters ===
      rho <- model$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      rho_signif <- get_signif_code(rho_p)
      
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # === Significance of the regression coefficients ===
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
      
      # === Goodness-of-fit metrics ===
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      # Number of parameters (regression coefficients + ρ + σ²)
      num_params <- length(coef(model)) + 1  #  add σ² only (the spatial parameter is included in coef)
      
      # === Four kinds of R² ===
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      # === Save the per-variable estimation results ===
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "MSAR",
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = rho_signif,
          lambda = NULL,
          lambda_se = NULL,
          lambda_z = NULL,
          lambda_p = NULL,
          lambda_signif = NULL,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        
        sigma2 = model$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      # === Detailed output ===
      if (verbose >= 1) {
        cat("  [Spatial correlation parameters]\n")
        cat(sprintf("    ρ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, rho, rho_se, rho_z, rho_p, rho_signif))
        
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "reject ρ=0", "fail to reject ρ=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R²(corr)       = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    Log-likelihood = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: lagsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed. Stopping execution.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial R matrix:\n")
    print(round(R_init, 4))
    cat("\nInitial Σ matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    R_init = R_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}

# 2. Extended initial estimation for MSEM

#' S1: initial estimation via per-variable MSEM models (extended version)
#' 
#' @param data_list output of prepare_data_extended()
#' @param verbose verbosity level
#' @return list(T_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sem_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation via per-variable MSEM models (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("the spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  T_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d:\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    }
    
    tryCatch({
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::errorsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::errorsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      T_init[i, i] <- model$lambda
      Sigma_init[i, i] <- model$s2
      
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      s <- summary(model)
      
      lambda <- model$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      lambda_signif <- get_signif_code(lambda_p)
      
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
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
      
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      num_params <- length(coef(model)) + 1  #  add σ² only (the spatial parameter is included in coef)
      
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "MSEM",
        
        spatial_params = list(
          rho = NULL,
          rho_se = NULL,
          rho_z = NULL,
          rho_p = NULL,
          rho_signif = NULL,
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = lambda_signif,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      if (verbose >= 1) {
        cat("  [Spatial correlation parameters]\n")
        cat(sprintf("    λ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, lambda, lambda_se, lambda_z, lambda_p, lambda_signif))
        
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "reject λ=0", "fail to reject λ=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R²(corr)       = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    Log-likelihood = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: errorsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial Λ matrix:\n")
    print(round(T_init, 4))
    cat("\nInitial Σ matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    T_init = T_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}

# 3. Extended initial estimation for MGNS

#' S1: initial estimation via per-variable MGNS models (extended version)
#' 
#' Uses spatialreg's sacsarlm()
#' 
#' @param data_list output of prepare_data_extended()
#' @param verbose verbosity level
#' @return list(R_init, T_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sdem_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation via per-variable MGNS models (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("the spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  R_init <- matrix(0, k, k)
  T_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d:\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    }
    
    tryCatch({
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::sacsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::sacsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      R_init[i, i] <- model$rho
      T_init[i, i] <- model$lambda
      Sigma_init[i, i] <- model$s2
      
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      s <- summary(model)
      
      rho <- model$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      rho_signif <- get_signif_code(rho_p)
      
      lambda <- model$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      lambda_signif <- get_signif_code(lambda_p)
      
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
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
      
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      num_params <- length(coef(model)) + 1  #  add σ² only (the spatial parameter is included in coef)
      
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "MGNS",
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = rho_signif,
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = lambda_signif,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      if (verbose >= 1) {
        cat("  [Spatial correlation parameters]\n")
        cat(sprintf("    ρ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, rho, rho_se, rho_z, rho_p, rho_signif))
        cat(sprintf("    λ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, lambda, lambda_se, lambda_z, lambda_p, lambda_signif))
        
        if (!is.na(LR_stat)) {
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f\n", LR_stat, LR_p))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R²(corr)       = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    Log-likelihood = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: sacsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial R matrix:\n")
    print(round(R_init, 4))
    cat("\nInitial Λ matrix:\n")
    print(round(T_init, 4))
    cat("\nInitial Σ matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    R_init = R_init,
    T_init = T_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}

# 4. Helper function

#' Extract the i-th yi and Xi from the data (with time-lag option)
#' 
#' @param data_list data list
#' @param i variable index
#' @param include_time_lag whether to include a time lag
#' @return list(yi, Xi, y_lag)
extract_data_for_yi_extended <- function(data_list, i, include_time_lag = TRUE) {
  n <- data_list$n
  k <- data_list$k
  
  yi_idx <- ((i-1)*n + 1):(i*n)
  yi <- data_list$y[yi_idx]
  
  x_vars <- data_list$data_info$x_vars
  include_intercept <- data_list$data_info$include_intercept
  
  x_counts <- sapply(x_vars, length)
  if (include_intercept) {
    x_counts <- x_counts + 1
  }
  
  if (i == 1) {
    x_col_start <- 1
  } else {
    x_col_start <- sum(x_counts[1:(i-1)]) + 1
  }
  x_col_end <- x_col_start + x_counts[i] - 1
  
  x_cols <- x_col_start:x_col_end
  Xi <- data_list$X[yi_idx, x_cols, drop = FALSE]
  
  y_lag <- NULL
  if (include_time_lag && !is.null(data_list$y_lag)) {
    y_lag <- data_list$y_lag
    colnames(y_lag) <- paste0("y", 1:k, "_lag")
  }
  
  return(list(
    yi = yi,
    Xi = Xi,
    y_lag = y_lag
  ))
}
