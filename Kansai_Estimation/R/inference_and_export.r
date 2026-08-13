# inference_and_export.r
# Export estimation results to CSV:
#   (1) Parameter estimates (long format: parameter, estimate, std_error, z_value,
#       p_value, signif; goodness-of-fit metrics appended)
#       -> export_results_csv(), export_multiple_models_csv()
#   (2) Fitted values / residuals (same long format as the input data:
#       region, time, <var>_obs/_pred/_resid)
#       -> export_fitted_residuals_csv(), export_fitted_residuals_all_models()

#' Compute the sandwich variance of the spatial parameters
#'
#' For penalized estimation with γ > 0:
#'   Var(θ̂_γ) = (H + γI)⁻¹ H (H + γI)⁻¹
#' For γ = 0:
#'   Var(θ̂) = H⁻¹  (ordinary inverse information matrix)
#'
#' @param H negative Hessian of the non-penalized profile likelihood (k₁ × k₁)
#' @param gamma penalty strength
#' @return variance-covariance matrix (k₁ × k₁)
compute_sandwich_vcov <- function(H, gamma) {
  k1 <- nrow(H)

  if (gamma == 0) {
    vcov <- tryCatch({
      solve(H)
    }, error = function(e) {
      warning("Error inverting the Hessian. Applying regularization: ", e$message)
      eigen_H <- eigen(H, symmetric = TRUE)
      ridge <- max(abs(min(eigen_H$values)), 1e-6)
      solve(H + diag(ridge, k1))
    })
  } else {
    # Sandwich variance: (H + γI)⁻¹ H (H + γI)⁻¹
    H_pen <- H + gamma * diag(k1)
    
    H_pen_inv <- tryCatch({
      solve(H_pen)
    }, error = function(e) {
      warning("Error inverting (H + γI). Applying regularization: ", e$message)
      solve(H_pen + diag(1e-6, k1))
    })
    
    vcov <- H_pen_inv %*% H %*% H_pen_inv
  }

  diag_vcov <- diag(vcov)
  if (any(diag_vcov < 0)) {
    warning("The variance-covariance matrix has negative diagonal elements")
  }
  
  return(vcov)
}


#' Compute the asymptotic standard errors of the error covariance matrix Σ
#'
#' For the normal distribution: Var(σ̂_ij) = (σ_ii σ_jj + σ_ij²) / n
#'
#' @param Sigma k × k estimated error covariance matrix
#' @param n number of regions
#' @return k × k standard-error matrix
compute_se_sigma <- function(Sigma, n) {
  k <- nrow(Sigma)
  se_Sigma <- matrix(NA, k, k)
  
  for (i in 1:k) {
    for (j in 1:k) {
      var_sigma_ij <- (Sigma[i, i] * Sigma[j, j] + Sigma[i, j]^2) / n
      se_Sigma[i, j] <- sqrt(max(var_sigma_ij, 0))
    }
  }
  
  return(se_Sigma)
}


# Main function: export_results_csv

#' Export the estimation result to a CSV file
#'
#' @param result object returned by build_result_object()
#' @param output_file output CSV file path
#' @param gamma penalty strength (used when it cannot be obtained from result)
#' @param compute_gic_flag whether to compute GIC
#' @param verbose verbose output
#' @return the exported data.frame (invisible)
export_results_csv <- function(
  result,
  output_file = "estimation_results.csv",
  gamma = NULL,
  compute_gic_flag = TRUE,
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== CSV export start ===\n")

  k <- result$data_info$k
  n <- result$data_info$n
  model_type <- result$model_type

  if (is.null(gamma)) {
    if (!is.null(result$penalty$gamma)) {
      gamma <- result$penalty$gamma
    } else {
      gamma <- 0
    }
  }
  
  include_time_lag <- result$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE
  
  y_vars <- result$data_info$y_vars
  x_vars <- result$data_info$x_vars

  # Standard errors of β (analytical — exact regardless of γ)
  if (verbose) cat("  Computing the standard errors of β...\n")

  se_beta <- NULL
  Lambda_for_vcov <- if (model_type %in% c("MSEM", "MGNS")) result$coefficients$Lambda else NULL
  W_for_vcov <- if (!is.null(Lambda_for_vcov)) result$model_data$W else NULL
  
  vcov_beta <- tryCatch({
    compute_vcov_beta(result$model_data$X, result$coefficients$Sigma, k, n,
                      Lambda_mat = Lambda_for_vcov, W = W_for_vcov)
  }, error = function(e) {
    warning("Error computing the variance-covariance matrix of β: ", e$message)
    NULL
  })
  
  if (!is.null(vcov_beta)) {
    se_beta <- sqrt(pmax(diag(vcov_beta), 0))
  }

  # Standard errors of the spatial parameters (sandwich variance)
  se_R <- NULL
  se_Lambda <- NULL
  
  if (model_type %in% c("MSAR", "MSEM", "MGNS")) {
    if (verbose) cat("  Computing the standard errors of the spatial parameters (sandwich variance)...\n")
    
    H_spatial <- tryCatch({
      if (model_type == "MSAR") {
        compute_profile_hessian_msar(
          R_hat = result$coefficients$R,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      } else if (model_type == "MSEM") {
        compute_profile_hessian_msem(
          Lambda_hat = result$coefficients$Lambda,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      } else {  # MGNS
        compute_profile_hessian_mgns(
          R_hat = result$coefficients$R,
          Lambda_hat = result$coefficients$Lambda,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      }
    }, error = function(e) {
      warning("Error computing the Hessian: ", e$message)
      NULL
    })
    
    if (!is.null(H_spatial)) {
      vcov_spatial <- tryCatch({
        compute_sandwich_vcov(H_spatial, gamma)
      }, error = function(e) {
        warning("Error computing the sandwich variance: ", e$message)
        NULL
      })
      
      if (!is.null(vcov_spatial)) {
        se_spatial <- sqrt(pmax(diag(vcov_spatial), 0))
        
        if (model_type == "MSAR") {
          se_R <- matrix(se_spatial, nrow = k, ncol = k, byrow = TRUE)
        } else if (model_type == "MSEM") {
          se_Lambda <- matrix(se_spatial, nrow = k, ncol = k, byrow = TRUE)
        } else {  # MGNS
          n_R <- k^2
          se_R <- matrix(se_spatial[1:n_R], nrow = k, ncol = k, byrow = TRUE)
          se_Lambda <- matrix(se_spatial[(n_R + 1):(2 * n_R)], nrow = k, ncol = k, byrow = TRUE)
        }
      }
    }
  }

  if (verbose) cat("  Computing the standard errors of Σ...\n")
  se_Sigma <- compute_se_sigma(result$coefficients$Sigma, n)

  if (verbose) cat("  Assembling the CSV table...\n")
  
  rows <- list()
  
  add_row <- function(param_name, estimate, se = NA) {
    z <- if (!is.na(se) && se > 0) estimate / se else NA
    p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
    sig <- if (!is.na(p)) get_signif_code(p) else ""
    
    # Estimate plus significance mark (e.g., "0.234567***")
    if (!is.na(estimate)) {
      est_sig <- paste0(sprintf("%.6f", estimate), sig)
    } else {
      est_sig <- ""
    }
    
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = param_name,
      estimate_signif = est_sig,
      estimate = estimate,
      std_error = se,
      z_value = z,
      p_value = p,
      signif = sig,
      stringsAsFactors = FALSE
    )
  }
  
  # For goodness-of-fit metrics (no SE/z/p, no significance mark)
  add_fit_row <- function(param_name, value) {
    est_sig <- if (!is.na(value)) sprintf("%.6f", value) else ""
    
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = param_name,
      estimate_signif = est_sig,
      estimate = value,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      signif = "",
      stringsAsFactors = FALSE
    )
  }
  
  # (a) R matrix — MSAR, MGNS only
  if (model_type %in% c("MSAR", "MGNS") && !is.null(result$coefficients$R)) {
    R_mat <- result$coefficients$R
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("R[%d,%d]", i, j)
        se_val <- if (!is.null(se_R)) se_R[i, j] else NA
        add_row(param_name, R_mat[i, j], se_val)
      }
    }
  } else {
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("R[%d,%d]", i, j), NA)
      }
    }
  }

  # (b) Lambda matrix — MSEM, MGNS only
  if (model_type %in% c("MSEM", "MGNS") && !is.null(result$coefficients$Lambda)) {
    Lambda_mat <- result$coefficients$Lambda
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("Lambda[%d,%d]", i, j)
        se_val <- if (!is.null(se_Lambda)) se_Lambda[i, j] else NA
        add_row(param_name, Lambda_mat[i, j], se_val)
      }
    }
  } else {
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("Lambda[%d,%d]", i, j), NA)
      }
    }
  }
  
  # (c) A matrix (time lag)
  if (include_time_lag && !is.null(result$coefficients$alpha)) {
    alpha <- result$coefficients$alpha

    # The SE of α is taken from the corresponding part of the β vector
    se_alpha <- NULL
    if (!is.null(se_beta) && !is.null(result$coefficients$beta)) {
      n_beta0 <- 0
      for (i in 1:k) {
        n_beta0 <- n_beta0 + length(x_vars[[i]])
        if (result$data_info$include_intercept) n_beta0 <- n_beta0 + 1
      }
      alpha_idx <- (n_beta0 + 1):(n_beta0 + k^2)
      if (max(alpha_idx) <= length(se_beta)) {
        se_alpha_vec <- se_beta[alpha_idx]
        se_alpha <- matrix(se_alpha_vec, nrow = k, ncol = k, byrow = TRUE)
      }
    }
    
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("A[%d,%d]", i, j)
        se_val <- if (!is.null(se_alpha)) se_alpha[i, j] else NA
        add_row(param_name, alpha[i, j], se_val)
      }
    }
  } else {
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("A[%d,%d]", i, j), NA)
      }
    }
  }
  
  # (d) Sigma matrix
  Sigma_mat <- result$coefficients$Sigma
  for (i in 1:k) {
    for (j in 1:k) {
      param_name <- sprintf("Sigma[%d,%d]", i, j)
      se_val <- se_Sigma[i, j]
      add_row(param_name, Sigma_mat[i, j], se_val)
    }
  }
  
  # (e) β regression coefficients (per variable)
  beta_vec <- result$coefficients$beta
  beta0 <- result$coefficients$beta0
  include_intercept <- result$data_info$include_intercept

  beta_idx <- 1

  for (i in 1:k) {
    var_name <- y_vars[i]

    if (include_intercept) {
      param_name <- sprintf("beta_intercept_%s", var_name)
      est <- beta_vec[beta_idx]
      se_val <- if (!is.null(se_beta) && beta_idx <= length(se_beta)) se_beta[beta_idx] else NA
      add_row(param_name, est, se_val)
      beta_idx <- beta_idx + 1
    }

    for (x_name in x_vars[[i]]) {
      param_name <- sprintf("beta_%s_%s", x_name, var_name)
      est <- beta_vec[beta_idx]
      se_val <- if (!is.null(se_beta) && beta_idx <= length(se_beta)) se_beta[beta_idx] else NA
      add_row(param_name, est, se_val)
      beta_idx <- beta_idx + 1
    }
  }
  
  # Time-lag coefficients: the tail of the β vector is already output as the
  # A matrix, so skip those β indices
  if (include_time_lag) {
    beta_idx <- beta_idx + k^2
  }

  # (f) Goodness-of-fit metrics
  add_fit_row("AIC", result$fit$AIC)
  add_fit_row("BIC", result$fit$BIC)

  # GIC (for full models when compute_gic_flag is set)
  if (compute_gic_flag && model_type %in% c("MSAR", "MSEM", "MGNS")) {
    if (!is.null(H_spatial)) {
      # Number of non-spatial parameters k₂
      n_beta0 <- 0
      for (i in 1:k) {
        n_beta0 <- n_beta0 + length(x_vars[[i]])
        if (include_intercept) n_beta0 <- n_beta0 + 1
      }
      if (include_time_lag) n_beta0 <- n_beta0 + k^2
      k2 <- n_beta0 + k * (k + 1) / 2  # β + vech(Σ)
      n_obs <- k * n
      
      gic_result <- tryCatch({
        compute_gic(H = H_spatial, gamma = gamma,
                    loglik = result$fit$loglik, k2 = k2, n_obs = n_obs)
      }, error = function(e) {
        warning("Error computing GIC: ", e$message)
        NULL
      })
      
      if (!is.null(gic_result)) {
        add_fit_row("GIC_AIC", gic_result$GIC_AIC)
        add_fit_row("GIC_BIC", gic_result$GIC_BIC)
        add_fit_row("df_eff", gic_result$df_eff)
      } else {
        add_fit_row("GIC_AIC", NA)
        add_fit_row("GIC_BIC", NA)
        add_fit_row("df_eff", NA)
      }
    } else {
      add_fit_row("GIC_AIC", NA)
      add_fit_row("GIC_BIC", NA)
      add_fit_row("df_eff", NA)
    }
  }
  
  # Pseudo R² (per variable + mean)
  # Unified definition: R²_pseudo,k = corr(y_k, ŷ_k)² with the trend prediction
  # ŷ (d = 0): MSAR/MGNS ŷ = (I−R̂⊗W)⁻¹Xβ̂, MSEM/VARX/OLS ŷ = Xβ̂
  # (computed in build_result_object; identical to pred in the fitted_residuals CSVs)
  r2_ind <- result$fit$R2_individual
  if (!is.null(r2_ind)) {
    for (i in 1:k) {
      vn <- y_vars[i]
      r2_val <- if (!is.null(names(r2_ind)) && vn %in% names(r2_ind)) r2_ind[[vn]] else r2_ind[i]
      add_fit_row(sprintf("R2_%s", vn), r2_val)
    }
  } else {
    for (i in 1:k) add_fit_row(sprintf("R2_%s", y_vars[i]), NA)
  }
  add_fit_row("R2", if (!is.null(result$fit$R2)) result$fit$R2 else NA)

  add_fit_row("n_param", result$fit$num_params)
  add_fit_row("loglik", result$fit$loglik)
  add_fit_row("n_obs", result$fit$num_obs)
  add_fit_row("n_region", n)
  add_fit_row("gamma", gamma)

  output_df <- do.call(rbind, rows)
  rownames(output_df) <- NULL

  write.csv(output_df, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s\n", output_file))
    cat(sprintf("  Number of rows: %d\n", nrow(output_df)))
    cat("=== CSV export complete ===\n\n")
  }
  
  invisible(output_df)
}


# Batch export of multiple models

#' Export the estimation results of multiple models to a single CSV
#'
#' Place each model's results column-wise
#'
#' @param results_list named list (model name -> result object)
#' @param output_file output CSV file path
#' @param gammas γ for each model (taken from result if NULL)
#' @param verbose verbose output
#' @return the exported data.frame (invisible)
export_multiple_models_csv <- function(
  results_list,
  output_file = "model_comparison.csv",
  gammas = NULL,
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== Batch CSV export of multiple models start ===\n")
  
  model_names <- names(results_list)
  n_models <- length(results_list)

  all_dfs <- list()
  
  for (i in seq_along(results_list)) {
    model_name <- model_names[i]
    result <- results_list[[i]]
    g <- if (!is.null(gammas)) gammas[i] else NULL
    
    if (verbose) cat(sprintf("  Processing %s...\n", model_name))
    
    df <- export_results_csv(
      result,
      output_file = tempfile(),  # do not write individual files
      gamma = g,
      verbose = FALSE
    )
    
    colnames(df)[2:7] <- paste0(colnames(df)[2:7], ".", model_name)

    all_dfs[[model_name]] <- df
  }

  merged <- all_dfs[[1]]
  if (n_models > 1) {
    for (i in 2:n_models) {
      merged <- merge(merged, all_dfs[[i]], by = "parameter",
                      all = TRUE, sort = FALSE)
    }
  }

  param_order <- all_dfs[[1]]$parameter
  merged <- merged[match(param_order, merged$parameter), ]
  rownames(merged) <- NULL
  
  write.csv(merged, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s\n", output_file))
    cat(sprintf("  Number of models: %d, number of parameter rows: %d\n", n_models, nrow(merged)))
    cat("=== Batch CSV export of multiple models complete ===\n\n")
  }
  
  invisible(merged)
}


# Export fitted values / residuals

#' Export fitted values / residuals to CSV in the same long format as the input data
#'
#' Output columns: region, time, <yvar>_obs, <yvar>_pred, <yvar>_resid,
#' <yvar>_resid_reduced (per dependent variable).
#' Same row structure as the slice of the input data (region, time, y1, y2, ...) at time = the used time point.
#'
#' Definitions (consistent with the shrinkage-function formulation, d = 0):
#'   - pred          = trend prediction ŷ (result$fitted$trend)
#'       OLS/VARX/MSEM: ŷ = Xβ̂ ; MSAR/MGNS: ŷ = (I−R̂⊗W)⁻¹Xβ̂
#'   - resid         = innovation residual ε̂ (result$residuals$raw), i.e. the
#'       residual of the likelihood / Σ estimation:
#'       OLS/VARX: ε̂ = y − Xβ̂ ; MSAR: ε̂ = (I−R̂⊗W)y − Xβ̂ ;
#'       MSEM: ε̂ = (I−Λ̂⊗W)(y−Xβ̂) ; MGNS: ε̂ = (I−Λ̂⊗W){(I−R̂⊗W)y − Xβ̂}
#'   - resid_reduced = obs − pred (reduced-form residual)
#' NOTE: obs = pred + resid holds only for OLS/VARX; for spatial models
#' obs = pred + resid_reduced instead (resid is the innovation ε̂ ≠ obs − pred).
#'
#' @param result object returned by build_result_object()
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return the exported data.frame (invisible). NULL if it cannot be exported.
export_fitted_residuals_csv <- function(
  result,
  output_file = "fitted_residuals.csv",
  verbose = FALSE
) {

  n <- result$data_info$n
  k <- result$data_info$k
  y_vars <- result$data_info$y_vars

  obs   <- result$model_data$y
  resid <- result$residuals$raw
  pred  <- result$fitted$trend

  if (is.null(obs) || is.null(resid)) {
    warning(sprintf("Cannot export fitted values / residuals (observed values or residuals are NULL): model_type=%s",
                    result$model_type))
    return(invisible(NULL))
  }
  if (is.null(pred)) {
    warning(sprintf("fitted$trend is missing; falling back to pred = obs − resid: model_type=%s",
                    result$model_type))
    pred <- as.numeric(obs) - as.numeric(resid)
  }
  if (length(obs) != k * n || length(resid) != k * n || length(pred) != k * n) {
    warning(sprintf("The length of observed/fitted/residual vectors does not match kn=%d (obs=%d, pred=%d, resid=%d). Skipping.",
                    k * n, length(obs), length(pred), length(resid)))
    return(invisible(NULL))
  }

  obs   <- as.numeric(obs)
  pred  <- as.numeric(pred)
  resid <- as.numeric(resid)

  # Region IDs (use the real IDs if stored in data_list, otherwise 1..n)
  regions <- if (!is.null(result$data_list$regions)) {
    result$data_list$regions
  } else {
    seq_len(n)
  }

  time_point <- result$data_info$time_point_used
  if (is.null(time_point)) time_point <- NA

  # Long format (row = region)
  out <- data.frame(region = regions, time = time_point, stringsAsFactors = FALSE)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    vn <- y_vars[i]
    out[[paste0(vn, "_obs")]]   <- obs[idx]
    out[[paste0(vn, "_pred")]]  <- pred[idx]
    out[[paste0(vn, "_resid")]] <- resid[idx]
    out[[paste0(vn, "_resid_reduced")]] <- obs[idx] - pred[idx]
  }

  write.csv(out, file = output_file, row.names = FALSE)

  if (verbose) {
    cat(sprintf("  Exported fitted values / residuals: %s (%d rows × %d cols)\n",
                output_file, nrow(out), ncol(out)))
  }

  invisible(out)
}


#' Batch-export fitted values / residuals for multiple models (one file per model)
#'
#' Output file name: <output_dir>/fitted_residuals_<model ID>.csv
#'
#' @param results_list named list (model ID -> result object)
#' @param output_dir output directory (created if absent)
#' @param verbose verbose output
#' @return vector of exported file paths (invisible)
export_fitted_residuals_all_models <- function(
  results_list,
  output_dir,
  verbose = FALSE
) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  if (verbose) {
    cat(sprintf("\n=== Batch export of fitted values / residuals (%d models) ===\n", length(results_list)))
  }

  written <- character(0)
  for (mid in names(results_list)) {
    result <- results_list[[mid]]
    if (is.null(result)) next

    f <- file.path(output_dir, sprintf("fitted_residuals_%s.csv", mid))
    df <- tryCatch(
      export_fitted_residuals_csv(result, output_file = f, verbose = verbose),
      error = function(e) {
        warning(sprintf("%s error exporting fitted values / residuals: %s", mid, e$message))
        NULL
      }
    )
    if (!is.null(df)) written <- c(written, f)
  }

  if (verbose) {
    cat(sprintf("=== Fitted values / residuals export complete (%d files -> %s) ===\n\n",
                length(written), output_dir))
  }

  invisible(written)
}
