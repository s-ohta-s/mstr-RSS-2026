# experiment_output_functions.r
# Table-generation functions for the numerical experiments: unified parameter
# extraction, comparison / model-selection / bias / SE-comparison tables.

# Unified parameter extraction

#' Extract all parameters from an estimation-result object in a unified format
#'
#' Supports all model types (full/diagonal/VARX/OLS).
#' Parameter names are output in a fixed order so rows line up across models.
#'
#' @param result multivar_spatial object
#' @param model_id model ID string (e.g., "1111", "d0dd")
#' @param gamma γ value used
#' @param pAIC pAIC value (equals AIC if NULL)
#' @param pBIC pBIC value (equals BIC if NULL)
#' @param d_eff effective degrees of freedom (num_params if NULL)
#' @return data.frame (parameter, estimate, se_psi, se_hessian, signif_psi, signif_hessian)
extract_params_uniform <- function(result, model_id, gamma = 0,
                                    pAIC = NULL, pBIC = NULL, d_eff = NULL) {
  
  k <- result$data_info$k
  n <- result$data_info$n
  model_type <- result$model_type
  y_vars <- result$data_info$y_vars
  x_vars <- result$data_info$x_vars
  include_time_lag <- result$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE
  include_intercept <- result$data_info$include_intercept
  if (is.null(include_intercept)) include_intercept <- TRUE
  
  rows <- list()
  
  safe_scalar <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    x <- x[1]
    if (is.na(x)) return(NA_real_)
    return(as.numeric(x))
  }
  
  safe_signif <- function(est, se) {
    est <- safe_scalar(est)
    se  <- safe_scalar(se)
    if (is.na(est) || is.na(se) || se <= 0) return("")
    p <- 2 * pnorm(-abs(est / se))
    return(get_signif_code(p))
  }
  
  add_row <- function(param, est, se_p = NA_real_, se_h = NA_real_) {
    est  <- safe_scalar(est)
    se_p <- safe_scalar(se_p)
    se_h <- safe_scalar(se_h)
    rows[[length(rows) + 1]] <<- data.frame(
      parameter      = as.character(param),
      estimate       = est,
      se_psi         = se_p,
      se_hessian     = se_h,
      signif_psi     = safe_signif(est, se_p),
      signif_hessian = safe_signif(est, se_h),
      stringsAsFactors = FALSE
    )
  }
  
  add_fit <- function(param, val) {
    val <- safe_scalar(val)
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = as.character(param), estimate = val,
      se_psi = NA_real_, se_hessian = NA_real_,
      signif_psi = "", signif_hessian = "",
      stringsAsFactors = FALSE
    )
  }
  
  is_full_model <- model_type %in% c("MSAR", "MSEM", "MGNS")
  is_diag_model <- model_type %in% c("IndSAR", "IndSEM", "IndGNS")

  se_beta_psi  <- result$std_errors$beta        # Ψ-based (add_inference)
  se_beta_hess <- result$std_errors$beta_full   # full Hessian (add_full_inference)
  se_R_hess <- result$std_errors$R_full
  se_Lambda_hess <- result$std_errors$Lambda_full

  # NOTE: se_psi / signif_psi are Ψ-based and defined only for the regression
  # coefficients β (incl. the time-lag block A). The R/Λ rows therefore leave
  # se_psi blank; their significance is judged with se_hessian (full Hessian).

  if (is_full_model && is.null(se_beta_psi) && !is.null(result$model_data$X)) {
    Lambda_for_vcov <- if (model_type %in% c("MSEM", "MGNS")) result$coefficients$Lambda else NULL
    W_for_vcov <- if (!is.null(Lambda_for_vcov) && !is.null(result$model_data)) result$model_data$W else NULL
    Psi <- tryCatch(
      compute_vcov_beta(result$model_data$X, result$coefficients$Sigma, k, n,
                        Lambda_mat = Lambda_for_vcov, W = W_for_vcov),
      error = function(e) NULL
    )
    if (!is.null(Psi)) se_beta_psi <- sqrt(pmax(diag(Psi), 0))
  }
  
  # (a) β regression coefficients
  beta_vec <- result$coefficients$beta  # full models (MSAR/MSEM/MGNS) only
  beta0    <- result$coefficients$beta0 # diagonal models / VARX / OLS

  if (!is.null(beta_vec)) {
    # Full model / VARX: extract from the unified β vector
    idx <- 1
    for (i in 1:k) {
      vn <- y_vars[i]
      if (include_intercept) {
        pname <- sprintf("beta_intercept_%s", vn)
        se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
        se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
        add_row(pname, beta_vec[idx], se_p, se_h)
        idx <- idx + 1
      }
      for (xn in x_vars[[i]]) {
        pname <- sprintf("beta_%s_%s", xn, vn)
        se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
        se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
        add_row(pname, beta_vec[idx], se_p, se_h)
        idx <- idx + 1
      }
    }
    
    # A matrix (time lag) — the tail of the β vector
    if (include_time_lag && idx <= length(beta_vec)) {
      for (i in 1:k) {
        for (j in 1:k) {
          pname <- sprintf("A[%d,%d]", i, j)
          se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
          se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
          add_row(pname, beta_vec[idx], se_p, se_h)
          idx <- idx + 1
        }
      }
    } else {
      # No time lag, or VARX (A is stored separately in coefficients$A)
      A_mat <- result$coefficients$A
      alpha_mat <- result$coefficients$alpha
      A_use <- if (!is.null(A_mat)) A_mat else alpha_mat
      for (i in 1:k) {
        for (j in 1:k) {
          pname <- sprintf("A[%d,%d]", i, j)
          est <- if (!is.null(A_use)) A_use[i, j] else NA
          add_row(pname, est, NA, NA)
        }
      }
    }
    
  } else if (!is.null(beta0)) {
    # Diagonal / OLS model: extract from the beta0 list (named numeric vectors).
    # Time-lag terms (e.g., "y1_lag") are output separately as the A matrix.
    lag_pattern <- "^y[0-9]+_lag$"

    for (i in 1:k) {
      vn <- y_vars[i]
      coef_i <- beta0[[vn]]
      if (is.null(coef_i)) coef_i <- beta0[[i]]
      if (is.null(coef_i)) {
        if (include_intercept) add_row(sprintf("beta_intercept_%s", vn), NA, NA, NA)
        for (xn in x_vars[[i]]) add_row(sprintf("beta_%s_%s", xn, vn), NA, NA, NA)
        next
      }
      
      cnames <- names(coef_i)
      im <- result$individual_models[[vn]]
      ct <- if (!is.null(im)) im$coefficients else NULL  # spatialreg coefficient table

      for (j in seq_along(coef_i)) {
        cn <- cnames[j]
        if (grepl(lag_pattern, cn)) next

        if (cn == "(Intercept)") {
          pname <- sprintf("beta_intercept_%s", vn)
        } else {
          pname <- sprintf("beta_%s_%s", cn, vn)
        }
        
        se_p <- NA
        if (!is.null(ct)) {
          match_row <- ct[ct$parameter == cn, ]
          if (nrow(match_row) > 0) se_p <- match_row$std_error[1]
        }
        
        add_row(pname, coef_i[j], se_p, NA)
      }
    }
    
    # A matrix: VARX -> coefficients$A/$alpha, diagonal models -> lag terms in beta0
    A_diag <- matrix(NA, k, k)
    A_se   <- matrix(NA, k, k)

    A_from_result <- result$coefficients$A
    if (is.null(A_from_result)) A_from_result <- result$coefficients$alpha
    
    if (!is.null(A_from_result)) {
      A_diag <- A_from_result
      if (!is.null(result$std_errors$A)) {
        A_se <- result$std_errors$A
      }
    } else {
      # Diagonal models: extract diagonal elements from lag terms in beta0
      for (i in 1:k) {
        vn <- y_vars[i]
        coef_i <- beta0[[vn]]
        if (is.null(coef_i)) coef_i <- beta0[[i]]
        if (is.null(coef_i)) next
        
        im <- result$individual_models[[vn]]
        ct <- if (!is.null(im)) im$coefficients else NULL
        
        for (j in 1:k) {
          lag_name <- sprintf("y%d_lag", j)
          if (lag_name %in% names(coef_i)) {
            A_diag[i, j] <- coef_i[[lag_name]]
            if (!is.null(ct)) {
              match_row <- ct[ct$parameter == lag_name, ]
              if (nrow(match_row) > 0) A_se[i, j] <- match_row$std_error[1]
            }
          }
        }
      }
    }
    
    for (i in 1:k) for (j in 1:k) {
      add_row(sprintf("A[%d,%d]", i, j), A_diag[i, j], A_se[i, j], NA)
    }
    
  } else {
    # Both beta and beta0 are NULL
    for (i in 1:k) {
      vn <- y_vars[i]
      if (include_intercept) add_row(sprintf("beta_intercept_%s", vn), NA, NA, NA)
      for (xn in x_vars[[i]]) add_row(sprintf("beta_%s_%s", xn, vn), NA, NA, NA)
    }
    for (i in 1:k) for (j in 1:k) add_row(sprintf("A[%d,%d]", i, j), NA, NA, NA)
  }
  
  # (b) R matrix
  R_mat <- result$coefficients$R
  for (i in 1:k) for (j in 1:k) {
    est  <- if (!is.null(R_mat)) R_mat[i, j] else NA
    se_h <- if (!is.null(se_R_hess)) se_R_hess[i, j] else NA
    add_row(sprintf("R[%d,%d]", i, j), est, NA, se_h)
  }
  
  # (c) Λ matrix
  Lambda_mat <- result$coefficients$Lambda
  for (i in 1:k) for (j in 1:k) {
    est  <- if (!is.null(Lambda_mat)) Lambda_mat[i, j] else NA
    se_h <- if (!is.null(se_Lambda_hess)) se_Lambda_hess[i, j] else NA
    add_row(sprintf("Lambda[%d,%d]", i, j), est, NA, se_h)
  }
  
  # (d) Σ matrix (upper triangle only)
  Sigma_mat <- result$coefficients$Sigma
  for (i in 1:k) for (j in i:k) {
    est <- if (!is.null(Sigma_mat)) Sigma_mat[i, j] else NA
    add_row(sprintf("Sigma[%d,%d]", i, j), est, NA, NA)
  }
  
  # (e) Goodness-of-fit metrics
  add_fit("loglik",     result$fit$loglik)
  add_fit("AIC",        result$fit$AIC)
  add_fit("BIC",        result$fit$BIC)
  add_fit("pAIC",       if (!is.null(pAIC)) pAIC else result$fit$AIC)
  add_fit("pBIC",       if (!is.null(pBIC)) pBIC else result$fit$BIC)
  add_fit("d_eff",      if (!is.null(d_eff)) d_eff else result$fit$num_params)
  add_fit("gamma",      gamma)
  add_fit("n_params",   result$fit$num_params)
  
  # pseudo_R2 / adj_R2: unified definition R²_pseudo,k = corr(y_k, ŷ_k)² with
  # the trend prediction ŷ (d = 0): MSAR/MGNS ŷ = (I−R̂⊗W)⁻¹Xβ̂,
  # MSEM/VARX/OLS ŷ = Xβ̂ — computed in build_result_object() for every model type
  add_fit("pseudo_R2",  result$fit$R2)
  add_fit("adj_R2",     result$fit$R2_adj)

  r2_ind <- result$fit$R2_individual
  if (!is.null(r2_ind)) {
    for (i in 1:k) {
      vn <- y_vars[i]
      r2_val <- if (!is.null(names(r2_ind)) && vn %in% names(r2_ind)) r2_ind[[vn]] else r2_ind[i]
      add_fit(sprintf("pseudo_R2_%s", vn), r2_val)
    }
  }
  
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  return(df)
}


# Cross-model parameter comparison table

#' Collect all models' estimation results into a wide comparison table
#'
#' @param results_list named list (model ID -> result object)
#' @param gamma_info  named list (model ID -> list(gamma, pAIC, pBIC, d_eff))
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return data.frame (invisible)
build_comparison_table <- function(
  results_list,
  gamma_info = NULL,
  output_file = "comparison_table.csv",
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== Building the cross-model parameter comparison table ===\n")
  
  model_ids <- names(results_list)
  n_models <- length(results_list)

  all_dfs <- list()
  for (mid in model_ids) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    if (verbose) cat(sprintf("  Processing %s...\n", mid))
    
    df <- extract_params_uniform(
      result, mid,
      gamma = gi$gamma,
      pAIC  = gi$pAIC,
      pBIC  = gi$pBIC,
      d_eff = gi$d_eff
    )
    
    colnames(df)[2:6] <- paste0(colnames(df)[2:6], ".", mid)
    all_dfs[[mid]] <- df
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
    cat(sprintf("  Output: %s (%d rows × %d models)\n", output_file, nrow(merged), n_models))
  }
  
  invisible(merged)
}


# Model-selection summary table

#' Build a summary table condensing model-selection metrics to one row per model
#'
#' @param results_list named list (model ID -> result object)
#' @param gamma_info_pAIC γ info for the pAIC criterion
#' @param gamma_info_pBIC γ info for the pBIC criterion
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return data.frame (invisible)
build_model_selection_table <- function(
  results_pAIC,
  results_pBIC = NULL,
  gamma_info_pAIC = NULL,
  gamma_info_pBIC = NULL,
  output_file = "model_selection_summary.csv",
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== Building the model-selection summary table ===\n")
  
  model_ids <- names(results_pAIC)
  
  rows <- list()
  for (mid in model_ids) {
    r_aic <- results_pAIC[[mid]]
    gi_aic <- if (!is.null(gamma_info_pAIC[[mid]])) gamma_info_pAIC[[mid]] else list(gamma = 0)
    gi_bic <- if (!is.null(gamma_info_pBIC) && !is.null(gamma_info_pBIC[[mid]])) gamma_info_pBIC[[mid]] else gi_aic
    
    r_bic <- if (!is.null(results_pBIC[[mid]])) results_pBIC[[mid]] else r_aic

    r2_val <- r_aic$fit$R2
    r2_adj_val <- r_aic$fit$R2_adj
    
    rows[[mid]] <- data.frame(
      model_id       = mid,
      n_params       = r_aic$fit$num_params,
      d_eff_pAIC     = if (!is.null(gi_aic$d_eff)) gi_aic$d_eff else r_aic$fit$num_params,
      d_eff_pBIC     = if (!is.null(gi_bic$d_eff)) gi_bic$d_eff else r_aic$fit$num_params,
      loglik_pAIC    = r_aic$fit$loglik,
      loglik_pBIC    = r_bic$fit$loglik,
      AIC            = r_aic$fit$AIC,
      pAIC           = if (!is.null(gi_aic$pAIC)) gi_aic$pAIC else r_aic$fit$AIC,
      BIC            = r_aic$fit$BIC,
      pBIC           = if (!is.null(gi_bic$pBIC)) gi_bic$pBIC else r_aic$fit$BIC,
      gamma_pAIC     = gi_aic$gamma,
      gamma_pBIC     = gi_bic$gamma,
      R2_mean        = if (!is.null(r2_val)) r2_val else NA,
      adj_R2_mean    = if (!is.null(r2_adj_val)) r2_adj_val else NA,
      stringsAsFactors = FALSE
    )
  }
  
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL

  tbl$delta_AIC  <- tbl$AIC  - min(tbl$AIC,  na.rm = TRUE)
  tbl$delta_pAIC <- tbl$pAIC - min(tbl$pAIC, na.rm = TRUE)
  tbl$delta_BIC  <- tbl$BIC  - min(tbl$BIC,  na.rm = TRUE)
  tbl$delta_pBIC <- tbl$pBIC - min(tbl$pBIC, na.rm = TRUE)

  tbl <- tbl[order(tbl$AIC), ]
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d models)\n", output_file, nrow(tbl)))
    cat(sprintf("  Best AIC:  %s (%.2f)\n", tbl$model_id[1], tbl$AIC[1]))
    cat(sprintf("  Best pAIC: %s (%.2f)\n",
                tbl$model_id[which.min(tbl$pAIC)], min(tbl$pAIC)))
    cat(sprintf("  Best BIC:  %s (%.2f)\n",
                tbl$model_id[which.min(tbl$BIC)], min(tbl$BIC)))
    cat(sprintf("  Best pBIC: %s (%.2f)\n",
                tbl$model_id[which.min(tbl$pBIC)], min(tbl$pBIC)))
  }
  
  invisible(tbl)
}


# Comparison table against true values (bias evaluation)

#' Build a table evaluating the bias between true and estimated values
#'
#' @param results_list named list (model ID -> result object)
#' @param true_params TRUE_PARAMS list
#' @param gamma_info γ info (passed to extract_params_uniform)
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return data.frame (invisible)
build_bias_table <- function(
  results_list,
  true_params,
  gamma_info = NULL,
  output_file = "bias_comparison.csv",
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== Building the comparison table against true values (bias evaluation) ===\n")
  
  # Vectorize the true values (parameter name -> value mapping)
  tv <- list()
  tv[["beta_intercept_y1"]] <- true_params$beta_intercept_y1
  tv[["beta_intercept_y2"]] <- true_params$beta_intercept_y2
  tv[["beta_x_common1_y1"]] <- true_params$beta_common1_y1
  tv[["beta_x_common1_y2"]] <- true_params$beta_common1_y2
  tv[["beta_x_common2_y1"]] <- true_params$beta_common2_y1
  tv[["beta_x_common2_y2"]] <- true_params$beta_common2_y2
  tv[["beta_x_specific1_1_y1"]] <- true_params$beta_specific1_1
  tv[["beta_x_specific2_1_y2"]] <- true_params$beta_specific2_1
  for (i in 1:2) for (j in 1:2) tv[[sprintf("R[%d,%d]", i, j)]] <- true_params$R[i, j]
  for (i in 1:2) for (j in 1:2) tv[[sprintf("Lambda[%d,%d]", i, j)]] <- true_params$Lambda[i, j]
  for (i in 1:2) for (j in 1:2) tv[[sprintf("A[%d,%d]", i, j)]] <- true_params$A[i, j]
  tv[["Sigma[1,1]"]] <- true_params$Sigma[1, 1]
  tv[["Sigma[1,2]"]] <- true_params$Sigma[1, 2]
  tv[["Sigma[2,2]"]] <- true_params$Sigma[2, 2]
  
  param_names <- names(tv)

  model_ids <- names(results_list)

  tbl <- data.frame(parameter = param_names, true = unlist(tv),
                    stringsAsFactors = FALSE)
  rownames(tbl) <- NULL
  
  for (mid in model_ids) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    df <- extract_params_uniform(result, mid, gamma = gi$gamma,
                                  pAIC = gi$pAIC, pBIC = gi$pBIC, d_eff = gi$d_eff)
    
    est_col <- rep(NA_real_, length(param_names))
    bias_col <- rep(NA_real_, length(param_names))
    
    for (p_idx in seq_along(param_names)) {
      pn <- param_names[p_idx]
      match_idx <- which(df$parameter == pn)
      if (length(match_idx) == 0) {
        # Handle naming variation (e.g., beta_x_common1_y1 vs beta_common1_y1)
        pn_alt <- gsub("^beta_x_", "beta_", pn)
        match_idx <- which(df$parameter == pn_alt)
      }
      if (length(match_idx) > 0) {
        val <- df$estimate[match_idx[1]]
        if (!is.na(val)) {
          est_col[p_idx] <- val
          bias_col[p_idx] <- val - tv[[pn]]
        }
      }
    }
    
    tbl[[paste0("est_", mid)]] <- est_col
    tbl[[paste0("bias_", mid)]] <- bias_col
  }
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d parameters × %d models)\n",
                output_file, length(param_names), length(model_ids)))
  }
  
  invisible(tbl)
}


# Ψ vs Hessian SE comparison table

#' Build a table comparing the Ψ-based and Hessian-based SE for full models
#'
#' @param results_list list of full-model estimation results (add_full_inference applied)
#' @param gamma_info γ info
#' @param criterion "pAIC" or "pBIC"
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return data.frame (invisible)
build_se_comparison_table <- function(
  results_list,
  gamma_info = NULL,
  criterion = "pAIC",
  output_file = "se_comparison.csv",
  verbose = FALSE
) {
  
  if (verbose) cat(sprintf("\n=== Building the Ψ vs Hessian SE comparison table (%s criterion) ===\n", criterion))
  
  all_rows <- list()
  
  for (mid in names(results_list)) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    df <- extract_params_uniform(result, mid, gamma = gi$gamma,
                                  pAIC = gi$pAIC, pBIC = gi$pBIC, d_eff = gi$d_eff)
    
    # Exclude goodness-of-fit rows
    fit_params <- c("loglik", "AIC", "BIC", "pAIC", "pBIC", "d_eff",
                    "gamma", "n_params", "pseudo_R2", "adj_R2",
                    "pseudo_R2_y1", "pseudo_R2_y2")
    df <- df[!df$parameter %in% fit_params, ]

    # SE_ratio is defined for β coefficients only
    df$SE_ratio <- ifelse(
      !is.na(df$se_psi) & !is.na(df$se_hessian) & df$se_psi > 0,
      df$se_hessian / df$se_psi,
      NA
    )
    
    df$model_id <- mid
    df$gamma <- gi$gamma
    df$criterion <- criterion
    
    all_rows[[mid]] <- df
  }
  
  tbl <- do.call(rbind, all_rows)
  rownames(tbl) <- NULL

  tbl <- tbl[, c("parameter", "model_id", "gamma", "criterion", "estimate",
                  "se_psi", "se_hessian", "signif_psi", "signif_hessian", "SE_ratio")]
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d rows)\n", output_file, nrow(tbl)))
  }
  
  invisible(tbl)
}


# CSV export of the γ-search trajectory table

#' Export the result of compare_gamma() to CSV
#'
#' @param gamma_result return value of compare_gamma()
#' @param model_id model ID string
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return summary data.frame (invisible)
export_gamma_search_csv <- function(gamma_result, model_id,
                                     output_file, verbose = FALSE) {
  
  summary_df <- gamma_result$summary
  results <- gamma_result$results

  k <- NULL
  extra_cols <- list()
  
  for (i in seq_len(nrow(summary_df))) {
    g_str <- as.character(summary_df$gamma[i])
    res <- results[[g_str]]
    
    if (is.null(res)) {
      extra_cols[[i]] <- data.frame(
        R11 = NA, R12 = NA, R21 = NA, R22 = NA,
        L11 = NA, L12 = NA, L21 = NA, L22 = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (is.null(k)) k <- res$data_info$k
    
    R_mat <- res$coefficients$R
    Lambda_mat <- res$coefficients$Lambda
    
    row_data <- data.frame(
      R11 = if (!is.null(R_mat)) R_mat[1, 1] else NA,
      R12 = if (!is.null(R_mat)) R_mat[1, 2] else NA,
      R21 = if (!is.null(R_mat)) R_mat[2, 1] else NA,
      R22 = if (!is.null(R_mat)) R_mat[2, 2] else NA,
      L11 = if (!is.null(Lambda_mat)) Lambda_mat[1, 1] else NA,
      L12 = if (!is.null(Lambda_mat)) Lambda_mat[1, 2] else NA,
      L21 = if (!is.null(Lambda_mat)) Lambda_mat[2, 1] else NA,
      L22 = if (!is.null(Lambda_mat)) Lambda_mat[2, 2] else NA,
      stringsAsFactors = FALSE
    )
    extra_cols[[i]] <- row_data
  }
  
  extra_df <- do.call(rbind, extra_cols)
  out_df <- cbind(summary_df[, c("gamma", "loglik", "penalty",
                                   "AIC", "BIC", "df_eff",
                                   "GIC_AIC", "GIC_BIC")],
                  extra_df)
  
  write.csv(out_df, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  γ-search trajectory: %s -> %s (%d rows)\n", model_id, output_file, nrow(out_df)))
  }
  
  invisible(out_df)
}


# Optimal-γ selection summary

#' Build a summary of the optimal γ from all γ-search results
#'
#' @param gamma_searches named list (model ID -> compare_gamma result)
#' @param output_file output CSV file path
#' @param verbose verbose output
#' @return data.frame (model_id, gamma_pAIC, pAIC, d_eff_pAIC, gamma_pBIC, pBIC, d_eff_pBIC)
build_gamma_optimal_summary <- function(
  gamma_searches,
  output_file = "gamma_optimal_summary.csv",
  verbose = FALSE
) {
  
  if (verbose) cat("\n=== Optimal-γ selection summary ===\n")
  
  rows <- list()
  for (mid in names(gamma_searches)) {
    gs <- gamma_searches[[mid]]
    s <- gs$summary

    # Use GIC_AIC if computed, otherwise AIC
    if ("GIC_AIC" %in% colnames(s) && any(!is.na(s$GIC_AIC))) {
      idx_aic <- which.min(s$GIC_AIC)
      paic_val <- s$GIC_AIC[idx_aic]
    } else {
      idx_aic <- which.min(s$AIC)
      paic_val <- s$AIC[idx_aic]
    }
    
    if ("GIC_BIC" %in% colnames(s) && any(!is.na(s$GIC_BIC))) {
      idx_bic <- which.min(s$GIC_BIC)
      pbic_val <- s$GIC_BIC[idx_bic]
    } else {
      idx_bic <- which.min(s$BIC)
      pbic_val <- s$BIC[idx_bic]
    }
    
    rows[[mid]] <- data.frame(
      model_id    = mid,
      gamma_pAIC  = s$gamma[idx_aic],
      pAIC        = paic_val,
      d_eff_pAIC  = if ("df_eff" %in% colnames(s)) s$df_eff[idx_aic] else NA,
      loglik_pAIC = s$loglik[idx_aic],
      gamma_pBIC  = s$gamma[idx_bic],
      pBIC        = pbic_val,
      d_eff_pBIC  = if ("df_eff" %in% colnames(s)) s$df_eff[idx_bic] else NA,
      loglik_pBIC = s$loglik[idx_bic],
      stringsAsFactors = FALSE
    )
    
    if (verbose) {
      cat(sprintf("  %s: γ*_pAIC=%.4g (pAIC=%.2f), γ*_pBIC=%.4g (pBIC=%.2f)\n",
                  mid, s$gamma[idx_aic], paic_val, s$gamma[idx_bic], pbic_val))
    }
  }
  
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) cat(sprintf("  Output: %s\n", output_file))
  
  invisible(tbl)
}
