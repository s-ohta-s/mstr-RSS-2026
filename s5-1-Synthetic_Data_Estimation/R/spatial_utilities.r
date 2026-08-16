# spatial_utilities.r
# Core functions shared by the multivariate spatial regression models (MSAR/MSEM/MGNS)
# Auxiliary estimation functions: data preparation (with time-lag support),
# variance-covariance and Hessian, coefficient of determination, and extraction
# of significance from spatialreg models.

# 1. Data preparation function (with time-lag option)

#' Data preparation (extended version with a time-lag option)
#'
#' @param data_file path to the panel-data CSV file
#' @param weight_file path to the spatial-weight-matrix CSV file
#' @param y_vars vector of dependent variable names (e.g., c("y1", "y2"))
#' @param x_vars list of regressor names for each yi
#' @param time_var time variable name
#' @param time_point time point to use
#' @param region_var region variable name
#' @param include_intercept whether to include an intercept (default: TRUE)
#' @param include_time_lag whether to include time-lag variables (default: TRUE)
#' @param verbose verbose output
#'
#' @return list(y, X, W, W_sp, W_listw, y_lag, eigen_W, n, k, regions, p0, data_info)
#'
#' Results are memoised: repeated calls with identical arguments (same files,
#' unchanged on disk) return the cached list instead of re-reading the CSVs and
#' recomputing eigen(W) / mat2listw. This makes per-γ refits in compare_gamma
#' nearly free of data-preparation cost.
.mgnst_prep_cache <- new.env(parent = emptyenv())

prepare_data_extended <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  verbose = FALSE
) {

  # --- memoisation: key on all arguments + file modification times ---
  cache_key <- paste(
    data_file, weight_file,
    paste(y_vars, collapse = ","),
    paste(names(x_vars), sapply(x_vars, paste, collapse = ","), collapse = ";"),
    time_var, ifelse(is.null(time_point), "NULL", time_point), region_var,
    include_intercept, include_time_lag,
    as.numeric(file.mtime(data_file)), as.numeric(file.mtime(weight_file)),
    sep = "|"
  )
  cached <- .mgnst_prep_cache[[cache_key]]
  if (!is.null(cached)) {
    if (verbose) cat("=== Data preparation: cache hit ===\n")
    return(cached)
  }

  if (verbose) cat("=== Data preparation start (extended) ===\n")

  if (verbose) cat("Loading panel data...\n")
  data <- read.csv(data_file, stringsAsFactors = FALSE)

  if (verbose) cat("Loading the spatial weight matrix...\n")
  W_raw <- read.csv(weight_file, header = TRUE, stringsAsFactors = FALSE)

  # Determine whether the first column is a region ID
  first_col <- W_raw[, 1]
  is_id_column <- !is.numeric(first_col) || all(first_col == 1:nrow(W_raw))

  if (is_id_column) {
    W <- as.matrix(W_raw[, -1])
    if (verbose) cat("Excluded the first column as a region ID\n")
  } else {
    W <- as.matrix(W_raw)
  }

  dimnames(W) <- NULL

  if (verbose) {
    cat(sprintf("Spatial weight matrix dimensions: %d×%d\n", nrow(W), ncol(W)))
  }

  if (!all(c(time_var, region_var, y_vars) %in% colnames(data))) {
    stop("Required columns not found")
  }

  k <- length(y_vars)
  regions <- sort(unique(data[[region_var]]))
  n <- length(regions)

  if (verbose) {
    cat(sprintf("Number of variables k = %d\n", k))
    cat(sprintf("Number of regions n = %d\n", n))
    cat(sprintf("Include time lag: %s\n", ifelse(include_time_lag, "yes", "no")))
  }

  if (nrow(W) != n || ncol(W) != n) {
    stop(sprintf("Spatial weight matrix size (%dx%d) does not match the number of regions (%d)",
                 nrow(W), ncol(W), n))
  }

  row_sums <- rowSums(W)
  if (!all(abs(row_sums - 1) < 1e-6)) {
    warning("The spatial weight matrix may not be row-standardized")
  }

  times <- sort(unique(data[[time_var]]))
  if (is.null(time_point)) {
    time_point <- max(times)
    if (verbose) cat(sprintf("time_point not specified; using the latest time point %d\n", time_point))
  }

  if (!time_point %in% times) {
    stop(sprintf("The specified time point %d does not exist in the data", time_point))
  }

  time_lag <- NULL
  y_lag <- NULL

  if (include_time_lag) {
    time_lag <- time_point - 1
    if (!time_lag %in% times) {
      stop(sprintf("The lag time point %d does not exist in the data. At least two time points are required.", time_lag))
    }

    if (verbose) cat(sprintf("Time point used: t=%d, lag time point: t=%d\n", time_point, time_lag))

    data_lag <- data[data[[time_var]] == time_lag, ]
    data_lag <- data_lag[order(data_lag[[region_var]]), ]

    if (nrow(data_lag) != n) {
      stop("The lag-time-point data does not contain all regions")
    }

    y_lag <- matrix(NA, nrow = n, ncol = k)
    for (i in 1:k) {
      y_lag[, i] <- data_lag[[y_vars[i]]]
    }
  } else {
    if (verbose) cat(sprintf("Time point used: t=%d (no time lag)\n", time_point))
  }

  data_current <- data[data[[time_var]] == time_point, ]
  data_current <- data_current[order(data_current[[region_var]]), ]

  if (nrow(data_current) != n) {
    stop("The current-time-point data does not contain all regions")
  }

  y <- numeric(k * n)
  for (i in 1:k) {
    y[((i-1)*n + 1):(i*n)] <- data_current[[y_vars[i]]]
  }

  if (verbose) {
    cat(sprintf("Dependent variable vector y dimension: %d×1\n", length(y)))
    if (include_time_lag) {
      cat(sprintf("Lag variable matrix y_lag dimension: %d×%d\n", nrow(y_lag), ncol(y_lag)))
    }
  }

  X <- build_design_matrix_extended(
    data = data_current,
    y_vars = y_vars,
    x_vars = x_vars,
    y_lag = y_lag,
    include_intercept = include_intercept,
    include_time_lag = include_time_lag,
    n = n,
    k = k,
    verbose = verbose
  )
  
  if (verbose) cat("Computing the eigenvalues of the spatial weight matrix...\n")
  eigen_W <- compute_eigen_W(W)

  if (verbose) cat("Converting to listw form...\n")
  W_listw <- spdep::mat2listw(W, style = "W")

  if (include_time_lag) {
    p0 <- ncol(X) - k^2
  } else {
    p0 <- ncol(X)
  }

  if (verbose) {
    cat(sprintf("Number of regressors p0 = %d (including intercept)\n", p0))
    cat(sprintf("Design matrix X dimension: %d×%d\n", nrow(X), ncol(X)))
    cat("=== Data preparation complete ===\n\n")
  }

  # Sparse copy of W for the hot multiplication paths (Queen contiguity is
  # ~8 neighbours per row, so sparse mat-vec is far cheaper at large n)
  W_sp <- if (requireNamespace("Matrix", quietly = TRUE)) {
    Matrix::Matrix(W, sparse = TRUE)
  } else {
    NULL
  }

  result <- list(
    y = y,
    X = X,
    W = W,
    W_sp = W_sp,
    W_listw = W_listw,
    y_lag = y_lag,
    eigen_W = eigen_W,
    n = n,
    k = k,
    regions = regions,
    p0 = p0,
    data_info = list(
      y_vars = y_vars,
      x_vars = x_vars,
      time_point = time_point,
      time_lag = time_lag,
      region_var = region_var,
      time_var = time_var,
      include_intercept = include_intercept,
      include_time_lag = include_time_lag
    )
  )

  .mgnst_prep_cache[[cache_key]] <- result

  return(result)
}

#' Build the design matrix X (with time-lag option)
#'
#' @param data data frame (current time point)
#' @param y_vars vector of dependent variable names
#' @param x_vars list of regressor names
#' @param y_lag n×k lag variable matrix (no lag if NULL)
#' @param include_intercept whether to include an intercept
#' @param include_time_lag whether to include a time lag
#' @param n number of regions
#' @param k number of variables
#' @param verbose verbose output
#'
#' @return design matrix
build_design_matrix_extended <- function(
  data,
  y_vars,
  x_vars,
  y_lag,
  include_intercept,
  include_time_lag,
  n,
  k,
  verbose = FALSE
) {
  
  if (verbose) cat("Building the design matrix X...\n")

  x_counts <- sapply(x_vars, length)
  if (include_intercept) {
    x_counts <- x_counts + 1
  }
  p0_total <- sum(x_counts)

  if (include_time_lag) {
    total_cols <- p0_total + k^2
  } else {
    total_cols <- p0_total
  }

  X <- matrix(0, nrow = k*n, ncol = total_cols)

  col_idx <- 1

  for (i in 1:k) {
    row_start <- (i-1)*n + 1
    row_end <- i*n

    Xi <- NULL

    if (include_intercept) {
      Xi <- cbind(Xi, rep(1, n))
    }

    if (length(x_vars[[i]]) > 0) {
      for (x_name in x_vars[[i]]) {
        if (!x_name %in% colnames(data)) {
          stop(sprintf("Regressor '%s' does not exist in the data for y%d", x_name, i))
        }
        Xi <- cbind(Xi, data[[x_name]])
      }
    }

    Xi_cols <- ncol(Xi)
    X[row_start:row_end, col_idx:(col_idx + Xi_cols - 1)] <- Xi
    col_idx <- col_idx + Xi_cols
  }

  if (include_time_lag && !is.null(y_lag)) {
    for (i in 1:k) {
      row_start <- (i-1)*n + 1
      row_end <- i*n

      for (j in 1:k) {
        X[row_start:row_end, col_idx] <- y_lag[, j]
        col_idx <- col_idx + 1
      }
    }
  }

  if (verbose) {
    cat(sprintf("Design matrix X: %d×%d\n", nrow(X), ncol(X)))
    cat(sprintf("  - regressor part: %d columns\n", p0_total))
    if (include_time_lag) {
      cat(sprintf("  - lag variable part: %d columns\n", k^2))
    }
  }

  return(X)
}

# 2. Variance-covariance matrix computation

#' Variance-covariance matrix of β (analytical computation)
#'
#' Ψ = {X'(Σ⁻¹⊗I)X}⁻¹
#'
#' @param X kn×p design matrix
#' @param Sigma k×k error covariance matrix
#' @param k number of variables
#' @param n number of regions
#' @return p×p variance-covariance matrix Ψ
compute_vcov_beta <- function(X, Sigma, k, n, Lambda_mat = NULL, W = NULL) {

  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    warning("Error inverting Σ. Applying ridge regularization.")
    solve(Sigma + diag(1e-6, k))
  })

  p <- ncol(X)

  # MSEM/MGNS: compute (I-Λ⊗W)X. MSAR/VARX/OLS: X as-is
  if (!is.null(Lambda_mat) && !is.null(W)) {
    AX <- X - compute_MW_times_M(Lambda_mat, W, X, k, n)
  } else {
    AX <- X
  }

  XtSigmaInvX <- matrix(0, p, p)

  # Compute AX'(Σ⁻¹⊗I)AX  (AX = (I-Λ⊗W)X or X)
  for (i in 1:k) {
    for (j in 1:k) {
      idx_i <- ((i-1)*n + 1):(i*n)
      idx_j <- ((j-1)*n + 1):(j*n)
      XtSigmaInvX <- XtSigmaInvX + Sigma_inv[i,j] * t(AX[idx_i, , drop=FALSE]) %*% AX[idx_j, , drop=FALSE]
    }
  }

  Psi <- tryCatch({
    solve(XtSigmaInvX)
  }, error = function(e) {
    warning("Error inverting X'(Σ⁻¹⊗I)X. Applying ridge regularization.")
    solve(XtSigmaInvX + diag(1e-6, p))
  })

  return(Psi)
}

#' Numerical Hessian computation for spatial parameters
#'
#' Approximation of the second derivative by central differences
#'
#' @param param_vec vector of spatial parameters
#' @param negative_loglik_fn negative log-likelihood function (takes param_vec as argument)
#' @param eps step size for numerical differentiation (default: 1e-5)
#' @return Hessian matrix
compute_hessian_numerical <- function(param_vec, negative_loglik_fn, eps = 1e-5) {

  n_params <- length(param_vec)
  H <- matrix(0, n_params, n_params)

  f0 <- negative_loglik_fn(param_vec)

  for (i in 1:n_params) {
    for (j in i:n_params) {
      ei <- ej <- rep(0, n_params)
      ei[i] <- eps
      ej[j] <- eps

      fpp <- negative_loglik_fn(param_vec + ei + ej)
      fpm <- negative_loglik_fn(param_vec + ei - ej)
      fmp <- negative_loglik_fn(param_vec - ei + ej)
      fmm <- negative_loglik_fn(param_vec - ei - ej)
      
      H[i, j] <- (fpp - fpm - fmp + fmm) / (4 * eps^2)
      H[j, i] <- H[i, j]
    }
  }
  
  return(H)
}

#' Compute the variance-covariance matrix from the Hessian
#'
#' @param hessian Hessian matrix (second derivative of the negative log-likelihood)
#' @return variance-covariance matrix
compute_vcov_from_hessian <- function(hessian, gamma = 0) {

if (gamma > 0) {
    # γ>0: sandwich estimator
    # Since this is the Hessian of the spatial parameters only, D = I
    # H_p = H + γI,  Avar = H_p⁻¹ H H_p⁻¹
    k1 <- nrow(hessian)
    H_pen <- hessian + gamma * diag(k1)

    H_pen_inv <- tryCatch({
      solve(H_pen)
    }, error = function(e) {
      warning("Error inverting (H + γI). Applying ridge regularization.")
      solve(H_pen + diag(1e-6, k1))
    })

    vcov <- H_pen_inv %*% hessian %*% H_pen_inv

  } else {
# γ=0: H⁻¹
    vcov <- tryCatch({
      solve(hessian)
    }, error = function(e) {
      warning("Error inverting the Hessian. Applying ridge regularization.")
      eigen_H <- eigen(hessian, symmetric = TRUE)
      min_eigenval <- min(eigen_H$values)
      if (min_eigenval <= 0) {
        ridge <- abs(min_eigenval) + 1e-6
        hessian <- hessian + diag(ridge, nrow(hessian))
      }
      solve(hessian)
    })
  }

  return(vcov)
}

# 3. Significance-test computation

#' Compute significance tests
#'
#' @param estimates vector of estimates
#' @param std_errors vector of standard errors
#' @param param_names vector of parameter names (optional)
#' @return data.frame (parameter, estimate, std_error, z_value, p_value, signif)
compute_inference <- function(estimates, std_errors, param_names = NULL) {

  n_params <- length(estimates)

  if (is.null(param_names)) {
    param_names <- paste0("param_", 1:n_params)
  }

  z_values <- estimates / std_errors

  p_values <- 2 * pnorm(-abs(z_values))

  signif_codes <- cut(p_values,
    breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    labels = c("***", "**", "*", ".", ""),
    right = FALSE,
    include.lowest = TRUE
  )
  
  result <- data.frame(
    parameter = param_names,
    estimate = estimates,
    std_error = std_errors,
    z_value = z_values,
    p_value = p_values,
    signif = as.character(signif_codes),
    stringsAsFactors = FALSE
  )
  
  return(result)
}

#' Get the significance-level symbol
#'
#' @param p_value p-value
#' @return significance-level symbol
get_signif_code <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  if (p_value < 0.1) return(".")
  return("")
}

# 4. Coefficient of determination computation

#' Compute four kinds of R²
#'
#' @param y observed values
#' @param fitted fitted values
#' @param residuals residuals
#' @param loglik model log-likelihood
#' @param num_params number of parameters (including spatial parameters)
#' @return list(R2, R2_adj, R2_cor, R2_pseudo)
compute_r_squared <- function(y, fitted, residuals, loglik = NULL, num_params = NULL) {

  n <- length(y)
  y_mean <- mean(y)

  # 1. SS-based R² (OLS-compatible)
  SST <- sum((y - y_mean)^2)
  SSE <- sum(residuals^2)
  R2 <- 1 - SSE/SST

  # 2. Adjusted R² (including spatial parameters)
  R2_adj <- 1 - (1 - R2) * (n - 1) / (n - num_params - 1)

  # 3. Correlation-based R²
  R2_cor <- cor(y, fitted)^2

  # 4. Pseudo R² (McFadden)
  # Log-likelihood of the null model (intercept only)
  null_loglik <- -n/2 * (log(2*pi) + 1 + log(var(y)))
  R2_pseudo <- 1 - (loglik / null_loglik)

  return(list(
    R2 = R2,
    R2_adj = R2_adj,
    R2_cor = R2_cor,
    R2_pseudo = R2_pseudo
  ))
}

# 4-1b. Fitted values ŷ of the multivariate model (for R² computation)

#' Compute the fitted values ŷ of the multivariate model
#'
#' Based on E[y] = (I - R⊗W)⁻¹ Xβ:
#'   MSAR/MGNS (R≠0): ŷ = (I - R⊗W)⁻¹ Xβ
#'   MSEM/VARX/OLS (R=0): ŷ = Xβ
#'
#' @param model_type model type string
#' @param R k×k spatial lag matrix (NULL allowed)
#' @param beta regression coefficient vector
#' @param X design matrix
#' @param W spatial weight matrix
#' @param k number of variables
#' @param n number of regions
#' @return Kn×1 fitted value vector
compute_fitted_multivar <- function(model_type, R, beta, X, W, k, n) {

  Xbeta <- as.numeric(X %*% beta)

  # MSAR/MGNS: ŷ = (I - R⊗W)⁻¹ Xβ
  if (model_type %in% c("MSAR", "MGNS", "IndSAR", "IndGNS") && !is.null(R)) {
    IKn_minus_RW <- diag(k * n) - kronecker(R, W)
    fitted_vals <- tryCatch({
      as.numeric(solve(IKn_minus_RW, Xbeta))
    }, error = function(e) {
      warning("Error solving the linear system (I - R⊗W). Using Xβ as a fallback: ", e$message)
      Xbeta
    })
  } else {
    # MSEM/VARX/OLS: ŷ = Xβ
    fitted_vals <- Xbeta
  }

  return(fitted_vals)
}

# 4-2. Multivariate pseudo R² computation

#' Compute the pseudo R² of the multivariate model
#'
#' R²_pseudo,k = corr(y_k, ŷ_k)²
#' R̄²_pseudo   = (1/K) Σ R²_pseudo,k
#'
#' ŷ is the fitted value computed by compute_fitted_multivar().
#'
#' @param y Kn×1 observed value vector
#' @param fitted Kn×1 fitted value vector (output of compute_fitted_multivar())
#' @param k number of variables
#' @param n number of regions
#' @return list(R2_pseudo_individual, R2_pseudo_mean)
compute_r_squared_multivar <- function(y, fitted, k, n) {

  R2_pseudo_vec <- numeric(k)

  for (i in 1:k) {
    idx <- ((i-1)*n + 1):(i*n)
    yi <- y[idx]
    fi <- fitted[idx]

    # pseudo R² = corr(y_k, ŷ_k)²
    R2_pseudo_vec[i] <- cor(yi, fi)^2
  }

  list(
    R2_pseudo_individual = R2_pseudo_vec,
    R2_pseudo_mean = mean(R2_pseudo_vec)   # R̄²_pseudo
  )
}

# 5. Extracting significance info from spatialreg models

#' Extract significance info from a spatialreg model (lagsarlm/errorsarlm)
#'
#' Function to obtain detailed info from the S1 per-variable estimation results
#'
#' @param model spatialreg model object
#' @param model_type "MSAR" (lagsarlm) or "MSEM" (errorsarlm)
#' @return list(coefficients, spatial_params, fit)
extract_inference_from_spatial_model <- function(model, model_type = "MSAR") {

  s <- summary(model)

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

  if (model_type == "MSAR") {
    rho <- model$rho
    rho_se <- s$rho.se
    rho_z <- rho / rho_se
    rho_p <- 2 * pnorm(-abs(rho_z))
    
    spatial_params <- list(
      rho = rho,
      rho_se = rho_se,
      rho_z = rho_z,
      rho_p = rho_p,
      rho_signif = get_signif_code(rho_p),
      lambda = NULL,
      lambda_se = NULL,
      lambda_z = NULL,
      lambda_p = NULL,
      lambda_signif = NULL,
      LR_statistic = if (!is.null(s$LR1)) s$LR1$statistic else NA,
      LR_p_value = if (!is.null(s$LR1)) s$LR1$p.value else NA
    )
    
  } else {
    lambda <- model$lambda
    lambda_se <- s$lambda.se
    lambda_z <- lambda / lambda_se
    lambda_p <- 2 * pnorm(-abs(lambda_z))
    
    spatial_params <- list(
      rho = NULL,
      rho_se = NULL,
      rho_z = NULL,
      rho_p = NULL,
      rho_signif = NULL,
      lambda = lambda,
      lambda_se = lambda_se,
      lambda_z = lambda_z,
      lambda_p = lambda_p,
      lambda_signif = get_signif_code(lambda_p),
      LR_statistic = if (!is.null(s$LR1)) s$LR1$statistic else NA,
      LR_p_value = if (!is.null(s$LR1)) s$LR1$p.value else NA
    )
  }
  
  y <- model$y
  fitted_vals <- fitted(model)
  residuals_vals <- residuals(model)
  loglik <- as.numeric(logLik(model))

  # Number of parameters (regression coefficients + spatial parameter + variance)
  num_params <- length(coef(model)) + 1 + 1
  
  r2_results <- compute_r_squared(
    y = y,
    fitted = fitted_vals,
    residuals = residuals_vals,
    loglik = loglik,
    num_params = num_params
  )
  
  fit <- list(
    loglik = loglik,
    AIC = AIC(model),
    BIC = BIC(model),
    R2 = r2_results$R2,
    R2_adj = r2_results$R2_adj,
    R2_cor = r2_results$R2_cor,
    R2_pseudo = r2_results$R2_pseudo,
    num_params = num_params,
    num_obs = length(y),
    sigma2 = model$s2
  )
  
  residuals_info <- list(
    raw = as.numeric(residuals_vals),
    fitted = as.numeric(fitted_vals)
  )
  
  return(list(
    coefficients = coef_table,
    spatial_params = spatial_params,
    fit = fit,
    residuals = residuals_info
  ))
}

# 6. Unified coefficient-table output

#' Formatted output of the coefficient table
#'
#' @param inference_table output of compute_inference()
#' @param title title
#' @param digits number of decimal places
print_inference_table <- function(inference_table, title = "Coefficient table", digits = 4) {

  cat(sprintf("\n%s:\n", title))
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("%-20s %10s %10s %10s %10s %5s\n",
              "Parameter", "Estimate", "Std.Error", "z-value", "p-value", ""))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  for (i in 1:nrow(inference_table)) {
    row <- inference_table[i, ]
    cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f %5s\n",
                row$parameter,
                row$estimate,
                row$std_error,
                row$z_value,
                row$p_value,
                row$signif))
  }
  
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("Significance levels: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")
}

# 7. Helper functions

#' Vectorize the spatial parameters
#'
#' Convert the R and Λ matrices into a vector
#'
#' @param R k×k R matrix (for MSAR, MGNS), NULL allowed
#' @param Lambda_mat k×k Λ matrix (for MSEM, MGNS), NULL allowed
#' @return parameter vector
spatial_params_to_vec <- function(R = NULL, Lambda_mat = NULL) {
  vec <- c()
  if (!is.null(R)) {
    vec <- c(vec, as.vector(t(R)))
  }
  if (!is.null(Lambda_mat)) {
    vec <- c(vec, as.vector(t(Lambda_mat)))
  }
  return(vec)
}

#' Convert a vector into spatial parameter matrices
#'
#' @param vec parameter vector
#' @param k number of variables
#' @param model_type "MSAR", "MSEM", "MGNS"
#' @return list(R, Lambda)
vec_to_spatial_params <- function(vec, k, model_type) {
  
  R <- NULL
  Lambda_mat <- NULL
  
  if (model_type == "MSAR") {
    R <- matrix(vec, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "MSEM") {
    Lambda_mat <- matrix(vec, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "MGNS") {
    n_R <- k^2
    R <- matrix(vec[1:n_R], nrow = k, ncol = k, byrow = TRUE)
    Lambda_mat <- matrix(vec[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
  }
  
  return(list(R = R, Lambda = Lambda_mat))
}
