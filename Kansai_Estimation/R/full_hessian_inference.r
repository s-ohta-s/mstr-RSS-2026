# full_hessian_inference.r
# Significance tests based on the joint Hessian over all parameters
# θ = (β', vec(R)', vec(Λ)')' (unlike compute_vcov_beta(), which conditions on
# R, Λ fixed at their estimates). Usage: result <- add_full_inference(result, gamma = 5)
#
# Theory:
#   Penalized log-likelihood: ℓ_p(θ) = ℓ(θ) - γ/2 * θ'Dθ
#   H(θ̂)  = -∂²ℓ/∂θ∂θ'|_{θ=θ̂}   (non-penalized Hessian)
#   H_p(θ̂) = H(θ̂) + γD           (penalized Hessian)
#   Avar(θ̂) ≈ H_p⁻¹ H H_p⁻¹      (sandwich estimator)
#   D = diag(0,...,0, 1,...,1)   (0 for β, 1 for the spatial parameters)

# Log-likelihoods with Σ profiled out (functions of θ = (β, spatial))

#' MSAR: Σ-profile log-likelihood — function of θ = (β, vec(R))
#'
#' From z = (I-R⊗W)y - Xβ, derive Σ̂ = (1/n)[z'_k z_l] analytically,
#' and return ℓ(θ) = -Kn/2 log(2πe) + log|I-R⊗W| - n/2 log|Σ̂|
#'
#' @param theta θ = c(β, vec(R)') — vector of length p + K²
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param eigen_W eigenvalue vector of W
#' @param k number of variables
#' @param n number of regions
#' @param p dimension of β
#' @return scalar (log-likelihood)
full_loglik_msar <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  R_vec <- theta[(p + 1):(p + k^2)]
  R <- matrix(R_vec, nrow = k, ncol = k, byrow = TRUE)

  log_det_R <- tryCatch(
    log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_R)) return(-1e10)

  # Residuals z = (I-R⊗W)y - Xβ
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta

  # Σ̂ = (1/n) [z'_k z_l]
  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }

  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  # ℓ = -Kn/2 log(2πe) + log|I-R⊗W| - n/2 log|Σ̂|
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_R - (n / 2) * log_det_Sigma
  return(loglik)
}


#' MSEM: Σ-profile log-likelihood — function of θ = (β, vec(Λ))
full_loglik_msem <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  Lambda_vec <- theta[(p + 1):(p + k^2)]
  Lambda_mat <- matrix(Lambda_vec, nrow = k, ncol = k, byrow = TRUE)

  log_det_Lambda <- tryCatch(
    log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_Lambda)) return(-1e10)
  
  # Residuals z = (I-Λ⊗W)(y - Xβ)
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }
  
  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_Lambda - (n / 2) * log_det_Sigma
  return(loglik)
}


#' MGNS: Σ-profile log-likelihood — function of θ = (β, vec(R), vec(Λ))
full_loglik_mgns <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  R_vec <- theta[(p + 1):(p + k^2)]
  Lambda_vec <- theta[(p + k^2 + 1):(p + 2 * k^2)]
  R <- matrix(R_vec, nrow = k, ncol = k, byrow = TRUE)
  Lambda_mat <- matrix(Lambda_vec, nrow = k, ncol = k, byrow = TRUE)

  log_det_R <- tryCatch(
    log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  log_det_Lambda <- tryCatch(
    log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_R) || !is.finite(log_det_Lambda)) return(-1e10)
  
  # Residuals z = (I-Λ⊗W){(I-R⊗W)y - Xβ}
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }
  
  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_R + log_det_Lambda - (n / 2) * log_det_Sigma
  return(loglik)
}


# Joint Hessian computation over all parameters

#' Numerical Hessian over all parameters θ = (β, spatial)
#'
#' Compute the second derivative of -ℓ(θ) using numDeriv::hessian
#'
#' @param object multivar_spatial object (already estimated)
#' @param verbose verbose output
#' @return list(H, theta_hat, p_beta, p_spatial, dim_theta) or NULL
compute_full_hessian <- function(object, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  beta <- object$coefficients$beta
  model_type <- object$model_type
  
  p <- length(beta)
  
  # Build θ̂ and select the model-specific likelihood function
  if (model_type == "MSAR") {
    R_vec <- as.vector(t(object$coefficients$R))
    theta_hat <- c(beta, R_vec)
    p_spatial <- k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_msar(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else if (model_type == "MSEM") {
    Lambda_vec <- as.vector(t(object$coefficients$Lambda))
    theta_hat <- c(beta, Lambda_vec)
    p_spatial <- k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_msem(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else if (model_type == "MGNS") {
    R_vec <- as.vector(t(object$coefficients$R))
    Lambda_vec <- as.vector(t(object$coefficients$Lambda))
    theta_hat <- c(beta, R_vec, Lambda_vec)
    p_spatial <- 2 * k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_mgns(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else {
    warning(sprintf("Model type '%s' is not supported for the full-parameter Hessian", model_type))
    return(NULL)
  }
  
  dim_theta <- length(theta_hat)
  
  if (verbose) {
    cat(sprintf("  Dimension of θ: %d (β: %d, spatial: %d)\n", dim_theta, p, p_spatial))
    cat("  Computing numDeriv::hessian...\n")
  }

  H <- tryCatch({
    numDeriv::hessian(neg_loglik, theta_hat)
  }, error = function(e) {
    warning("Error computing the full-parameter Hessian: ", e$message)
    NULL
  })
  
  if (is.null(H)) return(NULL)
  
  # Symmetrize (absorb numerical error)
  H <- (H + t(H)) / 2

  eig_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_H)
  
  if (verbose) {
    cat(sprintf("  Hessian eigenvalues: min = %.4e, max = %.4e\n",
                min_eig, max(eig_H)))
  }
  
  if (min_eig <= 0) {
    warning(sprintf("Hessian is not positive definite (min eigenvalue = %.2e). Applying ridge correction", min_eig))
    ridge <- abs(min_eig) + 1e-6
    H <- H + diag(ridge, dim_theta)
  }
  
  return(list(
    H         = H,
    theta_hat = theta_hat,
    p_beta    = p,
    p_spatial = p_spatial,
    dim_theta = dim_theta
  ))
}


# Sandwich variance estimation

#' Sandwich variance over all parameters
#'
#' D = diag(0,...,0, 1,...,1): no penalty on β, penalty on the spatial parameters
#' Avar(θ̂) = (H + γD)⁻¹ H (H + γD)⁻¹
#'
#' @param H non-penalized Hessian (dim_theta × dim_theta)
#' @param gamma penalty strength
#' @param p_beta dimension of β
#' @param p_spatial dimension of the spatial parameters
#' @return variance-covariance matrix (dim_theta × dim_theta)
compute_full_sandwich_vcov <- function(H, gamma, p_beta, p_spatial) {
  
  dim_theta <- nrow(H)
  
  if (gamma == 0) {
    # No penalty: ordinary inverse information matrix
    vcov <- tryCatch(solve(H), error = function(e) {
      warning("Error inverting the full Hessian. Applying regularization: ", e$message)
      eig <- eigen(H, symmetric = TRUE)
      ridge <- max(abs(min(eig$values)), 1e-6)
      solve(H + diag(ridge, dim_theta))
    })
  } else {
    # D matrix: 0 for β, 1 for the spatial parameters
    D_diag <- c(rep(0, p_beta), rep(1, p_spatial))
    D <- diag(D_diag)
    
    H_pen <- H + gamma * D
    
    H_pen_inv <- tryCatch(solve(H_pen), error = function(e) {
      warning("Error inverting (H + γD). Applying regularization: ", e$message)
      solve(H_pen + diag(1e-6, dim_theta))
    })
    
    # Sandwich: H_p⁻¹ H H_p⁻¹
    vcov <- H_pen_inv %*% H %*% H_pen_inv
  }
  
  return(vcov)
}


#' Compute the Z statistics for all parameters
#'
#' Z_j = θ̂_j / sqrt([H_p⁻¹ H H_p⁻¹]_jj)
#'
#' @param theta_hat θ̂ vector
#' @param vcov_full sandwich variance over all parameters
#' @return data.frame (index, estimate, se, z, p_value)
compute_full_z_statistics <- function(theta_hat, vcov_full) {
  
  se <- sqrt(pmax(diag(vcov_full), 0))
  z <- ifelse(se > 0, theta_hat / se, NA_real_)
  p_value <- ifelse(!is.na(z), 2 * pnorm(-abs(z)), NA_real_)
  
  data.frame(
    index    = seq_along(theta_hat),
    estimate = theta_hat,
    se       = se,
    z_value  = z,
    p_value  = p_value,
    stringsAsFactors = FALSE
  )
}


# Main function

#' Add full-parameter Hessian-based significance tests to the estimation result
#'
#' In addition to the Ψ-based SE computed by the existing add_inference(),
#' add the SE from the joint full-parameter Hessian -> sandwich estimator
#' as a "full" slot.
#'
#' @param object multivar_spatial object (MSAR/MSEM/MGNS)
#' @param gamma penalty strength (taken from object if NULL)
#' @param verbose verbose output
#' @return the object with full inference added
add_full_inference <- function(object, gamma = NULL, verbose = FALSE) {
  
  model_type <- object$model_type

  if (!(model_type %in% c("MSAR", "MSEM", "MGNS"))) {
    if (verbose) cat(sprintf("%s: the full-parameter Hessian only supports MSAR/MSEM/MGNS\n", model_type))
    return(object)
  }
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("Significance tests based on the joint full-parameter Hessian\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
  }

  if (is.null(gamma)) {
    if (!is.null(object$penalty$gamma)) {
      gamma <- object$penalty$gamma
    } else {
      gamma <- 0
    }
  }
  if (verbose) cat(sprintf("  γ = %g\n", gamma))
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars

  if (verbose) cat("\n[Step 1] Computing the full-parameter Hessian H(θ̂)\n")
  
  hess_result <- compute_full_hessian(object, verbose = verbose)
  
  if (is.null(hess_result)) {
    warning("Failed to compute the full-parameter Hessian")
    return(object)
  }
  
  H          <- hess_result$H
  theta_hat  <- hess_result$theta_hat
  p_beta     <- hess_result$p_beta
  p_spatial  <- hess_result$p_spatial
  dim_theta  <- hess_result$dim_theta

  if (verbose) cat("\n[Step 2] Sandwich variance Avar(θ̂) = Hₚ⁻¹ H Hₚ⁻¹\n")
  
  vcov_full <- tryCatch({
    compute_full_sandwich_vcov(H, gamma, p_beta, p_spatial)
  }, error = function(e) {
    warning("Error computing the sandwich variance: ", e$message)
    NULL
  })
  
  if (is.null(vcov_full)) return(object)

  if (verbose) cat("\n[Step 3] Constructing Z statistics\n")

  z_stats <- compute_full_z_statistics(theta_hat, vcov_full)

  # β block
  se_beta_full <- z_stats$se[1:p_beta]
  z_beta_full  <- z_stats$z_value[1:p_beta]
  p_beta_full  <- z_stats$p_value[1:p_beta]

  # Spatial-parameter block
  spatial_idx <- (p_beta + 1):dim_theta
  se_spatial_full <- z_stats$se[spatial_idx]
  z_spatial_full  <- z_stats$z_value[spatial_idx]
  p_spatial_full  <- z_stats$p_value[spatial_idx]

  if (verbose) {
    cat("\n[Step 4] Comparison of results\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")

    # Existing Ψ-based SE
    Lambda_for_vcov <- if (model_type %in% c("MSEM", "MGNS")) object$coefficients$Lambda else NULL
    W_for_vcov <- if (!is.null(Lambda_for_vcov)) object$model_data$W else NULL
    
    Psi <- tryCatch(
      compute_vcov_beta(object$model_data$X, object$coefficients$Sigma, k, n,
                        Lambda_mat = Lambda_for_vcov, W = W_for_vcov),
      error = function(e) NULL
    )
    se_beta_psi <- if (!is.null(Psi)) sqrt(pmax(diag(Psi), 0)) else rep(NA, p_beta)

    beta_names <- names(object$coefficients$beta)
    if (is.null(beta_names)) beta_names <- paste0("β[", 1:p_beta, "]")
    
    cat("\n  Comparison of β standard errors:\n")
    cat(sprintf("  %-30s %12s %12s %8s\n", "Parameter", "SE(Ψ)", "SE(Full)", "Ratio"))
    cat(paste0("  ", paste(rep("-", 62), collapse = ""), "\n"))
    
    for (i in 1:p_beta) {
      ratio <- if (!is.na(se_beta_psi[i]) && se_beta_psi[i] > 0) {
        se_beta_full[i] / se_beta_psi[i]
      } else {
        NA
      }
      sig_full <- get_signif_code(p_beta_full[i])
      cat(sprintf("  %-30s %12.6f %12.6f %7.3f %s\n",
                  beta_names[i], se_beta_psi[i], se_beta_full[i],
                  ifelse(is.na(ratio), NA, ratio), sig_full))
    }

    cat("\n  Spatial parameters:\n")
    cat(sprintf("  %-20s %12s %12s %12s %8s\n",
                "Parameter", "Estimate", "SE(Full)", "z", "Sig"))
    cat(paste0("  ", paste(rep("-", 64), collapse = ""), "\n"))
    
    sp_names <- character(0)
    if (model_type %in% c("MSAR", "MGNS")) {
      for (i in 1:k) for (j in 1:k) {
        sp_names <- c(sp_names, sprintf("R[%d,%d]", i, j))
      }
    }
    if (model_type %in% c("MSEM", "MGNS")) {
      for (i in 1:k) for (j in 1:k) {
        sp_names <- c(sp_names, sprintf("Λ[%d,%d]", i, j))
      }
    }
    
    for (i in seq_along(spatial_idx)) {
      sig <- get_signif_code(p_spatial_full[i])
      cat(sprintf("  %-20s %12.6f %12.6f %12.4f %s\n",
                  sp_names[i], theta_hat[spatial_idx[i]],
                  se_spatial_full[i], z_spatial_full[i], sig))
    }
  }

  # Store into the object; reshape the spatial-parameter SE into matrices
  se_R_full <- NULL
  se_Lambda_full <- NULL
  
  if (model_type == "MSAR") {
    se_R_full <- matrix(se_spatial_full, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "MSEM") {
    se_Lambda_full <- matrix(se_spatial_full, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "MGNS") {
    se_R_full <- matrix(se_spatial_full[1:(k^2)], nrow = k, ncol = k, byrow = TRUE)
    se_Lambda_full <- matrix(se_spatial_full[(k^2 + 1):(2 * k^2)], nrow = k, ncol = k, byrow = TRUE)
  }
  
  object$vcov$full <- vcov_full
  
  object$hessian$full <- list(
    matrix    = H,
    theta_hat = theta_hat,
    p_beta    = p_beta,
    p_spatial = p_spatial,
    method    = "numerical_full"
  )
  
  object$std_errors$beta_full <- se_beta_full
  object$std_errors$R_full    <- se_R_full
  object$std_errors$Lambda_full    <- se_Lambda_full
  
  object$inference$full <- list(
    z_statistics = z_stats,
    gamma        = gamma,
    method       = "full_hessian_sandwich"
  )
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("Added significance tests based on the full-parameter Hessian\n")
    cat("  Access:\n")
    cat("    object$std_errors$beta_full  — full-Hessian SE of β\n")
    cat("    object$std_errors$R_full     — SE matrix of R\n")
    cat("    object$std_errors$Lambda_full     — SE matrix of Λ\n")
    cat("    object$vcov$full             — full-parameter variance-covariance matrix\n")
    cat("    object$inference$full        — Z-statistics table\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
  }
  
  return(object)
}


# CSV-output helper

#' Generate a coefficient table including the full-Hessian SE
#'
#' @param object object that has had add_full_inference applied
#' @return data.frame (parameter, estimate, se_psi, se_full, z_psi, z_full, ...)
export_full_inference_table <- function(object) {
  
  if (is.null(object$inference$full)) {
    warning("No full-Hessian inference results. Run add_full_inference() first")
    return(NULL)
  }
  
  k <- object$data_info$k
  n <- object$data_info$n
  model_type <- object$model_type
  beta <- object$coefficients$beta
  p <- length(beta)

  # Ψ-based SE
  Lambda_for_vcov <- if (model_type %in% c("MSEM", "MGNS")) object$coefficients$Lambda else NULL
  W_for_vcov <- if (!is.null(Lambda_for_vcov)) object$model_data$W else NULL
  
  Psi <- tryCatch(
    compute_vcov_beta(object$model_data$X, object$coefficients$Sigma, k, n,
                      Lambda_mat = Lambda_for_vcov, W = W_for_vcov),
    error = function(e) NULL
  )
  se_psi <- if (!is.null(Psi)) sqrt(pmax(diag(Psi), 0)) else rep(NA, p)

  se_full <- object$std_errors$beta_full

  beta_names <- names(beta)
  if (is.null(beta_names)) beta_names <- paste0("beta[", 1:p, "]")

  z_psi <- ifelse(se_psi > 0, beta / se_psi, NA)
  z_full <- ifelse(se_full > 0, beta / se_full, NA)
  p_psi <- ifelse(!is.na(z_psi), 2 * pnorm(-abs(z_psi)), NA)
  p_full <- ifelse(!is.na(z_full), 2 * pnorm(-abs(z_full)), NA)
  
  tbl <- data.frame(
    parameter  = beta_names,
    estimate   = beta,
    se_psi     = se_psi,
    se_full    = se_full,
    z_psi      = z_psi,
    z_full     = z_full,
    p_psi      = p_psi,
    p_full     = p_full,
    signif_psi = sapply(p_psi, get_signif_code),
    signif_full = sapply(p_full, get_signif_code),
    stringsAsFactors = FALSE
  )

  # Add the spatial-parameter rows
  z_stats <- object$inference$full$z_statistics
  spatial_start <- p + 1
  spatial_end <- nrow(z_stats)
  
  sp_names <- character(0)
  if (model_type %in% c("MSAR", "MGNS")) {
    for (i in 1:k) for (j in 1:k) sp_names <- c(sp_names, sprintf("R[%d,%d]", i, j))
  }
  if (model_type %in% c("MSEM", "MGNS")) {
    for (i in 1:k) for (j in 1:k) sp_names <- c(sp_names, sprintf("Lambda[%d,%d]", i, j))
  }
  
  sp_rows <- z_stats[spatial_start:spatial_end, ]
  
  sp_tbl <- data.frame(
    parameter   = sp_names,
    estimate    = sp_rows$estimate,
    se_psi      = NA_real_,
    se_full     = sp_rows$se,
    z_psi       = NA_real_,
    z_full      = sp_rows$z_value,
    p_psi       = NA_real_,
    p_full      = sp_rows$p_value,
    signif_psi  = "",
    signif_full = sapply(sp_rows$p_value, get_signif_code),
    stringsAsFactors = FALSE
  )
  
  result_tbl <- rbind(tbl, sp_tbl)
  rownames(result_tbl) <- NULL
  
  return(result_tbl)
}
