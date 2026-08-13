# msar_likelihood.r
# Core likelihood and parameter-update functions for the MSAR model
# (log-likelihood, profile likelihood, AIC/BIC, β/Σ updates, residuals).

#' Compute the log-likelihood
#'
#' ℓ = -(kn/2)log(2π) + log|I-R⊗W| - (n/2)log|Σ| - Q/2
#'
#' @param R k×k spatial correlation matrix
#' @param beta (p0+k²)×1 regression coefficients
#' @param Sigma k×k error covariance matrix
#' @param y kn×1 dependent variable
#' @param X kn×(p0+k²) design matrix
#' @param W n×n spatial weight matrix
#' @param eigen_W eigenvalues of W
#' @param k, n dimensions
#' @param verbose verbosity level (0: none, 1: basic, 2: detailed)
#' @return value of the log-likelihood
compute_log_likelihood <- function(R, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {

  log_det_spatial_term <- log_det_spatial(R, eigen_W, k, n, verbose = (verbose >= 2))

  if (!is.finite(log_det_spatial_term)) {
    if (verbose >= 1) cat("Warning: log|I - R⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]

  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Σ| is not finite\n")
    return(-Inf)
  }

  # Residuals z = (I - R⊗W)y - Xβ
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta

  # Q = z'(Σ^{-1}⊗I)z, computed block by block
  Sigma_inv <- solve(Sigma)
  Q <- 0

  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }

  loglik <- -(k*n/2) * log(2*pi) + log_det_spatial_term - (n/2) * log_det_Sigma - Q/2

  if (verbose >= 2) {
    cat(sprintf("  log|I - R⊗W| = %.6f\n", log_det_spatial_term))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }

  return(loglik)
}


#' Compute the profile likelihood
#'
#' ℓ(R) = -(kn/2)log(2πe) + log|I-R⊗W| - (n/2)log|Σ̂(R)|
#'
#' @param R k×k spatial correlation matrix
#' @param beta_hat estimated β (depends on R)
#' @param Sigma_hat estimated Σ (depends on R)
#' @param y, X, W data
#' @param eigen_W eigenvalues of W
#' @param k, n dimensions
#' @param verbose verbosity level
#' @return profile log-likelihood
compute_profile_likelihood <- function(R, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {

  log_det_spatial_term <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = smooth)

  if (!is.finite(log_det_spatial_term)) {
    if (verbose >= 1) cat("Warning: log|I - R⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]

  if (!is.finite(log_det_Sigma_hat)) {
    if (verbose >= 1) cat("Warning: log|Σ̂| is not finite\n")
    return(-Inf)
  }

  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_spatial_term - (n/2) * log_det_Sigma_hat

  if (verbose >= 2) {
    cat(sprintf("  profile likelihood: %.6f\n", profile_loglik))
  }

  return(profile_loglik)
}


#' Compute AIC and BIC
#'
#' @param loglik log-likelihood
#' @param num_params number of parameters
#' @param num_obs number of observations
#' @return list(AIC, BIC)
compute_information_criteria <- function(loglik, num_params, num_obs) {
  AIC <- -2 * loglik + 2 * num_params
  BIC <- -2 * loglik + log(num_obs) * num_params
  
  return(list(AIC = AIC, BIC = BIC))
}


# Parameter-update functions

#' Update β
#'
#' β̂ = {X'(Σ⁻¹⊗I)X}⁻¹ X'(Σ⁻¹⊗I)(I-R⊗W)y
#'
#' @param R current R matrix
#' @param Sigma current Σ matrix
#' @param y, X, W data
#' @param k, n dimensions
#' @param ridge_eps ridge regularization parameter
#' @param verbose verbosity level
#' @return estimated β
update_beta <- function(R, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: error inverting Σ. Applying regularization\n")
    solve(Sigma + diag(ridge_eps, k))
  })

  # (I - R⊗W)y
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  y_transformed <- y - RWy

  # X'(Σ^{-1}⊗I)X and X'(Σ^{-1}⊗I)y via block operations
  XtSigmaInvX <- matrix(0, ncol(X), ncol(X))
  XtSigmaInvy <- numeric(ncol(X))

  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)

      Xi <- X[i_idx, , drop = FALSE]
      Xj <- X[j_idx, , drop = FALSE]

      XtSigmaInvX <- XtSigmaInvX + Sigma_inv[i,j] * t(Xi) %*% Xj
      XtSigmaInvy <- XtSigmaInvy + Sigma_inv[i,j] * t(Xi) %*% y_transformed[j_idx]
    }
  }

  beta_hat <- tryCatch({
    solve(XtSigmaInvX, XtSigmaInvy)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: singular matrix when estimating β. Applying ridge regularization\n")
    solve(XtSigmaInvX + diag(ridge_eps, ncol(XtSigmaInvX)), XtSigmaInvy)
  })

  if (verbose >= 2) {
    cat(sprintf("  β update done: ||β|| = %.6f\n", sqrt(sum(beta_hat^2))))
  }

  return(beta_hat)
}


#' Update Σ
#'
#' Σ̂ = (1/n) Z'Z (block structure, symmetric matrix)
#' where z = (I-R⊗W)y - Xβ
#'
#' @param R current R matrix
#' @param beta current β
#' @param y, X, W data
#' @param k, n dimensions
#' @param ridge_eps regularization parameter
#' @param verbose verbosity level
#' @return estimated Σ (symmetric matrix)
update_Sigma <- function(R, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  # Residuals z = (I - R⊗W)y - Xβ
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta

  Sigma_new <- matrix(0, k, k)

  for (i in 1:k) {
    for (j in i:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)

      Sigma_new[i,j] <- sum(z[zi_idx] * z[zj_idx]) / n

      if (i != j) {
        Sigma_new[j,i] <- Sigma_new[i,j]
      }
    }
  }

  # Check and fix positive definiteness
  eigen_vals <- eigen(Sigma_new, only.values = TRUE)$values
  min_eigen <- min(Re(eigen_vals))

  if (min_eigen <= 0) {
    if (verbose >= 1) {
      cat(sprintf("Warning: Σ is not positive definite (min eigenvalue = %.2e). Applying regularization\n", min_eigen))
    }
    Sigma_new <- Sigma_new + diag(ridge_eps, k)
  }

  if (verbose >= 2) {
    cat(sprintf("  Σ update done: log|Σ| = %.6f\n", determinant(Sigma_new, logarithm = TRUE)$modulus[1]))
  }

  return(Sigma_new)
}


#' Compute residuals
#'
#' z = (I-R⊗W)y - Xβ
#'
#' @param R, beta parameters
#' @param y, X, W data
#' @param k, n dimensions
#' @return kn×1 residual vector
compute_residuals <- function(R, beta, y, X, W, k, n) {
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta
  return(z)
}


#' Iterative update of β and Σ (the S2 loop)
#'
#' Fix R and alternately update β and Σ
#' Convergence: change in log|Σ| below the threshold
#'
#' @param R fixed R matrix
#' @param beta_init, Sigma_init initial values
#' @param y, X, W data
#' @param k, n dimensions
#' @param max_iter maximum number of iterations
#' @param tol convergence threshold
#' @param ridge_eps ridge regularization parameter
#' @param verbose verbosity level (0: none, 1: basic, 2: detailed)
#' @return list(beta, Sigma, log_det_Sigma_history, converged, iterations)
iterate_beta_sigma <- function(
  R, 
  beta_init, 
  Sigma_init, 
  y, X, W, 
  k, n, 
  max_iter = 100, 
  tol = 1e-6,
  ridge_eps = 1e-6,
  verbose = 0
) {
  
  beta <- beta_init
  Sigma <- Sigma_init
  log_det_Sigma_history <- numeric(0)

  if (verbose >= 1) {
    cat("  Inner loop start (updating β and Σ)\n")
  }

  for (iter in 1:max_iter) {

    log_det_Sigma_old <- determinant(Sigma, logarithm = TRUE)$modulus[1]

    # Save history every 10 iterations
    if (iter %% 10 == 1) {
      log_det_Sigma_history <- c(log_det_Sigma_history, log_det_Sigma_old)
    }

    beta <- update_beta(R, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    Sigma <- update_Sigma(R, beta, y, X, W, k, n, ridge_eps, verbose = verbose)

    log_det_Sigma_new <- determinant(Sigma, logarithm = TRUE)$modulus[1]

    change <- abs(log_det_Sigma_new - log_det_Sigma_old)

    if (verbose >= 2 && iter %% 10 == 0) {
      cat(sprintf("    inner iter %d: log|Σ| = %.6f, change = %.2e\n",
                  iter, log_det_Sigma_new, change))
    }

    if (change < tol) {
      if (verbose >= 1) {
        cat(sprintf("  Inner loop converged (%d iters): change = %.2e < %.2e\n", iter, change, tol))
      }
      return(list(
        beta = beta,
        Sigma = Sigma,
        log_det_Sigma_history = log_det_Sigma_history,
        converged = TRUE,
        iterations = iter
      ))
    }
  }

  if (verbose >= 1) {
    cat(sprintf("  Inner loop reached max iterations (%d)\n", max_iter))
  }
  
  return(list(
    beta = beta,
    Sigma = Sigma,
    log_det_Sigma_history = log_det_Sigma_history,
    converged = FALSE,
    iterations = max_iter
  ))
}
