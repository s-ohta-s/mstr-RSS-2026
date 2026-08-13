# msem_mgns_likelihood.r
# Core likelihood and parameter-update functions for the MSEM and MGNS models
# (log-likelihoods, profile likelihoods, β/Σ updates, residuals) plus block
# Kronecker-product utilities.

#' Efficient computation of (Λ⊗W)v (for MSEM)
#'
#' Implements (Λ⊗W)v with block operations instead of forming the Kronecker product explicitly
#'
#' @param Lambda k×k spatial error correlation matrix
#' @param W n×n spatial weight matrix
#' @param v kn×1 vector
#' @param k number of variables
#' @param n number of regions
#' @param verbose verbose output
#' @return kn×1 vector
compute_LambdaW_times_v <- function(Lambda_mat, W, v, k, n, verbose = FALSE) {
  v_matrix <- matrix(v, nrow = n, ncol = k)

  result <- numeric(k * n)

  # Block i = Σ_j λij * W * vj
  for (i in 1:k) {
    block_i <- numeric(n)

    for (j in 1:k) {
      block_i <- block_i + Lambda_mat[i, j] * (W %*% v_matrix[, j])
    }

    result[((i-1)*n + 1):(i*n)] <- block_i
  }

  if (verbose) {
    cat(sprintf("(Λ⊗W)v computation complete: result dimension = %d×1\n", length(result)))
  }

  return(result)
}


#' Compute (I - Λ⊗W)v
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param W n×n spatial weight matrix
#' @param v kn×1 vector
#' @param k number of variables
#' @param n number of regions
#' @return kn×1 vector
compute_I_minus_LambdaW_times_v <- function(Lambda_mat, W, v, k, n) {
  LambdaWv <- compute_LambdaW_times_v(Lambda_mat, W, v, k, n, verbose = FALSE)
  return(v - LambdaWv)
}


#' Compute (I - Λ'⊗W')v
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param W n×n spatial weight matrix
#' @param v kn×1 vector
#' @param k number of variables
#' @param n number of regions
#' @return kn×1 vector
compute_I_minus_LambdatWt_times_v <- function(Lambda_mat, W, v, k, n) {
  LambdatWtv <- compute_LambdaW_times_v(t(Lambda_mat), t(W), v, k, n, verbose = FALSE)
  return(v - LambdatWtv)
}


# MSEM (Spatial Error Model)

#' Compute the MSEM log-likelihood
#'
#' ℓ(Λ, β, Σ) = -(kn/2)log(2π) + log|I - Λ⊗W| - (n/2)log|Σ| - Q/2
#'
#' Residuals: z = (I - Λ⊗W)(y - Xβ)
#' Q = z'(Σ^{-1}⊗I)z
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param beta p×1 regression coefficients
#' @param Sigma k×k error covariance matrix
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param eigen_W eigenvalues of W
#' @param k, n dimensions
#' @param verbose verbose output
#' @return value of the log-likelihood
compute_log_likelihood_msem <- function(Lambda_mat, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {

  log_det_Lambda <- log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE)

  if (!is.finite(log_det_Lambda)) {
    if (verbose >= 1) cat("Warning: log|I - Λ⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]

  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Σ| is not finite\n")
    return(-Inf)
  }

  # Residuals z = (I - Λ⊗W)(y - Xβ)
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

  # Q = z'(Σ^{-1}⊗I)z
  Sigma_inv <- solve(Sigma)
  Q <- 0

  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }

  loglik <- -(k*n/2) * log(2*pi) + log_det_Lambda - (n/2) * log_det_Sigma - Q/2

  if (verbose >= 2) {
    cat(sprintf("  log|I - Λ⊗W| = %.6f\n", log_det_Lambda))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }

  return(loglik)
}


#' Compute the MSEM profile likelihood
#'
#' ℓ(Λ) = -(kn/2)log(2πe) + log|I - Λ⊗W| - (n/2)log|Σ̂(Λ)|
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param beta_hat estimated β
#' @param Sigma_hat estimated Σ
#' @param y, X, W data
#' @param eigen_W eigenvalues of W
#' @param k, n dimensions
#' @param verbose verbose output
#' @return profile log-likelihood
compute_profile_likelihood_msem <- function(Lambda_mat, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {

  log_det_Lambda <- log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE, smooth = smooth)

  if (!is.finite(log_det_Lambda)) {
    if (verbose >= 1) cat("Warning: log|I - Λ⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]

  if (!is.finite(log_det_Sigma_hat)) {
    if (verbose >= 1) cat("Warning: log|Σ̂| is not finite\n")
    return(-Inf)
  }

  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_Lambda - (n/2) * log_det_Sigma_hat

  if (verbose >= 2) {
    cat(sprintf("  profile likelihood: %.6f\n", profile_loglik))
  }

  return(profile_loglik)
}


#' Update β for MSEM
#'
#' β̂ = {X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)X}^{-1} X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)y
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param Sigma k×k error covariance matrix
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param k, n dimensions
#' @param ridge_eps ridge regularization parameter
#' @param verbose verbose output
#' @return estimated β
update_beta_msem <- function(Lambda_mat, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: error inverting Σ. Applying regularization\n")
    solve(Sigma + diag(ridge_eps, k))
  })

  p <- ncol(X)

  # (I-Λ⊗W)y and (I-Λ⊗W)X, then form the GLS normal equations blockwise
  ILambdaWy <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, y, k, n)

  ILambdaWX <- matrix(0, nrow = k*n, ncol = p)
  for (col in 1:p) {
    ILambdaWX[, col] <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, X[, col], k, n)
  }

  # X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)X = (ILambdaWX)'(Σ^{-1}⊗I)(ILambdaWX)
  XtAX <- matrix(0, p, p)
  XtAy <- numeric(p)

  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)

      XtAX <- XtAX + Sigma_inv[i,j] * t(ILambdaWX[i_idx, , drop = FALSE]) %*% ILambdaWX[j_idx, , drop = FALSE]
      XtAy <- XtAy + Sigma_inv[i,j] * t(ILambdaWX[i_idx, , drop = FALSE]) %*% ILambdaWy[j_idx]
    }
  }

  beta_hat <- tryCatch({
    solve(XtAX, XtAy)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: singular matrix when estimating β. Applying ridge regularization\n")
    solve(XtAX + diag(ridge_eps, p), XtAy)
  })

  if (verbose >= 2) {
    cat(sprintf("  β update done: ||β|| = %.6f\n", sqrt(sum(beta_hat^2))))
  }

  return(beta_hat)
}


#' Update Σ for MSEM
#'
#' Σ̂ = (1/n) * [ẑ'_i ẑ_j]
#' where ẑ = (I-Λ⊗W)(y - Xβ)
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param beta p×1 regression coefficients
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param k, n dimensions
#' @param ridge_eps regularization parameter
#' @param verbose verbose output
#' @return estimated Σ (symmetric matrix)
update_Sigma_msem <- function(Lambda_mat, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  # Residuals z = (I - Λ⊗W)(y - Xβ)
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

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


#' Compute MSEM residuals
#'
#' z = (I - Λ⊗W)(y - Xβ)
#'
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param beta p×1 regression coefficients
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param k, n dimensions
#' @return kn×1 residual vector
compute_residuals_msem <- function(Lambda_mat, beta, y, X, W, k, n) {
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)
  return(z)
}


#' Iterative update of β and Σ for MSEM
#'
#' Fix Λ and alternately update β and Σ
#'
#' @param Lambda_mat fixed Λ matrix
#' @param beta_init, Sigma_init initial values
#' @param y, X, W data
#' @param k, n dimensions
#' @param max_iter maximum number of iterations
#' @param tol convergence threshold
#' @param ridge_eps ridge regularization parameter
#' @param verbose verbosity level
#' @return list(beta, Sigma, converged, iterations)
iterate_beta_sigma_msem <- function(
  Lambda_mat, 
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

  if (verbose >= 1) {
    cat("  Inner loop start (updating β and Σ)\n")
  }

  for (iter in 1:max_iter) {

    log_det_Sigma_old <- determinant(Sigma, logarithm = TRUE)$modulus[1]

    beta <- update_beta_msem(Lambda_mat, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    Sigma <- update_Sigma_msem(Lambda_mat, beta, y, X, W, k, n, ridge_eps, verbose = verbose)

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
    converged = FALSE,
    iterations = max_iter
  ))
}

# MGNS (General Nesting Spatial Model)

#' Compute the MGNS log-likelihood
#'
#' ℓ(R, Λ, β, Σ) = -(kn/2)log(2π) + log|I-R⊗W| + log|I-Λ⊗W| - (n/2)log|Σ| - Q/2
#'
#' Residuals: z = (I-Λ⊗W)(y - (R⊗W)y - Xβ)
#' Q = z'(Σ^{-1}⊗I)z
#'
#' @param R k×k spatial correlation matrix of the response
#' @param Lambda_mat k×k spatial error correlation matrix
#' @param beta p×1 regression coefficients
#' @param Sigma k×k error covariance matrix
#' @param y kn×1 dependent variable
#' @param X kn×p design matrix
#' @param W n×n spatial weight matrix
#' @param eigen_W eigenvalues of W
#' @param k, n dimensions
#' @param verbose verbose output
#' @return value of the log-likelihood
compute_log_likelihood_mgns <- function(R, Lambda_mat, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {

  log_det_R <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE)

  if (!is.finite(log_det_R)) {
    if (verbose >= 1) cat("Warning: log|I - R⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Lambda <- log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE)

  if (!is.finite(log_det_Lambda)) {
    if (verbose >= 1) cat("Warning: log|I - Λ⊗W| is not finite\n")
    return(-Inf)
  }

  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]

  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Σ| is not finite\n")
    return(-Inf)
  }

  # Residuals z = (I - Λ⊗W)(y - (R⊗W)y - Xβ)
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

  # Q = z'(Σ^{-1}⊗I)z
  Sigma_inv <- solve(Sigma)
  Q <- 0

  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }

  loglik <- -(k*n/2) * log(2*pi) + log_det_R + log_det_Lambda - (n/2) * log_det_Sigma - Q/2

  if (verbose >= 2) {
    cat(sprintf("  log|I - R⊗W| = %.6f\n", log_det_R))
    cat(sprintf("  log|I - Λ⊗W| = %.6f\n", log_det_Lambda))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }

  return(loglik)
}


#' Compute the MGNS profile likelihood
#'
#' ℓ(R,Λ) = -(kn/2)log(2πe) + log|I-R⊗W| + log|I-Λ⊗W| - (n/2)log|Σ̂(R,Λ)|
compute_profile_likelihood_mgns <- function(R, Lambda_mat, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {
  
  log_det_R <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  log_det_Lambda <- log_det_spatial(Lambda_mat, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  
  if (!is.finite(log_det_R) || !is.finite(log_det_Lambda)) {
    return(-Inf)
  }
  
  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma_hat)) {
    return(-Inf)
  }
  
  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_R + log_det_Lambda - (n/2) * log_det_Sigma_hat
  
  return(profile_loglik)
}


#' Update β for MGNS
#'
#' β̂ = {X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)X}^{-1}
#'      X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)(I-R⊗W)y
update_beta_mgns <- function(R, Lambda_mat, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    solve(Sigma + diag(ridge_eps, k))
  })

  p <- ncol(X)

  # (I-R⊗W)y
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  IRWy <- y - RWy

  # (I-Λ⊗W)(I-R⊗W)y
  ILambdaW_IRWy <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, IRWy, k, n)

  # (I-Λ⊗W)X
  ILambdaWX <- matrix(0, nrow = k*n, ncol = p)
  for (col in 1:p) {
    ILambdaWX[, col] <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, X[, col], k, n)
  }

  # X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)X and X'(I-Λ'⊗W')(Σ^{-1}⊗I)(I-Λ⊗W)(I-R⊗W)y
  XtAX <- matrix(0, p, p)
  XtAy <- numeric(p)

  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)

      XtAX <- XtAX + Sigma_inv[i,j] * t(ILambdaWX[i_idx, , drop = FALSE]) %*% ILambdaWX[j_idx, , drop = FALSE]
      XtAy <- XtAy + Sigma_inv[i,j] * t(ILambdaWX[i_idx, , drop = FALSE]) %*% ILambdaW_IRWy[j_idx]
    }
  }

  beta_hat <- tryCatch({
    solve(XtAX, XtAy)
  }, error = function(e) {
    solve(XtAX + diag(ridge_eps, p), XtAy)
  })

  return(beta_hat)
}


#' Update Σ for MGNS
#'
#' z = (I-Λ⊗W)(y - (R⊗W)y - Xβ)
#' Σ̂ = (1/n) * [z'_i z_j]
update_Sigma_mgns <- function(R, Lambda_mat, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {

  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)

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

  eigen_vals <- eigen(Sigma_new, only.values = TRUE)$values
  min_eigen <- min(Re(eigen_vals))

  if (min_eigen <= 0) {
    Sigma_new <- Sigma_new + diag(ridge_eps, k)
  }

  return(Sigma_new)
}


#' Compute MGNS residuals
compute_residuals_mgns <- function(R, Lambda_mat, beta, y, X, W, k, n) {
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_LambdaW_times_v(Lambda_mat, W, residual_raw, k, n)
  return(z)
}


#' Iterative update of β and Σ for MGNS
iterate_beta_sigma_mgns <- function(
  R, Lambda_mat, 
  beta_init, Sigma_init, 
  y, X, W, 
  k, n, 
  max_iter = 100, 
  tol = 1e-6,
  ridge_eps = 1e-6,
  verbose = 0
) {
  
  beta <- beta_init
  Sigma <- Sigma_init
  
  for (iter in 1:max_iter) {
    log_det_Sigma_old <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    beta <- update_beta_mgns(R, Lambda_mat, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    Sigma <- update_Sigma_mgns(R, Lambda_mat, beta, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    log_det_Sigma_new <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    change <- abs(log_det_Sigma_new - log_det_Sigma_old)
    
    if (change < tol) {
      return(list(beta = beta, Sigma = Sigma, converged = TRUE, iterations = iter))
    }
  }
  
  return(list(beta = beta, Sigma = Sigma, converged = FALSE, iterations = max_iter))
}
