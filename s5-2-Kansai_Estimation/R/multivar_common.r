# multivar_common.r
#
# Foundational functions shared by all models: spatial-matrix utilities

# 1. Spatial-matrix utilities

#' Compute the eigenvalues of the spatial weight matrix
#'
#' Precomputed and stored for use in the determinant computation.
#'
#' @param W n×n spatial weight matrix
#' @param verbose verbose output
#' @return vector of eigenvalues (may be complex)
compute_eigen_W <- function(W, verbose = FALSE) {
  eigen_result <- eigen(W, only.values = TRUE)
  eigen_values <- eigen_result$values

  if (verbose) {
    cat(sprintf("Number of eigenvalues: %d\n", length(eigen_values)))
    cat(sprintf("Number of real eigenvalues: %d\n", sum(Im(eigen_values) == 0)))
    cat(sprintf("Number of complex eigenvalues: %d\n", sum(Im(eigen_values) != 0)))
  }

  return(eigen_values)
}


#' Compute the elementary symmetric polynomials
#'
#' ej(R) = elementary symmetric polynomial (symmetric function of the eigenvalues of R)
#' Computed from the traces r_j = tr(R^j)
#'
#' @param R k×k matrix
#' @param verbose verbose output
#' @return vector of length k+1 (e0, e1, ..., ek)
compute_elem_sym_poly <- function(R, verbose = FALSE) {
  k <- nrow(R)

  r <- numeric(k)
  R_power <- R
  for (j in 1:k) {
    r[j] <- sum(diag(R_power))
    if (j < k) {
      R_power <- R_power %*% R
    }
  }

  e <- numeric(k + 1)
  e[1] <- 1  # e0 = 1

  if (k >= 1) {
    e[2] <- r[1]  # e1 = r1
  }

  if (k >= 2) {
    e[3] <- (r[1]^2 - r[2]) / 2  # e2
  }

  if (k >= 3) {
    e[4] <- (r[1]^3 - 3*r[1]*r[2] + 2*r[3]) / 6  # e3
  }

  if (k >= 4) {
    e[5] <- (r[1]^4 - 6*r[1]^2*r[2] + 3*r[2]^2 + 8*r[1]*r[3] - 6*r[4]) / 24  # e4
  }

  if (k >= 5) {
    e[6] <- (r[1]^5 - 10*r[1]^3*r[2] + 15*r[1]*r[2]^2 + 20*r[1]^2*r[3] -
               20*r[2]*r[3] - 30*r[1]*r[4] + 24*r[5]) / 120  # e5
  }

  if (k >= 6) {
    e[7] <- (r[1]^6 - 15*r[1]^4*r[2] + 45*r[1]^2*r[2]^2 - 15*r[2]^3 +
               40*r[1]^3*r[3] - 120*r[1]*r[2]*r[3] + 40*r[3]^2 -
               90*r[1]^2*r[4] + 90*r[2]*r[4] + 144*r[1]*r[5] - 120*r[6]) / 720  # e6
  }

  if (k >= 7) {
    e[8] <- (r[1]^7 - 21*r[1]^5*r[2] + 105*r[1]^3*r[2]^2 - 105*r[1]*r[2]^3 +
               70*r[1]^4*r[3] - 420*r[1]^2*r[2]*r[3] + 280*r[2]^2*r[3] +
               210*r[1]*r[3]^2 - 210*r[1]^3*r[4] + 630*r[1]*r[2]*r[4] -
               420*r[3]*r[4] + 504*r[1]^2*r[5] - 504*r[2]*r[5] -
               840*r[1]*r[6] + 720*r[7]) / 5040  # e7
  }

  if (k > 7) {
    warning("Computation of elementary symmetric polynomials is not implemented for k > 7")
  }

  if (verbose) {
    cat("Elementary symmetric polynomials:\n")
    for (i in 0:min(k, 7)) {
      cat(sprintf("  e%d = %.6f\n", i, e[i+1]))
    }
  }

  return(e)
}


#' Efficient computation of log|Ikn - R⊗W|
#'
#' |Ikn - R⊗W| = ∏(i=1 to n) [∑(j=0 to k) (-1)^j ej(R) ωi^j]
#'
#' @param R k×k spatial correlation matrix
#' @param eigen_W eigenvalues of W (precomputed)
#' @param k number of variables
#' @param n number of regions
#' @param verbose verbose output
#' @param smooth if TRUE, use a C1-continuous linear extension near the
#'   stationarity boundary. Use TRUE during optimization and FALSE when
#'   reporting the final likelihood.
#' @return value of log|Ikn - R⊗W|
log_det_spatial <- function(R, eigen_W, k, n, verbose = FALSE, smooth = FALSE) {

  e <- compute_elem_sym_poly(R, verbose = FALSE)

  # Junction point of the smooth extension
  delta <- 1e-8

  log_det <- 0

  for (i in 1:n) {
    omega_i <- eigen_W[i]

    poly_value <- 0
    omega_power <- 1

    for (j in 0:k) {
      poly_value <- poly_value + ((-1)^j) * e[j+1] * omega_power
      omega_power <- omega_power * omega_i
    }

    if (is.complex(poly_value)) {
      poly_value <- Re(poly_value)
    }

    if (smooth) {
      if (poly_value > delta) {
        log_det <- log_det + log(poly_value)
      } else {
        # Out of domain / near boundary: extend along the tangent at log(delta)
        # value and first derivative match at x = delta (C1-continuous)
        log_det <- log_det + log(delta) + (poly_value - delta) / delta
      }
    } else {
      if (poly_value <= 0) {
        if (verbose) {
          cat(sprintf("Warning: poly_value = %.6e <= 0 at i=%d\n", i, poly_value))
        }
        return(-Inf)
      }
      log_det <- log_det + log(poly_value)
    }
  }

  if (verbose) {
    cat(sprintf("log|I - R⊗W| = %.6f\n", log_det))
  }

  return(log_det)
}


#' Stationarity-violation measure based on the spectral radius
#'
#' M is inside the stationarity region <=> return value < 1
#' For each eigenvalue ω_i, det(I_k - M·ω_i) > 0 is judged via
#' min_i det(I_k - M·ω_i)
#'
#' @param M k×k spatial parameter matrix
#' @param eigen_W vector of eigenvalues of W
#' @return stationarity-violation measure (< 1: stationary, >= 1: non-stationary)
spectral_radius <- function(M, eigen_W) {
  e <- compute_elem_sym_poly(M, verbose = FALSE)
  k <- nrow(M)
  min_poly <- Inf
  for (i in seq_along(eigen_W)) {
    omega_i <- eigen_W[i]
    poly_value <- 0; omega_power <- 1
    for (j in 0:k) {
      poly_value <- poly_value + ((-1)^j) * e[j+1] * omega_power
      omega_power <- omega_power * omega_i
    }
    if (is.complex(poly_value)) poly_value <- Re(poly_value)
    min_poly <- min(min_poly, poly_value)
  }
  if (min_poly > 0) {
    return(1 - min_poly / (1 + min_poly))  # 0 < result < 1
  } else {
    return(1 + abs(min_poly))  # result >= 1
  }
}


#' Efficient computation of (R⊗W)y
#'
#' Implemented with block operations instead of forming the Kronecker product
#' explicitly: reshape y into an n×k matrix and compute each block
#'
#' @param R k×k matrix
#' @param W n×n matrix
#' @param y kn×1 vector
#' @param k number of variables
#' @param n number of regions
#' @param verbose verbose output
#' @return kn×1 vector
compute_RW_times_y <- function(R, W, y, k, n, verbose = FALSE) {

  y_matrix <- matrix(y, nrow = n, ncol = k)

  result <- numeric(k * n)

  for (i in 1:k) {
    block_i <- numeric(n)

    for (j in 1:k) {
      # ρij * W * yj
      block_i <- block_i + R[i, j] * (W %*% y_matrix[, j])
    }

    result[((i-1)*n + 1):(i*n)] <- block_i
  }

  if (verbose) {
    cat(sprintf("(R⊗W)y computation complete: result dimension = %d×1\n", length(result)))
  }

  return(result)
}


#' Check the stationarity condition
#'
#' Simplified condition: |ρij| < 1 for all i,j
#'
#' @param R k×k matrix
#' @param verbose verbose output
#' @return logical (TRUE if the condition is satisfied)
check_stationarity <- function(R, verbose = FALSE) {

  max_abs <- max(abs(R))

  if (verbose) {
    cat(sprintf("Stationarity check: max|ρij| = %.6f\n", max_abs))
  }

  is_stationary <- max_abs < 1

  if (!is_stationary && verbose) {
    cat("Warning: stationarity condition not satisfied (|ρij| < 1)\n")
  }

  return(is_stationary)
}
