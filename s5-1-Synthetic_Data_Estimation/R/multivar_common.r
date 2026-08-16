# multivar_common.r
# Foundational functions shared by all models: spatial-matrix utilities
# Dependencies: none (foundational file loaded first in the source chain)

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

  # Compute the traces: r_j = tr(R^j)
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

  # Vectorized over eigenvalues: ∑(j=0 to k) (-1)^j ej(R) ω^j for all ω at once.
  # omega_power is built by successive multiplication (matches the scalar loop).
  poly <- .poly_values_spatial(e, eigen_W, k)

  if (smooth) {
    # Smooth extension: C1-continuous linear continuation below delta
    vals <- ifelse(poly > delta, log(pmax(poly, delta)), log(delta) + (poly - delta) / delta)
  } else {
    # Legacy behavior: return -Inf out of domain
    if (any(poly <= 0)) {
      if (verbose) {
        i_bad <- which(poly <= 0)[1]
        cat(sprintf("Warning: poly_value = %.6e <= 0 at i=%d\n", poly[i_bad], i_bad))
      }
      return(-Inf)
    }
    vals <- log(poly)
  }

  # Sequential accumulation (not sum()): R's sum() uses an extended-precision
  # accumulator, which perturbs the log-likelihood at ~1e-15 relative to the
  # historical per-eigenvalue loop. That noise is amplified ~1e9× by the
  # Richardson-extrapolated numerical Hessians (df_eff, se_hessian), so exact
  # reproduction of the pre-refactor sequence matters here.
  log_det <- 0
  for (i in seq_along(vals)) log_det <- log_det + vals[i]

  if (verbose) {
    cat(sprintf("log|I - R⊗W| = %.6f\n", log_det))
  }

  return(log_det)
}

#' Polynomial values ∑_j (-1)^j e_j ω_i^j for all eigenvalues at once (internal)
#'
#' Vectorized counterpart of the per-eigenvalue scalar loop. Powers of ω are
#' accumulated by successive multiplication so the floating-point sequence is
#' identical to the original implementation.
#'
#' @param e vector (e0, ..., ek) from compute_elem_sym_poly
#' @param eigen_W eigenvalues of W (possibly complex)
#' @param k number of variables
#' @return real vector of polynomial values (length = length(eigen_W))
.poly_values_spatial <- function(e, eigen_W, k) {
  poly <- rep.int(1, length(eigen_W)) * e[1]   # j = 0 term (e0 = 1)
  omega_power <- rep.int(1, length(eigen_W))
  if (is.complex(eigen_W)) {
    poly <- as.complex(poly)
    omega_power <- as.complex(omega_power)
  }
  for (j in 1:k) {
    omega_power <- omega_power * eigen_W
    poly <- poly + ((-1)^j) * e[j + 1] * omega_power
  }
  Re(poly)
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
  min_poly <- min(.poly_values_spatial(e, eigen_W, k))
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

  # Reshape y into an n×k matrix, apply W once to all columns, then combine
  # blocks: block_i = Σ_j R[i,j]·(W y_j)  <=>  (W·Y) %*% t(R), stacked by column.
  # Works for both dense and Matrix::sparse W.
  y_matrix <- matrix(y, nrow = n, ncol = k)
  WY <- as.matrix(W %*% y_matrix)
  result <- as.numeric(WY %*% t(R))

  if (verbose) {
    cat(sprintf("(R⊗W)y computation complete: result dimension = %d×1\n", length(result)))
  }

  return(result)
}

#' Compute (M⊗W)V for a kn×p matrix V (all columns at once)
#'
#' Matrix version of compute_RW_times_y / compute_LambdaW_times_v: applies the
#' Kronecker operator to every column of V with k sparse/dense W-multiplications
#' instead of k² per column.
#'
#' @param M k×k coefficient matrix
#' @param W n×n spatial weight matrix (dense or Matrix::sparse)
#' @param V kn×p matrix
#' @param k number of variables
#' @param n number of regions
#' @return kn×p matrix (M⊗W)V
compute_MW_times_M <- function(M, W, V, k, n) {
  p <- ncol(V)
  WV <- vector("list", k)
  for (j in 1:k) {
    idx_j <- ((j - 1) * n + 1):(j * n)
    WV[[j]] <- as.matrix(W %*% V[idx_j, , drop = FALSE])
  }
  out <- matrix(0, nrow = k * n, ncol = p)
  for (i in 1:k) {
    acc <- M[i, 1] * WV[[1]]
    if (k >= 2) {
      for (j in 2:k) acc <- acc + M[i, j] * WV[[j]]
    }
    out[((i - 1) * n + 1):(i * n), ] <- acc
  }
  out
}

#' Precompute GLS cross blocks for the inner β update (internal)
#'
#' For fixed spatial parameters the transformed design SX and target Sy do not
#' change across the inner β/Σ iterations; only Σ⁻¹ does. Precomputing the k²
#' cross blocks t(SX_i)·SX_j and t(SX_i)·Sy_j removes the dominant per-iteration
#' cost. The assembled XtAX/XtAy are numerically identical to the direct loops.
#'
#' @param SX kn×p transformed design matrix
#' @param Sy kn×1 transformed target vector
#' @param k number of variables
#' @param n number of regions
#' @return list(B = k×k list of p×p blocks, by = k×k list of p-vectors)
precompute_gls_blocks <- function(SX, Sy, k, n) {
  B  <- vector("list", k * k)
  by <- vector("list", k * k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    SXi_t <- t(SX[i_idx, , drop = FALSE])
    for (j in 1:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      B[[(i - 1) * k + j]]  <- SXi_t %*% SX[j_idx, , drop = FALSE]
      by[[(i - 1) * k + j]] <- SXi_t %*% Sy[j_idx]
    }
  }
  list(B = B, by = by)
}

#' Assemble XtAX / XtAy from precomputed blocks and Σ⁻¹ (internal)
#'
#' Accumulation order (i outer, j inner) matches the original double loops.
#'
#' @param blocks output of precompute_gls_blocks
#' @param Sigma_inv k×k inverse error covariance
#' @param k number of variables
#' @return list(XtAX = p×p, XtAy = p-vector)
assemble_gls_from_blocks <- function(blocks, Sigma_inv, k) {
  p <- ncol(blocks$B[[1]])
  XtAX <- matrix(0, p, p)
  XtAy <- numeric(p)
  for (i in 1:k) {
    for (j in 1:k) {
      idx <- (i - 1) * k + j
      XtAX <- XtAX + Sigma_inv[i, j] * blocks$B[[idx]]
      XtAy <- XtAy + Sigma_inv[i, j] * blocks$by[[idx]]
    }
  }
  list(XtAX = XtAX, XtAy = XtAy)
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
