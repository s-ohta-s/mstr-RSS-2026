# identifiability_diag.r
# Numerical diagnostics for the identifiability conditions of §3.2 of the paper
# Both propositions give SUFFICIENT conditions;
# these diagnostics quantify whether — and by how much margin — they hold in
# each simulated cell.
#   Proposition 1 (R, β):  rank[G_µ0, X] = K² + p  (full column rank), where
#     µ0 = (I − R⊗W)^{-1} X β is the true conditional mean and G_µ0 stacks the
#     spatially lagged means (Eq. (10) of the paper).
#   Proposition 2 (Λ, Σ):  {I_n, W, W', W'W} linearly independent.

#' Build G_µ0 (Kn × K²) from the stacked mean vector µ0 (paper Eq. (10))
#'
#' Row block i (i = 1..K) contains [W µ_1, ..., W µ_K] placed in columns
#' ((i−1)K+1):(iK); all other entries are zero.
#'
#' @param mu0 Kn×1 stacked mean vector
#' @param W n×n spatial weight matrix (dense or sparse)
#' @param k number of responses
#' @param n number of regions
#' @return Kn × K² matrix
build_G_mu0 <- function(mu0, W, k, n) {
  mu_mat <- matrix(mu0, nrow = n, ncol = k)
  Wmu <- as.matrix(W %*% mu_mat)              # n × k: [Wµ_1, ..., Wµ_K]
  G <- matrix(0, nrow = k * n, ncol = k * k)
  for (i in 1:k) {
    rows <- ((i - 1) * n + 1):(i * n)
    cols <- ((i - 1) * k + 1):(i * k)
    G[rows, cols] <- Wmu
  }
  G
}

#' Proposition 1 margin: rank and scaled minimum singular value of [G_µ0, X]
#'
#' µ0 is computed from the TRUE parameters (available in simulation), so this
#' verifies the population-level rank condition on the realized design.
#' Columns are normalized to unit length before the SVD so the margin is
#' invariant to the (large) scale differences between covariates.
#'
#' @param R_true k×k true spatial lag matrix (zero matrix for 0111 DGPs)
#' @param beta_true p-vector of true coefficients (ordering must match X)
#' @param X Kn×p realized design matrix (including intercepts and y_lag columns)
#' @param W n×n spatial weight matrix
#' @param k, n dimensions
#' @return list(rank_ok, rank, required, min_sv, cond)
prop1_rank_margin <- function(R_true, beta_true, X, W, k, n) {

  Xbeta <- as.numeric(X %*% beta_true)

  # µ0 = (I − R⊗W)^{-1} Xβ  (identity solve when R = 0)
  if (all(R_true == 0)) {
    mu0 <- Xbeta
  } else {
    S <- diag(k * n) - kronecker(R_true, as.matrix(W))
    mu0 <- as.numeric(solve(S, Xbeta))
  }

  G <- build_G_mu0(mu0, W, k, n)
  M <- cbind(G, X)

  # Column normalization (scale-invariant margin)
  norms <- sqrt(colSums(M^2))
  norms[norms == 0] <- 1     # an all-zero column will show up as rank deficiency
  Mn <- sweep(M, 2, norms, "/")

  sv <- svd(Mn, nu = 0, nv = 0)$d
  required <- ncol(M)
  rank_val <- sum(sv > max(dim(Mn)) * .Machine$double.eps * max(sv))

  list(
    rank_ok  = (rank_val == required),
    rank     = rank_val,
    required = required,
    min_sv   = min(sv),
    cond     = max(sv) / min(sv)
  )
}

#' Proposition 2 diagnostic: linear independence of {I_n, W, W', W'W}
#'
#' Computed once per W. Uses the normalized Gram matrix of the vectorized
#' matrices; a positive minimum eigenvalue certifies linear independence, and
#' its magnitude is the margin (reported because it shrinks with n).
#'
#' @param W n×n spatial weight matrix (dense)
#' @return list(min_eig, eigvals, W_symmetric)
prop2_gram_diag <- function(W) {
  W <- as.matrix(W)
  mats <- list(diag(nrow(W)), W, t(W), crossprod(W))   # I, W, W', W'W
  G <- matrix(0, 4, 4)
  for (i in 1:4) {
    for (j in i:4) {
      G[i, j] <- sum(mats[[i]] * mats[[j]])
      G[j, i] <- G[i, j]
    }
  }
  Gn <- G / sqrt(outer(diag(G), diag(G)))
  ev <- eigen(Gn, symmetric = TRUE, only.values = TRUE)$values
  list(
    min_eig = min(ev),
    eigvals = ev,
    W_symmetric = isTRUE(all.equal(W, t(W), tolerance = 1e-12))
  )
}

#' Minimum eigenvalue of the unpenalized profile Hessian at the estimate
#'
#' Local-identification diagnostic: positive definiteness of the observed
#' information for the spatial parameters at γ = 0. Reuses the existing
#' compute_profile_hessian_* machinery (penalized_spatial.r).
#'
#' @param result build_result_object() output
#' @param data_list data_list used for the fit
#' @return smallest eigenvalue of the profile Hessian (NA on failure)
profile_hessian_min_eig <- function(result, data_list) {
  H <- tryCatch({
    if (result$model_type == "MSAR") {
      compute_profile_hessian_msar(
        R_hat = result$coefficients$R,
        beta_hat = result$coefficients$beta,
        Sigma_hat = result$coefficients$Sigma,
        data_list = data_list)
    } else if (result$model_type == "MSEM") {
      compute_profile_hessian_msem(
        Lambda_hat = result$coefficients$Lambda,
        beta_hat = result$coefficients$beta,
        Sigma_hat = result$coefficients$Sigma,
        data_list = data_list)
    } else {
      compute_profile_hessian_mgns(
        R_hat = result$coefficients$R,
        Lambda_hat = result$coefficients$Lambda,
        beta_hat = result$coefficients$beta,
        Sigma_hat = result$coefficients$Sigma,
        data_list = data_list)
    }
  }, error = function(e) NULL)
  if (is.null(H)) return(NA_real_)
  min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
}
