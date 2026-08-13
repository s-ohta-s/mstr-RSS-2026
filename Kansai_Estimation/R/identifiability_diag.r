# identifiability_diag.r
# Numerical diagnostics for the identifiability conditions (Propositions 1 and 2,
# paper S3.2). Both give SUFFICIENT conditions; these diagnostics quantify
# whether -- and by how much margin -- they hold.
#   Proposition 1 (R, beta):  rank[G_mu, X] = K^2 + p  (full column rank), where
#     G_mu stacks the spatially lagged means (paper Eq. (10)).
#     * synthetic data: mu = mu0 = (I - R x W)^{-1} X beta at the TRUE values
#     * empirical data: mu = mu_hat, the FITTED mean vector
#   Proposition 2 (Lambda, Sigma):  {I_n, W, W', W'W} linearly independent.


#' Build G_mu (Kn x K^2) from the stacked mean vector mu (paper Eq. (10))
#'
#' Row block i (i = 1..K) contains [W mu_1, ..., W mu_K] placed in columns
#' ((i-1)K+1):(iK); all other entries are zero.
#'
#' @param mu Kn x 1 stacked mean vector
#' @param W n x n spatial weight matrix
#' @param k number of responses
#' @param n number of regions
#' @return Kn x K^2 matrix
build_G_mu <- function(mu, W, k, n) {
  mu_mat <- matrix(mu, nrow = n, ncol = k)
  Wmu <- as.matrix(W %*% mu_mat)              # n x k: [W mu_1, ..., W mu_K]
  G <- matrix(0, nrow = k * n, ncol = k * k)
  for (i in 1:k) {
    rows <- ((i - 1) * n + 1):(i * n)
    cols <- ((i - 1) * k + 1):(i * k)
    G[rows, cols] <- Wmu
  }
  G
}


#' Mean vector mu = (I - R x W)^{-1} X beta
#'
#' Falls back to X beta when R is NULL or the zero matrix (MSEM / VARX / OLS),
#' which is both the correct value and a cheaper path.
#'
#' @param R k x k spatial lag matrix (NULL or zero matrix if absent)
#' @param beta p-vector of coefficients (ordering must match the columns of X)
#' @param X Kn x p design matrix
#' @param W n x n spatial weight matrix
#' @param k,n dimensions
#' @return Kn-vector
compute_mu <- function(R, beta, X, W, k, n) {
  Xbeta <- as.numeric(X %*% beta)
  if (is.null(R) || all(R == 0)) return(Xbeta)
  S <- diag(k * n) - kronecker(R, as.matrix(W))
  as.numeric(solve(S, Xbeta))
}


#' Proposition 1 margin: rank and scaled minimum singular value of [G_mu, X]
#'
#' Columns are normalized to unit length before the SVD so the margin is
#' invariant to the (large) scale differences between covariates.
#'
#' @param mu Kn x 1 mean vector (mu0 for synthetic data, mu_hat for empirical)
#' @param X Kn x p realized design matrix (including intercepts and y_lag columns)
#' @param W n x n spatial weight matrix
#' @param k,n dimensions
#' @return list(rank_ok, rank, required, min_sv, max_sv, cond)
prop1_rank_margin_from_mu <- function(mu, X, W, k, n) {

  G <- build_G_mu(mu, W, k, n)
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
    max_sv   = max(sv),
    cond     = max(sv) / min(sv)
  )
}


#' Proposition 1 margin from (R, beta): builds mu internally
#'
#' Convenience wrapper matching the synthetic-data usage, where mu0 is computed
#' from the true parameters.
#'
#' @param R k x k spatial lag matrix (zero matrix / NULL when absent)
#' @param beta p-vector of coefficients (ordering must match X)
#' @param X Kn x p realized design matrix
#' @param W n x n spatial weight matrix
#' @param k,n dimensions
#' @return list(rank_ok, rank, required, min_sv, max_sv, cond)
prop1_rank_margin <- function(R, beta, X, W, k, n) {
  mu <- compute_mu(R, beta, X, W, k, n)
  prop1_rank_margin_from_mu(mu, X, W, k, n)
}


#' Proposition 2 diagnostic: linear independence of {I_n, W, W', W'W}
#'
#' Computed once per W. Uses the normalized Gram matrix of the vectorized
#' matrices; a positive minimum eigenvalue certifies linear independence, and
#' its magnitude is the margin (reported because it shrinks with n).
#'
#' min_sv is the smallest singular value of the column-normalized matrix
#' [vec(I_n), vec(W), vec(W'), vec(W'W)] and equals sqrt(min_eig); it is
#' reported as well so that either wording can be used in the manuscript.
#'
#' @param W n x n spatial weight matrix (dense)
#' @return list(min_eig, min_sv, eigvals, W_symmetric)
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
    min_sv  = sqrt(max(min(ev), 0)),
    eigvals = ev,
    W_symmetric = isTRUE(all.equal(W, t(W), tolerance = 1e-12))
  )
}
