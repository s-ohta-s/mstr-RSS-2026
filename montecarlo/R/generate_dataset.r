# generate_dataset.r
# In-memory data generators for the Monte Carlo experiments
# Both generators return list(df_current, df_lag): one row per region, columns
# region, time, y1, y2, x_common1, x_common2, x_specific1_1, x_specific2_1 —
# directly consumable by attach_data() (data_prep_fast.r). Nothing is written
# to disk.

#' Queen-contiguity row-standardized weight matrix for a square lattice
#'
#'
#' @param nrow_grid,ncol_grid lattice dimensions
#' @return n×n row-standardized dense matrix
create_queen_weight_matrix <- function(nrow_grid, ncol_grid) {
  n <- nrow_grid * ncol_grid
  W <- matrix(0, n, n)
  for (i in 1:n) {
    row_i <- (i - 1) %/% ncol_grid + 1
    col_i <- (i - 1) %% ncol_grid + 1
    for (j in 1:n) {
      if (i != j) {
        row_j <- (j - 1) %/% ncol_grid + 1
        col_j <- (j - 1) %% ncol_grid + 1
        if (abs(row_i - row_j) <= 1 && abs(col_i - col_j) <= 1) W[i, j] <- 1
      }
    }
  }
  W / rowSums(W)
}

#' Load the per-n weight matrix CSV, creating and saving it on first use
#'
#' n must be a perfect square (10×10 / 20×20 / 30×30 lattices).
#'
#' @param n number of regions (100 / 400 / 900)
#' @param dir directory holding spatial_weights_n{N}.csv
#' @return dense W matrix
get_or_create_W <- function(n, dir = file.path("data", "simulated")) {
  path <- file.path(dir, sprintf("spatial_weights_n%d.csv", n))
  if (file.exists(path)) {
    W <- as.matrix(read.csv(path, header = TRUE))
    dimnames(W) <- NULL
    return(W)
  }
  m <- sqrt(n)
  if (m != round(m)) stop(sprintf("n = %d is not a square lattice size", n))
  W <- create_queen_weight_matrix(m, m)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  write.csv(W, path, row.names = FALSE)
  W
}

#' Cached sparse solvers for (I − R⊗W) and (I − Λ⊗W)
#'
#' The spatial systems are fixed within a (paramset, dgp, n) cell, so the LU
#' factorizations are computed once and reused across periods and replications.
#'
#' @param params list with R, Lambda (k×k)
#' @param static prepare_static() output
#' @return list(SR_solve = function(v) ..., SLam_solve = function(v) ...)
make_spatial_solvers <- function(params, static) {
  k <- nrow(params$R)
  n <- static$n
  W_sp <- if (!is.null(static$W_sp)) static$W_sp else Matrix::Matrix(static$W, sparse = TRUE)
  I_kn <- Matrix::Diagonal(k * n)

  make_one <- function(M) {
    if (all(M == 0)) return(identity)
    # Cache the sparse system matrix; solve(<sparseLU>, .) is unimplemented in
    # this Matrix version, so let solve(<dgCMatrix>, v) factorize internally
    # (kn ≤ 1800: milliseconds per call)
    S <- I_kn - Matrix::kronecker(M, W_sp)
    function(v) as.numeric(Matrix::solve(S, v))
  }

  list(
    SR_solve   = make_one(params$R),
    SLam_solve = make_one(params$Lambda)
  )
}

#' Draw the covariate columns for one period
#'
#' weak/mid/strong: i.i.d. N(0,1). realistic: N(mean, sd²) matched to the real
#' covariates (params$x_moments).
#'
#' @param n number of regions
#' @param x_moments NULL (standard normal) or the realistic x_moments list
#' @return data.frame with the four covariate columns
draw_covariates <- function(n, x_moments = NULL) {
  if (is.null(x_moments)) {
    data.frame(
      x_common1     = rnorm(n),
      x_common2     = rnorm(n),
      x_specific1_1 = rnorm(n),
      x_specific2_1 = rnorm(n)
    )
  } else {
    data.frame(
      x_common1     = rnorm(n, x_moments$x_common1["mean"],     x_moments$x_common1["sd"]),
      x_common2     = rnorm(n, x_moments$x_common2["mean"],     x_moments$x_common2["sd"]),
      x_specific1_1 = rnorm(n, x_moments$x_specific1_1["mean"], x_moments$x_specific1_1["sd"]),
      x_specific2_1 = rnorm(n, x_moments$x_specific2_1["mean"], x_moments$x_specific2_1["sd"])
    )
  }
}

#' Xβ for one period (K = 2, stacked as c(block y1, block y2))
.xbeta_stacked <- function(Xdf, p) {
  c(p$beta_intercept[1] +
      p$beta_common1[1] * Xdf$x_common1 +
      p$beta_common2[1] * Xdf$x_common2 +
      p$beta_specific1_1 * Xdf$x_specific1_1,
    p$beta_intercept[2] +
      p$beta_common1[2] * Xdf$x_common1 +
      p$beta_common2[2] * Xdf$x_common2 +
      p$beta_specific2_1 * Xdf$x_specific2_1)
}

#' One model transition: y_t = SR^{-1}{Xβ + (A⊗I)y_{t-1} + SΛ^{-1}ε}
.one_transition <- function(Xdf, y_prev_mat, p, solvers, chol_Sigma, n) {
  eps <- matrix(rnorm(2 * n), n, 2) %*% chol_Sigma        # rows ~ N(0, Σ)
  u <- solvers$SLam_solve(as.numeric(eps))                # (I−Λ⊗W)^{-1}ε
  a_term <- c(p$A[1, 1] * y_prev_mat[, 1] + p$A[1, 2] * y_prev_mat[, 2],
              p$A[2, 1] * y_prev_mat[, 1] + p$A[2, 2] * y_prev_mat[, 2])
  rhs <- .xbeta_stacked(Xdf, p) + a_term + u
  y <- solvers$SR_solve(rhs)                              # (I−R⊗W)^{-1}rhs
  matrix(y, n, 2)
}

#' Burn-in generator (stationary parameter sets: weak / mid / strong)
#'
#' Adaptive burn-in: burnin = max(burnin_floor, ceil(log(eps)/log ρ(M_dyn))).
#' Refuses to generate when ρ(M_dyn) ≥ 1 (the recursion would diverge) — use
#' generate_dataset_conditional() for such parameter sets.
#'
#' @param params get_dgp_params() output (generation = "burnin")
#' @param static prepare_static() output
#' @param seed RNG seed for this replication
#' @param burnin "auto" or an integer
#' @param eps residual initial-condition influence tolerated after burn-in
#' @param burnin_floor minimum burn-in length
#' @param solvers optional make_spatial_solvers() output (cached per cell)
#' @return list(df_current (time = 2), df_lag (time = 1), burnin, rho_dyn)
generate_dataset <- function(params, static, seed,
                             burnin = "auto", eps = 0.005, burnin_floor = 50,
                             solvers = NULL) {

  n <- static$n
  set.seed(seed)

  rho_M <- temporal_spectral_radius(params$R, params$A, static$eigen_W)
  if (rho_M >= 1) {
    stop(sprintf(
      "generate_dataset: rho(M_dyn) = %.4f >= 1 (explosive); use generate_dataset_conditional()",
      rho_M))
  }
  if (identical(burnin, "auto")) {
    burnin <- if (rho_M <= 0) burnin_floor
              else max(burnin_floor, ceiling(log(eps) / log(rho_M)))
  }

  if (is.null(solvers)) solvers <- make_spatial_solvers(params, static)
  chol_Sigma <- chol(params$Sigma)

  T_keep <- 2
  T_total <- burnin + T_keep

  y_prev <- matrix(0, n, 2)          # start at zero; burn-in washes it out
  kept <- vector("list", T_keep)

  for (t in 1:T_total) {
    Xdf <- draw_covariates(n)
    y_now <- .one_transition(Xdf, y_prev, params, solvers, chol_Sigma, n)
    if (t > burnin) {
      tk <- t - burnin
      kept[[tk]] <- data.frame(
        region = 1:n, time = tk,
        y1 = y_now[, 1], y2 = y_now[, 2],
        Xdf
      )
    }
    y_prev <- y_now
  }

  list(df_current = kept[[2]], df_lag = kept[[1]],
       burnin = burnin, rho_dyn = rho_M)
}

#' Conditional single-transition generator (realistic parameter set)
#'
#' y_{t-1} ~ N₂(m, V) i.i.d. over regions (V from the real-data moments,
#' including the y1–y2 cross correlation), X from the real covariate moments,
#' then ONE model transition. Well-defined even for explosive A; matches the
#' paper's conditional inference f(y_t | y_{t-1}, X_t) exactly.
#'
#' @param params get_dgp_params("realistic", dgp) output (generation = "conditional")
#' @param static prepare_static() output
#' @param seed RNG seed for this replication
#' @param solvers optional make_spatial_solvers() output (cached per cell)
#' @return list(df_current (time = 2), df_lag (time = 1), rho_dyn)
generate_dataset_conditional <- function(params, static, seed, solvers = NULL) {

  n <- static$n
  set.seed(seed)

  if (is.null(solvers)) solvers <- make_spatial_solvers(params, static)
  chol_Sigma <- chol(params$Sigma)

  # y_{t-1} ~ N2(m, V) with the real-data cross correlation
  m  <- params$y_prev_moments$mean
  sd <- params$y_prev_moments$sd
  r  <- params$y_prev_moments$cor12
  V  <- matrix(c(sd[1]^2, r * sd[1] * sd[2],
                 r * sd[1] * sd[2], sd[2]^2), 2, 2)
  y_prev <- matrix(rnorm(2 * n), n, 2) %*% chol(V)
  y_prev[, 1] <- y_prev[, 1] + m[1]
  y_prev[, 2] <- y_prev[, 2] + m[2]

  # Covariates: lag-period draw only fills the format (estimation uses X_t);
  # current-period draw enters the transition
  X_lag <- draw_covariates(n, params$x_moments)
  X_cur <- draw_covariates(n, params$x_moments)

  y_now <- .one_transition(X_cur, y_prev, params, solvers, chol_Sigma, n)

  rho_M <- temporal_spectral_radius(params$R, params$A, static$eigen_W)

  list(
    df_current = data.frame(region = 1:n, time = 2,
                            y1 = y_now[, 1], y2 = y_now[, 2], X_cur),
    df_lag     = data.frame(region = 1:n, time = 1,
                            y1 = y_prev[, 1], y2 = y_prev[, 2], X_lag),
    rho_dyn = rho_M
  )
}

#' Dispatch on the parameter set's generation mode
#'
#' @param params get_dgp_params() output
#' @param static prepare_static() output
#' @param seed RNG seed
#' @param solvers optional cached solvers
#' @return list(df_current, df_lag, ...) from the appropriate generator
generate_cell_dataset <- function(params, static, seed, solvers = NULL, ...) {
  if (identical(params$generation, "conditional")) {
    generate_dataset_conditional(params, static, seed, solvers = solvers)
  } else {
    generate_dataset(params, static, seed, solvers = solvers, ...)
  }
}
