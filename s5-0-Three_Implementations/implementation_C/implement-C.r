# ==============================================================================
# Complete Implementation based on mstr.pdf
# (Implement-C)

# ==============================================================================
# 0.0 Timer Initialization
# ==============================================================================
# Record the system time when the script starts to measure total elapsed time
script_start_time <- Sys.time()
cat("Execution started at:", format(script_start_time, "%Y-%m-%d %H:%M:%S"), "\n\n")

# ==============================================================================
# 0. Load Required Packages
# ==============================================================================
packages <- c("Matrix", "parallel", "doParallel", "foreach", "openxlsx", "numDeriv")
new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if (length(new_packages) > 0) install.packages(new_packages)

library(Matrix)
library(parallel)
library(doParallel)
library(foreach)
library(openxlsx)
library(numDeriv)

# ==============================================================================
# 1. User Settings and Global Variables
# ==============================================================================
data_file <- "../simulated_data_1111_n400_T5.csv"
w_file    <- "../spatial_weights_n400.csv"

region_var <- "region"
time_var   <- "time"
time_now   <- 2
time_lag   <- 1

# K denotes the number of response variables (K-variate GNR models)
K <- 2 
model_ids <- c("0011", "1011", "0111", "1111", "1001", "0101", "1101", "000d", "d0dd", "0ddd", "dddd")

# Initial coarse grid for the penalty parameter gamma
coarse_grid <- c(0, 10^seq(-2, 4, by = 0.5))

# TRUE: Gamma search is performed only for models with full spatial matrices.
use_gamma_only_for_full_spatial <- TRUE

# ==============================================================================
# 2. Data Loading and Preprocessing
# ==============================================================================
if (!file.exists(data_file)) stop(sprintf("Data file not found: %s", data_file))
if (!file.exists(w_file))    stop(sprintf("Spatial weights file not found: %s", w_file))

df <- read.csv(data_file, header = TRUE, stringsAsFactors = FALSE)
required_cols <- c(region_var, time_var, "y1", "y2", "x_common1", "x_common2", "x_specific1_1", "x_specific2_1")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns in data: %s", paste(missing_cols, collapse = ", ")))
}

regions <- sort(unique(df[[region_var]]))
n <- length(regions)
if (n <= 0) stop("Number of regions (n) is 0.")

tmp_W <- read.csv(w_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Clean W matrix by ignoring the ID column if present
if (ncol(tmp_W) == n + 1) {
  W_mat <- as.matrix(data.frame(lapply(tmp_W[, -1, drop = FALSE], as.numeric), check.names = FALSE))
} else if (ncol(tmp_W) > n) {
  W_mat <- as.matrix(data.frame(lapply(tmp_W[, (ncol(tmp_W) - n + 1):ncol(tmp_W), drop = FALSE], as.numeric), check.names = FALSE))
} else {
  W_mat <- as.matrix(data.frame(lapply(tmp_W, as.numeric), check.names = FALSE))
}

dimnames(W_mat) <- NULL
if (nrow(W_mat) != n || ncol(W_mat) != n) {
  stop(sprintf("Dimensions of W (%d x %d) do not match the number of regions n=%d.", nrow(W_mat), ncol(W_mat), n))
}

# Row-normalize the spatial weight matrix W (Section 1)
row_sums <- rowSums(W_mat)
if (max(abs(row_sums - 1)) > 1e-8) {
  W_mat <- sweep(W_mat, 1, row_sums, "/")
}
eigen_W <- eigen(W_mat, only.values = TRUE)$values

# Extract data for current time (t=2) and lagged time (t=1)
df_t2 <- df[df[[time_var]] == time_now, , drop = FALSE]
df_t1 <- df[df[[time_var]] == time_lag, , drop = FALSE]

df_t2 <- df_t2[order(df_t2[[region_var]]), , drop = FALSE]
df_t1 <- df_t1[order(df_t1[[region_var]]), , drop = FALSE]

# Construct response vector Y of size (Kn x 1) (Eqn 4)
Y <- c(df_t2$y1, df_t2$y2)

base_data_list <- list(
  K = K, n = n, W_mat = W_mat, Y = Y, eigen_W = eigen_W,
  y1 = df_t2$y1, y2 = df_t2$y2, regions = regions, df_t2 = df_t2, df_t1 = df_t1
)

# ==============================================================================
# 3. Model Specifications and Constraints (Table 1)
# ==============================================================================
get_R_code      <- function(model_id) substr(model_id, 1, 1)
get_Lambda_code <- function(model_id) substr(model_id, 2, 2)
get_A_code      <- function(model_id) substr(model_id, 3, 3)
get_Sigma_code  <- function(model_id) substr(model_id, 4, 4)

# Get the count of spatial parameters (R and Lambda)
get_num_spatial_params <- function(model_id, K = 2) {
  r_code <- get_R_code(model_id)
  l_code <- get_Lambda_code(model_id)
  n_R <- if (r_code == "1") K^2 else if (r_code == "d") K else if (r_code == "0") 0
  n_L <- if (l_code == "1") K^2 else if (l_code == "d") K else if (l_code == "0") 0
  n_R + n_L
}

# Get the count of Sigma (covariance) parameters
get_num_sigma_params <- function(model_id, K = 2) {
  s_code <- get_Sigma_code(model_id)
  if (s_code == "1") return(K * (K + 1) / 2) # Full symmetric covariance
  if (s_code == "d") return(K)               # Diagonal covariance
}

is_sigma_diagonal <- function(model_id) get_Sigma_code(model_id) == "d"

# Build R and Lambda matrices depending on model codes
build_matrices <- function(theta_spatial, model_id, K) {
  R <- matrix(0, K, K)
  Lambda <- matrix(0, K, K)
  idx <- 1

  r_code <- get_R_code(model_id)
  if (r_code == "1") {
    R <- matrix(theta_spatial[idx:(idx + K^2 - 1)], K, K, byrow = TRUE)
    idx <- idx + K^2
  } else if (r_code == "d") {
    R <- diag(theta_spatial[idx:(idx + K - 1)], K, K)
    idx <- idx + K
  }

  l_code <- get_Lambda_code(model_id)
  if (l_code == "1") {
    Lambda <- matrix(theta_spatial[idx:(idx + K^2 - 1)], K, K, byrow = TRUE)
  } else if (l_code == "d") {
    Lambda <- diag(theta_spatial[idx:(idx + K - 1)], K, K)
  }

  list(R = R, Lambda = Lambda)
}

# Construct the design matrix X_block according to the A (time-AR) constraint
make_X_block <- function(df_t2, df_t1, model_id) {
  A_code <- get_A_code(model_id)
  n_local <- nrow(df_t2)

  X1_base <- cbind(
    y1_intercept     = rep(1, n_local),
    y1_x_common1    = df_t2$x_common1,
    y1_x_common2    = df_t2$x_common2,
    y1_x_specific1_1 = df_t2$x_specific1_1
  )
  X2_base <- cbind(
    y2_intercept     = rep(1, n_local),
    y2_x_common1    = df_t2$x_common1,
    y2_x_common2    = df_t2$x_common2,
    y2_x_specific2_1 = df_t2$x_specific2_1
  )

  # Include time-AR terms if constrained by model_id
  if (A_code == "1") {
    X1 <- cbind(X1_base, y1_alpha_y1_lag = df_t1$y1, y1_alpha_y2_lag = df_t1$y2)
    X2 <- cbind(X2_base, y2_alpha_y1_lag = df_t1$y1, y2_alpha_y2_lag = df_t1$y2)
  } else if (A_code == "d") {
    X1 <- cbind(X1_base, y1_alpha_y1_lag = df_t1$y1)
    X2 <- cbind(X2_base, y2_alpha_y2_lag = df_t1$y2)
  } else if (A_code == "0") {
    X1 <- X1_base
    X2 <- X2_base
  }

  X_block <- matrix(0, nrow = 2 * n_local, ncol = ncol(X1) + ncol(X2))
  X_block[1:n_local, 1:ncol(X1)] <- X1
  X_block[(n_local + 1):(2 * n_local), (ncol(X1) + 1):(ncol(X1) + ncol(X2))] <- X2
  colnames(X_block) <- c(colnames(X1), colnames(X2))
  X_block
}

make_data_list_for_model <- function(base_dl, model_id) {
  dl <- base_dl
  dl$X_block <- make_X_block(base_dl$df_t2, base_dl$df_t1, model_id)
  dl
}

# ==============================================================================
# 4. Numerical Computation Functions for Estimation
# ==============================================================================

# Ensure Sigma remains a positive definite symmetric matrix
regularize_symmetric <- function(S, eps = 1e-8) {
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  min_ev <- min(ev)
  if (!is.finite(min_ev)) return(S + diag(eps, nrow(S)))
  if (min_ev <= eps) S <- S + diag(eps - min_ev, nrow(S))
  S
}

apply_sigma_constraint <- function(Sigma, model_id) {
  Sigma <- (Sigma + t(Sigma)) / 2
  if (is_sigma_diagonal(model_id)) {
    Sigma <- diag(pmax(diag(Sigma), 1e-10), nrow = nrow(Sigma), ncol = ncol(Sigma))
  }
  regularize_symmetric(Sigma)
}

safe_solve <- function(A, b = NULL, ridge = 1e-8) {
  if (is.null(b)) {
    tryCatch(solve(A), error = function(e) solve(A + diag(ridge, nrow(A))))
  } else {
    tryCatch(solve(A, b), error = function(e) solve(A + diag(ridge, nrow(A)), b))
  }
}

# Computes log |I_{Kn} - M %x% W| using eigenvalues (Eqn 18)
# Includes strict check for stationarity (real parts of elements must be strictly positive)
log_det_kron_spatial <- function(M, eigen_W, eps = 1e-12) {
  eig_M <- eigen(M, only.values = TRUE)$values
  factors <- as.vector(outer(eig_M, eigen_W, function(a, w) 1 - a * w))
  
  # Strict Stationarity Check: Reject regions where the determinant approaches 0 or goes negative.
  real_factors <- factors[Im(factors) == 0]
  if (any(Re(real_factors) <= eps)) return(-Inf)
  if (any(!is.finite(Mod(factors)))) return(-Inf)
  
  sum(log(Mod(factors)))
}

# Computes the profile log-likelihood l_c(theta_1) for a given set of spatial parameters (Eqn 35)
profile_loglik <- function(theta_spatial, model_id, dl) {
  mats <- build_matrices(theta_spatial, model_id, dl$K)
  I_Kn <- diag(dl$K * dl$n)
  A_R <- I_Kn - kronecker(mats$R, dl$W_mat)
  A_L <- I_Kn - kronecker(mats$Lambda, dl$W_mat)

  z_temp <- A_L %*% (A_R %*% dl$Y)
  X_temp <- A_L %*% dl$X_block

  # Step 1: Initial OLS Estimation
  beta_hat <- tryCatch(
    as.vector(qr.solve(X_temp, z_temp)),
    error = function(e) rep(0, ncol(dl$X_block))
  )
  names(beta_hat) <- colnames(dl$X_block)

  res_vec <- as.vector(z_temp - X_temp %*% beta_hat)
  res_mat <- matrix(res_vec, nrow = dl$n, ncol = dl$K)
  Sigma_hat <- apply_sigma_constraint(t(res_mat) %*% res_mat / dl$n, model_id)

  # Step 2: Iterative SUR/FGLS Updates for beta hat (Eqns 13, 14) and Sigma hat (Eqn 15)
  max_iter <- 50
  tol <- 1e-6
  for (iter in seq_len(max_iter)) {
    Sigma_inv <- safe_solve(Sigma_hat)
    Sigma_inv_I <- kronecker(Sigma_inv, diag(dl$n))

    XtWX <- t(X_temp) %*% Sigma_inv_I %*% X_temp
    XtWz <- t(X_temp) %*% Sigma_inv_I %*% z_temp

    beta_new <- tryCatch(
      as.vector(safe_solve(XtWX, XtWz)),
      error = function(e) beta_hat
    )
    names(beta_new) <- colnames(dl$X_block)

    res_vec_new <- as.vector(z_temp - X_temp %*% beta_new)
    res_mat_new <- matrix(res_vec_new, nrow = dl$n, ncol = dl$K)
    Sigma_new <- apply_sigma_constraint(t(res_mat_new) %*% res_mat_new / dl$n, model_id)

    change <- max(abs(beta_new - beta_hat), abs(Sigma_new - Sigma_hat))
    beta_hat <- beta_new
    Sigma_hat <- Sigma_new

    if (change < tol) break
  }

  # Step 3: Calculate Log-Likelihood l(R, Lambda, beta, Sigma) (Eqn 9)
  det_R <- log_det_kron_spatial(mats$R, dl$eigen_W)
  det_L <- log_det_kron_spatial(mats$Lambda, dl$eigen_W)
  if (!is.finite(det_R) || !is.finite(det_L)) {
    return(list(ll = -1e10, beta = beta_hat, Sigma = Sigma_hat, R = mats$R, Lambda = mats$Lambda, iterations = iter))
  }

  det_S <- determinant(Sigma_hat, logarithm = TRUE)
  log_det_S <- as.numeric(det_S$modulus[1])
  if (!is.finite(log_det_S) || det_S$sign <= 0) {
    return(list(ll = -1e10, beta = beta_hat, Sigma = Sigma_hat, R = mats$R, Lambda = mats$Lambda, iterations = iter))
  }

  # Calculate the exact Q value (Eqn 10): Q = z'(Sigma^{-1} %x% I_n)z
  Sigma_inv_exact <- safe_solve(Sigma_hat)
  Q <- 0
  for (i in 1:dl$K) {
    for (j in 1:dl$K) {
      Q <- Q + Sigma_inv_exact[i, j] * sum(res_mat_new[, i] * res_mat_new[, j])
    }
  }

  # Final Log-Likelihood l(R, Lambda, beta, Sigma) according to Eqn (9)
  ll <- - (dl$K * dl$n / 2) * log(2 * pi) + det_R + det_L - (dl$n / 2) * log_det_S - Q / 2

  list(
    ll = as.numeric(ll),
    beta = beta_hat,
    Sigma = Sigma_hat,
    R = mats$R,
    Lambda = mats$Lambda,
    iterations = iter
  )
}

# Penalized partial log-likelihood l_cp(theta_1) (Eqn 36)
neg_penalized_loglik <- function(theta_spatial, model_id, dl, gamma) {
  res <- profile_loglik(theta_spatial, model_id, dl)$ll
  if (!is.finite(res) || res <= -1e10) return(1e10)
  Re(-res) + (gamma / 2) * sum(theta_spatial^2)
}

# Calculates effective degrees of freedom (d_eff) used for penalized Information Criteria
compute_effective_df <- function(theta_spatial, model_id, dl, gamma) {
  n_spatial <- length(theta_spatial)
  
  # beta_count = p (exogenous vars) + K^2 (time-AR). As noted in Table 1, 
  # parameters excluding time-AR are denoted by 'p'. Therefore, beta_count correctly 
  # represents all unpenalized regression parameters.
  beta_count <- ncol(dl$X_block)
  sigma_count <- get_num_sigma_params(model_id, dl$K)

  if (n_spatial > 0) {
    if (gamma == 0) {
      d_spatial <- n_spatial
    } else {
      H <- tryCatch(
        numDeriv::hessian(function(t) -profile_loglik(t, model_id, dl)$ll, theta_spatial),
        error = function(e) NULL
      )
      if (!is.null(H) && all(is.finite(H))) {
        # Eqn (24): d_eff for regularized spatial parameters = trace[H %*% (H + gamma*I)^{-1}]
        d_spatial <- tryCatch(
          sum(diag(H %*% safe_solve(H + gamma * diag(n_spatial)))),
          error = function(e) n_spatial
        )
        if (!is.finite(d_spatial)) d_spatial <- n_spatial
      } else {
        d_spatial <- n_spatial
      }
    }
  } else {
    d_spatial <- 0
  }

  list(
    d_spatial = as.numeric(d_spatial),
    beta_count = beta_count,
    sigma_count = sigma_count,
    d_eff = as.numeric(d_spatial + beta_count + sigma_count), # For pAIC / pBIC (Eqn 25)
    d_eff_unpenalized = as.numeric(n_spatial + beta_count + sigma_count) # For standard AIC / BIC
  )
}

search_gamma <- function(model_id, grid, dl, init_theta) {
  best_pAIC <- Inf
  best_pBIC <- Inf
  best_gAIC <- NA
  best_gBIC <- NA
  best_mAIC <- NULL
  best_mBIC <- NULL
  best_dfAIC <- NULL
  best_dfBIC <- NULL
  cur_theta <- init_theta
  n_spatial <- length(init_theta)

  for (g in grid) {
    if (n_spatial > 0) {
      opt <- optim(
        par = cur_theta,
        fn = neg_penalized_loglik,
        model_id = model_id,
        dl = dl,
        gamma = g,
        method = "L-BFGS-B",
        lower = rep(-0.99, n_spatial),
        upper = rep(0.99, n_spatial),
        control = list(maxit = 500)
      )
      cur_theta <- opt$par
    }

    m <- profile_loglik(cur_theta, model_id, dl)
    df_info <- compute_effective_df(cur_theta, model_id, dl, g)

    # Penalized Information Criteria pAIC and pBIC (Eqn 25)
    pa <- -2 * m$ll + 2 * df_info$d_eff
    pb <- -2 * m$ll + log(dl$K * dl$n) * df_info$d_eff

    if (pa < best_pAIC) {
      best_pAIC <- pa; best_gAIC <- g; best_mAIC <- m; best_dfAIC <- df_info
    }
    if (pb < best_pBIC) {
      best_pBIC <- pb; best_gBIC <- g; best_mBIC <- m; best_dfBIC <- df_info
    }
  }

  list(
    pAIC = best_pAIC, gAIC = best_gAIC, mAIC = best_mAIC, dfAIC = best_dfAIC,
    pBIC = best_pBIC, gBIC = best_gBIC, mBIC = best_mBIC, dfBIC = best_dfBIC,
    last_t = cur_theta
  )
}

# ==============================================================================
# 5. Result Formatting
# ==============================================================================
fmt_coef <- function(b, p_vals, name, digits = 4) {
  if (is.null(names(b)) || !(name %in% names(b))) return(NA_character_)
  idx <- match(name, names(b))
  stars <- function(p) {
    ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))
  }
  paste0(round(b[idx], digits), stars(p_vals[idx]))
}

make_result_row <- function(m, g_opt, id, crit_name, crit_val, dl, df_info) {
  b <- m$beta
  S <- m$Sigma
  R <- m$R
  L <- m$Lambda

  X_tilde <- (diag(dl$K * dl$n) - kronecker(L, dl$W_mat)) %*% dl$X_block
  Sigma_inv <- safe_solve(S)
  V <- tryCatch(
    safe_solve(t(X_tilde) %*% kronecker(Sigma_inv, diag(dl$n)) %*% X_tilde),
    error = function(e) matrix(NA_real_, length(b), length(b))
  )

  if (any(is.na(V))) {
    p_vals <- rep(NA_real_, length(b))
  } else {
    se <- sqrt(pmax(diag(V), 0))
    z_val <- ifelse(se > 0, b / se, NA_real_)
    p_vals <- 2 * (1 - pnorm(abs(z_val)))
  }
  names(p_vals) <- names(b)

  # pseudo R2 (Eqn 26): Computed via expected values E[y_t] = (I_{Kn} - R %x% W)^{-1} X_t beta (Note 5)
  I_Kn <- diag(dl$K * dl$n)
  A_R <- I_Kn - kronecker(R, dl$W_mat)
  xb <- as.vector(dl$X_block %*% b)
  fitted_y <- tryCatch(as.vector(solve(A_R, xb)), error = function(e) xb)
  r2_1 <- suppressWarnings(cor(dl$y1, fitted_y[1:dl$n])^2)
  r2_2 <- suppressWarnings(cor(dl$y2, fitted_y[(dl$n + 1):(2 * dl$n)])^2)

  data.frame(
    Model_ID = paste0(id, "_", crit_name),
    Optimized_For = crit_name,
    Optimal_gamma = g_opt,
    Criterion_Value = crit_val, # Holds the specific value for pAIC or pBIC

    free_spatial_param_count = get_num_spatial_params(id, dl$K),
    beta_count = df_info$beta_count,
    sigma_count = df_info$sigma_count,
    d_spatial_eff = df_info$d_spatial,
    d_eff_used = df_info$d_eff,
    d_eff_unpenalized = df_info$d_eff_unpenalized,

    rho11 = R[1, 1], rho12 = R[1, 2], rho21 = R[2, 1], rho22 = R[2, 2],
    lambda11 = L[1, 1], lambda12 = L[1, 2], lambda21 = L[2, 1], lambda22 = L[2, 2],

    alpha11 = fmt_coef(b, p_vals, "y1_alpha_y1_lag"),
    alpha12 = fmt_coef(b, p_vals, "y1_alpha_y2_lag"),
    alpha21 = fmt_coef(b, p_vals, "y2_alpha_y1_lag"),
    alpha22 = fmt_coef(b, p_vals, "y2_alpha_y2_lag"),

    sigma11 = S[1, 1], sigma12 = S[1, 2], sigma21 = S[2, 1], sigma22 = S[2, 2],

    beta_intercept_y1 = fmt_coef(b, p_vals, "y1_intercept"),
    beta_intercept_y2 = fmt_coef(b, p_vals, "y2_intercept"),
    beta_common1_y1 = fmt_coef(b, p_vals, "y1_x_common1"),
    beta_common1_y2 = fmt_coef(b, p_vals, "y2_x_common1"),
    beta_common2_y1 = fmt_coef(b, p_vals, "y1_x_common2"),
    beta_common2_y2 = fmt_coef(b, p_vals, "y2_x_common2"),
    beta_specific1_1 = fmt_coef(b, p_vals, "y1_x_specific1_1"),
    beta_specific2_1 = fmt_coef(b, p_vals, "y2_x_specific2_1"),

    logLik = m$ll,
    beta_sigma_iterations = m$iterations,
    pseudo_R2_avg = mean(c(r2_1, r2_2), na.rm = TRUE),

    # Standard Unpenalized Criteria Output
    AIC = -2 * m$ll + 2 * df_info$d_eff_unpenalized,
    BIC = -2 * m$ll + log(dl$K * dl$n) * df_info$d_eff_unpenalized,
    
    # Penalized Information Criteria Output (Explicit representation)
    pAIC = -2 * m$ll + 2 * df_info$d_eff,
    pBIC = -2 * m$ll + log(dl$K * dl$n) * df_info$d_eff,

    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 6. Parallel Execution
# ==============================================================================
num_cores <- min(length(model_ids), max(1, detectCores() - 1))
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("Starting parallel estimation. Cores in use: %d\n", num_cores))

results_combined <- foreach(
  id = model_ids,
  .packages = c("Matrix", "numDeriv"),
  .combine = rbind,
  .export = c(
    "get_R_code", "get_Lambda_code", "get_A_code", "get_Sigma_code",
    "get_num_spatial_params", "get_num_sigma_params", "is_sigma_diagonal",
    "build_matrices", "make_X_block", "make_data_list_for_model",
    "regularize_symmetric", "apply_sigma_constraint", "safe_solve",
    "log_det_kron_spatial", "profile_loglik", "neg_penalized_loglik",
    "compute_effective_df", "search_gamma", "fmt_coef", "make_result_row",
    "coarse_grid", "use_gamma_only_for_full_spatial", "base_data_list"
  )
) %dopar% {

  dl_id <- make_data_list_for_model(base_data_list, id)

  n_spatial <- get_num_spatial_params(id, dl_id$K)
  init_theta <- if (n_spatial > 0) rep(0.01, n_spatial) else numeric(0)

  has_full_spatial <- (get_R_code(id) == "1" || get_Lambda_code(id) == "1")
  use_gamma <- if (use_gamma_only_for_full_spatial) has_full_spatial else (n_spatial > 0)

  if (!use_gamma) {
    # Non-spatial models and diagonal-only models are treated as standard AIC/BIC with gamma = 0
    f_res <- search_gamma(id, c(0), dl_id, init_theta)
  } else {
    # 1st Pass: Coarse Grid Search
    c_res <- search_gamma(id, coarse_grid, dl_id, init_theta)

    # 2nd Pass: Fine Grid Search around the best gamma
    best_gs <- unique(c(c_res$gAIC, c_res$gBIC))
    fine_grid <- c()

    for (g_star in best_gs) {
      if (is.finite(g_star) && g_star > 0 && g_star < max(coarse_grid)) {
        pos <- which(abs(coarse_grid - g_star) < 1e-10)
        if (length(pos) > 0) {
          pos <- pos[1]
          g_low <- if (pos > 2) coarse_grid[pos - 1] else g_star / 10
          g_high <- if (pos < length(coarse_grid)) coarse_grid[pos + 1] else max(coarse_grid)
          fine_grid <- c(fine_grid, 10^seq(log10(g_low), log10(g_high), by = 0.05))
        }
      }
    }
    fine_grid <- sort(unique(fine_grid))

    if (length(fine_grid) > 0) {
      f_res <- search_gamma(id, fine_grid, dl_id, c_res$last_t)
      # Retain the better result between 1st and 2nd pass
      if (c_res$pAIC < f_res$pAIC) {
        f_res$pAIC <- c_res$pAIC; f_res$gAIC <- c_res$gAIC
        f_res$mAIC <- c_res$mAIC; f_res$dfAIC <- c_res$dfAIC
      }
      if (c_res$pBIC < f_res$pBIC) {
        f_res$pBIC <- c_res$pBIC; f_res$gBIC <- c_res$gBIC
        f_res$mBIC <- c_res$mBIC; f_res$dfBIC <- c_res$dfBIC
      }
    } else {
      f_res <- c_res
    }
  }

  row_aic <- make_result_row(f_res$mAIC, f_res$gAIC, id, "pAIC", f_res$pAIC, dl_id, f_res$dfAIC)
  row_bic <- make_result_row(f_res$mBIC, f_res$gBIC, id, "pBIC", f_res$pBIC, dl_id, f_res$dfBIC)

  row_aic$Group <- "pAIC"
  row_bic$Group <- "pBIC"

  rbind(row_aic, row_bic)
}

stopCluster(cl)

# ==============================================================================
# 7. Export Results separated by pAIC / pBIC
# ==============================================================================
results_pAIC <- results_combined[results_combined$Group == "pAIC", ]
results_pBIC <- results_combined[results_combined$Group == "pBIC", ]

transpose_and_format <- function(df) {
  df$Group <- NULL
  rownames(df) <- df$Model_ID
  t_df <- as.data.frame(t(df[, -1, drop = FALSE]))
  t_df <- cbind(Parameter = rownames(t_df), t_df)
  rownames(t_df) <- NULL
  t_df
}

transposed_pAIC <- transpose_and_format(results_pAIC)
transposed_pBIC <- transpose_and_format(results_pBIC)

write.xlsx(transposed_pAIC, "MSTR_Final_Transposed_pAIC.xlsx", rowNames = FALSE)
write.xlsx(transposed_pBIC, "MSTR_Final_Transposed_pBIC.xlsx", rowNames = FALSE)
write.csv(results_combined, "MSTR_Final_Long_Format.csv", row.names = FALSE)

cat("Done. Please check the following files:\n")
cat("  - MSTR_Final_Transposed_pAIC.xlsx\n")
cat("  - MSTR_Final_Transposed_pBIC.xlsx\n")
cat("  - MSTR_Final_Long_Format.csv\n")

# ==============================================================================
# 8. Execution Time Measurement
# ==============================================================================
# Record the system time when the script finishes
script_end_time <- Sys.time()

# Calculate the total elapsed time in seconds
elapsed_time <- difftime(script_end_time, script_start_time, units = "secs")
elapsed_secs <- as.numeric(elapsed_time)
elapsed_mins <- elapsed_secs / 60

# Create a clean summary message array
time_summary <- c(
  "========================================",
  sprintf("Execution started at:  %s", format(script_start_time, "%Y-%m-%d %H:%M:%S")),
  sprintf("Execution finished at: %s", format(script_end_time, "%Y-%m-%d %H:%M:%S")),
  sprintf("Total elapsed time:    %.2f seconds (%.2f minutes)", elapsed_secs, elapsed_mins),
  "========================================"
)

# Print the timing summary to the console
cat("\n")
cat(paste(time_summary, collapse = "\n"), "\n")

# Save the execution time summary to a text file in the working directory
writeLines(time_summary, "MSTR_Final_Execution_Time.txt")
cat("\nTiming results have been saved to 'MSTR_Final_Execution_Time.txt'.\n")