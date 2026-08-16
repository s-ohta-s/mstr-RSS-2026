# penalized_spatial.r
# Penalized-likelihood estimation of spatial parameters
#   When the sample size is small (e.g., n=46 prefectures), the spatial parameters (ρ, λ)
#   tend to stick to the boundary (±0.99); an L2 penalty is introduced to mitigate this.
# Penalized likelihood:
#   ℓ_pen(θ) = ℓ(θ) - γ/2 * ||θ_spatial||²
#   MSAR:  ℓ_pen = ℓ - γ/2 * ||R||²_F
#   MSEM:  ℓ_pen = ℓ - γ/2 * ||Λ||²_F
#   MGNS: ℓ_pen = ℓ - γ/2 * (||R||²_F + ||Λ||²_F)
#   where ||·||_F is the Frobenius norm
#   - The penalty is used only during optimization; the final AIC/BIC use the non-penalized likelihood
#   - γ=0 means no penalty (identical to ordinary estimation)
#   - The larger γ is, the more the spatial parameters shrink toward 0 (similar to ridge regression)
#   - γ can be selected via cross-validation or BIC-like criteria

# 1. Penalty function

#' Compute the L2 penalty (squared Frobenius norm)
#' 
#' @param ... one or more matrices
#' @param gamma penalty strength
#' @return γ/2 * Σ||M_i||²_F
compute_penalty <- function(..., gamma) {
  matrices <- list(...)
  total <- 0
  for (M in matrices) {
    total <- total + sum(M^2)
  }
  return(gamma / 2 * total)
}

# 1b. Hessian and effective-degrees-of-freedom computation for GIC (generalized information criterion)
#   GIC(γ) = -2ℓ(θ̂_γ) + 2·df_eff(γ)
#   df_eff(γ) = tr[H(H + 2γI)^{-1}] + k₂
#   H: negative Hessian of the non-penalized profile log-likelihood (observed information matrix)
#   k₂: number of non-spatial parameters (elements of β, Σ)

#' Numerical Hessian of the non-penalized profile likelihood (MSAR model)
#'
#' @param R_hat estimated R matrix
#' @param beta_hat estimated β vector
#' @param Sigma_hat estimated Σ matrix
#' @param data_list return value of prepare_data_extended
#' @param max_iter_inner maximum number of inner iterations
#' @param tol_inner convergence threshold of the inner iteration
#' @return k₁×k₁ Hessian matrix (second derivative of the negative log-likelihood)
compute_profile_hessian_msar <- function(
  R_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_R <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  neg_profile_loglik <- function(param) {
    R <- param_to_R(param)
    
# Do not use check_stationarity (a smooth objective improves Hessian accuracy)
    
    inner <- tryCatch({
      iterate_beta_sigma(R = R, beta_init = beta_ref, Sigma_init = Sigma_ref,
                         y = y, X = X, W = W, k = k, n = n,
                         max_iter = max_iter_inner, tol = tol_inner,
                         ridge_eps = 1e-6, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(inner)) {
      rho_R <- spectral_radius(R, eigen_W)
      return(1e6 * max(rho_R, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood(R = R, beta_hat = inner$beta, Sigma_hat = inner$Sigma,
                                  y = y, X = X, W = W, eigen_W = eigen_W,
                                  k = k, n = n, verbose = 0, smooth = TRUE)
    }, error = function(e) -Inf)
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(R, eigen_W)
      return(1e6 * max(rho_R, 1))
    }
    
    return(-prof_lik)
  }
  
  theta1_hat <- as.vector(t(R_hat))
  H <- numDeriv::hessian(neg_profile_loglik, theta1_hat)
  return(H)
}

#' Numerical Hessian of the non-penalized profile likelihood (MSEM model)
compute_profile_hessian_msem <- function(
  Lambda_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_T <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  neg_profile_loglik <- function(param) {
    Lambda_mat <- param_to_T(param)
    
# Do not use check_stationarity (a smooth objective improves Hessian accuracy)
    
    inner <- tryCatch({
      iterate_beta_sigma_msem(Lambda_mat = Lambda_mat, beta_init = beta_ref, Sigma_init = Sigma_ref,
                              y = y, X = X, W = W, k = k, n = n,
                              max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(inner)) {
      rho_Lambda <- spectral_radius(Lambda_mat, eigen_W)
      return(1e6 * max(rho_Lambda, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_msem(Lambda_mat = Lambda_mat, beta_hat = inner$beta, Sigma_hat = inner$Sigma,
                                      y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
                                      smooth = TRUE)
    }, error = function(e) -Inf)
    if (!is.finite(prof_lik)) {
      rho_Lambda <- spectral_radius(Lambda_mat, eigen_W)
      return(1e6 * max(rho_Lambda, 1))
    }
    
    return(-prof_lik)
  }
  
  theta1_hat <- as.vector(t(Lambda_hat))
  H <- numDeriv::hessian(neg_profile_loglik, theta1_hat)
  return(H)
}

#' Numerical Hessian of the non-penalized profile likelihood (MGNS model)
compute_profile_hessian_mgns <- function(
  R_hat, Lambda_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_RT <- function(param) {
    R <- matrix(param[1:(k*k)], nrow = k, ncol = k, byrow = TRUE)
    Lambda_mat <- matrix(param[(k*k + 1):(2*k*k)], nrow = k, ncol = k, byrow = TRUE)
    list(R = R, Lambda = Lambda_mat)
  }
  
  neg_profile_loglik <- function(param) {
    RT <- param_to_RT(param)
    
# A smooth objective improves Hessian accuracy
    
    inner <- tryCatch({
      iterate_beta_sigma_mgns(R = RT$R, Lambda_mat = RT$Lambda,
                               beta_init = beta_ref, Sigma_init = Sigma_ref,
                               y = y, X = X, W = W, k = k, n = n,
                               max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(inner)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_Lambda <- spectral_radius(RT$Lambda, eigen_W)
      return(1e6 * max(rho_R, rho_Lambda, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_mgns(R = RT$R, Lambda_mat = RT$Lambda,
                                       beta_hat = inner$beta, Sigma_hat = inner$Sigma,
                                       y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
                                       smooth = TRUE)
    }, error = function(e) -Inf)
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_Lambda <- spectral_radius(RT$Lambda, eigen_W)
      return(1e6 * max(rho_R, rho_Lambda, 1))
    }
    
    return(-prof_lik)
  }
  
  theta1_hat <- c(as.vector(t(R_hat)), as.vector(t(Lambda_hat)))
  H <- numDeriv::hessian(neg_profile_loglik, theta1_hat)
  return(H)
}

#' Compute the GIC (generalized information criterion)
#'
#' @param H k₁×k₁ negative Hessian of the non-penalized profile likelihood (observed information matrix)
#' @param gamma penalty strength
#' @param loglik non-penalized log-likelihood
#' @param k2 number of non-spatial parameters
#' @param n_obs number of observations (k*n)
#' @return list containing df_eff, df_spatial, GIC_AIC, GIC_BIC
compute_gic <- function(H, gamma, loglik, k2, n_obs) {
  k1 <- nrow(H)
  
  if (gamma == 0) {
    # γ=0: df_eff = k1 + k2 (matches ordinary AIC/BIC)
    df_spatial <- k1
  } else {
    # df_spatial = tr[H(H + γD)^{-1}]
    H_reg <- H + gamma * diag(k1)
    
    H_reg_cond <- tryCatch(kappa(H_reg), error = function(e) Inf)
    if (H_reg_cond > 1e12) {
      warning(sprintf("GIC: large condition number of H + γI (κ=%.2e, γ=%.4g)", H_reg_cond, gamma))
    }
    
    S <- H %*% solve(H_reg)
    df_spatial <- sum(diag(S))
  }
  
  df_eff <- df_spatial + k2
  
  GIC_AIC <- -2 * loglik + 2 * df_eff
  GIC_BIC <- -2 * loglik + log(n_obs) * df_eff
  
  return(list(
    df_eff     = df_eff,
    df_spatial = df_spatial,
    GIC_AIC    = GIC_AIC,
    GIC_BIC    = GIC_BIC,
    H          = H
  ))
}

# 2. MSAR penalized BFGS

#' Optimize the R matrix via penalized BFGS (MSAR model)
#'
#' @param gamma penalty strength (0 = no penalty)
optimize_R_bfgs_penalized <- function(
  R_init,
  beta_current,
  Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  R_to_param <- function(R) as.vector(t(R))
  param_to_R <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  # Penalized objective function (smooth version)
  objective <- function(param) {
    R <- param_to_R(param)
    
    # Do not use check_stationarity (rely on the natural barrier of the log-likelihood)
    
    inner_result <- tryCatch({
      iterate_beta_sigma(
        R = R, beta_init = beta_ref, Sigma_init = Sigma_ref,
        y = y, X = X, W = W, k = k, n = n,
        max_iter = max_iter_inner, tol = tol_inner,
        ridge_eps = 1e-6, verbose = 0)
    }, error = function(e) NULL)
    
    if (is.null(inner_result)) {
      rho_R <- spectral_radius(R, eigen_W)
      return(1e6 * max(rho_R, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood(
        R = R, beta_hat = inner_result$beta, Sigma_hat = inner_result$Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n, verbose = 0,
        smooth = TRUE)
    }, error = function(e) -Inf)
    
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(R, eigen_W)
      return(1e6 * max(rho_R, 1))
    }
    
    penalty <- compute_penalty(R, gamma = gamma)
    
    return(-prof_lik + penalty)  # minimization, so negative likelihood + penalty
  }
  
  param_init <- R_to_param(R_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalized BFGS optimization start (MSAR) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Number of parameters: %d, γ = %.4f\n", k*k, gamma))
    init_obj <- objective(param_init)
    cat(sprintf("  Initial penalized objective: %.4f\n", init_obj))
  }
  
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS did not converge -> retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  R_final <- param_to_R(result$par)
  
  final_inner <- iterate_beta_sigma(
    R = R_final, beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner,
    ridge_eps = 1e-6, verbose = 0)
  
  # Non-penalized likelihood (for AIC/BIC computation)
  final_loglik <- compute_log_likelihood(
    R = R_final, beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  # Penalized likelihood (for reference)
  penalized_loglik <- final_loglik - compute_penalty(R_final, gamma = gamma)
  
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s, iterations: %d\n", ifelse(converged, "yes", "no"), result$counts["function"]))
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalized log-likelihood: %.4f\n", penalized_loglik))
    cat(sprintf("  Penalty amount: %.4f\n", compute_penalty(R_final, gamma = gamma)))
    cat("  Final R:\n"); print(round(R_final, 4))
    cat("=== Penalized BFGS done (MSAR) ===\n\n")
  }
  
  return(list(
    R = R_final,
    beta = final_inner$beta,
    Sigma = final_inner$Sigma,
    loglik = final_loglik,           # non-penalized (for AIC/BIC)
    penalized_loglik = penalized_loglik,  # penalized (for reference)
    gamma = gamma,
    penalty = compute_penalty(R_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}

# 3. MSEM penalized BFGS

optimize_Lambda_bfgs_penalized <- function(
  T_init,
  beta_current, Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  T_to_param <- function(Lambda_mat) as.vector(t(Lambda_mat))
  param_to_T <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  objective <- function(param) {
    Lambda_mat <- param_to_T(param)
    
    # Do not use check_stationarity (rely on the natural barrier of the log-likelihood)
    
    inner_result <- tryCatch({
      iterate_beta_sigma_msem(
        Lambda_mat = Lambda_mat, beta_init = beta_ref, Sigma_init = Sigma_ref,
        y = y, X = X, W = W, k = k, n = n,
        max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    
    if (is.null(inner_result)) {
      rho_Lambda <- spectral_radius(Lambda_mat, eigen_W)
      return(1e6 * max(rho_Lambda, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_msem(
        Lambda_mat = Lambda_mat, beta_hat = inner_result$beta, Sigma_hat = inner_result$Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
        smooth = TRUE)
    }, error = function(e) -Inf)
    
    if (!is.finite(prof_lik)) {
      rho_Lambda <- spectral_radius(Lambda_mat, eigen_W)
      return(1e6 * max(rho_Lambda, 1))
    }
    
    penalty <- compute_penalty(Lambda_mat, gamma = gamma)
    
    return(-prof_lik + penalty)
  }
  
  param_init <- T_to_param(T_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalized BFGS optimization start (MSEM) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Number of parameters: %d, γ = %.4f\n", k*k, gamma))
  }
  
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS did not converge -> retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  Lambda_final <- param_to_T(result$par)
  
  final_inner <- iterate_beta_sigma_msem(
    Lambda_mat = Lambda_final, beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
  
  # Non-penalized likelihood
  final_loglik <- compute_log_likelihood_msem(
    Lambda_mat = Lambda_final, beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  penalized_loglik <- final_loglik - compute_penalty(Lambda_final, gamma = gamma)
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s\n", ifelse(converged, "yes", "no")))
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalty amount: %.4f\n", compute_penalty(Lambda_final, gamma = gamma)))
    cat("  Final Λ:\n"); print(round(Lambda_final, 4))
    cat("=== Penalized BFGS done (MSEM) ===\n\n")
  }
  
  return(list(
    Lambda = Lambda_final,
    beta = final_inner$beta,
    Sigma = final_inner$Sigma,
    loglik = final_loglik,
    penalized_loglik = penalized_loglik,
    gamma = gamma,
    penalty = compute_penalty(Lambda_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}

# 4. MGNS penalized BFGS

optimize_RLambda_bfgs_penalized <- function(
  R_init, T_init,
  beta_current, Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W
  eigen_W <- data_list$eigen_W
  
  param_to_RT <- function(param) {
    R <- matrix(param[1:(k*k)], nrow = k, ncol = k, byrow = TRUE)
    Lambda_mat <- matrix(param[(k*k + 1):(2*k*k)], nrow = k, ncol = k, byrow = TRUE)
    list(R = R, Lambda = Lambda_mat)
  }
  RT_to_param <- function(R, Lambda_mat) c(as.vector(t(R)), as.vector(t(Lambda_mat)))
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  objective <- function(param) {
    RT <- param_to_RT(param)
    
    # Do not use check_stationarity (rely on the natural barrier of the log-likelihood)
    
    inner_result <- tryCatch({
      iterate_beta_sigma_mgns(
        R = RT$R, Lambda_mat = RT$Lambda,
        beta_init = beta_ref, Sigma_init = Sigma_ref,
        y = y, X = X, W = W, k = k, n = n,
        max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    
    if (is.null(inner_result)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_Lambda <- spectral_radius(RT$Lambda, eigen_W)
      return(1e6 * max(rho_R, rho_Lambda, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_mgns(
        R = RT$R, Lambda_mat = RT$Lambda,
        beta_hat = inner_result$beta, Sigma_hat = inner_result$Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
        smooth = TRUE)
    }, error = function(e) -Inf)
    
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_Lambda <- spectral_radius(RT$Lambda, eigen_W)
      return(1e6 * max(rho_R, rho_Lambda, 1))
    }
    
    penalty <- compute_penalty(RT$R, RT$Lambda, gamma = gamma)
    
    return(-prof_lik + penalty)
  }
  
  param_init <- RT_to_param(R_init, T_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalized BFGS optimization start (MGNS) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Number of parameters: %d (R: %d, Λ: %d), γ = %.4f\n", 2*k*k, k*k, k*k, gamma))
  }
  
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS did not converge -> retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  RT_opt <- param_to_RT(result$par)
  R_final <- RT_opt$R; Lambda_final <- RT_opt$Lambda
  
  final_inner <- iterate_beta_sigma_mgns(
    R = R_final, Lambda_mat = Lambda_final,
    beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
  
  final_loglik <- compute_log_likelihood_mgns(
    R = R_final, Lambda_mat = Lambda_final,
    beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  penalized_loglik <- final_loglik - compute_penalty(R_final, Lambda_final, gamma = gamma)
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s\n", ifelse(converged, "yes", "no")))
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalty amount: %.4f\n", compute_penalty(R_final, Lambda_final, gamma = gamma)))
    cat("  Final R:\n"); print(round(R_final, 4))
    cat("  Final Λ:\n"); print(round(Lambda_final, 4))
    cat("=== Penalized BFGS done (MGNS) ===\n\n")
  }
  
  return(list(
    R = R_final, Lambda = Lambda_final,
    beta = final_inner$beta, Sigma = final_inner$Sigma,
    loglik = final_loglik,
    penalized_loglik = penalized_loglik,
    gamma = gamma,
    penalty = compute_penalty(R_final, Lambda_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}

# 5. Integrated fit functions (penalty-aware version)

#' Penalized multivariate MSAR estimation
#'
#' gamma > 0: penalized estimation. gamma = 0: ordinary estimation.
#' AIC/BIC are computed from the non-penalized likelihood.
fit_msar_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  R_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = FALSE,
  data_list = NULL
) {

  start_time <- Sys.time()

  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalized multivariate MSAR estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }

  # S0: data preparation (skipped when a prebuilt data_list is supplied,
  # e.g. from attach_data() in the Monte Carlo driver)
  if (is.null(data_list)) {
    data_list <- prepare_data_extended(
      data_file = data_file, weight_file = weight_file,
      y_vars = y_vars, x_vars = x_vars,
      time_var = time_var, time_point = time_point,
      region_var = region_var, include_intercept = include_intercept,
      include_time_lag = include_time_lag, verbose = FALSE)
  }
  
  k <- data_list$k; n <- data_list$n
  actual_include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(actual_include_time_lag)) actual_include_time_lag <- TRUE
  
  if (verbose) cat(sprintf("  Data: k=%d, n=%d\n", k, n))
  
  # S1: initial estimation
  if (is.null(R_init)) {
    if (verbose) cat("S1: initial estimation via per-variable MSAR...\n")
    init_result <- initial_estimation_sly_extended(data_list, verbose = ifelse(verbose, 1, 0))
    R <- init_result$R_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    R <- R_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  # S2-S4: penalized optimization
  opt_result <- optimize_R_bfgs_penalized(
    R_init = R, beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  R_final <- opt_result$R
  beta_final <- opt_result$beta
  Sigma_final <- opt_result$Sigma
  loglik <- opt_result$loglik  # non-penalized
  
  # Information criteria (based on the non-penalized likelihood)
  num_params <- ncol(data_list$X) + k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  γ = %.4f, penalty amount = %.4f\n", gamma, opt_result$penalty))
    cat(sprintf("  R diagonal: [%s]\n", paste(sprintf("%.4f", diag(R_final)), collapse=", ")))
    cat(sprintf("  Execution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "MSAR",
    R                = R_final,
    beta             = beta_final,
    Sigma            = Sigma_final,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_bfgs",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}

#' Penalized multivariate MSEM estimation
fit_msem_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  T_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = FALSE,
  data_list = NULL
) {

  start_time <- Sys.time()

  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalized multivariate MSEM estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }

  if (is.null(data_list)) {
    data_list <- prepare_data_extended(
      data_file = data_file, weight_file = weight_file,
      y_vars = y_vars, x_vars = x_vars,
      time_var = time_var, time_point = time_point,
      region_var = region_var, include_intercept = include_intercept,
      include_time_lag = include_time_lag, verbose = FALSE)
  }
  
  k <- data_list$k; n <- data_list$n
  if (verbose) cat(sprintf("  Data: k=%d, n=%d\n", k, n))
  
  if (is.null(T_init)) {
    if (verbose) cat("S1: initial estimation via per-variable MSEM...\n")
    init_result <- initial_estimation_sem_extended(data_list, verbose = ifelse(verbose, 1, 0))
    Lambda_mat <- init_result$T_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    Lambda_mat <- T_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  opt_result <- optimize_Lambda_bfgs_penalized(
    T_init = Lambda_mat, beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  Lambda_final <- opt_result$Lambda
  loglik <- opt_result$loglik
  
  num_params <- ncol(data_list$X) + k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  γ = %.4f, penalty amount = %.4f\n", gamma, opt_result$penalty))
    cat(sprintf("  Λ diagonal: [%s]\n", paste(sprintf("%.4f", diag(Lambda_final)), collapse=", ")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "MSEM",
    Lambda_mat            = Lambda_final,
    beta             = opt_result$beta,
    Sigma            = opt_result$Sigma,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_bfgs",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}

#' Penalized multivariate MGNS estimation
fit_mgns_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  R_init = NULL, T_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = FALSE,
  data_list = NULL
) {

  start_time <- Sys.time()

  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalized multivariate MGNS estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }

  if (is.null(data_list)) {
    data_list <- prepare_data_extended(
      data_file = data_file, weight_file = weight_file,
      y_vars = y_vars, x_vars = x_vars,
      time_var = time_var, time_point = time_point,
      region_var = region_var, include_intercept = include_intercept,
      include_time_lag = include_time_lag, verbose = FALSE)
  }
  
  k <- data_list$k; n <- data_list$n
  if (verbose) cat(sprintf("  Data: k=%d, n=%d\n", k, n))
  
  if (is.null(R_init) || is.null(T_init)) {
    if (verbose) cat("S1: initial estimation via per-variable MGNS...\n")
    init_result <- initial_estimation_sdem_extended(data_list, verbose = ifelse(verbose, 1, 0))
    R <- if(is.null(R_init)) init_result$R_init else R_init
    Lambda_mat <- if(is.null(T_init)) init_result$T_init else T_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    R <- R_init; Lambda_mat <- T_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  opt_result <- optimize_RLambda_bfgs_penalized(
    R_init = R, T_init = Lambda_mat,
    beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  R_final <- opt_result$R; Lambda_final <- opt_result$Lambda
  loglik <- opt_result$loglik
  
  num_params <- ncol(data_list$X) + 2*k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Non-penalized log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  γ = %.4f, penalty amount = %.4f\n", gamma, opt_result$penalty))
    cat("  R:\n"); print(round(R_final, 4))
    cat("  Λ:\n"); print(round(Lambda_final, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "MGNS",
    R                = R_final,
    Lambda_mat            = Lambda_final,
    beta             = opt_result$beta,
    Sigma            = opt_result$Sigma,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_bfgs",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}

# 6. γ-search utilities

#' Run penalized estimation for multiple γ values and select the optimal γ via GIC
#'
#' @param fit_func one of fit_msar_penalized, fit_msem_penalized, fit_mgns_penalized
#' @param gammas vector of γ values to try (e.g., c(0, 10^seq(-6, 2, by=0.5)))
#' @param compute_gic_flag whether to compute GIC (numDeriv required if TRUE)
#' @param ... arguments passed to fit_func
#' @return a data frame of results and a list of all results
compare_gamma <- function(
  fit_func,
  gammas = c(0, 0.5, 1, 2, 5, 10, 20, 50),
  compute_gic_flag = TRUE,
  ...
) {
  
  cat("\n", paste(rep("#", 70), collapse=""), "\n")
  cat("### γ search: comparison of penalty strengths (GIC-aware version) ###\n")
  cat(paste(rep("#", 70), collapse=""), "\n\n")
  
  if (compute_gic_flag && !requireNamespace("numDeriv", quietly = TRUE)) {
    warning("The numDeriv package is not installed. Skipping GIC computation.")
    compute_gic_flag <- FALSE
  }
  
  results <- list()
  summary_df <- data.frame(
    gamma = numeric(),
    loglik = numeric(),
    AIC = numeric(),
    BIC = numeric(),
    df_eff = numeric(),
    GIC_AIC = numeric(),
    GIC_BIC = numeric(),
    penalty = numeric(),
    spatial_params = character(),
    stringsAsFactors = FALSE
  )
  
  # Phase 1: run estimation for all γ
  cat("--- Phase 1: estimation ---\n")
  for (i in seq_along(gammas)) {
    g <- gammas[i]
    cat(sprintf("  γ = %.4g (%d/%d) ... ", g, i, length(gammas)))
    
    res <- tryCatch({
      fit_func(gamma = g, verbose = FALSE, ...)
    }, error = function(e) {
      cat(sprintf("Error: %s\n", e$message))
      NULL
    })
    
    if (!is.null(res)) {
      results[[as.character(g)]] <- res
      cat(sprintf("loglik = %.4f\n", res$fit$loglik))
    }
  }
  
  # Phase 2: Hessian computation and GIC calculation
  if (compute_gic_flag && length(results) > 0) {
    cat("\n--- Phase 2: Hessian computation / GIC calculation ---\n")
    
    first_res <- results[[1]]
    model_type <- first_res$model_type
    k <- first_res$data_info$k
    n_obs <- first_res$fit$num_obs
    
    # k1: number of spatial parameters, k2: number of non-spatial parameters
    if (model_type == "MSAR") {
      k1 <- k^2
    } else if (model_type == "MSEM") {
      k1 <- k^2
    } else {
      k1 <- 2 * k^2
    }
    k2 <- first_res$fit$num_params - k1
    
    cat(sprintf("  Model: %s, k₁=%d (spatial), k₂=%d (non-spatial), n_obs=%d\n",
                model_type, k1, k2, n_obs))
  }
  
  for (g_str in names(results)) {
    g <- as.numeric(g_str)
    res <- results[[g_str]]
    
    format_matrix_params <- function(M, name) {
      k <- nrow(M)
      diag_str <- paste(sprintf("%.3f", diag(M)), collapse=",")
      if (k >= 2) {
        offdiag <- M[row(M) != col(M)]
        offdiag_str <- paste(sprintf("%.3f", offdiag), collapse=",")
        sprintf("%s_diag=[%s] %s_offdiag=[%s]", name, diag_str, name, offdiag_str)
      } else {
        sprintf("%s=[%s]", name, diag_str)
      }
    }
    
    if (res$model_type == "MSAR") {
      sp <- format_matrix_params(res$coefficients$R, "R")
    } else if (res$model_type == "MSEM") {
      sp <- format_matrix_params(res$coefficients$Lambda, "Λ")
    } else {
      sp <- paste(format_matrix_params(res$coefficients$R, "R"),
                  format_matrix_params(res$coefficients$Lambda, "Λ"))
    }
    
    pen_val <- if (!is.null(res$penalty)) res$penalty$value else 0
    
    df_eff_val <- res$fit$num_params  # default (num_params when GIC is not computed)
    gic_aic_val <- res$fit$AIC
    gic_bic_val <- res$fit$BIC
    
    if (compute_gic_flag) {
      cat(sprintf("  γ = %.4g: computing Hessian...", g))
      
      hessian_result <- tryCatch({
        if (res$model_type == "MSAR") {
          compute_profile_hessian_msar(
            R_hat = res$coefficients$R,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        } else if (res$model_type == "MSEM") {
          compute_profile_hessian_msem(
            Lambda_hat = res$coefficients$Lambda,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        } else {
          compute_profile_hessian_mgns(
            R_hat = res$coefficients$R,
            Lambda_hat = res$coefficients$Lambda,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        }
      }, error = function(e) {
        cat(sprintf(" Error: %s", e$message))
        NULL
      })
      
      if (!is.null(hessian_result)) {
        gic_result <- tryCatch({
          compute_gic(H = hessian_result, gamma = g,
                      loglik = res$fit$loglik, k2 = k2, n_obs = n_obs)
        }, error = function(e) {
          cat(sprintf(" GICError: %s", e$message))
          NULL
        })
        
        if (!is.null(gic_result)) {
          df_eff_val <- gic_result$df_eff
          gic_aic_val <- gic_result$GIC_AIC
          gic_bic_val <- gic_result$GIC_BIC
          cat(sprintf(" df_eff=%.2f, GIC_AIC=%.4f\n", df_eff_val, gic_aic_val))
        } else {
          cat(" GIC computation failed\n")
        }
      } else {
        cat(" Hessian computation failed\n")
      }
    }
    
    summary_df <- rbind(summary_df, data.frame(
      gamma = g,
      loglik = res$fit$loglik,
      AIC = res$fit$AIC,
      BIC = res$fit$BIC,
      df_eff = df_eff_val,
      GIC_AIC = gic_aic_val,
      GIC_BIC = gic_bic_val,
      penalty = pen_val,
      spatial_params = sp,
      stringsAsFactors = FALSE
    ))
  }
  
  summary_df <- summary_df[order(summary_df$gamma), ]
  rownames(summary_df) <- NULL
  
  cat("\n", paste(rep("=", 90), collapse=""), "\n")
  cat("### γ-search result summary ###\n")
  cat(paste(rep("=", 90), collapse=""), "\n\n")
  
  if (compute_gic_flag) {
    summary_numeric <- summary_df[, c("gamma", "loglik", "AIC", "BIC",
                                       "df_eff", "GIC_AIC", "GIC_BIC", "penalty")]
  } else {
    summary_numeric <- summary_df[, c("gamma", "loglik", "AIC", "BIC", "penalty")]
  }
  print(summary_numeric, row.names = FALSE, right = FALSE)
  
  cat("\n--- Evolution of the spatial parameters ---\n")
  for (i in seq_len(nrow(summary_df))) {
    g <- summary_df$gamma[i]
    res <- results[[as.character(g)]]
    if (is.null(res)) next
    cat(sprintf("\nγ = %.4g:\n", g))
    if (res$model_type == "MSAR" || res$model_type == "MGNS") {
      if (nrow(res$coefficients$R) == 2) {
        cat(sprintf("  R = [%7.4f %7.4f]\n", res$coefficients$R[1,1], res$coefficients$R[1,2]))
        cat(sprintf("      [%7.4f %7.4f]\n", res$coefficients$R[2,1], res$coefficients$R[2,2]))
      }
    }
    if (res$model_type == "MSEM" || res$model_type == "MGNS") {
      if (nrow(res$coefficients$Lambda) == 2) {
        cat(sprintf("  Λ = [%7.4f %7.4f]\n", res$coefficients$Lambda[1,1], res$coefficients$Lambda[1,2]))
        cat(sprintf("      [%7.4f %7.4f]\n", res$coefficients$Lambda[2,1], res$coefficients$Lambda[2,2]))
      }
    }
  }
  
  if (nrow(summary_df) > 0) {
    cat("\n--- Selection of the best γ ---\n")
    
    best_aic_idx <- which.min(summary_df$AIC)
    best_bic_idx <- which.min(summary_df$BIC)
    cat(sprintf("  Best AIC:     γ = %.4g (AIC = %.4f)\n",
                summary_df$gamma[best_aic_idx], summary_df$AIC[best_aic_idx]))
    cat(sprintf("  Best BIC:     γ = %.4g (BIC = %.4f)\n",
                summary_df$gamma[best_bic_idx], summary_df$BIC[best_bic_idx]))
    
    if (compute_gic_flag) {
      best_gic_aic_idx <- which.min(summary_df$GIC_AIC)
      best_gic_bic_idx <- which.min(summary_df$GIC_BIC)
      cat(sprintf("  Best GIC_AIC: γ = %.4g (GIC_AIC = %.4f, df_eff = %.2f)\n",
                  summary_df$gamma[best_gic_aic_idx], summary_df$GIC_AIC[best_gic_aic_idx],
                  summary_df$df_eff[best_gic_aic_idx]))
      cat(sprintf("  Best GIC_BIC: γ = %.4g (GIC_BIC = %.4f, df_eff = %.2f)\n",
                  summary_df$gamma[best_gic_bic_idx], summary_df$GIC_BIC[best_gic_bic_idx],
                  summary_df$df_eff[best_gic_bic_idx]))
    }
  }
  
  cat(paste(rep("=", 90), collapse=""), "\n")
  
  return(list(summary = summary_df, results = results))
}

# 7. Two-stage γ search

#' Two-stage γ search: coarse search -> fine search around the optimal γ
#'
#' Stage 1: search coarsely over gammas_coarse and identify the best γ by pAIC/pBIC
#' Stage 2: refine the neighborhood of the best γ in log10 steps of 0.05
#'
#' Skip conditions (when Stage 2 is not performed):
#'   - γ* == 0 (no penalty is best)
#'   - γ* == max(gammas_coarse) (monotone decreasing -> boundary solution)
#'
#' @param fit_func one of fit_msar_penalized, fit_msem_penalized, fit_mgns_penalized
#' @param gammas_coarse coarse γ grid for Stage 1 (e.g., c(0, 10^seq(-2, 4, by=0.5)))
#' @param fine_log10_step log10 step width for Stage 2 (default: 0.05)
#' @param compute_gic_flag whether to compute GIC
#' @param ... arguments passed to fit_func
#' @return same format as compare_gamma: list(summary, results)
compare_gamma_twostage <- function(
  fit_func,
  gammas_coarse = c(0, 10^seq(-2, 4, by = 0.5)),
  fine_log10_step = 0.05,
  compute_gic_flag = TRUE,
  ...
) {
  
  cat("\n", paste(rep("#", 70), collapse=""), "\n")
  cat("### Two-stage γ search ###\n")
  cat(paste(rep("#", 70), collapse=""), "\n\n")
  
  # Stage 1: coarse search
  cat("========== Stage 1: coarse search ==========\n")
  cat(sprintf("  Grid: %d points (%.4g ~ %.4g)\n",
              length(gammas_coarse), min(gammas_coarse), max(gammas_coarse)))
  
  stage1 <- compare_gamma(
    fit_func = fit_func,
    gammas = gammas_coarse,
    compute_gic_flag = compute_gic_flag,
    ...
  )
  
  s1 <- stage1$summary
  
  if (nrow(s1) == 0) {
    warning("Stage 1 produced no valid results")
    return(stage1)
  }
  
  # Identify the best γ by pAIC/pBIC from Stage 1
  gamma_max <- max(gammas_coarse)
  
  if ("GIC_AIC" %in% colnames(s1) && any(!is.na(s1$GIC_AIC))) {
    idx_aic <- which.min(s1$GIC_AIC)
    g_star_aic <- s1$gamma[idx_aic]
    criterion_aic <- "GIC_AIC"
  } else {
    idx_aic <- which.min(s1$AIC)
    g_star_aic <- s1$gamma[idx_aic]
    criterion_aic <- "AIC"
  }
  
  if ("GIC_BIC" %in% colnames(s1) && any(!is.na(s1$GIC_BIC))) {
    idx_bic <- which.min(s1$GIC_BIC)
    g_star_bic <- s1$gamma[idx_bic]
    criterion_bic <- "GIC_BIC"
  } else {
    idx_bic <- which.min(s1$BIC)
    g_star_bic <- s1$gamma[idx_bic]
    criterion_bic <- "BIC"
  }
  
  cat(sprintf("\n  Stage 1 result: best %s γ* = %.4g, best %s γ* = %.4g\n",
              criterion_aic, g_star_aic, criterion_bic, g_star_bic))
  
  # Stage 2: build the refined grid
  # Extract only positive values from the coarse grid and sort (for the log10 transform)
  gammas_positive <- sort(gammas_coarse[gammas_coarse > 0])
  
  build_fine_grid <- function(g_star, criterion_name) {
    # Skip condition 1: γ* == 0
    if (g_star == 0) {
      cat(sprintf("  %s: γ* = 0 -> skip Stage 2 (no penalty is best)\n",
                  criterion_name))
      return(numeric(0))
    }
    # Skip condition 2: γ* == max(gammas_coarse)
    if (g_star >= gamma_max) {
      cat(sprintf("  %s: γ* = %.4g = max(grid) -> skip Stage 2 (boundary solution)\n",
                  criterion_name, g_star))
      return(numeric(0))
    }
    
    pos <- which(abs(gammas_positive - g_star) < 1e-12)
    if (length(pos) == 0) {
      # When g_star is not found in gammas_positive (g_star==0 already handled above)
      cat(sprintf("  %s: γ* = %.4g not in the grid -> skip Stage 2\n",
                  criterion_name, g_star))
      return(numeric(0))
    }
    pos <- pos[1]
    
    # Lower bound: the previous positive γ, or g_star / 10 if none
    if (pos > 1) {
      g_low <- gammas_positive[pos - 1]
    } else {
      g_low <- g_star / 10
    }
    # Upper bound: the next γ, or gamma_max if none
    if (pos < length(gammas_positive)) {
      g_high <- gammas_positive[pos + 1]
    } else {
      g_high <- gamma_max
    }
    
    log10_low  <- log10(g_low)
    log10_high <- log10(g_high)
    fine_grid <- 10^seq(log10_low, log10_high, by = fine_log10_step)
    
    cat(sprintf("  %s: γ* = %.4g -> refine [%.4g, %.4g] (%d points)\n",
                criterion_name, g_star, g_low, g_high, length(fine_grid)))
    
    return(fine_grid)
  }
  
  cat("\n========== Stage 2: preparing the refined search ==========\n")
  fine_aic <- build_fine_grid(g_star_aic, criterion_aic)
  fine_bic <- build_fine_grid(g_star_bic, criterion_bic)
  
  fine_all <- sort(unique(c(fine_aic, fine_bic)))
  
  # Exclude γ already computed in Stage 1 (accounting for floating-point rounding error)
  is_already_computed <- sapply(fine_all, function(g) {
    any(abs(s1$gamma - g) / max(g, 1e-10) < 1e-6)
  })
  gammas_stage2 <- fine_all[!is_already_computed]
  
  if (length(gammas_stage2) == 0) {
    cat("  No additional search needed -> returning the Stage 1 results as-is\n")
    return(stage1)
  }
  
  cat(sprintf("  Stage 2 additional γ: %d points\n", length(gammas_stage2)))
  
  # Stage 2: run the refined search
  cat("\n========== Stage 2: running the refined search ==========\n")
  
  stage2 <- compare_gamma(
    fit_func = fit_func,
    gammas = gammas_stage2,
    compute_gic_flag = compute_gic_flag,
    ...
  )
  
  # Merge Stage 1 + Stage 2
  cat("\n========== Merging results ==========\n")
  
  combined_summary <- rbind(stage1$summary, stage2$summary)
  combined_summary <- combined_summary[order(combined_summary$gamma), ]
  rownames(combined_summary) <- NULL
  
  combined_results <- c(stage1$results, stage2$results)
  
  cat("\n", paste(rep("=", 90), collapse=""), "\n")
  cat("### Two-stage γ search final result summary ###\n")
  cat(sprintf("  Stage 1: %d points, Stage 2: %d points, total: %d points\n",
              nrow(stage1$summary), nrow(stage2$summary), nrow(combined_summary)))
  cat(paste(rep("=", 90), collapse=""), "\n\n")
  
  if (compute_gic_flag && "GIC_AIC" %in% colnames(combined_summary)) {
    summary_numeric <- combined_summary[, c("gamma", "loglik", "AIC", "BIC",
                                             "df_eff", "GIC_AIC", "GIC_BIC", "penalty")]
  } else {
    summary_numeric <- combined_summary[, c("gamma", "loglik", "AIC", "BIC", "penalty")]
  }
  print(summary_numeric, row.names = FALSE, right = FALSE)
  
  if (nrow(combined_summary) > 0) {
    cat("\n--- Selection of the best γ (after two-stage search) ---\n")
    
    best_aic_idx <- which.min(combined_summary$AIC)
    best_bic_idx <- which.min(combined_summary$BIC)
    cat(sprintf("  Best AIC:     γ = %.4g (AIC = %.4f)\n",
                combined_summary$gamma[best_aic_idx], combined_summary$AIC[best_aic_idx]))
    cat(sprintf("  Best BIC:     γ = %.4g (BIC = %.4f)\n",
                combined_summary$gamma[best_bic_idx], combined_summary$BIC[best_bic_idx]))
    
    if (compute_gic_flag && "GIC_AIC" %in% colnames(combined_summary)) {
      best_gic_aic_idx <- which.min(combined_summary$GIC_AIC)
      best_gic_bic_idx <- which.min(combined_summary$GIC_BIC)
      cat(sprintf("  Best GIC_AIC: γ = %.4g (GIC_AIC = %.4f, df_eff = %.2f)\n",
                  combined_summary$gamma[best_gic_aic_idx],
                  combined_summary$GIC_AIC[best_gic_aic_idx],
                  combined_summary$df_eff[best_gic_aic_idx]))
      cat(sprintf("  Best GIC_BIC: γ = %.4g (GIC_BIC = %.4f, df_eff = %.2f)\n",
                  combined_summary$gamma[best_gic_bic_idx],
                  combined_summary$GIC_BIC[best_gic_bic_idx],
                  combined_summary$df_eff[best_gic_bic_idx]))
    }
  }
  
  cat(paste(rep("=", 90), collapse=""), "\n")
  
  return(list(summary = combined_summary, results = combined_results))
}
