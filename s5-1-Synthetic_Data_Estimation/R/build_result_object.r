# build_result_object.r
# Helper functions to unify the output format of all fit functions

# 0. Helper functions: extracting structure from the beta vector

#' Extract the time-lag matrix A from the beta vector
#'
#' The last k^2 entries of beta are the elements of A (row-major)
#' A[i,j] = beta[p0_total + (i-1)*k + j]
#'
#' @param beta p×1 regression coefficient vector
#' @param data_list output of prepare_data_extended()
#' @param k number of variables
#' @return k×k A matrix
extract_alpha <- function(beta, data_list, k) {
  x_vars <- data_list$data_info$x_vars
  include_intercept <- data_list$data_info$include_intercept

  # p0_total: number of columns in the regressor part (intercept + x_vars, all variables)
  p0_total <- 0
  for (i in 1:k) {
    p0_total <- p0_total + length(x_vars[[i]])
    if (include_intercept) p0_total <- p0_total + 1
  }

  alpha_idx <- (p0_total + 1):(p0_total + k^2)
  if (max(alpha_idx) > length(beta)) {
    stop("beta vector too short: cannot extract the time-lag part")
  }

  alpha_vec <- beta[alpha_idx]
  A <- matrix(alpha_vec, nrow = k, ncol = k, byrow = TRUE)

  y_vars <- data_list$data_info$y_vars
  rownames(A) <- y_vars
  colnames(A) <- paste0(y_vars, "_lag")
  
  return(A)
}

#' Structure the beta vector into a named list by variable
#'
#' @param beta p×1 regression coefficient vector
#' @param data_list output of prepare_data_extended()
#' @param k number of variables
#' @return named list keyed by variable name
structure_beta <- function(beta, data_list, k) {
  x_vars <- data_list$data_info$x_vars
  y_vars <- data_list$data_info$y_vars
  include_intercept <- data_list$data_info$include_intercept

  beta0 <- list()
  idx <- 1

  for (i in 1:k) {
    coef_names <- c()
    coef_vals <- c()

    if (include_intercept) {
      coef_names <- c(coef_names, "(Intercept)")
      coef_vals <- c(coef_vals, beta[idx])
      idx <- idx + 1
    }

    for (x_name in x_vars[[i]]) {
      coef_names <- c(coef_names, x_name)
      coef_vals <- c(coef_vals, beta[idx])
      idx <- idx + 1
    }
    
    names(coef_vals) <- coef_names
    beta0[[y_vars[i]]] <- coef_vals
  }
  
  return(beta0)
}

# 1. Main function: build_result_object

#' Build an object with the unified output structure
#'
#' @param model_type model type string
#'   "MSAR", "MSEM", "MGNS",
#'   "IndSAR", "IndSEM", "IndGNS",
#'   "VARX", "IndReg"
#' @param R spatial lag matrix (NULL if not applicable to the model)
#' @param Lambda_mat spatial error matrix (NULL if not applicable to the model)
#' @param beta regression coefficient vector (omitted if NULL)
#' @param Sigma error covariance matrix
#' @param loglik non-penalized log-likelihood
#' @param num_params number of parameters
#' @param converged convergence flag
#' @param method name of the estimation method
#' @param iterations number of iterations (NA allowed)
#' @param data_list data list returned by prepare_data_extended
#' @param gamma penalty strength (0 or NULL = no penalty)
#' @param penalty_value penalty amount
#' @param penalized_loglik penalized log-likelihood
#' @param execution_time execution time (difftime)
#' @param individual_models list of per-variable estimation results for the diagonal model (NULL allowed)
#' @param individual_estimates list of per-variable results from the S1 initial estimation (NULL allowed)
#' @param std_errors_R standard-error matrix of the R matrix (NULL allowed)
#' @param std_errors_Lambda standard-error matrix of the Lambda matrix (NULL allowed)
#' @param beta0 per-variable beta (list form; auto-built from beta if NULL)
#' @param alpha time-lag matrix A (auto-extracted from beta if NULL)
#' @param residuals_raw residual vector (auto-computed if NULL)
#' @param residuals_std standardized residual vector (auto-computed if NULL)
#' @return object of class multivar_spatial
#'
build_result_object <- function(
  model_type,
  R = NULL,
  Lambda_mat = NULL,
  beta = NULL,
  Sigma,
  loglik,
  num_params,
  converged,
  method,
  iterations = NA,
  data_list,
  gamma = NULL,
  penalty_value = NULL,
  penalized_loglik = NULL,
  execution_time = NULL,
  individual_models = NULL,
  individual_estimates = NULL,
  std_errors_R = NULL,
  std_errors_Lambda = NULL,
  beta0 = NULL,
  alpha = NULL,
  residuals_raw = NULL,
  residuals_std = NULL
) {
  
  k <- data_list$k
  n <- data_list$n
  
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE

  # Information criteria
  ic <- compute_information_criteria(loglik, num_params, k * n)

  # Build beta0 (per-variable beta)
  if (is.null(beta0) && !is.null(beta) && !is.null(data_list$data_info$x_vars)) {
    beta0 <- tryCatch({
      structure_beta(beta, data_list, k)
    }, error = function(e) NULL)
  }

  # Extract alpha (time-lag matrix)
  if (is.null(alpha) && !is.null(beta) && include_time_lag) {
    alpha <- tryCatch({
      extract_alpha(beta, data_list, k)
    }, error = function(e) NULL)
  }

  # Compute residuals
  if (is.null(residuals_raw) && !is.null(beta)) {
    residuals_raw <- tryCatch({
      if (model_type %in% c("MSAR", "IndSAR")) {
        compute_residuals(R, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else if (model_type %in% c("MSEM", "IndSEM")) {
        compute_residuals_msem(Lambda_mat, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else if (model_type %in% c("MGNS", "IndGNS")) {
        compute_residuals_mgns(R, Lambda_mat, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else {
        # OLS / VARX: y - Xβ
        data_list$y - data_list$X %*% beta
      }
    }, error = function(e) NULL)
  }

  if (is.null(residuals_std) && !is.null(residuals_raw)) {
    residuals_std <- as.numeric(residuals_raw)
    for (i in 1:k) {
      idx <- ((i - 1) * n + 1):(i * n)
      sigma_i <- sqrt(Sigma[i, i])
      if (sigma_i > 0) {
        residuals_std[idx] <- residuals_std[idx] / sigma_i
      }
    }
  }
  
  # Compute pseudo R²
  # From ŷ = (I - R⊗W)⁻¹ Xβ:  R²_pseudo,k = corr(y_k, ŷ_k)²
  R2 <- NULL
  R2_adj <- NULL
  if (!is.null(beta)) {
    fitted_vals <- tryCatch({
      compute_fitted_multivar(
        model_type = model_type, R = R, beta = beta,
        X = data_list$X, W = data_list$W, k = k, n = n
      )
    }, error = function(e) NULL)
    
    if (!is.null(fitted_vals)) {
      r2_result <- tryCatch({
        compute_r_squared_multivar(
          y = data_list$y, fitted = fitted_vals,
          k = k, n = n
        )
      }, error = function(e) NULL)
      
      if (!is.null(r2_result)) {
        R2 <- r2_result$R2_pseudo_mean
        R2_adj <- NA  # no adjusted version for pseudo R²
      }
    }
  }

  # penalty section
  penalty_section <- NULL
  if (!is.null(gamma) && gamma > 0) {
    penalty_section <- list(
      gamma = gamma,
      value = if (!is.null(penalty_value)) penalty_value else 0,
      penalized_loglik = if (!is.null(penalized_loglik)) penalized_loglik else loglik
    )
  }
  
  # Determine the class name
  class_map <- list(
    "MSAR"           = c("multivar_msar", "multivar_spatial"),
    "MSEM"           = c("multivar_msem", "multivar_spatial"),
    "MGNS"          = c("multivar_mgns", "multivar_spatial"),
    "IndSAR"  = c("multivar_indsar", "multivar_msar", "multivar_spatial"),
    "IndSEM"  = c("multivar_indsem", "multivar_msem", "multivar_spatial"),
    "IndGNS" = c("multivar_indgns", "multivar_mgns", "multivar_spatial"),
    "VARX"          = c("multivar_varx", "multivar_spatial"),
    "IndReg"  = c("multivar_indreg", "multivar_indreg", "multivar_spatial")
  )
  
  obj_class <- class_map[[model_type]]
  if (is.null(obj_class)) {
    obj_class <- c("multivar_spatial")
    warning(sprintf("Unknown model_type: %s", model_type))
  }

  # Build the object
  result <- structure(
    list(
      model_type = model_type,
      model_description = get_model_description(model_type),

      coefficients = list(
        R     = R,
        Lambda     = Lambda_mat,
        beta  = beta,
        beta0 = beta0,
        alpha = alpha,
        Sigma = Sigma
      ),

      fit = list(
        loglik     = loglik,
        AIC        = ic$AIC,
        BIC        = ic$BIC,
        num_params = num_params,
        num_obs    = k * n,
        R2         = R2,
        R2_adj     = R2_adj
      ),

      convergence = list(
        converged  = converged,
        iterations = iterations,
        method     = method,
        message    = ifelse(converged, "converged", "max iterations reached or warning")
      ),

      # Penalty information (full models only, when γ>0)
      penalty = penalty_section,

      data_info = list(
        n                 = n,
        k                 = k,
        y_vars            = data_list$data_info$y_vars,
        x_vars            = data_list$data_info$x_vars,
        time_point_used   = data_list$data_info$time_point,
        time_lag_used     = data_list$data_info$time_lag,
        include_time_lag  = include_time_lag,
        region_var        = data_list$data_info$region_var,
        time_var          = data_list$data_info$time_var,
        include_intercept = data_list$data_info$include_intercept
      ),

      residuals = list(
        raw          = if (!is.null(residuals_raw)) as.numeric(residuals_raw) else NULL,
        standardized = if (!is.null(residuals_std)) as.numeric(residuals_std) else NULL
      ),

      model_data = list(
        y       = data_list$y,
        X       = data_list$X,
        W       = data_list$W,
        W_listw = data_list$W_listw,
        y_lag   = data_list$y_lag,
        eigen_W = data_list$eigen_W
      ),

      # Initial values (S1 estimation results of the full model)
      initial_values = if (!is.null(individual_estimates)) {
        list(individual_estimates = individual_estimates)
      } else {
        NULL
      },

      # Standard errors (already obtained from spatialreg for diagonal models)
      std_errors = list(
        beta = NULL,
        R    = std_errors_R,
        Lambda    = std_errors_Lambda
      ),

      # Inference results (added later by add_inference)
      inference = NULL,
      vcov      = NULL,
      hessian   = NULL,

      individual_models = individual_models,

      execution = list(
        time      = execution_time,
        call      = sys.call(-1),
        R_version = R.version.string
      ),

      # The original data_list (needed by downstream processing)
      data_list = data_list
    ),
    class = obj_class
  )
  
  return(result)
}

# 2. Model-name helper

#' Get the description string from a model type
#'
#' @param model_type model type string
#' @return description string
get_model_description <- function(model_type) {
  descriptions <- list(
    "MSAR"           = "Multivariate Spatial Lag of Y Model (full)",
    "MSEM"           = "Multivariate Spatial Error Model (full)",
    "MGNS"          = "Multivariate General Nesting Spatial Model (full)",
    "IndSAR"  = "Multivariate Spatial Lag of Y Model (diagonal)",
    "IndSEM"  = "Multivariate Spatial Error Model (diagonal)",
    "IndGNS" = "Multivariate General Nesting Spatial Model (diagonal)",
    "VARX"          = "Vector Autoregressive with Exogenous Variables Model",
    "IndReg"  = "Independent Regression Model (no spatial dependence)"
  )
  desc <- descriptions[[model_type]]
  if (is.null(desc)) return("Unknown Model")
  return(desc)
}

# 3. For diagonal models: minimal data_list construction helper

#' Build a minimal data_list for diagonal models
#'
#' Diagonal models do not use prepare_data_extended, so this creates the
#' minimal data_list structure to pass to build_result_object.
#'
#' @param y y vector (k*n × 1)
#' @param X X matrix (NULL allowed; may be unnecessary since spatialreg handles it internally for diagonal models)
#' @param W spatial weight matrix (n × n)
#' @param W_listw listw object
#' @param k number of variables
#' @param n number of regions
#' @param y_vars dependent variable names
#' @param x_vars list of regressor names
#' @param time_point time point used
#' @param include_time_lag whether a time lag is included
#' @param include_intercept whether an intercept is included
#' @param region_var region variable name
#' @param time_var time variable name
#' @return data_list-compatible structure
#'
build_data_list_from_parts <- function(
  y, X = NULL, W = NULL, W_listw = NULL,
  k, n, y_vars, x_vars,
  time_point = NULL,
  include_time_lag = FALSE,
  include_intercept = TRUE,
  region_var = "region",
  time_var = "time"
) {
  
  eigen_W <- NULL
  if (!is.null(W)) {
    eigen_W <- tryCatch({
      eigen(W, only.values = TRUE)$values
    }, error = function(e) NULL)
  }
  
  list(
    y = y,
    X = X,
    W = W,
    W_listw = W_listw,
    y_lag = NULL,
    eigen_W = eigen_W,
    k = k,
    n = n,
    data_info = list(
      y_vars = y_vars,
      x_vars = x_vars,
      time_point = time_point,
      time_lag = NULL,
      include_time_lag = include_time_lag,
      include_intercept = include_intercept,
      region_var = region_var,
      time_var = time_var
    )
  )
}
