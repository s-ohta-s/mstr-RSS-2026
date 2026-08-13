# sim_parameters.r
# Parameter sets for the Monte Carlo redesign
# Four dependence levels: weak / mid / strong (common A, Σ, β; stationary,
# burn-in generation) and realistic (real-data estimates from
# Kansai_Estimation/output/comparison_table_pAIC.csv; EXPLOSIVE in time ->
# conditional single-transition generation).
# DGP-level zeroing (1111 / 0111 / 1011) is applied by get_dgp_params().

## ---------------------------------------------------------------- common β
## weak/mid/strong share the historical base_params (X ~ N(0,1))
.SIM_BETA_COMMON <- list(
  beta_intercept   = c(0.2, -0.1),
  beta_common1     = c(1.2,  0.8),
  beta_common2     = c(-0.6, -0.3),
  beta_specific1_1 = 0.8,
  beta_specific2_1 = 0.5
)

## Common A (spectral ≈ 0.32; keeps the strong level jointly stationary) and common Σ
.SIM_A_COMMON <- matrix(c(0.30, 0.03,
                          0.05, 0.25), 2, 2, byrow = TRUE)
.SIM_SIGMA_COMMON <- matrix(c(0.10, 0.03,
                              0.03, 0.10), 2, 2, byrow = TRUE)

## ---------------------------------------------------------------- level sets
SIM_PARAMSETS <- list(

  weak = c(list(
    R      = matrix(c(0.12, 0.04, 0.03, 0.10), 2, 2, byrow = TRUE),
    Lambda = matrix(c(0.10, 0.03, 0.02, 0.08), 2, 2, byrow = TRUE),
    A      = .SIM_A_COMMON,
    Sigma  = .SIM_SIGMA_COMMON,
    generation = "burnin"
  ), .SIM_BETA_COMMON),

  mid = c(list(
    R      = matrix(c(0.30, 0.10, 0.08, 0.25), 2, 2, byrow = TRUE),
    Lambda = matrix(c(0.25, 0.08, 0.05, 0.20), 2, 2, byrow = TRUE),
    A      = .SIM_A_COMMON,
    Sigma  = .SIM_SIGMA_COMMON,
    generation = "burnin"
  ), .SIM_BETA_COMMON),

  strong = c(list(
    R      = matrix(c(0.55, 0.15, 0.12, 0.45), 2, 2, byrow = TRUE),
    Lambda = matrix(c(0.50, 0.12, 0.10, 0.40), 2, 2, byrow = TRUE),
    A      = .SIM_A_COMMON,
    Sigma  = .SIM_SIGMA_COMMON,
    generation = "burnin"
  ), .SIM_BETA_COMMON),

  ## Realistic set: per-DGP true values = the real-data estimates
  ## (Kansai_Estimation/output/comparison_table_pAIC.csv, estimate.{1111,0111,1011}).
  ## All three are temporally explosive (ρ(M_dyn) ≈ 1.07 / 1.026 / 1.023)
  ## -> generation = "conditional" (single transition).
  realistic = list(
    generation = "conditional",

    dgp = list(
      "1111" = list(
        R      = matrix(c(-0.057669572, 0.012324534,
                          -0.044217939, 0.033401891), 2, 2, byrow = TRUE),
        Lambda = matrix(c( 0.516369645, 0.110591912,
                          -0.025634856, 0.302714362), 2, 2, byrow = TRUE),
        A      = matrix(c( 1.037044530, -0.021314174,
                          -0.028430180,  0.916785162), 2, 2, byrow = TRUE),
        Sigma  = matrix(c( 0.021996041, -0.001528895,
                          -0.001528895,  0.019506258), 2, 2, byrow = TRUE),
        beta_intercept   = c(2.605430033, 3.652736966),
        beta_common1     = c(-0.275504915, -0.490846230),
        beta_common2     = c(-0.012884830,  0.020352708),
        beta_specific1_1 = -0.150158559,
        beta_specific2_1 =  0.000889674
      ),
      "0111" = list(
        R      = matrix(0, 2, 2),
        Lambda = matrix(c( 0.479376505, -0.136536809,
                           0.150606682,  0.374415400), 2, 2, byrow = TRUE),
        A      = matrix(c( 1.024122793, -0.004195401,
                          -0.046656678,  0.937071517), 2, 2, byrow = TRUE),
        Sigma  = matrix(c( 0.022897354, -0.001260798,
                          -0.001260798,  0.019808341), 2, 2, byrow = TRUE),
        beta_intercept   = c(1.982256443, 3.155625299),
        beta_common1     = c(-0.198605728, -0.431142004),
        beta_common2     = c(-0.013328957,  0.021172079),
        beta_specific1_1 = -0.149969544,
        beta_specific2_1 =  0.000971082
      ),
      "1011" = list(
        R      = matrix(c(-0.002356956, 0.015464768,
                          -0.026124980, 0.031793231), 2, 2, byrow = TRUE),
        Lambda = matrix(0, 2, 2),
        A      = matrix(c( 1.016852035, -0.012828768,
                          -0.042452876,  0.923650580), 2, 2, byrow = TRUE),
        Sigma  = matrix(c( 0.027198906, -0.001151820,
                          -0.001151820,  0.020875362), 2, 2, byrow = TRUE),
        beta_intercept   = c(1.806478977, 3.808447930),
        beta_common1     = c(-0.184540514, -0.505670106),
        beta_common2     = c(-0.005007925,  0.013953814),
        beta_specific1_1 = -0.156635951,
        beta_specific2_1 =  0.000963422
      )
    ),

    ## Sample moments of the real data (Data_I/transformed_data198_for_R.csv,
    ## n = 198). y_{t-1} from time = 1 (with cross-correlation); X from time = 2.
    y_prev_moments = list(
      mean = c(y1 = -0.175212, y2 = 0.058638),
      sd   = c(y1 =  0.938637, y2 = 1.005752),
      cor12 = -0.252493
    ),
    x_moments = list(
      x_common1     = c(mean = 8.019864,  sd = 0.135370),
      x_common2     = c(mean = 6.831375,  sd = 0.823563),
      x_specific1_1 = c(mean = -0.404646, sd = 0.604984),
      x_specific2_1 = c(mean = 38.490404, sd = 22.595992)
    )
  )
)

#' Resolve the (paramset, dgp) cell into a concrete parameter list
#'
#' For weak/mid/strong the DGP id zeroes R (0111) or Λ (1011). For realistic
#' the per-DGP real-data estimates are returned as-is.
#'
#' @param paramset_id "weak" | "mid" | "strong" | "realistic"
#' @param dgp_id "1111" | "0111" | "1011"
#' @return list(R, Lambda, A, Sigma, beta_*..., generation, [y_prev_moments,
#'         x_moments for realistic])
get_dgp_params <- function(paramset_id, dgp_id) {
  ps <- SIM_PARAMSETS[[paramset_id]]
  if (is.null(ps)) stop(sprintf("unknown paramset '%s'", paramset_id))

  if (paramset_id == "realistic") {
    p <- ps$dgp[[dgp_id]]
    if (is.null(p)) stop(sprintf("realistic set has no DGP '%s'", dgp_id))
    p$generation <- ps$generation
    p$y_prev_moments <- ps$y_prev_moments
    p$x_moments <- ps$x_moments
    return(p)
  }

  p <- ps
  chars <- strsplit(dgp_id, "")[[1]]
  if (chars[1] == "0") p$R <- matrix(0, 2, 2)
  if (chars[2] == "0") p$Lambda <- matrix(0, 2, 2)
  if (chars[3] == "0") p$A <- matrix(0, 2, 2)
  # Σ is never zeroed in the A/B experiments (all DGPs are *1**-type on Σ)
  p
}

#' Flatten a parameter list into a true-value table
#'
#' Parameter names follow extract_params_uniform() (experiment_output_functions.r)
#' so the Monte Carlo driver can join estimates and true values directly:
#' beta_intercept_y1, beta_x_common1_y1, ..., A[i,j], R[i,j], Lambda[i,j],
#' Sigma[i,j].
params_to_true_df <- function(p) {
  mat_rows <- function(M, name) {
    data.frame(
      parameter = sprintf("%s[%d,%d]", name,
                          rep(1:2, each = 2), rep(1:2, times = 2)),
      true = as.vector(t(M)),
      stringsAsFactors = FALSE
    )
  }
  rbind(
    data.frame(
      parameter = c("beta_intercept_y1", "beta_x_common1_y1",
                    "beta_x_common2_y1", "beta_x_specific1_1_y1",
                    "beta_intercept_y2", "beta_x_common1_y2",
                    "beta_x_common2_y2", "beta_x_specific2_1_y2"),
      true = c(p$beta_intercept[1], p$beta_common1[1], p$beta_common2[1],
               p$beta_specific1_1,
               p$beta_intercept[2], p$beta_common1[2], p$beta_common2[2],
               p$beta_specific2_1),
      stringsAsFactors = FALSE
    ),
    mat_rows(p$A, "A"),
    mat_rows(p$R, "R"),
    mat_rows(p$Lambda, "Lambda"),
    mat_rows(p$Sigma, "Sigma")
  )
}

#' True β vector in the design-matrix column order (incl. the y_lag block)
#'
#' Order: [intercept, x_common1, x_common2, x_specific1_1 | intercept,
#' x_common1, x_common2, x_specific2_1 | α11, α12, α21, α22] — matches
#' build_design_matrix_extended() with the standard x_vars and time lag.
params_to_beta_vector <- function(p, include_time_lag = TRUE) {
  b <- c(p$beta_intercept[1], p$beta_common1[1], p$beta_common2[1], p$beta_specific1_1,
         p$beta_intercept[2], p$beta_common1[2], p$beta_common2[2], p$beta_specific2_1)
  if (include_time_lag) b <- c(b, as.vector(t(p$A)))
  b
}
