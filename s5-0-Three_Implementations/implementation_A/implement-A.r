################################################################################
# implement-A.r
#
# Integrated implementation of the Multivariate General Nesting Spatial (MGNS)
# model within the Multivariate Spatio-Temporal Regression (MSTR) framework.
#
# Theory reference: mstr.pdf
#   Section 3.1 — MGNS model definition (Eqs. 5–6, 8)
#   Section 3.2 — Maximum likelihood estimation (Eqs. 9–16)
#   Section 3.3 — Estimation algorithm S1–S5 (Steps S1–S5)
#   Section 3.4 — Efficient determinant computation (Eqs. 18–23)
#   Section 3.5 — Model hierarchy / Table 1 (11 model IDs)
#   Section 3.6 — Penalised information criteria pAIC / pBIC (Eqs. 24–26)
#   Appendix A  — Penalised likelihood formulation (Eq. 36)
#
# Model specification (Table 1):
#   ID    R      Λ      A      Σ        Name
#   0011  0      0      full   full     VARX
#   1011  0      full   full   full     MSAR
#   0111  0      full   full   full     MSEM
#   1111  full   full   full   full     MGNS (most general)
#   1001  full   0      0      full     MSAR without time-AR
#   0101  0      full   0      full     MSEM without time-AR
#   1101  full   full   0      full     MGNS without time-AR
#   000d  0      0      0      diag     Independent OLS
#   d0dd  diag   0      diag   diag     Diagonal SAR
#   0ddd  0      diag   diag   diag     Diagonal SEM
#   dddd  diag   diag   diag   diag     Diagonal GNS
#
# Estimation workflow:
#   Phase 1  Data preparation and spatial matrix utilities (Eqs. 7, 18–23)
#   Phase 2  Model fitting: initial estimates (S1), inner loop (S2),
#            profile likelihood (S3), BFGS optimisation of (R, Λ) (S4–S5)
#   Phase 3  Post-estimation inference: standard errors via Ψ (Eq. 16)
#            and full numerical Hessian sandwich estimator (Appendix A)
#   Phase 4  Output tables: pAIC/pBIC model comparison, bias assessment
#   Phase 5  Save results to RDS and CSV
#
# Required packages: spdep, spatialreg, systemfit, numDeriv
################################################################################

rm(list = ls())

################################################################################
# Program-wide elapsed-time measurement (START marker)
# ---------------------------------------------------------------------------
# We record a timestamp at the very beginning of execution — placed *after*
# `rm(list = ls())` so this variable is not itself wiped by that call. Because
# it sits before all `library()` calls and function definitions, the measured
# span covers essentially the entire program: package loading, source setup,
# and all five estimation phases.
#
# `.program_start_time` is wall-clock time (Sys.time) used for the headline
# "elapsed" figure. `.program_start_proc` additionally captures CPU time via
# proc.time(), letting us also report user/system CPU usage at the end.
# The leading dot keeps these names out of the way of the analysis variables.
################################################################################
.program_start_time <- Sys.time()
.program_start_proc <- proc.time()
cat(sprintf("[timing] Program started at %s\n",
            format(.program_start_time, "%Y-%m-%d %H:%M:%S")))

################################################################################
# Parallel infrastructure (multi-core support)
# ---------------------------------------------------------------------------
# The estimation ALGORITHM is unchanged. We only parallelise computations that
# are embarrassingly parallel (mutually independent), namely:
#   (1) compare_gamma()  Phase 1: the fit at each gamma value on the grid
#   (2) compare_gamma()  Phase 2: the profile-likelihood Hessian at each gamma
#   (3) compute_hessian_numerical(): the (i,j) entries of the central-difference
#       Hessian (each entry needs 4 independent likelihood evaluations)
#   (4) run_experiment Phase 3: add_inference()/add_full_inference() per model
# Each parallel task computes exactly the same quantities, in the same way, as
# the original sequential loop; only the order of evaluation changes, and all
# computations are deterministic (no RNG), so the results are identical.
#
# Configuration:
#   - Number of worker cores: set the environment variable MSTR_CORES, or
#     options(mstr.cores = N), before sourcing this script. Default:
#     (number of physical cores) - 1, at least 1.
#   - Setting MSTR_CORES=1 restores fully sequential execution.
#   - On Unix/macOS, fork-based parallel::mclapply is used (zero-copy, fast).
#     On Windows, a PSOCK cluster is created lazily on first use; the global
#     environment (all functions/config defined by this script) is exported to
#     the workers at that time.
#
# Note: when running many R workers, multi-threaded BLAS (OpenBLAS/MKL) can
# oversubscribe the CPU. We therefore pin BLAS to 1 thread per process below
# (best effort; harmless if BLAS is single-threaded).
################################################################################

suppressPackageStartupMessages(library(parallel))

MSTR_CORES <- local({
  env_val <- suppressWarnings(as.integer(Sys.getenv("MSTR_CORES", "")))
  if (!is.na(env_val) && env_val >= 1L) return(env_val)
  opt_val <- suppressWarnings(as.integer(getOption("mstr.cores", NA_integer_)))
  if (!is.na(opt_val) && opt_val >= 1L) return(opt_val)
  n_cores <- tryCatch(parallel::detectCores(logical = FALSE),
                      error = function(e) NA_integer_)
  if (is.na(n_cores)) n_cores <- tryCatch(parallel::detectCores(),
                                          error = function(e) 1L)
  if (is.na(n_cores)) n_cores <- 1L
  max(1L, n_cores - 1L)
})

# Best-effort: avoid BLAS-thread oversubscription when running many workers
if (MSTR_CORES > 1L) {
  Sys.setenv(OMP_NUM_THREADS = "1",
             OPENBLAS_NUM_THREADS = "1",
             MKL_NUM_THREADS = "1")
}

.mstr_use_fork <- (.Platform$OS.type == "unix")
.mstr_cluster  <- NULL
.mstr_exported <- character(0)  # symbols already exported to PSOCK workers

#' Lazily create (and cache) a PSOCK cluster for Windows.
#' All objects in the global environment are exported so that the workers can
#' call every function defined in this script.
mstr_get_cluster <- function() {
  if (!is.null(.mstr_cluster)) return(.mstr_cluster)
  cl <- parallel::makeCluster(MSTR_CORES)
  parallel::clusterCall(cl, function() {
    Sys.setenv(MSTR_IN_WORKER = "1",
               OMP_NUM_THREADS = "1",
               OPENBLAS_NUM_THREADS = "1",
               MKL_NUM_THREADS = "1")
    suppressPackageStartupMessages({
      library(spdep)
      library(spatialreg)
      requireNamespace("numDeriv", quietly = TRUE)
      requireNamespace("Matrix",   quietly = TRUE)
    })
    invisible(NULL)
  })
  syms <- ls(globalenv())
  parallel::clusterExport(cl, varlist = syms, envir = globalenv())
  .mstr_exported <<- syms
  .mstr_cluster  <<- cl
  cl
}

#' Export to the PSOCK workers any FUNCTION defined in the global environment
#' after the cluster was created (large data objects are deliberately not
#' re-exported; data should be passed through the job list instead).
mstr_sync_cluster_functions <- function(cl) {
  new_syms <- setdiff(ls(globalenv()), .mstr_exported)
  if (length(new_syms) == 0) return(invisible(NULL))
  new_funs <- Filter(function(nm) {
    is.function(get(nm, envir = globalenv()))
  }, new_syms)
  if (length(new_funs) > 0) {
    parallel::clusterExport(cl, varlist = new_funs, envir = globalenv())
    .mstr_exported <<- union(.mstr_exported, new_funs)
  }
  invisible(NULL)
}

#' Stop the cached PSOCK cluster (no-op if none exists).
mstr_stop_cluster <- function() {
  if (!is.null(.mstr_cluster)) {
    try(parallel::stopCluster(.mstr_cluster), silent = TRUE)
    .mstr_cluster <<- NULL
  }
  invisible(NULL)
}

#' TRUE when executing inside a parallel worker (used to forbid nested
#' parallelism, which would oversubscribe the CPU).
mstr_in_worker <- function() {
  identical(Sys.getenv("MSTR_IN_WORKER"), "1")
}

#' Parallel lapply wrapper.
#' Falls back to plain lapply when: only one job, one core configured, or
#' already running inside a worker (no nested parallelism). Uses fork-based
#' mclapply on Unix (load-balanced, mc.preschedule = FALSE because the cost
#' per gamma/model varies a lot), and a PSOCK cluster with parLapplyLB on
#' Windows.
mstr_parallel_lapply <- function(X, FUN) {
  n_jobs <- length(X)
  if (n_jobs <= 1L || MSTR_CORES <= 1L || mstr_in_worker()) {
    return(lapply(X, FUN))
  }
  if (.mstr_use_fork) {
    wrapped <- function(x) {
      Sys.setenv(MSTR_IN_WORKER = "1")
      FUN(x)
    }
    parallel::mclapply(X, wrapped,
                       mc.cores = min(MSTR_CORES, n_jobs),
                       mc.preschedule = FALSE)
  } else {
    cl <- mstr_get_cluster()
    mstr_sync_cluster_functions(cl)
    parallel::parLapplyLB(cl, X, FUN)
  }
}

cat(sprintf("[parallel] Using up to %d worker core(s) (%s backend)\n",
            MSTR_CORES,
            if (.mstr_use_fork) "fork/mclapply" else "PSOCK cluster"))

################################################################################
# START OF FILE: multivar_sly_phase1_1.R
# Phase 1: Foundation layer for the Multivariate Spatio-Temporal Regression (MSTR)
# framework described in mstr.pdf.
#
# This file provides:
#   - prepare_data()            : Load panel CSV + spatial weight matrix; build y, X, W
#   - build_design_matrix()     : Construct block-diagonal design matrix X_t (Eq. 7)
#   - compute_eigen_W()         : Pre-compute eigenvalues of W for determinant formulae
#   - compute_elem_sym_poly()   : Elementary symmetric polynomials e_j(R) (Eqs. 19–23)
#   - log_det_spatial()         : log|I_{Kn} - R⊗W| via eigenvalue expansion (Eq. 18)
#   - spectral_radius()         : Check stationarity condition |I - M·ω_i| > 0
#   - compute_RW_times_y()      : Efficient computation of (R⊗W)y without Kronecker fill
#   - check_stationarity()      : Simple stationarity check max|ρ_{ij}| < 1
################################################################################


################################################################################
# Phase 1: Foundation implementation
# Multivariate Spatial Regressive (MGNS) model — data preparation and
# spatial matrix utilities corresponding to Sections 3.1–3.4 of mstr.pdf.
################################################################################

################################################################################
# File: data_preparation.R
# Data loading, reshaping, and design matrix construction.
# Implements the composite design matrix X_t = [X*_t, Y_{t-1}] from Eq. (4),
# which combines exogenous regressors and AR(1) lag terms for each response.
################################################################################

#' Load panel data from CSV and reshape for MGNS/MSTR model estimation.
#'
#' Constructs the stacked response vector y_t (Kn×1) from K response variables
#' observed at n regions and time t, together with the lag matrix Y_{t-1} (n×K)
#' used as the AR(1) component in the composite design matrix X_t (Eq. 7).
#'
#' @param data_file   Path to panel data CSV (columns: region, time, y_vars, x_vars)
#' @param weight_file Path to row-normalised spatial weight matrix CSV
#' @param y_vars      Character vector of K response variable names
#' @param x_vars      Named list of exogenous covariate names for each response
#' @param time_var    Name of the time index column (default: "time")
#' @param time_point  Time period to use (default: maximum available)
#' @param region_var  Name of the region identifier column (default: "region")
#' @param include_intercept Logical; include an intercept column (default: TRUE)
#' @param verbose     Logical; print progress messages
#'
#' @return A list with elements: y (Kn×1), X (Kn×p), W (n×n), W_listw,
#'         y_lag (n×K), eigen_W, n, k, p0, data_info
prepare_data <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  verbose = FALSE
) {
  
  if (verbose) cat("=== Data preparation start ===\n")
  
  # Load panel data from CSV
  if (verbose) cat("Loading panel data...\n")
  data <- read.csv(data_file, stringsAsFactors = FALSE)
  
  # Robust CSV reading: first try with header=TRUE to detect ID column
  if (verbose) cat("Loading spatial weight matrix...\n")
  
  # Read with header so we can inspect the first column
  W_raw <- read.csv(weight_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Detect whether the first column is a region ID column
  first_col <- W_raw[, 1]
  is_id_column <- !is.numeric(first_col) || all(first_col == 1:nrow(W_raw))
  
  if (is_id_column) {
    # Exclude the first column (region ID), keep only the numeric weight entries
    W <- as.matrix(W_raw[, -1])
    if (verbose) cat("First column treated as region ID and excluded.\n")
  } else {
    # All columns are numeric weight entries
    W <- as.matrix(W_raw)
  }
  
  # Drop row/column names from W to ensure clean numeric matrix
  dimnames(W) <- NULL
  
  if (verbose) {
    cat(sprintf("Spatial weight matrix dimension: %d x %d\n", nrow(W), ncol(W)))
  }
  
  # Basic validation: check that required columns exist in the data frame
  if (!all(c(time_var, region_var, y_vars) %in% colnames(data))) {
    stop("Required columns not found")
  }
  
  # Determine the number of response variables K and regions n
  k <- length(y_vars)
  regions <- sort(unique(data[[region_var]]))
  n <- length(regions)
  
  if (verbose) {
    cat(sprintf("Number of variables k = %d\n", k))
    cat(sprintf("Number of regions n = %d\n", n))
  }
  
  # Verify that W has dimensions n×n matching the number of regions
  if (nrow(W) != n || ncol(W) != n) {
    stop(sprintf("Spatial weight matrix size (%dx%d) does not match the number of regions (%d)", 
                 nrow(W), ncol(W), n))
  }
  
  # Verify that W is row-normalised (each row sums to 1)
  row_sums <- rowSums(W)
  if (!all(abs(row_sums - 1) < 1e-6)) {
    warning("The spatial weight matrix may not be row-normalised")
  }
  
  # Identify available time periods and select the analysis time point
  times <- sort(unique(data[[time_var]]))
  if (is.null(time_point)) {
    time_point <- max(times)
    if (verbose) cat(sprintf("time_point not specified; using the latest time point %d\n", time_point))
  }
  
  if (!time_point %in% times) {
    stop(sprintf("The specified time point %d does not exist in the data", time_point))
  }
  
  # Verify that the lag period t-1 exists in the data (needed for AR(1) terms)
  time_lag <- time_point - 1
  if (!time_lag %in% times) {
    stop(sprintf("Lag time point %d does not exist in the data. At least two time points are required.", time_lag))
  }
  
  if (verbose) cat(sprintf("Time point used: t=%d, lag time point: t=%d\n", time_point, time_lag))
  
  # Extract data rows for the current period t and the lag period t-1
  data_current <- data[data[[time_var]] == time_point, ]
  data_lag <- data[data[[time_var]] == time_lag, ]
  
  # Sort both datasets by region identifier to ensure alignment
  data_current <- data_current[order(data_current[[region_var]]), ]
  data_lag <- data_lag[order(data_lag[[region_var]]), ]
  
  # Confirm that both periods contain all n regions
  if (nrow(data_current) != n || nrow(data_lag) != n) {
    stop("The data for each time point does not cover all regions")
  }
  
  # Build the stacked response vector y_t (Kn×1): [y_{1,t}; ...; y_{K,t}]
  y <- numeric(k * n)
  for (i in 1:k) {
    y[((i-1)*n + 1):(i*n)] <- data_current[[y_vars[i]]]
  }
  
  # Build the n×K lag response matrix Y_{t-1} for the AR(1) component (Eq. 7)
  y_lag <- matrix(NA, nrow = n, ncol = k)
  for (i in 1:k) {
    y_lag[, i] <- data_lag[[y_vars[i]]]
  }
  
  if (verbose) {
    cat(sprintf("Response vector y dimension: %d x 1\n", length(y)))
    cat(sprintf("Lag-variable matrix y_lag dimension: %d x %d\n", nrow(y_lag), ncol(y_lag)))
  }
  
  # Build the composite block-diagonal design matrix X_t (Eq. 7 of mstr.pdf)
  X <- build_design_matrix(
    data = data_current,
    y_vars = y_vars,
    x_vars = x_vars,
    y_lag = y_lag,
    include_intercept = include_intercept,
    n = n,
    k = k,
    verbose = verbose
  )
  
  # Pre-compute eigenvalues of W for efficient determinant evaluation (Eq. 18)
  if (verbose) cat("Computing eigenvalues of W...\n")
  eigen_W <- compute_eigen_W(W)
  
  # Convert W to spdep listw format (required by spatialreg::lagsarlm)
  if (verbose) cat("Converting W to listw format...\n")
  W_listw <- spdep::mat2listw(W, style = "W")
  
  # Compute p0: total number of exogenous regressors (excluding AR(1) lag terms)
  p0 <- ncol(X) - k^2
  
  if (verbose) {
    cat(sprintf("Number of explanatory variables p0 = %d (intercept included)\n", p0))
    cat(sprintf("Design matrix X dimension: %d x %d\n", nrow(X), ncol(X)))
    cat("=== Data preparation complete ===\n\n")
  }
  
  # Return all prepared objects as a named list
  result <- list(
    y = y,                    # kn x 1 response vector
    X = X,                    # kn x (p0+k^2) design matrix
    W = W,                    # n x n spatial weight matrix
    W_listw = W_listw,        # listw format
    y_lag = y_lag,            # n x k lag-variable matrix
    eigen_W = eigen_W,        # eigenvalues of W
    n = n,                    # number of regions
    k = k,                    # number of variables
    p0 = p0,                  # number of explanatory variables (with intercept)
    data_info = list(
      y_vars = y_vars,
      x_vars = x_vars,
      time_point = time_point,
      time_lag = time_lag,
      region_var = region_var,
      time_var = time_var,
      include_intercept = include_intercept
    )
  )
  
  return(result)
}


#' Construct the block-diagonal composite design matrix X_t (Eq. 7 of mstr.pdf).
#'
#' The design matrix has a block structure where the i-th block (rows for y_i)
#' contains [X*_{i,t} | Y_{t-1}] — exogenous covariates followed by all K lagged
#' responses as AR(1) terms.  Off-diagonal blocks are zero, giving the overall
#' structure:
#'   X_t = blockdiag( [X*_{1,t}, Y_{t-1}], ..., [X*_{K,t}, Y_{t-1}] )
#'
#' @param data            Data frame for the current time period
#' @param y_vars          Character vector of response variable names
#' @param x_vars          List of exogenous covariate names per response
#' @param y_lag           n×K matrix of lagged response values Y_{t-1}
#' @param include_intercept Logical; prepend a column of ones
#' @param n               Number of regions
#' @param k               Number of response variables (K)
#' @param verbose         Logical; print progress
#'
#' @return Kn × (p0 + K²) block-diagonal design matrix
build_design_matrix <- function(
  data,
  y_vars,
  x_vars,
  y_lag,
  include_intercept,
  n,
  k,
  verbose = FALSE
) {
  
  if (verbose) cat("Building block-diagonal design matrix X_t (Eq. 7)...\n")
  
  # Count the number of exogenous columns per response variable (+ intercept)
  x_counts <- sapply(x_vars, length)
  if (include_intercept) {
    x_counts <- x_counts + 1  # Add the intercept
  }
  p0_total <- sum(x_counts)
  
  # Total number of columns = exogenous regressors + K² AR(1) lag coefficients
  total_cols <- p0_total + k^2
  
  # Initialise the Kn × (p0 + K²) design matrix with zeros
  X <- matrix(0, nrow = k*n, ncol = total_cols)
  
  # Track the current column index while filling in blocks
  col_idx <- 1
  
  # Fill in the i-th row block of X_t: [X*_{i,t}]  (exogenous part)
  for (i in 1:k) {
    row_start <- (i-1)*n + 1
    row_end <- i*n
    
    # Construct the sub-matrix X*_{i,t} for response y_{i,t}
    Xi <- NULL
    
    # Prepend an intercept column (column of ones) if requested
    if (include_intercept) {
      Xi <- cbind(Xi, rep(1, n))
    }
    
    # Append the exogenous covariate columns for this response
    if (length(x_vars[[i]]) > 0) {
      for (x_name in x_vars[[i]]) {
        if (!x_name %in% colnames(data)) {
          stop(sprintf("Explanatory variable '%s' does not exist in the data for y%d", x_name, i))
        }
        Xi <- cbind(Xi, data[[x_name]])
      }
    }
    
    # Place X*_{i,t} in the correct block of X
    Xi_cols <- ncol(Xi)
    X[row_start:row_end, col_idx:(col_idx + Xi_cols - 1)] <- Xi
    col_idx <- col_idx + Xi_cols
  }
  
  # Fill in the AR(1) lag block: for each response y_i, append all K lagged
  # response columns Y_{t-1} = [y_{1,t-1}, ..., y_{K,t-1}] (Eq. 7)
  for (i in 1:k) {
    row_start <- (i-1)*n + 1
    row_end <- i*n
    
    for (j in 1:k) {
      X[row_start:row_end, col_idx] <- y_lag[, j]
      col_idx <- col_idx + 1
    }
  }
  
  if (verbose) {
    cat(sprintf("design matrix X: %d×%d\n", nrow(X), ncol(X)))
    cat(sprintf("  - explanatory-variable part: %d columns\n", p0_total))
    cat(sprintf("  - lag-variable part: %d columns\n", k^2))
  }
  
  return(X)
}


################################################################################
# File: spatial_utils.R
# Spatial matrix calculation utilities.
#
# Key formula implemented here (Section 3.4 of mstr.pdf):
#   |I_{Kn} - R⊗W| = ∏_{i=1}^{n} { Σ_{j=0}^{K} (-1)^j e_j(R) ω_i^j }   (Eq. 18)
# where ω_i are eigenvalues of W and e_j(R) are elementary symmetric polynomials
# of R computed from traces r_j = tr(R^j) via Newton's identities (Eqs. 19–23).
################################################################################

#' Pre-compute eigenvalues of the spatial weight matrix W.
#'
#' The eigenvalues ω_i (i = 1, ..., n) are stored once and reused in every
#' evaluation of log|I_{Kn} - R⊗W| via the product formula (Eq. 18), avoiding
#' repeated full determinant computations during optimisation.
#'
#' @param W       n×n row-normalised spatial weight matrix
#' @param verbose Logical; print diagnostics
#' @return Complex-valued eigenvalue vector of length n
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


#' Compute elementary symmetric polynomials e_j(R) (Eqs. 19–23, mstr.pdf).
#'
#' Newton's identity expresses e_j in terms of power-sum traces r_j = tr(R^j):
#'   e_0 = 1
#'   e_1 = r_1
#'   e_2 = (r_1^2 - r_2) / 2
#'   e_3 = (r_1^3 - 3 r_1 r_2 + 2 r_3) / 6
#'   ... (up to e_7 for K ≤ 7)
#'
#' These polynomials appear in the determinant expansion (Eq. 18):
#'   |I_{Kn} - R⊗W| = ∏_i Σ_j (-1)^j e_j(R) ω_i^j
#'
#' @param R       K×K spatial lag coefficient matrix
#' @param verbose Logical; print computed values
#' @return Numeric vector of length K+1: (e_0, e_1, ..., e_K)
compute_elem_sym_poly <- function(R, verbose = FALSE) {
  k <- nrow(R)
  
  # Compute power-sum traces r_j = tr(R^j) for j = 1, ..., K
  r <- numeric(k)
  R_power <- R
  for (j in 1:k) {
    r[j] <- sum(diag(R_power))
    if (j < k) {
      R_power <- R_power %*% R
    }
  }
  
  # Compute elementary symmetric polynomials e_j from traces (Eqs. 19–23)
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
    warning("For k > 7, computation of elementary symmetric polynomials is not implemented")
  }
  
  if (verbose) {
    cat("Elementary symmetric polynomials:\n")
    for (i in 0:min(k, 7)) {
      cat(sprintf("  e%d = %.6f\n", i, e[i+1]))
    }
  }
  
  return(e)
}


#' Compute log|I_{Kn} - R⊗W| efficiently using the eigenvalue expansion (Eq. 18).
#'
#' Instead of filling the full Kn×Kn Kronecker product, uses:
#'   log|I_{Kn} - R⊗W| = Σ_{i=1}^{n} log{ Σ_{j=0}^{K} (-1)^j e_j(R) ω_i^j }
#' where ω_i are eigenvalues of W (pre-computed) and e_j(R) are elementary
#' symmetric polynomials from compute_elem_sym_poly().
#'
#' The same formula applies to |I_{Kn} - Λ⊗W| by substituting Λ for R.
#'
#' @param R        K×K spatial lag coefficient matrix (or Λ for error model)
#' @param eigen_W  Pre-computed eigenvalues of W (length n, possibly complex)
#' @param k        Number of response variables K
#' @param n        Number of regions
#' @param verbose  Logical; print intermediate values
#' @param smooth   If TRUE, apply C1-continuous linear extension near the
#'                 stationarity boundary (use TRUE during optimisation to
#'                 maintain differentiability; FALSE for final likelihood)
#' @return Scalar: log|I_{Kn} - R⊗W|, or -Inf if the matrix is singular
log_det_spatial <- function(R, eigen_W, k, n, verbose = FALSE, smooth = FALSE) {
  
  # Compute elementary symmetric polynomials e_j from traces (Eqs. 19–23)
  e <- compute_elem_sym_poly(R, verbose = FALSE)
  
  # Threshold δ for C1-continuous boundary extension
  delta <- 1e-8
  
  # Evaluate the characteristic polynomial for each eigenvalue ω_i (Eq. 18)
  log_det <- 0
  
  for (i in 1:n) {
    omega_i <- eigen_W[i]
    
    # Evaluate the characteristic polynomial Σ_j (-1)^j e_j(R) ω_i^j  (Eq. 18)
    poly_value <- 0
    omega_power <- 1
    
    for (j in 0:k) {  # Accumulate Σ_{j=0}^{K} (-1)^j e_j(R) ω_i^j
      poly_value <- poly_value + ((-1)^j) * e[j+1] * omega_power
      omega_power <- omega_power * omega_i
    }
    
    # W may have complex eigenvalues; the determinant is real so take Re(.)
    if (is.complex(poly_value)) {
      poly_value <- Re(poly_value)
    }
    
    if (smooth) {
      # Smooth extension: C1-continuous linear extrapolation at the boundary
      if (poly_value > delta) {
        # Normal case: polynomial value is positive, accumulate log
        log_det <- log_det + log(poly_value)
      } else {
        # Near/outside boundary: linear tangent extension at x = δ
        # log(delta) + (poly_value - delta)/delta
        # C1 continuity: value and derivative match at x = δ
        log_det <- log_det + log(delta) + (poly_value - delta) / delta
      }
    } else {
      # Standard behaviour: return -Inf if outside the stationarity region
      if (poly_value <= 0) {
        if (verbose) {
          cat(sprintf("Warning: at i=%d, poly_value = %.6e <= 0\n", i, poly_value))
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


#' Compute the stationarity violation degree based on spectral radius.
#' 
#' The model is stationary ⟺ returned value < 1.
#' Tests min_i det(I_K - M·ω_i) > 0 for all eigenvalues ω_i of W.
#' If this minimum is positive the spatial process is stable.
#' 
#' @param M     K×K spatial parameter matrix (R or Λ)
#' @param eigen_W Eigenvalues of W (pre-computed)
#' @return Violation degree: < 1 → stationary; ≥ 1 → non-stationary
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


#' Efficient computation of (R⊗W)y without forming the Kn×Kn Kronecker product.
#'
#' Using the block-matrix identity (R⊗W)y = vec(W Y R'), where Y is the n×K
#' matrix obtained by reshaping y, the function computes for each block i:
#'   [(R⊗W)y]_i = Σ_j ρ_{ij} · W y_j
#' This avoids the O((Kn)^2) memory cost of the explicit Kronecker product.
#'
#' @param R K×K spatial lag coefficient matrix
#' @param W n×n spatial weight matrix
#' @param y Kn×1 stacked response vector
#' @param k Number of response variables K
#' @param n Number of regions
#' @param verbose Logical; print diagnostics
#' @return Kn×1 vector equal to (R⊗W)y
compute_RW_times_y <- function(R, W, y, k, n, verbose = FALSE) {
  
  # Reshape y (Kn×1) into n×K matrix (columns = K responses)
  y_matrix <- matrix(y, nrow = n, ncol = k)
  
  # Initialise result vector (Kn×1) to zero
  result <- numeric(k * n)
  
  # Compute each i-th block of (R⊗W)y
  for (i in 1:k) {
    block_i <- numeric(n)
    
    for (j in 1:k) {
      # Add contribution ρ_{ij} · (W y_j) to block i
      block_i <- block_i + R[i, j] * (W %*% y_matrix[, j])
    }
    
    # Store computed block in the result vector
    result[((i-1)*n + 1):(i*n)] <- block_i
  }
  
  if (verbose) {
    cat(sprintf("(R(x)W)y computation complete: result dimension = %d x 1\n", length(result)))
  }
  
  return(result)
}


#' Simple stationarity check used as a feasibility guard in Step S4.
#' 
#' Sufficient (not necessary) condition: |ρ_{ij}| < 1 for all i, j.
#' 
#' @param R       K×K spatial parameter matrix
#' @param verbose Logical; print diagnostic messages
#' @return Logical TRUE if max|R_{ij}| < 1 (stationary); FALSE otherwise
check_stationarity <- function(R, verbose = FALSE) {
  
  max_abs <- max(abs(R))
  
  if (verbose) {
    cat(sprintf("Stationarity check: max|rho_ij| = %.6f\n", max_abs))
  }
  
  is_stationary <- max_abs < 1
  
  if (!is_stationary && verbose) {
    cat("Warning: stationarity condition not satisfied (|rho_ij| < 1)\n")
  }
  
  return(is_stationary)
}


################################################################################
# File: test_phase1.R
# Unit tests for Phase 1 functions
################################################################################

#' Run unit tests for Phase 1 functions
#' 
#' @param data_file   Path to panel data CSV
#' @param weight_file Path to spatial weight matrix CSV
test_phase1 <- function(
  data_file = "../simulatedData-data1.csv",
  weight_file = "../spatial_weights.csv"
) {
  
  cat("\n")
  cat(strrep("=", 80), "\n", sep="")
  cat("Phase 1: unit tests start\n")
  cat(strrep("=", 80), "\n", sep="")
  cat("\n")
  
  # Test 1: Data loading and reshaping
  cat("--- Test 1: data loading and reshaping ---\n")
  tryCatch({
    data_list <- prepare_data(
      data_file = data_file,
      weight_file = weight_file,
      y_vars = c("y1", "y2"),
      x_vars = list(
        y1 = c("x_common1", "x_common2", "x_specific1_1"),
        y2 = c("x_common1", "x_common2", "x_specific2_1")
      ),
      time_var = "time",
      time_point = 2,
      region_var = "region",
      include_intercept = TRUE,
      verbose = TRUE
    )
    
    cat("OK: data preparation succeeded\n")
    cat(sprintf("  - y dimension: %d x 1\n", length(data_list$y)))
    cat(sprintf("  - X dimension: %d x %d\n", nrow(data_list$X), ncol(data_list$X)))
    cat(sprintf("  - W dimension: %d x %d\n", nrow(data_list$W), ncol(data_list$W)))
    cat(sprintf("  - number of regions n = %d\n", data_list$n))
    cat(sprintf("  - number of variables k = %d\n", data_list$k))
    cat(sprintf("  - number of explanatory variables p0 = %d\n", data_list$p0))
    cat("\n")
    
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
    return(NULL)
  })
  
  # Test 2: Eigenvalue computation of W
  cat("--- Test 2: eigenvalue computation ---\n")
  tryCatch({
    eigen_vals <- compute_eigen_W(data_list$W, verbose = TRUE)
    cat("OK: eigenvalue computation succeeded\n")
    cat(sprintf("  - eigenvalue range: [%.4f, %.4f]\n", 
                min(Re(eigen_vals)), max(Re(eigen_vals))))
    cat("\n")
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
  })
  
  # Test 3: Elementary symmetric polynomials e_j(R)
  cat("--- Test 3: elementary symmetric polynomial computation ---\n")
  tryCatch({
    # Simple 2×2 test matrix for K = 2
    R_test <- matrix(c(0.3, 0.1, 0.1, 0.4), 2, 2)
    cat("Test matrix R:\n")
    print(R_test)
    
    e <- compute_elem_sym_poly(R_test, verbose = TRUE)
    
    # Verify against analytical values (K=2 case): e1 = tr(R), e2 = det(R)
    e1_expected <- sum(diag(R_test))  # tr(R) = 0.7
    e2_expected <- det(R_test)         # det(R) = 0.11
    
    cat(sprintf("Check: e1 = %.6f (expected: %.6f)\n", e[2], e1_expected))
    cat(sprintf("Check: e2 = %.6f (expected: %.6f)\n", e[3], e2_expected))
    
    if (abs(e[2] - e1_expected) < 1e-10 && abs(e[3] - e2_expected) < 1e-10) {
      cat("OK: elementary symmetric polynomial computation is accurate\n")
    } else {
      cat("FAIL: elementary symmetric polynomial computation has errors\n")
    }
    cat("\n")
    
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
  })
  
  # Test 4: log|I_{Kn} - R⊗W| via eigenvalue expansion (Eq. 18)
  cat("--- Test 4: computation of log|I - R(x)W| ---\n")
  tryCatch({
    R_test <- matrix(c(0.3, 0.1, 0.1, 0.4), 2, 2)
    log_det_val <- log_det_spatial(
      R = R_test, 
      eigen_W = data_list$eigen_W, 
      k = 2, 
      n = data_list$n,
      verbose = TRUE
    )
    
    cat(sprintf("log|I - R⊗W| = %.6f\n", log_det_val))
    
    if (is.finite(log_det_val)) {
      cat("OK: determinant computation succeeded\n")
    } else {
      cat("FAIL: determinant is -Inf\n")
    }
    cat("\n")
    
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
  })
  
  # Test 5: Efficient Kronecker product (R⊗W)y computation
  cat("--- Test 5: computation of (R(x)W)y ---\n")
  tryCatch({
    R_test <- matrix(c(0.3, 0.1, 0.1, 0.4), 2, 2)
    result <- compute_RW_times_y(
      R = R_test,
      W = data_list$W,
      y = data_list$y,
      k = data_list$k,
      n = data_list$n,
      verbose = TRUE
    )
    
    cat(sprintf("Result dimension: %d x 1\n", length(result)))
    cat(sprintf("Result range: [%.4f, %.4f]\n", min(result), max(result)))
    cat("OK: Kronecker-product operation succeeded\n")
    cat("\n")
    
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
  })
  
  # Test 6: Stationarity guard check_stationarity()
  cat("--- Test 6: stationarity check ---\n")
  tryCatch({
    R_stationary <- matrix(c(0.3, 0.1, 0.1, 0.4), 2, 2)
    R_nonstationary <- matrix(c(1.2, 0.1, 0.1, 0.4), 2, 2)
    
    cat("Stationary matrix:\n")
    print(R_stationary)
    is_stat1 <- check_stationarity(R_stationary, verbose = TRUE)
    cat(sprintf("Result: %s\n\n", ifelse(is_stat1, "stationary", "non-stationary")))
    
    cat("Non-stationary matrix:\n")
    print(R_nonstationary)
    is_stat2 <- check_stationarity(R_nonstationary, verbose = TRUE)
    cat(sprintf("Result: %s\n\n", ifelse(is_stat2, "stationary", "non-stationary")))
    
    if (is_stat1 && !is_stat2) {
      cat("OK: stationarity check passed\n")
    } else {
      cat("FAIL: stationarity check has a problem\n")
    }
    cat("\n")
    
  }, error = function(e) {
    cat("FAIL - error:", e$message, "\n\n")
  })
  
  cat(strrep("=", 80), "\n", sep="")
  cat("Phase 1: unit tests complete\n")
  cat(strrep("=", 80), "\n", sep="")
  cat("\n")
  
  return(data_list)
}


################################################################################
# File: main_phase1.R
# Integration test / run script for Phase 1
################################################################################

#' Run the Phase 1 integration test with real data
#' 
#' Tests the full pipeline from data loading to spatial matrix operations.
main_phase1 <- function() {
  
  cat("\n")
  cat(strrep("=", 80), "\n", sep="")
  cat("Phase 1: integration test\n")
  cat(strrep("=", 80), "\n", sep="")
  cat("\n")
  
  # Load required packages (spdep)
  if (!require("spdep", quietly = TRUE)) {
    stop("spdep package is required: install.packages('spdep')")
  }
  
  # Data preparation (prepare_data)
  cat("Step 1: Data preparation\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="") 
  
  data_list <- prepare_data(
    data_file = "simulatedData-data1.csv",
    weight_file = "spatial_weights.csv",
    y_vars = c("y1", "y2"),
    x_vars = list(
      y1 = c("x_common1", "x_common2", "x_specific1_1"),
      y2 = c("x_common1", "x_common2", "x_specific2_1")
    ),
    time_var = "time",
    time_point = 2,
    region_var = "region",
    include_intercept = TRUE,
    verbose = TRUE
  )
  
  cat("\n")
  cat("Step 2: Basic statistics\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="") 
  
  cat("Response vector y statistics:\n")
  cat(sprintf("  Mean: %.4f\n", mean(data_list$y)))
  cat(sprintf("  Std dev: %.4f\n", sd(data_list$y)))
  cat(sprintf("  Min: %.4f\n", min(data_list$y)))
  cat(sprintf("  Max: %.4f\n", max(data_list$y)))
  
  cat("\nDesign matrix X diagnostics:\n")
  cat(sprintf("  Condition number: %.2e\n", kappa(data_list$X)))
  cat(sprintf("  Rank: %d / %d\n", qr(data_list$X)$rank, ncol(data_list$X)))
  
  cat("\n")
  cat("Step 3: Spatial matrix computation test\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="") 
  
  # Define a test R matrix for Step 3 verification
  R_test <- matrix(c(0.3, 0.1, 0.1, 0.4), 2, 2)
  cat("Test R matrix:\n")
  print(R_test)
  cat("\n")
  
  # Compute log|I - R⊗W| using eigenvalue expansion
  log_det <- log_det_spatial(R_test, data_list$eigen_W, data_list$k, data_list$n)
  cat(sprintf("log|I - R⊗W| = %.6f\n", log_det))
  
  # Compute the spatially filtered response (I - R⊗W)y
  RWy <- compute_RW_times_y(R_test, data_list$W, data_list$y, 
                             data_list$k, data_list$n)
  transformed_y <- data_list$y - RWy
  cat(sprintf("Statistics of (I - R(x)W)y:\n"))
  cat(sprintf("  Mean: %.4f\n", mean(transformed_y)))
  cat(sprintf("  Std dev: %.4f\n", sd(transformed_y)))
  
  cat("\n")
  cat("Step 4: Design matrix structure check\n")
  cat(paste(rep("-", 80), collapse=""), "\n", sep="") 
  
  # Inspect block structure of the design matrix X
  n <- data_list$n
  k <- data_list$k
  
  cat("Row block for y1 (first 5 rows × first 10 cols):\n")
  print(data_list$X[1:5, 1:10])
  
  cat("\nRow block for y2 (next 5 rows × first 10 cols):\n")
  print(data_list$X[(n+1):(n+5), 1:10])
  
  cat("\n")
  cat(paste(rep("=", 80), collapse=""), "\n", sep="") 
  cat("Phase 1: integration test complete\n")
  cat("All basic functions worked correctly.\n")
  cat("You can proceed to implementing Phase 2 (estimation logic).\n")
  cat(paste(rep("=", 80), collapse=""), "\n", sep="") 
  cat("\n")
  
  return(data_list)
}


################################################################################
# Usage examples
################################################################################

# Usage examples
# 
# library(spdep)
# test_results <- test_phase1()

# Usage examples
# 
# library(spdep)
# data_list <- main_phase1()

################################################################################
# START OF FILE: phase2_implementation.r
################################################################################

################################################################################
# phase2_implementation.r
#
# Core likelihood and parameter update functions for the MGNS/MSTR model
# described in mstr.pdf (Sections 3.2–3.3).
#
# This file implements the maximum penalized likelihood estimation procedure:
#   S1  Initial estimation: individual GNS models (Section 3.3, Step S1)
#   S2  Inner loop: iterative update of β̂ and Σ̂ given (R, Λ) (Step S2, Eqs. 13, 15)
#   S3  Profile likelihood evaluation (Step S3, Section 3.3)
#   S4  BFGS quasi-Newton update of (R, Λ) (Step S4, Section 3.3)
#   S5  Outer convergence loop (Step S5)
#
# Functions provided:
#   - compute_log_likelihood()      : Log-likelihood computation (Eq. 9)
#   - compute_profile_likelihood()  : Profile log-likelihood (Step S3)
#   - compute_information_criteria() : AIC / BIC computation
#   - update_beta()                 : β̂ update — GLS estimator (Eq. 13)
#   - update_Sigma()                : Σ̂ update (Eq. 15)
#   - compute_residuals()           : Residual vector z = (I-R⊗W)y - Xβ
#   - iterate_beta_sigma()          : Inner loop: iterate β̂ and Σ̂ updates (Step S2)
#
# Dependencies:
#
# Called from:
#   penalized_spatial.r  →  optimize_R_lbfgsb_penalized()
#   build_output.r       →  build_result_object()
#
# Changelog:
#   Phase 2 refactoring: removed legacy functions
# Removed: (legacy functions removed during refactoring)
# Removed: (legacy functions removed during refactoring)
# Removed: (legacy functions removed during refactoring)
################################################################################

################################################################################
# File: likelihood_functions.R
# Log-likelihood and profile log-likelihood computation for the MGNS model.
#
# Log-likelihood (Eq. 9):
#   ℓ(R, Λ, β, Σ) = -(Kn/2)log(2π) + log|I-R⊗W| + log|I-Λ⊗W|
#                   - (n/2)log|Σ| - (1/2)Q
# where Q = z'(Σ^{-1}⊗I_n)z  and  z = (I-Λ⊗W){y - (R⊗W)y - Xβ}  (Eq. 10)
#
# Profile log-likelihood (Step S3):
#   ℓ(R,Λ) = -(Kn/2)log(2πe) + log|I-R⊗W| + log|I-Λ⊗W| - (n/2)log|Σ̂(R,Λ)|
################################################################################

#' Compute the MGNS log-likelihood (Eq. 9 of mstr.pdf).
#'
#' ℓ(R, Λ, β, Σ) = -(Kn/2)log(2π) + log|I-R⊗W| + log|I-Λ⊗W|
#'                  - (n/2)log|Σ| - (1/2) Q
#'
#' where Q = z'(Σ^{-1}⊗I_n)z  and  z = (I-Λ⊗W){y - (R⊗W)y - Xβ}  (Eq. 10).
#'
#' Note: For SLY models (no Λ) pass the identity for the Λ term; for SEM
#' models (no R) pass the identity for the R term.  The SLY implementation
#' here uses z = (I-R⊗W)y - Xβ  (Λ = 0 case, i.e. no spatial error filter).
#'
#' @param R       K×K spatial lag matrix R
#' @param beta    p×1 regression coefficient vector
#' @param Sigma   K×K error covariance matrix Σ
#' @param y       Kn×1 stacked response vector
#' @param X       Kn×p composite design matrix
#' @param W       n×n spatial weight matrix
#' @param eigen_W Pre-computed eigenvalues of W
#' @param k       Number of responses K
#' @param n       Number of regions
#' @param verbose Integer verbosity level (0 = silent, 1 = basic, 2 = detailed)
#' @return Scalar log-likelihood value
compute_log_likelihood <- function(R, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {
  
  # Evaluate log|I_{Kn} - R⊗W| via eigenvalue expansion (Eq. 18)
  log_det_spatial_term <- log_det_spatial(R, eigen_W, k, n, verbose = (verbose >= 2))
  
  if (!is.finite(log_det_spatial_term)) {
    if (verbose >= 1) cat("Warning: log|I - R(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|Σ| using R's built-in determinant function
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Sigma| is not finite\n")
    return(-Inf)
  }
  
  # Compute residual vector z = (I-R⊗W)y - Xβ  (SLY case, Λ = 0)
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta
  
  # Compute Q = z'(Σ^{-1}⊗I_n)z block-by-block (Eq. 10)
  Sigma_inv <- solve(Sigma)
  Q <- 0
  
  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }
  
  # Assemble the full log-likelihood (Eq. 9)
  loglik <- -(k*n/2) * log(2*pi) + log_det_spatial_term - (n/2) * log_det_Sigma - Q/2
  
  if (verbose >= 2) {
    cat(sprintf("  log|I - R⊗W| = %.6f\n", log_det_spatial_term))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }
  
  return(loglik)
}


#' Compute the SLY profile log-likelihood (Step S3 of mstr.pdf, Section 3.3).
#'
#' With β̂ and Σ̂ concentrated out, the profile likelihood is:
#'   ℓ(R) = -(Kn/2)log(2πe) + log|I-R⊗W| - (n/2)log|Σ̂(R)|
#'
#' This is maximised over R in the outer BFGS loop (Step S4).
#'
#' @param R         K×K spatial lag matrix
#' @param beta_hat  Current β̂ estimate
#' @param Sigma_hat Current Σ̂ estimate
#' @param y,X,W     Data objects
#' @param eigen_W   Pre-computed eigenvalues of W
#' @param k,n       Dimensions K, n
#' @param verbose   Integer verbosity
#' @param smooth    If TRUE, use C1 extension at boundary (for BFGS stability)
#' @return Scalar profile log-likelihood
compute_profile_likelihood <- function(R, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {
  
  # Evaluate log|I_{Kn} - R⊗W| via eigenvalue expansion (Eq. 18)
  log_det_spatial_term <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  
  if (!is.finite(log_det_spatial_term)) {
    if (verbose >= 1) cat("Warning: log|I - R(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|Σ̂| (log-determinant of estimated covariance)
  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma_hat)) {
    if (verbose >= 1) cat("Warning: log|hat{Sigma}| is not finite\n")
    return(-Inf)
  }
  
  # Assemble the profile log-likelihood
  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_spatial_term - (n/2) * log_det_Sigma_hat
  
  if (verbose >= 2) {
    cat(sprintf("  profile log-likelihood: %.6f\n", profile_loglik))
  }
  
  return(profile_loglik)
}


#' Compute AIC and BIC from log-likelihood and number of parameters
#' 
#' @param loglik     Maximised log-likelihood
#' @param num_params  Number of estimated parameters
#' @param num_obs    Number of observations (effective sample size)
#' @return list(AIC, BIC)
compute_information_criteria <- function(loglik, num_params, num_obs) {
  AIC <- -2 * loglik + 2 * num_params
  BIC <- -2 * loglik + log(num_obs) * num_params
  
  return(list(AIC = AIC, BIC = BIC))
}


################################################################################
# File: parameter_update.R
# Closed-form parameter update functions for the inner EM-like loop (Step S2).
#
# Given fixed (R, Λ), the MLE equations (Eqs. 11–12) yield closed-form solutions:
#
#   β̂(R, Λ, Σ)  =  [X'(I-Λ⊗W)'(Σ^{-1}⊗I)(I-Λ⊗W)X]^{-1}
#                    × X'(I-Λ⊗W)'(Σ^{-1}⊗I)(I-Λ⊗W)(I-R⊗W)y       (Eq. 13)
#
#   Σ̂_{kl}(R, Λ, β)  =  (1/n) z_k' z_l                             (Eq. 15)
#
# where z = (I-Λ⊗W){y - (R⊗W)y - Xβ} is the transformed residual vector.
################################################################################

#' Update β̂ given fixed (R, Σ) — closed-form GLS estimator (Eq. 13).
#'
#' For SLY (Λ = 0) the estimator simplifies to:
#'   β̂ = [X'(Σ^{-1}⊗I)X]^{-1} · X'(Σ^{-1}⊗I)(I-R⊗W)y
#'
#' This corresponds to Step (c) of S2 in mstr.pdf, Section 3.3.
#' Block-by-block computation avoids forming the full Kn×Kn matrix.
#'
#' @param R         K×K current spatial lag matrix
#' @param Sigma     K×K current error covariance matrix
#' @param y,X,W     Data objects
#' @param k,n       Dimensions
#' @param ridge_eps Ridge regularisation ε added to diagonal if XtSX is singular
#' @param verbose   Integer verbosity
#' @return p×1 updated regression coefficient vector β̂
update_beta <- function(R, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  # Compute Σ^{-1} (invert the K×K error covariance matrix)
  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: error inverting Sigma. Applying regularisation\n")
    solve(Sigma + diag(ridge_eps, k))
  })
  
  # Compute the spatially filtered response (I - R⊗W)y
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  y_transformed <- y - RWy
  
  # Compute X'(Σ^{-1}⊗I_n)X and X'(Σ^{-1}⊗I_n)(I-R⊗W)y block-by-block
  # (Eq. 13): avoid forming the full Kn×Kn Kronecker product
  XtSigmaInvX <- matrix(0, ncol(X), ncol(X))
  XtSigmaInvy <- numeric(ncol(X))
  
  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)
      
      Xi <- X[i_idx, , drop = FALSE]
      Xj <- X[j_idx, , drop = FALSE]
      
      XtSigmaInvX <- XtSigmaInvX + Sigma_inv[i,j] * t(Xi) %*% Xj
      XtSigmaInvy <- XtSigmaInvy + Sigma_inv[i,j] * t(Xi) %*% y_transformed[j_idx]
    }
  }
  
  # Solve for β̂ = [X'(Σ^{-1}⊗I)X]^{-1} X'(Σ^{-1}⊗I)(I-R⊗W)y (Eq. 13)
  beta_hat <- tryCatch({
    solve(XtSigmaInvX, XtSigmaInvy)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: singular matrix in beta estimation. Applying ridge regularisation\n")
    solve(XtSigmaInvX + diag(ridge_eps, ncol(XtSigmaInvX)), XtSigmaInvy)
  })
  
  if (verbose >= 2) {
    cat(sprintf("  beta update complete: ||beta|| = %.6f\n", sqrt(sum(beta_hat^2))))
  }
  
  return(beta_hat)
}


#' Update Σ̂ given fixed (R, β) — closed-form MLE (Eq. 15).
#'
#' The (k, l) element of the updated covariance matrix is:
#'   σ̂_{kl} = (1/n) z_k' z_l
#'
#' where z = (I-R⊗W)y - Xβ is the residual vector (SLY case; Λ = 0).
#' This corresponds to Step (c) of S2 in mstr.pdf, Section 3.3.
#'
#' @param R         K×K current spatial lag matrix
#' @param beta      Current coefficient vector
#' @param y,X,W     Data objects
#' @param k,n       Dimensions
#' @param ridge_eps Small positive value added to diagonal if Σ̂ is not PD
#' @param verbose   Integer verbosity
#' @return K×K symmetric positive-definite covariance matrix Σ̂
update_Sigma <- function(R, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  # Compute residual vector z = (I-R⊗W)y - Xβ  (SLY case, Λ = 0)
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta
  
  # Compute Σ̂_{kl} = (1/n) z_k' z_l for all k, l  (Eq. 15)
  Sigma_new <- matrix(0, k, k)
  
  for (i in 1:k) {
    for (j in i:k) {  # j >= i compute only (using symmetry)
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      
      Sigma_new[i,j] <- sum(z[zi_idx] * z[zj_idx]) / n
      
      if (i != j) {
        Sigma_new[j,i] <- Sigma_new[i,j]  # Enforce symmetry Σ = Σ'
      }
    }
  }
  
  # Ensure positive-definiteness; add ridge ε·I if smallest eigenvalue ≤ 0
  eigen_vals <- eigen(Sigma_new, only.values = TRUE)$values
  min_eigen <- min(Re(eigen_vals))
  
  if (min_eigen <= 0) {
    if (verbose >= 1) {
      cat(sprintf("Warning: Sigma is not positive definite (min eigenvalue = %.2e). Applying regularisation\n", min_eigen))
    }
    Sigma_new <- Sigma_new + diag(ridge_eps, k)
  }
  
  if (verbose >= 2) {
    cat(sprintf("  Sigma update complete: log|Sigma| = %.6f\n", determinant(Sigma_new, logarithm = TRUE)$modulus[1]))
  }
  
  return(Sigma_new)
}


#' Compute the residual vector z = (I-R⊗W)y - Xβ
#' 
#' z = (I-R⊗W)y - Xβ
#' 
#' @param R,beta  Current spatial lag matrix and coefficient vector
#' @param y,X,W  Data objects (response, design matrix, weight matrix)
#' @param k,n  Dimensions (K = responses, n = regions)
#' @return  Numeric vector
compute_residuals <- function(R, beta, y, X, W, k, n) {
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta
  return(z)
}


#' Inner iterative loop: alternate updates of β̂ and Σ̂ with R fixed (Step S2).
#'
#' Implements the inner coordinate-ascent cycle of Section 3.3:
#'   (c)  β̂ ← update_beta(R, Σ)  [Eq. 13]
#'   (c)  Σ̂ ← update_Sigma(R, β̂) [Eq. 15]
#'   (d)  Repeat until convergence.
#'
#' Convergence is assessed by |log|Σ̂_new| - log|Σ̂_old|| < tol.
#'
#' @param R           K×K fixed spatial lag matrix
#' @param beta_init   Initial coefficient vector
#' @param Sigma_init  Initial covariance matrix
#' @param y,X,W       Data objects
#' @param k,n         Dimensions
#' @param max_iter    Maximum number of inner iterations
#' @param tol         Convergence tolerance on log|Σ|
#' @param ridge_eps   Ridge regularisation parameter
#' @param verbose     Integer verbosity (0 = silent)
#' @return list(beta, Sigma, log_det_Sigma_history, converged, iterations)
iterate_beta_sigma <- function(
  R, 
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
  log_det_Sigma_history <- numeric(0)
  
  if (verbose >= 1) {
    cat("  Inner loop start: iterating beta and Sigma updates\n")
  }
  
  for (iter in 1:max_iter) {
    
    # Record current log|Σ| for convergence monitoring
    log_det_Sigma_old <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    # Store history every 10 iterations
    if (iter %% 10 == 1) {
      log_det_Sigma_history <- c(log_det_Sigma_history, log_det_Sigma_old)
    }
    
    # Update β̂ (Eq. 13) — Step (c) of the inner loop
    beta <- update_beta(R, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    # Update Σ̂ (Eq. 15) — Step (c) of the inner loop
    Sigma <- update_Sigma(R, beta, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    # Compute new log|Σ̂| after the update
    log_det_Sigma_new <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    # Convergence check: |log|Σ̂_new| - log|Σ̂_old|| < tol
    change <- abs(log_det_Sigma_new - log_det_Sigma_old)
    
    if (verbose >= 2 && iter %% 10 == 0) {
      cat(sprintf("    Inner iteration %d: log|Sigma| = %.6f, change = %.2e\n", 
                  iter, log_det_Sigma_new, change))
    }
    
    if (change < tol) {
      if (verbose >= 1) {
        cat(sprintf("  Inner loop converged after %d iterations: change = %.2e < tol = %.2e\n", iter, change, tol))
      }
      return(list(
        beta = beta,
        Sigma = Sigma,
        log_det_Sigma_history = log_det_Sigma_history,
        converged = TRUE,
        iterations = iter
      ))
    }
  }
  
  if (verbose >= 1) {
    cat(sprintf("  Inner loop: maximum iterations (%d) reached without convergence\n", max_iter))
  }
  
  return(list(
    beta = beta,
    Sigma = Sigma,
    log_det_Sigma_history = log_det_Sigma_history,
    converged = FALSE,
    iterations = max_iter
  ))
}

################################################################################
# START OF FILE: spatial_core_functions.r
################################################################################

################################################################################
# spatial_core_functions.r
# 
# Shared core functions for multivariate spatial regression (SLY / SEM / SDEM)
# 
# Usage:
#
# Functions added in this module:
#   - prepare_data_extended()      : Data preparation with optional temporal lag
#   - build_design_matrix_extended(): Design matrix with optional temporal lag
#   - compute_vcov_beta()          : Analytic variance-covariance matrix of β̂ (Eq. 16)
#   - compute_hessian_numerical()  : Numerical Hessian computation
# Compute significance tests (SE, z-value, p-value)
# Compute R² statistics
#   - extract_inference_from_spatial_model(): Extract inference from spatialreg objects
# Compute R² statistics
#
################################################################################

cat("spatial_core_functions.r loaded
")

################################################################################
# 1. Data preparation (with time-lag option)
################################################################################

#' Extended data preparation with optional AR(1) temporal lag.
#' 
#' @param data_file   Path to panel data CSV
#' @param weight_file Path to spatial weight matrix CSV
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param verbose Logical; print diagnostic messages
#' 
#' @return list(y, X, W, W_listw, y_lag, eigen_W, n, k, p0, data_info)
prepare_data_extended <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  verbose = FALSE
) {
  
  if (verbose) cat("=== Data preparation start (extended) ===\n")
  
  # Load panel data from CSV
  if (verbose) cat("Loading panel data...\n")
  data <- read.csv(data_file, stringsAsFactors = FALSE)
  
  # Load spatial weight matrix from CSV
  if (verbose) cat("Loading spatial weight matrix...\n")
  W_raw <- read.csv(weight_file, header = TRUE, stringsAsFactors = FALSE)
  
  # Detect whether the first column is a region ID column
  first_col <- W_raw[, 1]
  is_id_column <- !is.numeric(first_col) || all(first_col == 1:nrow(W_raw))
  
  if (is_id_column) {
    W <- as.matrix(W_raw[, -1])
    if (verbose) cat("The first column was treated as the region ID and excluded\n")
  } else {
    W <- as.matrix(W_raw)
  }
  
  dimnames(W) <- NULL
  
  if (verbose) {
    cat(sprintf("Spatial weight matrix dimension: %d x %d\n", nrow(W), ncol(W)))
  }
  
  # Basic validation: check that required columns exist in the data frame
  if (!all(c(time_var, region_var, y_vars) %in% colnames(data))) {
    stop("Required columns not found")
  }
  
  # Determine the number of response variables K and regions n
  k <- length(y_vars)
  regions <- sort(unique(data[[region_var]]))
  n <- length(regions)
  
  if (verbose) {
    cat(sprintf("Number of variables k = %d\n", k))
    cat(sprintf("Number of regions n = %d\n", n))
    cat(sprintf("Include temporal lag: %s\n", ifelse(include_time_lag, "Yes", "No")))
  }
  
  # Verify that W has dimensions n×n matching the number of regions
  if (nrow(W) != n || ncol(W) != n) {
    stop(sprintf("Spatial weight matrix size (%dx%d) does not match the number of regions (%d)", 
                 nrow(W), ncol(W), n))
  }
  
  # Verify that W is row-normalised (each row sums to 1)
  row_sums <- rowSums(W)
  if (!all(abs(row_sums - 1) < 1e-6)) {
    warning("The spatial weight matrix may not be row-normalised")
  }
  
  # Identify available time periods and select the analysis time point
  times <- sort(unique(data[[time_var]]))
  if (is.null(time_point)) {
    time_point <- max(times)
    if (verbose) cat(sprintf("time_point not specified; using the latest time point %d\n", time_point))
  }
  
  if (!time_point %in% times) {
    stop(sprintf("The specified time point %d does not exist in the data", time_point))
  }
  
  # Temporal lag processing
  time_lag <- NULL
  y_lag <- NULL
  
  if (include_time_lag) {
    # Verify that the lag period t-1 exists in the data (needed for AR(1) terms)
    time_lag <- time_point - 1
    if (!time_lag %in% times) {
      stop(sprintf("Lag time point %d does not exist in the data. At least two time points are required.", time_lag))
    }
    
    if (verbose) cat(sprintf("Time point used: t=%d, lag time point: t=%d\n", time_point, time_lag))
    
    # Extract data for the lag period t-1
    data_lag <- data[data[[time_var]] == time_lag, ]
    data_lag <- data_lag[order(data_lag[[region_var]]), ]
    
    if (nrow(data_lag) != n) {
      stop("The lag-time-point data does not cover all regions")
    }
    
    # Build the n×K lag response matrix Y_{t-1} for the AR(1) component (Eq. 7)
    y_lag <- matrix(NA, nrow = n, ncol = k)
    for (i in 1:k) {
      y_lag[, i] <- data_lag[[y_vars[i]]]
    }
  } else {
    if (verbose) cat(sprintf("Time point used: t=%d (no temporal lag)\n", time_point))
  }
  
  # Extract data for the current period t
  data_current <- data[data[[time_var]] == time_point, ]
  data_current <- data_current[order(data_current[[region_var]]), ]
  
  if (nrow(data_current) != n) {
    stop("The current-time-point data does not cover all regions")
  }
  
  # Build the stacked response vector y_t (Kn×1): [y_{1,t}; ...; y_{K,t}]
  y <- numeric(k * n)
  for (i in 1:k) {
    y[((i-1)*n + 1):(i*n)] <- data_current[[y_vars[i]]]
  }
  
  if (verbose) {
    cat(sprintf("Response vector y dimension: %d x 1\n", length(y)))
    if (include_time_lag) {
      cat(sprintf("Lag-variable matrix y_lag dimension: %d x %d\n", nrow(y_lag), ncol(y_lag)))
    }
  }
  
  # Build the composite block-diagonal design matrix X_t (Eq. 7 of mstr.pdf)
  X <- build_design_matrix_extended(
    data = data_current,
    y_vars = y_vars,
    x_vars = x_vars,
    y_lag = y_lag,
    include_intercept = include_intercept,
    include_time_lag = include_time_lag,
    n = n,
    k = k,
    verbose = verbose
  )
  
  # Pre-compute eigenvalues of W for efficient determinant evaluation (Eq. 18)
  if (verbose) cat("Computing eigenvalues of W...\n")
  eigen_W <- compute_eigen_W(W)
  
  # Convert W to spdep listw format (required by spatialreg::lagsarlm)
  if (verbose) cat("Converting W to listw format...\n")
  W_listw <- spdep::mat2listw(W, style = "W")
  
  # Compute p0: total number of exogenous regressors (excluding AR(1) lag terms)
  if (include_time_lag) {
    p0 <- ncol(X) - k^2
  } else {
    p0 <- ncol(X)
  }
  
  if (verbose) {
    cat(sprintf("Number of explanatory variables p0 = %d (intercept included)\n", p0))
    cat(sprintf("Design matrix X dimension: %d x %d\n", nrow(X), ncol(X)))
    cat("=== Data preparation complete ===\n\n")
  }
  
  # Return all prepared objects as a named list
  result <- list(
    y = y,
    X = X,
    W = W,
    W_listw = W_listw,
    y_lag = y_lag,
    eigen_W = eigen_W,
    n = n,
    k = k,
    p0 = p0,
    data_info = list(
      y_vars = y_vars,
      x_vars = x_vars,
      time_point = time_point,
      time_lag = time_lag,
      region_var = region_var,
      time_var = time_var,
      include_intercept = include_intercept,
      include_time_lag = include_time_lag
    )
  )
  
  return(result)
}


#' Build the composite design matrix X_t with optional AR(1) lag columns (Eq. 7).
#' 
#' @param data  Data frame for the current time period
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param y_lag  n×K lag response matrix Y_{t-1} (NULL if include_time_lag = FALSE)
#' @param include_intercept  Logical; include intercept column
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param n  Number of regions
#' @param k  Number of response variables K
#' @param verbose Logical; print diagnostic messages
#' 
#' @return  Numeric matrix
build_design_matrix_extended <- function(
  data,
  y_vars,
  x_vars,
  y_lag,
  include_intercept,
  include_time_lag,
  n,
  k,
  verbose = FALSE
) {
  
  if (verbose) cat("Building block-diagonal design matrix X_t (Eq. 7)...\n")
  
  # Count the number of exogenous columns per response variable (+ intercept)
  x_counts <- sapply(x_vars, length)
  if (include_intercept) {
    x_counts <- x_counts + 1
  }
  p0_total <- sum(x_counts)
  
  # Total number of columns = exogenous regressors + K² AR(1) lag coefficients
  if (include_time_lag) {
    total_cols <- p0_total + k^2
  } else {
    total_cols <- p0_total
  }
  
  # Initialise the Kn × (p0 + K²) design matrix with zeros
  X <- matrix(0, nrow = k*n, ncol = total_cols)
  
  # Track the current column index while filling in blocks
  col_idx <- 1
  
  # Fill in the i-th row block of X_t: [X*_{i,t}]  (exogenous part)
  for (i in 1:k) {
    row_start <- (i-1)*n + 1
    row_end <- i*n
    
    # Construct the sub-matrix X*_{i,t} for response y_{i,t}
    Xi <- NULL
    
    # Prepend an intercept column (column of ones) if requested
    if (include_intercept) {
      Xi <- cbind(Xi, rep(1, n))
    }
    
    # Append the exogenous covariate columns for this response
    if (length(x_vars[[i]]) > 0) {
      for (x_name in x_vars[[i]]) {
        if (!x_name %in% colnames(data)) {
          stop(sprintf("Explanatory variable '%s' does not exist in the data for y%d", x_name, i))
        }
        Xi <- cbind(Xi, data[[x_name]])
      }
    }
    
    # Place X*_{i,t} in the correct block of X
    Xi_cols <- ncol(Xi)
    X[row_start:row_end, col_idx:(col_idx + Xi_cols - 1)] <- Xi
    col_idx <- col_idx + Xi_cols
  }
  
  # Append the K² AR(1) lag columns to X_t (only when include_time_lag = TRUE)
  if (include_time_lag && !is.null(y_lag)) {
    for (i in 1:k) {
      row_start <- (i-1)*n + 1
      row_end <- i*n
      
      for (j in 1:k) {
        X[row_start:row_end, col_idx] <- y_lag[, j]
        col_idx <- col_idx + 1
      }
    }
  }
  
  if (verbose) {
    cat(sprintf("design matrix X: %d×%d\n", nrow(X), ncol(X)))
    cat(sprintf("  - explanatory-variable part: %d columns\n", p0_total))
    if (include_time_lag) {
      cat(sprintf("  - lag-variable part: %d columns\n", k^2))
    }
  }
  
  return(X)
}


################################################################################
# 2. Variance-covariance matrix computation
################################################################################

#' Analytic variance-covariance matrix of β̂ (Eq. 16, mstr.pdf).
#'
#' For SLY/VARX/OLS (Λ = 0):
#'   Ψ = [X'(Σ^{-1}⊗I_n)X]^{-1}
#'
#' For SEM/SDEM (Λ ≠ 0), replace X with (I-Λ⊗W)X:
#'   Ψ = [(I-Λ⊗W)X]'(Σ^{-1}⊗I)[(I-Λ⊗W)X]]^{-1}
#'
#' This matrix is used for significance testing of regression coefficients β.
#'
#' @param X       Kn×p composite design matrix
#' @param Sigma   K×K estimated error covariance matrix Σ̂
#' @param k       Number of responses K
#' @param n       Number of regions
#' @param T_mat   K×K spatial error matrix Λ (NULL for SLY/VARX/OLS)
#' @param W       n×n spatial weight matrix (NULL if T_mat is NULL)
#' @return p×p variance-covariance matrix Ψ of β̂
compute_vcov_beta <- function(X, Sigma, k, n, T_mat = NULL, W = NULL) {
  
  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    warning("Error inverting Sigma. Applying ridge regularisation.")
    solve(Sigma + diag(1e-6, k))
  })
  
  p <- ncol(X)
  
  # For SEM/SDEM: pre-multiply X by (I-Λ⊗W) to form the filtered design matrix.
  # For SLY/VARX/OLS: use X directly (Λ = 0).
  if (!is.null(T_mat) && !is.null(W)) {
    AX <- matrix(0, nrow = k*n, ncol = p)
    for (col in 1:p) {
      AX[, col] <- compute_I_minus_TW_times_v(T_mat, W, X[, col], k, n)
    }
  } else {
    AX <- X
  }
  
  XtSigmaInvX <- matrix(0, p, p)
  
  # Compute [(I-Λ⊗W)X]'(Σ^{-1}⊗I)[(I-Λ⊗W)X] block-by-block
  for (i in 1:k) {
    for (j in 1:k) {
      idx_i <- ((i-1)*n + 1):(i*n)
      idx_j <- ((j-1)*n + 1):(j*n)
      XtSigmaInvX <- XtSigmaInvX + Sigma_inv[i,j] * t(AX[idx_i, , drop=FALSE]) %*% AX[idx_j, , drop=FALSE]
    }
  }
  
  # Invert the information matrix to get Psi = [X'(Sigma^{-1}⊗I)X]^{-1}
  Psi <- tryCatch({
    solve(XtSigmaInvX)
  }, error = function(e) {
    warning("Error inverting X'(Sigma^{-1}(x)I)X. Applying ridge regularisation.")
    solve(XtSigmaInvX + diag(1e-6, p))
  })
  
  return(Psi)
}

#' Compute numerical Hessian of the negative log-likelihood w.r.t. spatial parameters.
#' 
#' Numerical Hessian by central differences.
#' 
#' @param param_vec        Vector of spatial parameters (e.g. vec(R) or vec(Λ))
#' @param negative_loglik_fn Function returning the negative log-likelihood given param_vec
#' @param eps  Step size for central-difference approximation (default: 1e-5)
#' @return  Numerical Hessian matrix (negative second derivative of log-likelihood)
compute_hessian_numerical <- function(param_vec, negative_loglik_fn, eps = 1e-5) {
  
  n_params <- length(param_vec)
  H <- matrix(0, n_params, n_params)
  
  f0 <- negative_loglik_fn(param_vec)
  
  # Enumerate the upper-triangular (i, j) index pairs. Each Hessian entry is
  # computed from 4 independent likelihood evaluations, so the pairs can be
  # evaluated in parallel; the central-difference formula itself is unchanged.
  pair_list <- list()
  for (i in 1:n_params) {
    for (j in i:n_params) {
      pair_list[[length(pair_list) + 1L]] <- c(i, j)
    }
  }
  
  pair_values <- mstr_parallel_lapply(pair_list, function(ij) {
    i <- ij[1L]
    j <- ij[2L]
    ei <- ej <- rep(0, n_params)
    ei[i] <- eps
    ej[j] <- eps
    
    # Central-difference approximation of the (i,j) Hessian entry:
    # [f(x+e_i+e_j) - f(x+e_i-e_j) - f(x-e_i+e_j) + f(x-e_i-e_j)] / (4ε²)
    fpp <- negative_loglik_fn(param_vec + ei + ej)
    fpm <- negative_loglik_fn(param_vec + ei - ej)
    fmp <- negative_loglik_fn(param_vec - ei + ej)
    fmm <- negative_loglik_fn(param_vec - ei - ej)
    
    (fpp - fpm - fmp + fmm) / (4 * eps^2)
  })
  
  for (idx in seq_along(pair_list)) {
    i <- pair_list[[idx]][1L]
    j <- pair_list[[idx]][2L]
    H[i, j] <- pair_values[[idx]]
    H[j, i] <- H[i, j]
  }
  
  return(H)
}


#' Compute variance-covariance matrix from a Hessian matrix.
#' 
#' @param hessian  Hessian matrix (negative second derivative of log-likelihood)
#' @return  Numeric matrix
compute_vcov_from_hessian <- function(hessian, gamma = 0) {
  
if (gamma > 0) {
    # γ > 0: sandwich variance estimator (Appendix A, mstr.pdf)
    # Penalty matrix D = I (spatial params only), so:
    # H_p = H + γI,  Avar(θ̂_γ) = H_p^{-1} H H_p^{-1}
    k1 <- nrow(hessian)
    H_pen <- hessian + gamma * diag(k1)
    
    H_pen_inv <- tryCatch({
      solve(H_pen)
    }, error = function(e) {
      warning("Error inverting (H + gamma I). Applying ridge regularisation.")
      solve(H_pen + diag(1e-6, k1))
    })
    
    vcov <- H_pen_inv %*% hessian %*% H_pen_inv
    
  } else {
    # γ = 0: standard inverse information matrix H^{-1}
    vcov <- tryCatch({
      solve(hessian)
    }, error = function(e) {
      warning("Error inverting the Hessian. Applying ridge regularisation.")
      eigen_H <- eigen(hessian, symmetric = TRUE)
      min_eigenval <- min(eigen_H$values)
      if (min_eigenval <= 0) {
        ridge <- abs(min_eigenval) + 1e-6
        hessian <- hessian + diag(ridge, nrow(hessian))
      }
      solve(hessian)
    })
  }
  
  return(vcov)
}


################################################################################
# 3. Significance testing
################################################################################

#' Compute significance tests (SE, z-value, p-value)
#' 
#' @param estimates   Numeric vector of parameter estimates
#' @param std_errors  Numeric vector of standard errors
#' @param param_names Character vector of parameter names (optional)
#' @return data.frame (parameter, estimate, std_error, z_value, p_value, signif)
compute_inference <- function(estimates, std_errors, param_names = NULL) {
  
  n_params <- length(estimates)
  
  if (is.null(param_names)) {
    param_names <- paste0("param_", 1:n_params)
  }
  
  # Compute z-statistics: z = estimate / SE
  z_values <- estimates / std_errors
  
  # Compute two-sided p-values: p = 2·Φ(-|z|)
  p_values <- 2 * pnorm(-abs(z_values))
  
  # Significance codes: *** ** * . (p < 0.001, 0.01, 0.05, 0.1)
  signif_codes <- cut(p_values,
    breaks = c(0, 0.001, 0.01, 0.05, 0.1, 1),
    labels = c("***", "**", "*", ".", ""),
    right = FALSE,
    include.lowest = TRUE
  )
  
  result <- data.frame(
    parameter = param_names,
    estimate = estimates,
    std_error = std_errors,
    z_value = z_values,
    p_value = p_values,
    signif = as.character(signif_codes),
    stringsAsFactors = FALSE
  )
  
  return(result)
}


#' Return significance code string for a given p-value.
#' 
#' @param p_value  Numeric p-value
#' @return  Character string: '***', '**', '*', '.', or ''
get_signif_code <- function(p_value) {
  if (is.na(p_value)) return("")
  if (p_value < 0.001) return("***")
  if (p_value < 0.01) return("**")
  if (p_value < 0.05) return("*")
  if (p_value < 0.1) return(".")
  return("")
}


################################################################################
# 4. Coefficient of determination
################################################################################

#' Compute R² statistics
#' 
#' @param y          Observed values
#' @param fitted     Fitted (predicted) values
#' @param residuals  Residual vector
#' @param loglik     Model log-likelihood
#' @param num_params  Number of estimated parameters (including spatial parameters)
#' @return list(R2, R2_adj, R2_cor, R2_pseudo)
compute_r_squared <- function(y, fitted, residuals, loglik, num_params) {
  
  n <- length(y)
  y_mean <- mean(y)
  
  # 1. SS-based R² (OLS-compatible, not generally meaningful for spatial models)
  SST <- sum((y - y_mean)^2)
  SSE <- sum(residuals^2)
  R2 <- 1 - SSE/SST
  
  # 2. Adjusted R² (counts all parameters including spatial)
  R2_adj <- 1 - (1 - R2) * (n - 1) / (n - num_params - 1)
  
  # 3. Correlation-based R² = corr(y, ŷ)²
  R2_cor <- cor(y, fitted)^2
  
  # 4. Pseudo R² (likelihood-ratio based, McFadden's definition)
  # Log-likelihood of the null model (intercept only) for pseudo R²
  null_loglik <- -n/2 * (log(2*pi) + 1 + log(var(y)))
  R2_pseudo <- 1 - (loglik / null_loglik)
  
  return(list(
    R2 = R2,
    R2_adj = R2_adj,
    R2_cor = R2_cor,
    R2_pseudo = R2_pseudo
  ))
}

################################################################################
# 4-1b. Compute fitted values ŷ for the multivariate model (used for R²)
################################################################################

#' Compute fitted values ŷ for the multivariate spatial model.
#'
#' Based on Note 5 (mstr.pdf): E[y] = (I - R⊗W)^{-1} Xβ
#'   SLY/SDEM (R≠0): ŷ = (I - R⊗W)⁻¹ Xβ
#'   SEM/VARX/OLS (R=0): ŷ = Xβ
#'
#' @param model_type  Model type string (e.g. 'SLY', 'SEM', 'SDEM')
#' @param R  K×K spatial lag matrix
#' @param beta  Regression coefficient vector β
#' @param X  Kn×p composite design matrix
#' @param W  n×n spatial weight matrix
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @return  Numeric vector
compute_fitted_multivar <- function(model_type, R, beta, X, W, k, n) {
  
  Xbeta <- as.numeric(X %*% beta)
  
  # SLY/SDEM: ŷ = (I - R⊗W)⁻¹ Xβ
  if (model_type %in% c("SLY", "SDEM", "SLY_diagonal", "SDEM_diagonal") && !is.null(R)) {
    IKn_minus_RW <- diag(k * n) - kronecker(R, W)
    fitted_vals <- tryCatch({
      as.numeric(solve(IKn_minus_RW, Xbeta))
    }, error = function(e) {
      warning("Error solving the (I - R(x)W) linear system. Using X*beta as fallback: ", e$message)
      Xbeta
    })
  } else {
    # SEM/VARX/OLS: ŷ = Xβ
    fitted_vals <- Xbeta
  }
  
  return(fitted_vals)
}


################################################################################
# Compute R² statistics
################################################################################

#' Compute averaged pseudo R² for the multivariate model (Eq. 26, mstr.pdf).
#'
#' Eq. 26: R²_{k,pseudo} = corr(y_k, ŷ_k)²  (averaged over K responses)
#'              R̄²_pseudo   = (1/K) Σ R²_pseudo,k
#'
#' ŷ is computed by compute_fitted_multivar().
#'
#' @param y       Kn×1 stacked observed response vector
#' @param fitted  Kn×1 fitted value vector from compute_fitted_multivar()
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @return list(R2_pseudo_individual, R2_pseudo_mean)
compute_r_squared_multivar <- function(y, fitted, k, n) {
  
  R2_pseudo_vec <- numeric(k)
  
  for (i in 1:k) {
    idx <- ((i-1)*n + 1):(i*n)
    yi <- y[idx]
    fi <- fitted[idx]
    
    # pseudo R² = corr(y_k, ŷ_k)²
    R2_pseudo_vec[i] <- cor(yi, fi)^2
  }
  
  list(
    R2_pseudo_individual = R2_pseudo_vec,
    R2_pseudo_mean = mean(R2_pseudo_vec)   # Eq. 26: mean pseudo R^2
  )
}

################################################################################
# 5. Extract inference information from spatialreg model objects
################################################################################

#' Extract significance information from a spatialreg model (lagsarlm/errorsarlm).
#' 
#' Used to obtain detailed inference from the Step S1 individual model estimates.
#' 
#' @param model  A spatialreg model object (lagsarlm or errorsarlm)
#' @param model_type "SLY" (lagsarlm) or "SEM" (errorsarlm)
#' @return list(coefficients, spatial_params, fit)
extract_inference_from_spatial_model <- function(model, model_type = "SLY") {
  
  s <- summary(model)
  
  # Regression coefficient significance
  coef_mat <- s$Coef
  coef_table <- data.frame(
    parameter = rownames(coef_mat),
    estimate = coef_mat[, 1],
    std_error = coef_mat[, 2],
    z_value = coef_mat[, 3],
    p_value = coef_mat[, 4],
    stringsAsFactors = FALSE
  )
  coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
  rownames(coef_table) <- NULL
  
  # Spatial parameter significance test
  if (model_type == "SLY") {
    rho <- model$rho
    rho_se <- s$rho.se
    rho_z <- rho / rho_se
    rho_p <- 2 * pnorm(-abs(rho_z))
    
    spatial_params <- list(
      rho = rho,
      rho_se = rho_se,
      rho_z = rho_z,
      rho_p = rho_p,
      rho_signif = get_signif_code(rho_p),
      lambda = NULL,
      lambda_se = NULL,
      lambda_z = NULL,
      lambda_p = NULL,
      lambda_signif = NULL,
      LR_statistic = if (!is.null(s$LR1)) s$LR1$statistic else NA,
      LR_p_value = if (!is.null(s$LR1)) s$LR1$p.value else NA
    )
    
  } else {  # SEM
    lambda <- model$lambda
    lambda_se <- s$lambda.se
    lambda_z <- lambda / lambda_se
    lambda_p <- 2 * pnorm(-abs(lambda_z))
    
    spatial_params <- list(
      rho = NULL,
      rho_se = NULL,
      rho_z = NULL,
      rho_p = NULL,
      rho_signif = NULL,
      lambda = lambda,
      lambda_se = lambda_se,
      lambda_z = lambda_z,
      lambda_p = lambda_p,
      lambda_signif = get_signif_code(lambda_p),
      LR_statistic = if (!is.null(s$LR1)) s$LR1$statistic else NA,
      LR_p_value = if (!is.null(s$LR1)) s$LR1$p.value else NA
    )
  }
  
  # Goodness-of-fit information
  y <- model$y
  fitted_vals <- fitted(model)
  residuals_vals <- residuals(model)
  loglik <- as.numeric(logLik(model))
  
  # Total parameter count: regression coefficients + spatial parameters + variance
  num_params <- length(coef(model)) + 1 + 1
  
  r2_results <- compute_r_squared(
    y = y,
    fitted = fitted_vals,
    residuals = residuals_vals,
    loglik = loglik,
    num_params = num_params
  )
  
  fit <- list(
    loglik = loglik,
    AIC = AIC(model),
    BIC = BIC(model),
    R2 = r2_results$R2,
    R2_adj = r2_results$R2_adj,
    R2_cor = r2_results$R2_cor,
    R2_pseudo = r2_results$R2_pseudo,
    num_params = num_params,
    num_obs = length(y),
    sigma2 = model$s2
  )
  
  # Residual information
  residuals_info <- list(
    raw = as.numeric(residuals_vals),
    fitted = as.numeric(fitted_vals)
  )
  
  return(list(
    coefficients = coef_table,
    spatial_params = spatial_params,
    fit = fit,
    residuals = residuals_info
  ))
}


################################################################################
# 6. Output the unified coefficient table
################################################################################

#' Format and print the coefficient table with significance codes
#' 
#' @param inference_table  Output of compute_inference()
#' @param title title
#' @param digits  Number of decimal places for display
print_inference_table <- function(inference_table, title = "coefficient table", digits = 4) {
  
  cat(sprintf("\n%s:\n", title))
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("%-20s %10s %10s %10s %10s %5s\n",
              "Parameter", "Estimate", "Std.Error", "z-value", "p-value", ""))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  for (i in 1:nrow(inference_table)) {
    row <- inference_table[i, ]
    cat(sprintf("%-20s %10.4f %10.4f %10.4f %10.4f %5s\n",
                row$parameter,
                row$estimate,
                row$std_error,
                row$z_value,
                row$p_value,
                row$signif))
  }
  
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("significance level: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")
}


################################################################################
# 7. helper function
################################################################################

#' Vectorise the spatial parameters
#' 
#' Convert the R and T matrices into a vector
#' 
#' @param R k x k R matrix (SLY, SDEM case); NULL allowed
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @return  Numeric vector
spatial_params_to_vec <- function(R = NULL, T_mat = NULL) {
  vec <- c()
  if (!is.null(R)) {
    vec <- c(vec, as.vector(t(R)))
  }
  if (!is.null(T_mat)) {
    vec <- c(vec, as.vector(t(T_mat)))
  }
  return(vec)
}


#' Convert a vector into the spatial-parameter matrices
#' 
#' @param vec parameter vector
#' @param k  Number of response variables K
#' @param model_type "SLY", "SEM", "SDEM"
#' @return list(R, T)
vec_to_spatial_params <- function(vec, k, model_type) {
  
  R <- NULL
  T_mat <- NULL
  
  if (model_type == "SLY") {
    R <- matrix(vec, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "SEM") {
    T_mat <- matrix(vec, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "SDEM") {
    n_R <- k^2
    R <- matrix(vec[1:n_R], nrow = k, ncol = k, byrow = TRUE)
    T_mat <- matrix(vec[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
  }
  
  return(list(R = R, T = T_mat))
}


################################################################################
# Usage examples
################################################################################

# cat("
# ================================================================================
# Usage examples
# ================================================================================

# [Data preparation (no temporal lag)]
#   data_list <- prepare_data_extended(
#     data_file = 'data.csv',
#     weight_file = 'W.csv',
#     y_vars = c('y1', 'y2'),
#     x_vars = list(y1 = c('x1', 'x2'), y2 = c('x1', 'x3')),
#     include_time_lag = FALSE,    # No temporal lag
#     verbose = TRUE
#   )

# [Variance-covariance matrix of beta]
#   Psi <- compute_vcov_beta(X, Sigma, k, n)
#   se_beta <- sqrt(diag(Psi))

# [significance test]
#   inference <- compute_inference(beta_hat, se_beta, param_names)
#   print_inference_table(inference)

# [Information extraction from S1]
#   s1_info <- extract_inference_from_spatial_model(model, model_type = 'SLY')

# ================================================================================
# ")

################################################################################
# START OF FILE: initial_estimation_extended.r
################################################################################

################################################################################
# initial_estimation_extended.r
# 
# extended version of S1 initial estimation
# - significance information (SE, z-value, p-value)
# Compute R² statistics
# - detailed verbose output
#
# Usage:
#
################################################################################

cat("Loaded initial_estimation_extended.r\n")

################################################################################
# 1. Extended initial estimation for SLY
################################################################################

#' S1: Initial estimation (individual models per response variable)
#' 
#' Compute R² statistics
#' 
#' @param data_list  Output of prepare_data_extended()
#' @param verbose  Integer verbosity level
#' @return list(R_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sly_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: Individual SAR initial estimation (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Package availability check
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  # matrix to store initial values
  R_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  # estimate each variable individually
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d...\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    # Extract the data for the i-th variable
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    # Construct data frame for spatialreg::lagsarlm / errorsarlm
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    # Clean up column names in the data frame
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    } else {
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[2:(1+length(xi_names))] <- xi_names
      }
    }
    
    # Estimate the model with lagsarlm
    tryCatch({
      if (verbose >= 2) {
        cat("  Calling lagsarlm()...\n")
      }
      
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::lagsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::lagsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      # estimation result extraction
      R_init[i, i] <- model$rho
      Sigma_init[i, i] <- model$s2
      
      # Save the regression coefficients
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      # detailed-information extraction
      s <- summary(model)
      
      # === significance of the spatial parameters ===
      rho <- model$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      rho_signif <- get_signif_code(rho_p)
      
      # Likelihood-ratio test (H_0: spatial parameter = 0)
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # === significance of the regression coefficients ===
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # === goodness-of-fitmeasures ===
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      # number of parameters (regression coefficients + rho + sigma^2)
      num_params <- length(coef(model)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      # Compute R² statistics
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      # Store individual model estimation results
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "SLY",
        
        # Spatial parameters
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = rho_signif,
          lambda = NULL,
          lambda_se = NULL,
          lambda_z = NULL,
          lambda_p = NULL,
          lambda_signif = NULL,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        # Regression coefficients
        coefficients = coef_table,
        
        # Error variance
        sigma2 = model$s2,
        
        # Goodness of fit
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        # Residuals
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      # === verbose output ===
      if (verbose >= 1) {
        # spatial correlation parameter
        cat("  [Spatial lag parameters]\n")
        cat(sprintf("    ρ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, rho, rho_se, rho_z, rho_p, rho_signif))
        
        # Likelihood-ratio test (H_0: spatial parameter = 0)
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "Reject rho=0", "Fail to reject rho=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        # Error variance
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        # Regression coefficients
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        # Goodness of fit
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R^2(correlation) = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    log-likelihood   = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: lagsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed. Aborting.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial R matrix:\n")
    print(round(R_init, 4))
    cat("\nInitial Sigma matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: Initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    R_init = R_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}


################################################################################
# 2. Extended initial estimation for SEM
################################################################################

#' S1: Initial estimation (individual models per response variable)
#' 
#' @param data_list  Output of prepare_data_extended()
#' @param verbose  Integer verbosity level
#' @return list(T_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sem_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: Individual SEM initial estimation (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  T_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d...\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    }
    
    tryCatch({
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::errorsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::errorsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      T_init[i, i] <- model$lambda
      Sigma_init[i, i] <- model$s2
      
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      s <- summary(model)
      
      # Spatial parameter significance test
      lambda <- model$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      lambda_signif <- get_signif_code(lambda_p)
      
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # Regression coefficient significance
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Goodness of fit
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      num_params <- length(coef(model)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "SEM",
        
        spatial_params = list(
          rho = NULL,
          rho_se = NULL,
          rho_z = NULL,
          rho_p = NULL,
          rho_signif = NULL,
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = lambda_signif,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      if (verbose >= 1) {
        cat("  [Spatial lag parameters]\n")
        cat(sprintf("    λ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, lambda, lambda_se, lambda_z, lambda_p, lambda_signif))
        
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "Reject lambda=0", "Fail to reject lambda=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R^2(correlation) = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    log-likelihood   = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: errorsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial T matrix:\n")
    print(round(T_init, 4))
    cat("\nInitial Sigma matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: Initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    T_init = T_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}


################################################################################
# 3. extended initial estimation for SDEM
################################################################################

#' S1: Initial estimation (individual models per response variable)
#' 
#' Uses sacsarlm() from spatialreg
#' 
#' @param data_list  Output of prepare_data_extended()
#' @param verbose  Integer verbosity level
#' @return list(R_init, T_init, Sigma_init, beta_list, individual_estimates)
initial_estimation_sdem_extended <- function(data_list, verbose = 0) {
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: initial estimation by individual SDEM models (extended) ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  
  k <- data_list$k
  n <- data_list$n
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- TRUE
  
  R_init <- matrix(0, k, k)
  T_init <- matrix(0, k, k)
  Sigma_init <- matrix(0, k, k)
  beta_list <- list()
  individual_estimates <- list()
  
  for (i in 1:k) {
    if (verbose >= 1) {
      cat(sprintf("\nEstimating variable y%d...\n", i))
      cat(paste(rep("-", 50), collapse=""), "\n")
    }
    
    yi_data <- extract_data_for_yi_extended(data_list, i, include_time_lag)
    
    if (include_time_lag) {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi,
        yi_data$y_lag
      )
    } else {
      df_i <- data.frame(
        y = yi_data$yi,
        yi_data$Xi
      )
    }
    
    if (data_list$data_info$include_intercept) {
      names(df_i)[2] <- "intercept"
      xi_names <- data_list$data_info$x_vars[[i]]
      if (length(xi_names) > 0) {
        names(df_i)[3:(2+length(xi_names))] <- xi_names
      }
    }
    
    tryCatch({
      if (data_list$data_info$include_intercept) {
        df_i_no_intercept <- df_i[, -2, drop = FALSE]
        
        model <- spatialreg::sacsarlm(
          formula = y ~ .,
          data = df_i_no_intercept,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      } else {
        model <- spatialreg::sacsarlm(
          formula = y ~ . - 1,
          data = df_i,
          listw = data_list$W_listw,
          method = "eigen",
          quiet = (verbose < 2)
        )
      }
      
      R_init[i, i] <- model$rho
      T_init[i, i] <- model$lambda
      Sigma_init[i, i] <- model$s2
      
      beta_coef <- coef(model)
      beta_list[[i]] <- beta_coef
      
      s <- summary(model)
      
      # Spatial parameter significance test (ρ)
      rho <- model$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      rho_signif <- get_signif_code(rho_p)
      
      # Spatial parameter significance test (λ)
      lambda <- model$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      lambda_signif <- get_signif_code(lambda_p)
      
      # Likelihood-ratio test (H_0: spatial parameter = 0)
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # Regression coefficients
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Goodness of fit
      y <- model$y
      fitted_vals <- fitted(model)
      residuals_vals <- residuals(model)
      loglik_i <- as.numeric(logLik(model))
      aic_i <- AIC(model)
      bic_i <- BIC(model)
      
      num_params <- length(coef(model)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      r2_results <- compute_r_squared(
        y = y,
        fitted = fitted_vals,
        residuals = residuals_vals,
        loglik = loglik_i,
        num_params = num_params
      )
      
      individual_estimates[[paste0("y", i)]] <- list(
        model_type = "SDEM",
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = rho_signif,
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = lambda_signif,
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = aic_i,
          BIC = bic_i,
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params,
          num_obs = length(y)
        ),
        
        residuals = list(
          raw = as.numeric(residuals_vals),
          fitted = as.numeric(fitted_vals)
        )
      )
      
      if (verbose >= 1) {
        cat("  [Spatial lag parameters]\n")
        cat(sprintf("    ρ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, rho, rho_se, rho_z, rho_p, rho_signif))
        cat(sprintf("    λ%d%d = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    i, i, lambda, lambda_se, lambda_z, lambda_p, lambda_signif))
        
        if (!is.na(LR_stat)) {
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f\n", LR_stat, LR_p))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ%d%d = %.6f\n", i, i, model$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-18s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²             = %10.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²         = %10.4f\n", r2_results$R2_adj))
        cat(sprintf("    R^2(correlation) = %10.4f\n", r2_results$R2_cor))
        cat(sprintf("    Pseudo R²      = %10.4f\n", r2_results$R2_pseudo))
        cat(sprintf("    log-likelihood   = %10.4f\n", loglik_i))
        cat(sprintf("    AIC            = %10.4f\n", aic_i))
        cat(sprintf("    BIC            = %10.4f\n", bic_i))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: sacsarlm() failed when estimating variable y%d\n", i))
      cat("Error message:", e$message, "\n")
      stop("Initial estimation (S1) failed.")
    })
  }
  
  if (verbose >= 1) {
    cat("\n")
    cat(paste(rep("-", 50), collapse=""), "\n")
    cat("Initial R matrix:\n")
    print(round(R_init, 4))
    cat("\nInitial T matrix:\n")
    print(round(T_init, 4))
    cat("\nInitial Sigma matrix:\n")
    print(round(Sigma_init, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== S1: Initial estimation complete ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n\n")
  }
  
  return(list(
    R_init = R_init,
    T_init = T_init,
    Sigma_init = Sigma_init,
    beta_list = beta_list,
    individual_estimates = individual_estimates
  ))
}


################################################################################
# 4. helper function
################################################################################

#' Extract the i-th yi and Xi from the data (supports the temporal-lag option)
#' 
#' @param data_list  Output list from prepare_data_extended()
#' @param i variable index
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @return list(yi, Xi, y_lag)
extract_data_for_yi_extended <- function(data_list, i, include_time_lag = TRUE) {
  n <- data_list$n
  k <- data_list$k
  
  # Extract yi 
  yi_idx <- ((i-1)*n + 1):(i*n)
  yi <- data_list$y[yi_idx]
  
  # Extract Xi 
  x_vars <- data_list$data_info$x_vars
  include_intercept <- data_list$data_info$include_intercept
  
  x_counts <- sapply(x_vars, length)
  if (include_intercept) {
    x_counts <- x_counts + 1
  }
  
  if (i == 1) {
    x_col_start <- 1
  } else {
    x_col_start <- sum(x_counts[1:(i-1)]) + 1
  }
  x_col_end <- x_col_start + x_counts[i] - 1
  
  x_cols <- x_col_start:x_col_end
  Xi <- data_list$X[yi_idx, x_cols, drop = FALSE]
  
  # lag variables
  y_lag <- NULL
  if (include_time_lag && !is.null(data_list$y_lag)) {
    y_lag <- data_list$y_lag
    colnames(y_lag) <- paste0("y", 1:k, "_lag")
  }
  
  return(list(
    yi = yi,
    Xi = Xi,
    y_lag = y_lag
  ))
}


################################################################################
# Usage examples
################################################################################

# cat("
# ================================================================================
# Usage examples
# ================================================================================

# [SLY use]
#   init_sly <- initial_estimation_sly_extended(data_list, verbose = 1)
#   # Result: init_sly$individual_estimates$y1$spatial_params$rho_p etc.

# [SEM use]
#   init_sem <- initial_estimation_sem_extended(data_list, verbose = 1)
#   # Result: init_sem$individual_estimates$y1$spatial_params$lambda_p etc.

# [SDEM use]
#   init_sdem <- initial_estimation_sdem_extended(data_list, verbose = 1)
#   # Result: init_sdem$individual_estimates$y1$spatial_params etc.

# Compute R² statistics
#   r2 <- init_sly$individual_estimates$y1$fit$R2
#   r2_adj <- init_sly$individual_estimates$y1$fit$R2_adj

# [Regression-coefficient table retrieval]
#   coef_table <- init_sly$individual_estimates$y1$coefficients
#   print(coef_table)

# ================================================================================
# ")

################################################################################
# START OF FILE: multivar_sem_sdem_v2.r
################################################################################

################################################################################
# multivar_sem_sdem_v2.r
#
# Core likelihood / parameter-update functions for the SEM/SDEM models
#
# Functions provided:
#   [Shared utilities]
#     compute_TW_times_v()            : Block computation of (Λ⊗W)v
#     compute_I_minus_TW_times_v()    : (I-T⊗W)v
#     compute_I_minus_TtWt_times_v()  : (I-T'⊗W')v
#
#   [SEM (Spatial Error Model)]
#     compute_log_likelihood_sem()    : SEM log-likelihood (Eq. 9)
#     compute_profile_likelihood_sem(): SEM profile log-likelihood
#     update_beta_sem()               : SEM beta update (Eq. 13, with Λ)
#     update_Sigma_sem()              : SEM Sigma update (Eq. 15, with Λ)
#     compute_residuals_sem()         : SEM residual vector z
#     iterate_beta_sigma_sem()        : SEM inner loop — Step S2
#
#   [SDEM (Spatial Durbin Error Model = MGNS)]
#     compute_log_likelihood_sdem()   : SDEM log-likelihood (Eq. 9)
#     compute_profile_likelihood_sdem(): SDEM profile log-likelihood (Step S3)
#     update_beta_sdem()              : SDEM beta update (Eq. 13, with R and Λ)
#     update_Sigma_sdem()             : SDEM Sigma update (Eq. 15, with R and Λ)
#     compute_residuals_sdem()        : SDEM residual vector z
#     iterate_beta_sigma_sdem()       : SDEM inner loop — Step S2
#
# Dependencies:
#
# Called from:
#   penalized_spatial.r  →  optimize_T_lbfgsb_penalized() (SEM)
#   penalized_spatial.r  →  optimize_RT_lbfgsb_penalized() (SDEM)
#   build_output.r       →  build_result_object()
#
# Changelog:
#   Phase 2B refactoring: removed legacy functions
# Removed: (legacy functions removed during refactoring)
#           (manual gradient ascent replaced by BFGS in penalized_spatial.r)
# Removed: (legacy functions removed during refactoring)
#           (replaced by: initial_estimation_extended.r)
# Removed: (legacy functions removed during refactoring)
#           (replaced by: fit_sem/sdem_penalized in penalized_spatial.r)
# Removed: (legacy functions removed during refactoring)
#           (replaced by: build_output.r)
# Removed: (legacy functions removed during refactoring)
#           (replaced by: spatial_output_functions.r)
################################################################################

################################################################################
# 1. Common utility functions
################################################################################

#' Efficient computation of (T(x)W)v (for SEM)
#' 
#' Implement (T(x)W)v via block computation without explicitly forming the Kronecker product
#' 
#' @param T k x k spatial-correlation matrix of the errors
#' @param W    n×n row-normalised spatial weight matrix
#' @param v  Kn×1 vector
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @param verbose Logical; print diagnostic messages
#' @return  Numeric vector
compute_TW_times_v <- function(T_mat, W, v, k, n, verbose = FALSE) {
  # Reshape v into an n x k matrix
  v_matrix <- matrix(v, nrow = n, ncol = k)
  
  # Initialise result vector (Kn×1) to zero
  result <- numeric(k * n)
  
  # Compute each block
  for (i in 1:k) {
    block_i <- numeric(n)
    
    for (j in 1:k) {
      # Compute λij * W * vj 
      block_i <- block_i + T_mat[i, j] * (W %*% v_matrix[, j])
    }
    
    # Store computed block in the result vector
    result[((i-1)*n + 1):(i*n)] <- block_i
  }
  
  if (verbose) {
    cat(sprintf("(T(x)W)v computation complete: result dimension = %d x 1\n", length(result)))
  }
  
  return(result)
}


#' (I - T⊗W)v  computation
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param W    n×n row-normalised spatial weight matrix
#' @param v  Kn×1 vector
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @return  Numeric vector
compute_I_minus_TW_times_v <- function(T_mat, W, v, k, n) {
  TWv <- compute_TW_times_v(T_mat, W, v, k, n, verbose = FALSE)
  return(v - TWv)
}


#' (I - T'⊗W')v  computation
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param W    n×n row-normalised spatial weight matrix
#' @param v  Kn×1 vector
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @return  Numeric vector
compute_I_minus_TtWt_times_v <- function(T_mat, W, v, k, n) {
  # Use T' (x) W'
  TtWtv <- compute_TW_times_v(t(T_mat), t(W), v, k, n, verbose = FALSE)
  return(v - TtWtv)
}


################################################################################
# 2. SEM (Spatial Error Model)  implementation
################################################################################

#' SEM log-likelihood computation
#' 
#' ℓ(T, β, Σ) = -(kn/2)log(2π) + log|I - T⊗W| - (n/2)log|Σ| - Q/2
#' 
#' residuals: z = (I - T⊗W)(y - Xβ)
#' Q = z'(Σ^{-1}⊗I)z
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta  Regression coefficient vector β
#' @param Sigma  K×K error covariance matrix Σ
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param eigen_W  Pre-computed eigenvalues of W
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param verbose Logical; print diagnostic messages
#' @return  Scalar log-likelihood value
compute_log_likelihood_sem <- function(T_mat, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {
  
  # Evaluate log|I_{Kn} - Λ⊗W| via eigenvalue expansion (Eq. 18)
  log_det_T <- log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE)
  
  if (!is.finite(log_det_T)) {
    if (verbose >= 1) cat("Warning: log|I - T(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|Σ| using R's built-in determinant function
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Sigma| is not finite\n")
    return(-Inf)
  }
  
  # Compute residual z = (I-Λ⊗W)(y - Xβ)  for the SEM model
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Compute Q = z'(Σ^{-1}⊗I_n)z (quadratic form in Eq. 9)
  Sigma_inv <- solve(Sigma)
  Q <- 0
  
  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }
  
  # Assemble the full log-likelihood (Eq. 9)
  loglik <- -(k*n/2) * log(2*pi) + log_det_T - (n/2) * log_det_Sigma - Q/2
  
  if (verbose >= 2) {
    cat(sprintf("  log|I - T⊗W| = %.6f\n", log_det_T))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }
  
  return(loglik)
}


#' SEM profile log-likelihood computation
#' 
#' ℓ(T) = -(kn/2)log(2πe) + log|I - T⊗W| - (n/2)log|Σ̂(T)|
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta_hat βestimated value
#' @param Sigma_hat Σestimated value
#' @param y,X,W  Data objects (response, design matrix, weight matrix)
#' @param eigen_W  Pre-computed eigenvalues of W
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param verbose Logical; print diagnostic messages
#' @return  Scalar profile log-likelihood value
compute_profile_likelihood_sem <- function(T_mat, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {
  
  # Evaluate log|I_{Kn} - Λ⊗W| via eigenvalue expansion (Eq. 18)
  log_det_T <- log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  
  if (!is.finite(log_det_T)) {
    if (verbose >= 1) cat("Warning: log|I - T(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|Σ̂| (log-determinant of estimated covariance)
  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma_hat)) {
    if (verbose >= 1) cat("Warning: log|hat{Sigma}| is not finite\n")
    return(-Inf)
  }
  
  # Assemble the profile log-likelihood
  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_T - (n/2) * log_det_Sigma_hat
  
  if (verbose >= 2) {
    cat(sprintf("  profile log-likelihood: %.6f\n", profile_loglik))
  }
  
  return(profile_loglik)
}


#' SEM beta update
#' 
#' Eqs. 13-14 (general form, here for SEM): hat{beta} = {X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)X}^{-1} X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)y
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param Sigma  K×K error covariance matrix Σ
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param ridge_eps  Ridge regularisation parameter ε
#' @param verbose Logical; print diagnostic messages
#' @return  List of estimated model parameters
update_beta_sem <- function(T_mat, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  # Compute Σ^{-1} (invert the K×K error covariance matrix)
  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: error inverting Sigma. Applying regularisation\n")
    solve(Sigma + diag(ridge_eps, k))
  })
  
  p <- ncol(X)
  
  # Computation of X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)X and X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)y
  # 
  # For efficient computation:
  # 1. Compute (I-T(x)W)X and (I-T(x)W)y
  # 2. Then apply (I-T'(x)W')
  
  # Compute (I-T⊗W)y 
  ITWy <- compute_I_minus_TW_times_v(T_mat, W, y, k, n)
  
  # Computation of (I-T(x)W)X (for each column)
  ITWX <- matrix(0, nrow = k*n, ncol = p)
  for (col in 1:p) {
    ITWX[, col] <- compute_I_minus_TW_times_v(T_mat, W, X[, col], k, n)
  }
  
  # Compute X'(I-T'⊗W')(Σ^{-1}⊗I)(I-T⊗W)X 
  # = (ITWX)'(Σ^{-1}⊗I)(ITWX)
  XtAX <- matrix(0, p, p)
  XtAy <- numeric(p)
  
  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)
      
      XtAX <- XtAX + Sigma_inv[i,j] * t(ITWX[i_idx, , drop = FALSE]) %*% ITWX[j_idx, , drop = FALSE]
      XtAy <- XtAy + Sigma_inv[i,j] * t(ITWX[i_idx, , drop = FALSE]) %*% ITWy[j_idx]
    }
  }
  
  # Estimate beta (with regularisation)
  beta_hat <- tryCatch({
    solve(XtAX, XtAy)
  }, error = function(e) {
    if (verbose >= 1) cat("Warning: singular matrix in beta estimation. Applying ridge regularisation\n")
    solve(XtAX + diag(ridge_eps, p), XtAy)
  })
  
  if (verbose >= 2) {
    cat(sprintf("  beta update complete: ||beta|| = %.6f\n", sqrt(sum(beta_hat^2))))
  }
  
  return(beta_hat)
}


#' SEM Sigma update
#' 
#' Eq. 15: hat{Sigma} = (1/n) * [hat{z}'_i hat{z}_j]
#' where hat{z} = (I-T(x)W)(y - X*beta)
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta  Regression coefficient vector β
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param ridge_eps  Ridge regularisation parameter ε
#' @param verbose Logical; print diagnostic messages
#' @return  List of estimated model parameters
update_Sigma_sem <- function(T_mat, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  # Compute residual z = (I-Λ⊗W)(y - Xβ)  for the SEM model
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Compute the elements of Sigma (block structure)
  Sigma_new <- matrix(0, k, k)
  
  for (i in 1:k) {
    for (j in i:k) {  # j >= i compute only (using symmetry)
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      
      Sigma_new[i,j] <- sum(z[zi_idx] * z[zj_idx]) / n
      
      if (i != j) {
        Sigma_new[j,i] <- Sigma_new[i,j]  # Enforce symmetry Σ = Σ'
      }
    }
  }
  
  # Ensure positive-definiteness; add ridge ε·I if smallest eigenvalue ≤ 0
  eigen_vals <- eigen(Sigma_new, only.values = TRUE)$values
  min_eigen <- min(Re(eigen_vals))
  
  if (min_eigen <= 0) {
    if (verbose >= 1) {
      cat(sprintf("Warning: Sigma is not positive definite (min eigenvalue = %.2e). Applying regularisation\n", min_eigen))
    }
    Sigma_new <- Sigma_new + diag(ridge_eps, k)
  }
  
  if (verbose >= 2) {
    cat(sprintf("  Sigma update complete: log|Sigma| = %.6f\n", determinant(Sigma_new, logarithm = TRUE)$modulus[1]))
  }
  
  return(Sigma_new)
}


#' SEM residual computation
#' 
#' z = (I - T⊗W)(y - Xβ)
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta  Regression coefficient vector β
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param k,n  Dimensions (K = responses, n = regions)
#' @return  Numeric vector
compute_residuals_sem <- function(T_mat, beta, y, X, W, k, n) {
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  return(z)
}


#' Iterative update of beta and Sigma for SEM
#' 
#' With T fixed, update beta and Sigma alternately
#' 
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta_init,Sigma_init  Initial values for β and Σ
#' @param y,X,W  Data objects (response, design matrix, weight matrix)
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param max_iter  Maximum number of iterations
#' @param tol  Convergence tolerance
#' @param ridge_eps  Ridge regularisation parameter ε
#' @param verbose  Integer verbosity level
#' @return  list(beta, Sigma, converged, iterations)
iterate_beta_sigma_sem <- function(
  T_mat, 
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
    cat("  Inner loop start (updating beta and Sigma)\n")
  }
  
  for (iter in 1:max_iter) {
    
    # Record current log|Σ| for convergence monitoring
    log_det_Sigma_old <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    # Update β
    beta <- update_beta_sem(T_mat, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    # Update Σ
    Sigma <- update_Sigma_sem(T_mat, beta, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    # Compute new log|Σ̂| after the update
    log_det_Sigma_new <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    
    # Convergence check: |log|Σ̂_new| - log|Σ̂_old|| < tol
    change <- abs(log_det_Sigma_new - log_det_Sigma_old)
    
    if (verbose >= 2 && iter %% 10 == 0) {
      cat(sprintf("    Inner iteration %d: log|Sigma| = %.6f, change = %.2e\n", 
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
    cat(sprintf("  Inner loop reached maximum iterations (%d)\n", max_iter))
  }
  
  return(list(
    beta = beta,
    Sigma = Sigma,
    converged = FALSE,
    iterations = max_iter
  ))
}

################################################################################
# 3. SDEM (Spatial Durbin Error Model)  implementation
################################################################################

#' SDEM log-likelihood computation
#' 
#' ℓ(R, T, β, Σ) = -(kn/2)log(2π) + log|I-R⊗W| + log|I-T⊗W| - (n/2)log|Σ| - Q/2
#' 
#' residuals: z = (I-T⊗W)(y - (R⊗W)y - Xβ)
#' Q = z'(Σ^{-1}⊗I)z
#' 
#' @param R  K×K spatial lag coefficient matrix
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta  Regression coefficient vector β
#' @param Sigma  K×K error covariance matrix Σ
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param eigen_W  Pre-computed eigenvalues of W
#' @param k,n  Dimensions (K = responses, n = regions)
#' @param verbose Logical; print diagnostic messages
#' @return  Scalar log-likelihood value
compute_log_likelihood_sdem <- function(R, T_mat, beta, Sigma, y, X, W, eigen_W, k, n, verbose = 0) {
  
  # Evaluate log|I_{Kn} - R⊗W| via eigenvalue expansion (Eq. 18)
  log_det_R <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE)
  
  if (!is.finite(log_det_R)) {
    if (verbose >= 1) cat("Warning: log|I - R(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|I_{Kn} - Λ⊗W| via eigenvalue expansion (Eq. 18)
  log_det_T <- log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE)
  
  if (!is.finite(log_det_T)) {
    if (verbose >= 1) cat("Warning: log|I - T(x)W| is not finite\n")
    return(-Inf)
  }
  
  # Evaluate log|Σ| using R's built-in determinant function
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma)) {
    if (verbose >= 1) cat("Warning: log|Sigma| is not finite\n")
    return(-Inf)
  }
  
  # residuals z = (I - T⊗W)(y - (R⊗W)y - Xβ)  computation
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Compute Q = z'(Σ^{-1}⊗I_n)z (quadratic form in Eq. 9)
  Sigma_inv <- solve(Sigma)
  Q <- 0
  
  for (i in 1:k) {
    for (j in 1:k) {
      zi_idx <- ((i-1)*n + 1):(i*n)
      zj_idx <- ((j-1)*n + 1):(j*n)
      Q <- Q + Sigma_inv[i,j] * sum(z[zi_idx] * z[zj_idx])
    }
  }
  
  # Assemble the full log-likelihood (Eq. 9)
  loglik <- -(k*n/2) * log(2*pi) + log_det_R + log_det_T - (n/2) * log_det_Sigma - Q/2
  
  if (verbose >= 2) {
    cat(sprintf("  log|I - R⊗W| = %.6f\n", log_det_R))
    cat(sprintf("  log|I - T⊗W| = %.6f\n", log_det_T))
    cat(sprintf("  log|Σ| = %.6f\n", log_det_Sigma))
    cat(sprintf("  Q = %.6f\n", Q))
    cat(sprintf("  log-likelihood = %.6f\n", loglik))
  }
  
  return(loglik)
}


#' SDEM profile log-likelihood computation
#' 
#' ℓ(R,T) = -(kn/2)log(2πe) + log|I-R⊗W| + log|I-T⊗W| - (n/2)log|Σ̂(R,T)|
compute_profile_likelihood_sdem <- function(R, T_mat, beta_hat, Sigma_hat, y, X, W, eigen_W, k, n, verbose = 0, smooth = FALSE) {
  
  log_det_R <- log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  log_det_T <- log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE, smooth = smooth)
  
  if (!is.finite(log_det_R) || !is.finite(log_det_T)) {
    return(-Inf)
  }
  
  log_det_Sigma_hat <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  
  if (!is.finite(log_det_Sigma_hat)) {
    return(-Inf)
  }
  
  profile_loglik <- -(k*n/2) * log(2*pi*exp(1)) + log_det_R + log_det_T - (n/2) * log_det_Sigma_hat
  
  return(profile_loglik)
}


#' SDEM beta update
#' 
#' Eqs. 13-14 (general form, here for SDEM): hat{beta} = {X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)X}^{-1} 
#'                  X'(I-T'⊗W')(Σ^{-1}⊗I)(I-T⊗W)(I-R⊗W)y
update_beta_sdem <- function(R, T_mat, Sigma, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  Sigma_inv <- tryCatch({
    solve(Sigma)
  }, error = function(e) {
    solve(Sigma + diag(ridge_eps, k))
  })
  
  p <- ncol(X)
  
  # Compute (I-R⊗W)y 
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  IRWy <- y - RWy
  
  # Compute (I-T⊗W)(I-R⊗W)y 
  ITW_IRWy <- compute_I_minus_TW_times_v(T_mat, W, IRWy, k, n)
  
  # Compute (I-T⊗W)X 
  ITWX <- matrix(0, nrow = k*n, ncol = p)
  for (col in 1:p) {
    ITWX[, col] <- compute_I_minus_TW_times_v(T_mat, W, X[, col], k, n)
  }
  
  # X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)X and X'(I-T'(x)W')(Sigma^{-1}(x)I)(I-T(x)W)(I-R(x)W)y
  XtAX <- matrix(0, p, p)
  XtAy <- numeric(p)
  
  for (i in 1:k) {
    for (j in 1:k) {
      i_idx <- ((i-1)*n + 1):(i*n)
      j_idx <- ((j-1)*n + 1):(j*n)
      
      XtAX <- XtAX + Sigma_inv[i,j] * t(ITWX[i_idx, , drop = FALSE]) %*% ITWX[j_idx, , drop = FALSE]
      XtAy <- XtAy + Sigma_inv[i,j] * t(ITWX[i_idx, , drop = FALSE]) %*% ITW_IRWy[j_idx]
    }
  }
  
  beta_hat <- tryCatch({
    solve(XtAX, XtAy)
  }, error = function(e) {
    solve(XtAX + diag(ridge_eps, p), XtAy)
  })
  
  return(beta_hat)
}


#' SDEM Sigma update
#' 
#' z = (I-T⊗W)(y - (R⊗W)y - Xβ)
#' Σ̂ = (1/n) * [z'_i z_j]
update_Sigma_sdem <- function(R, T_mat, beta, y, X, W, k, n, ridge_eps = 1e-6, verbose = 0) {
  
  # residualscomputation
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Compute Σ
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
  
  # positive-definiteness check
  eigen_vals <- eigen(Sigma_new, only.values = TRUE)$values
  min_eigen <- min(Re(eigen_vals))
  
  if (min_eigen <= 0) {
    Sigma_new <- Sigma_new + diag(ridge_eps, k)
  }
  
  return(Sigma_new)
}


#' SDEM residual computation
compute_residuals_sdem <- function(R, T_mat, beta, y, X, W, k, n) {
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  return(z)
}


#' Iterative update of beta and Sigma for SDEM
iterate_beta_sigma_sdem <- function(
  R, T_mat, 
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
    
    beta <- update_beta_sdem(R, T_mat, Sigma, y, X, W, k, n, ridge_eps, verbose = verbose)
    Sigma <- update_Sigma_sdem(R, T_mat, beta, y, X, W, k, n, ridge_eps, verbose = verbose)
    
    log_det_Sigma_new <- determinant(Sigma, logarithm = TRUE)$modulus[1]
    change <- abs(log_det_Sigma_new - log_det_Sigma_old)
    
    if (change < tol) {
      return(list(beta = beta, Sigma = Sigma, converged = TRUE, iterations = iter))
    }
  }
  
  return(list(beta = beta, Sigma = Sigma, converged = FALSE, iterations = max_iter))
}

################################################################################
# START OF FILE: spatial_output_functions.r
################################################################################

################################################################################
# spatial_output_functions.r
# 
# Unified output function for the multivariate spatial regression model (SLY/SEM/SDEM)
# 
# Usage:
#
# Functions provided:
#   - print.multivar_spatial()     : Basic result display
#   - summary.multivar_spatial()   : Detailed summary with significance tests
#   - coef.multivar_spatial()      : Extract coefficients
#   - vcov.multivar_spatial()      : Extract variance-covariance matrix
#   - fitted.multivar_spatial()    : Extract fitted values
#   - residuals.multivar_spatial() : Extract residuals
#   - logLik.multivar_spatial()    : Extract log-likelihood
#   - AIC.multivar_spatial()       : Extract AIC
#   - BIC.multivar_spatial()       : Extract BIC
#   - confint.multivar_spatial()   : Compute confidence intervals
#
################################################################################

cat("spatial_output_functions.r loaded
")

################################################################################
# 1. print method
################################################################################

#' Basic display of a multivar_spatial object
#' 
#' @param x multivar_spatial object
#' @param digits  Number of decimal places for display
#' @param ...  Additional arguments (passed to methods)
#' @export
print.multivar_spatial <- function(x, digits = 4, ...) {
  
  model_name <- switch(x$model_type,
    "SLY" = "Spatial Lag of Y Model (full)",
    "SEM" = "Spatial Error Model (full)",
    "SDEM" = "Spatial Durbin Error Model (full)",
    "SLY_diagonal" = "Spatial Lag of Y Model (diagonal)",
    "SEM_diagonal" = "Spatial Error Model (diagonal)",
    "SDEM_diagonal" = "Spatial Durbin Error Model (diagonal)",
    "VARX" = "VARX Model",
    "OLS_diagonal" = "OLS Model (diagonal)",
    x$model_type  # fallback
  )
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat(sprintf("Multivariate spatial regression model (%s)\n", model_name))
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  # Model specification
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  
  # Presence of temporal lag
  if (!is.null(x$data_info$include_time_lag)) {
    cat(sprintf("  Temporal lag: %s\n", 
                ifelse(x$data_info$include_time_lag, "with", "without")))
  }
  
  # Convergence status
  cat("\n[Convergence status]\n")
  cat(sprintf("  Converged: %s\n", ifelse(x$convergence$converged, "Yes", "No")))
  cat(sprintf("  Iterations: %d\n", x$convergence$iterations))
  cat(sprintf("  Optimisation method: %s\n", x$convergence$method))
  
  # Spatial parameters
  if (x$model_type %in% c("SLY", "SDEM", "SLY_diagonal", "SDEM_diagonal") && !is.null(x$coefficients$R)) {
    cat("\n[Spatial-lag matrix R]\n")
    cat(sprintf("  Diagonal elements: %s\n", 
                paste(sprintf("%.4f", diag(x$coefficients$R)), collapse=", ")))
    
    # Report off-diagonal elements if any are substantial
    R <- x$coefficients$R
    offdiag <- R[row(R) != col(R)]
    if (any(abs(offdiag) > 0.001)) {
      cat(sprintf("  Max off-diagonal absolute value: %.4f\n", max(abs(offdiag))))
    }
  }
  
  if (x$model_type %in% c("SEM", "SDEM", "SEM_diagonal", "SDEM_diagonal") && !is.null(x$coefficients$T)) {
    cat("\n[Spatial-error matrix T]\n")
    cat(sprintf("  Diagonal elements: %s\n", 
                paste(sprintf("%.4f", diag(x$coefficients$T)), collapse=", ")))
    
    T_mat <- x$coefficients$T
    offdiag <- T_mat[row(T_mat) != col(T_mat)]
    if (any(abs(offdiag) > 0.001)) {
      cat(sprintf("  Max off-diagonal absolute value: %.4f\n", max(abs(offdiag))))
    }
  }
  
  # Goodness of fit
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  
  cat("\n")
  cat("Use summary() for details\n")
  cat("\n")
  
  invisible(x)
}


# Aliases for SLY, SEM, and SDEM subclasses
print.multivar_sly <- function(x, ...) print.multivar_spatial(x, ...)
print.multivar_sem <- function(x, ...) print.multivar_spatial(x, ...)
print.multivar_sdem <- function(x, ...) print.multivar_spatial(x, ...)


################################################################################
# 2. summary method
################################################################################

#' Detailed display of a multivar_spatial object
#' 
#' @param object  A 'multivar_spatial' model object
#' @param show_s1 whether to display the S1 initial-estimation results
#' @param show_inference whether to display the significance-test results
#' @param digits  Number of decimal places for display
#' @param ...  Additional arguments (passed to methods)
#' @export
summary.multivar_spatial <- function(object, show_s1 = FALSE, show_inference = TRUE, 
                                      digits = 4, ...) {
  
  x <- object
  
  model_name <- switch(x$model_type,
    "SLY" = "Spatial Lag of Y Model (SLY, full)",
    "SEM" = "Spatial Error Model (SEM, full)",
    "SDEM" = "Spatial Durbin Error Model (SDEM, full)",
    "SLY_diagonal" = "Spatial Lag of Y Model (SLY, diagonal)",
    "SEM_diagonal" = "Spatial Error Model (SEM, diagonal)",
    "SDEM_diagonal" = "Spatial Durbin Error Model (SDEM, diagonal)",
    "VARX" = "VARX Model",
    "OLS_diagonal" = "OLS Model (diagonal)",
    x$model_type  # fallback
  )
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat(sprintf("Multivariate spatial regression model - detailed results\n"))
  cat(sprintf("Model type: %s\n", model_name))
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # ============================================================
  # basic information
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  
  if (!is.null(x$data_info$time_point_used)) {
    cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  }
  
  if (!is.null(x$data_info$include_time_lag)) {
    cat(sprintf("  Temporal lag: %s\n", 
                ifelse(x$data_info$include_time_lag, "with (y_{t-1})", "without")))
  }
  
  # ============================================================
  # Spatial parameters
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[2. Spatial parameters]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  # Rmatrix (SLY, SDEM)
  if (x$model_type %in% c("SLY", "SDEM", "SLY_diagonal", "SDEM_diagonal") && !is.null(x$coefficients$R)) {
    cat("\nSpatial-lag matrix R:\n")
    print(round(x$coefficients$R, digits))
    
    # Show significance information if available
    if (show_inference && !is.null(x$std_errors$R)) {
      cat("\nStandard errors of the R matrix:\n")
      print(round(x$std_errors$R, digits))
    }
  }
  
  # Tmatrix (SEM, SDEM)
  if (x$model_type %in% c("SEM", "SDEM", "SEM_diagonal", "SDEM_diagonal") && !is.null(x$coefficients$T)) {
    cat("\nSpatial-error matrix T:\n")
    print(round(x$coefficients$T, digits))
    
    if (show_inference && !is.null(x$std_errors$T)) {
      cat("\nStandard errors of the T matrix:\n")
      print(round(x$std_errors$T, digits))
    }
  }
  
  # ============================================================
  # Regression coefficients
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[3. Regression coefficients]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  # If a significance-test table is available
  if (show_inference && !is.null(x$inference$coefficients_table)) {
    cat("\nCoefficient table (with significance tests):\n")
    print_coef_table(x$inference$coefficients_table, digits)
  } else {
    # Standard coefficient display
    for (i in 1:x$data_info$k) {
      var_name <- x$data_info$y_vars[i]
      cat(sprintf("\nCoefficients of %s:\n", var_name))
      
      # beta0 — retrieve
      if (!is.null(x$coefficients$beta0)) {
        coef_i <- x$coefficients$beta0[[var_name]]
        if (is.null(coef_i)) {
          coef_i <- x$coefficients$beta0[[paste0("y", i)]]
        }
        
        if (!is.null(coef_i)) {
          for (j in seq_along(coef_i)) {
            cat(sprintf("  %-20s: %10.4f\n", names(coef_i)[j], coef_i[j]))
          }
        }
      }
    }
  }
  
  # Temporal AR(1) coefficient matrix A
  if (!is.null(x$coefficients$alpha) && 
      !is.null(x$data_info$include_time_lag) && 
      x$data_info$include_time_lag) {
    cat("\nTemporal-lag coefficient matrix alpha:\n")
    print(round(x$coefficients$alpha, digits))
  }
  
  # VARX: display of the A matrix
  if (!is.null(x$coefficients$A)) {
    cat("\nTemporal-lag coefficient matrix A:\n")
    print(round(x$coefficients$A, digits))
  }
  
  # ============================================================
  # Error covariance matrix Σ
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[4. Error covariance matrix Sigma]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\n")
  print(round(x$coefficients$Sigma, digits))
  
  # Also display the correlation matrix
  Sigma <- x$coefficients$Sigma
  if (nrow(Sigma) > 1) {
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    cat("\nError correlation matrix:\n")
    print(round(corr_mat, digits))
  }
  
  # ============================================================
  # Goodness of fit
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[5. Goodness-of-fit statistics]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  log-likelihood:        %12.4f\n", x$fit$loglik))
  
  if (!is.null(x$fit$profile_loglik)) {
    cat(sprintf("  profile log-likelihood:%12.4f\n", x$fit$profile_loglik))
  }
  
  cat(sprintf("  AIC:                %12.4f\n", x$fit$AIC))
  cat(sprintf("  BIC:                %12.4f\n", x$fit$BIC))
  
  # Compute R² statistics
  if (!is.null(x$fit$R2)) {
    cat(sprintf("  R²:                 %12.4f\n", x$fit$R2))
  }
  if (!is.null(x$fit$R2_adj)) {
    cat(sprintf("  adjusted R^2:          %12.4f\n", x$fit$R2_adj))
  }
  
  # ============================================================
  # Convergence diagnostics
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[6. Convergence diagnostics]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  Convergence status: %s\n", 
              ifelse(x$convergence$converged, "Converged", "Not converged")))
  cat(sprintf("  Iterations: %d\n", x$convergence$iterations))
  cat(sprintf("  Optimisation method: %s\n", x$convergence$method))
  
  if (!is.null(x$convergence$final_gradient_norm)) {
    cat(sprintf("  Final gradient norm: %.2e\n", x$convergence$final_gradient_norm))
  }
  if (!is.null(x$convergence$final_param_change)) {
    cat(sprintf("  Final parameter change: %.2e\n", x$convergence$final_param_change))
  }
  if (!is.null(x$convergence$final_lik_change)) {
    cat(sprintf("  Final likelihood change: %.2e\n", x$convergence$final_lik_change))
  }
  
  cat(sprintf("  Message: %s\n", x$convergence$message))
  
  # ============================================================
  # Residual summary
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[7. Residual summary]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (!is.null(x$residuals$standardized)) {
    cat("\nDistribution of standardised residuals:\n")
    res_summary <- summary(x$residuals$standardized)
    print(res_summary)
  } else if (!is.null(x$residuals$raw)) {
    cat("\nDistribution of residuals:\n")
    res_summary <- summary(x$residuals$raw)
    print(res_summary)
  }
  
  # ============================================================
  # S1 initial-estimation results (optional)
  # ============================================================
  
  if (show_s1 && !is.null(x$initial_values$individual_estimates)) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("[8. S1 initial-estimation results]\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    for (var_name in names(x$initial_values$individual_estimates)) {
      s1 <- x$initial_values$individual_estimates[[var_name]]
      cat(sprintf("\n%s:\n", var_name))
      
      # Spatial parameters
      if (!is.null(s1$spatial_params$rho)) {
        cat(sprintf("  ρ = %.4f (SE: %.4f, p: %.4f %s)\n",
                    s1$spatial_params$rho,
                    s1$spatial_params$rho_se,
                    s1$spatial_params$rho_p,
                    s1$spatial_params$rho_signif))
      }
      if (!is.null(s1$spatial_params$lambda)) {
        cat(sprintf("  λ = %.4f (SE: %.4f, p: %.4f %s)\n",
                    s1$spatial_params$lambda,
                    s1$spatial_params$lambda_se,
                    s1$spatial_params$lambda_p,
                    s1$spatial_params$lambda_signif))
      }
      
      # Goodness of fit
      if (!is.null(s1$fit)) {
        cat(sprintf("  R² = %.4f, Adj.R² = %.4f\n", 
                    s1$fit$R2, s1$fit$R2_adj))
        cat(sprintf("  log-likelihood = %.4f, AIC = %.4f\n",
                    s1$fit$loglik, s1$fit$AIC))
      }
    }
  }
  
  # ============================================================
  # Execution metadata
  # ============================================================
  
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[Execution information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  if (!is.null(x$execution$time)) {
    cat(sprintf("\n  Execution time: %.2f s\n", as.numeric(x$execution$time)))
  }
  if (!is.null(x$execution$R_version)) {
    cat(sprintf("  R version: %s\n", x$execution$R_version))
  }
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("significance level: *** p<0.001, ** p<0.01, * p<0.05, . p<0.1\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("\n")
  
  invisible(x)
}


# Aliases for SLY, SEM, and SDEM subclasses
summary.multivar_sly <- function(object, ...) summary.multivar_spatial(object, ...)
summary.multivar_sem <- function(object, ...) summary.multivar_spatial(object, ...)
summary.multivar_sdem <- function(object, ...) summary.multivar_spatial(object, ...)


################################################################################
# 3. Formatted output of the coefficient table
################################################################################

#' Format and print the coefficient table with significance codes
#' 
#' @param coef_table data.frame (parameter, estimate, std_error, z_value, p_value, signif)
#' @param digits  Number of decimal places for display
print_coef_table <- function(coef_table, digits = 4) {
  
  cat(paste(rep("-", 75), collapse=""), "\n")
  cat(sprintf("%-25s %10s %10s %10s %10s %5s\n",
              "Parameter", "Estimate", "Std.Error", "z-value", "p-value", ""))
  cat(paste(rep("-", 75), collapse=""), "\n")
  
  for (i in 1:nrow(coef_table)) {
    row <- coef_table[i, ]
    
    # Handle NA
    est <- if (is.na(row$estimate)) "NA" else sprintf("%10.4f", row$estimate)
    se <- if (is.na(row$std_error)) "NA" else sprintf("%10.4f", row$std_error)
    z <- if (is.na(row$z_value)) "NA" else sprintf("%10.4f", row$z_value)
    p <- if (is.na(row$p_value)) "NA" else sprintf("%10.4f", row$p_value)
    sig <- if (is.na(row$signif)) "" else row$signif
    
    cat(sprintf("%-25s %10s %10s %10s %10s %5s\n",
                row$parameter, est, se, z, p, sig))
  }
  
  cat(paste(rep("-", 75), collapse=""), "\n")
}


################################################################################
# 4. Standard method functions
################################################################################

#' Coefficient extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param type extraction type ("all", "R", "T", "beta", "alpha", "Sigma")
#' @param ...  Additional arguments (passed to methods)
#' @export
coef.multivar_spatial <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(object$coefficients$R)
  } else if (type == "T") {
    return(object$coefficients$T)
  } else if (type == "beta") {
    return(object$coefficients$beta)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else if (type == "alpha") {
    return(object$coefficients$alpha)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

coef.multivar_sly <- function(object, ...) coef.multivar_spatial(object, ...)
coef.multivar_sem <- function(object, ...) coef.multivar_spatial(object, ...)
coef.multivar_sdem <- function(object, ...) coef.multivar_spatial(object, ...)


#' variance-covariance matrix extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param type extraction type ("beta", "spatial", "full")
#' @param ...  Additional arguments (passed to methods)
#' @export
vcov.multivar_spatial <- function(object, type = "beta", ...) {
  
  if (is.null(object$vcov)) {
    warning("variance-covariance matrix has not been computed")
    return(NULL)
  }
  
  if (type == "beta") {
    return(object$vcov$beta)
  } else if (type == "spatial") {
    return(object$vcov$spatial)
  } else if (type == "full") {
    return(object$vcov$full)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

vcov.multivar_sly <- function(object, ...) vcov.multivar_spatial(object, ...)
vcov.multivar_sem <- function(object, ...) vcov.multivar_spatial(object, ...)
vcov.multivar_sdem <- function(object, ...) vcov.multivar_spatial(object, ...)


#' Fitted-value extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param ...  Additional arguments (passed to methods)
#' @export
fitted.multivar_spatial <- function(object, ...) {
  
  if (!is.null(object$residuals$fitted)) {
    return(object$residuals$fitted)
  }
  
  # Compute from residuals
  if (!is.null(object$model_data$y) && !is.null(object$residuals$raw)) {
    return(object$model_data$y - object$residuals$raw)
  }
  
  warning("Cannot compute fitted values")
  return(NULL)
}

fitted.multivar_sly <- function(object, ...) fitted.multivar_spatial(object, ...)
fitted.multivar_sem <- function(object, ...) fitted.multivar_spatial(object, ...)
fitted.multivar_sdem <- function(object, ...) fitted.multivar_spatial(object, ...)


#' Residual extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param type residual type ("raw", "standardized")
#' @param ...  Additional arguments (passed to methods)
#' @export
residuals.multivar_spatial <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop(sprintf("Unknown type: %s", type))
  }
}

residuals.multivar_sly <- function(object, ...) residuals.multivar_spatial(object, ...)
residuals.multivar_sem <- function(object, ...) residuals.multivar_spatial(object, ...)
residuals.multivar_sdem <- function(object, ...) residuals.multivar_spatial(object, ...)


#' log-likelihood extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param ...  Additional arguments (passed to methods)
#' @export
logLik.multivar_spatial <- function(object, ...) {
  
  ll <- object$fit$loglik
  attr(ll, "df") <- object$fit$num_params
  attr(ll, "nobs") <- object$fit$num_obs
  class(ll) <- "logLik"
  
  return(ll)
}

logLik.multivar_sly <- function(object, ...) logLik.multivar_spatial(object, ...)
logLik.multivar_sem <- function(object, ...) logLik.multivar_spatial(object, ...)
logLik.multivar_sdem <- function(object, ...) logLik.multivar_spatial(object, ...)


#' AIC extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param ...  Additional arguments (passed to methods)
#' @param k penalty (default: 2)
#' @export
AIC.multivar_spatial <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}

AIC.multivar_sly <- function(object, ...) AIC.multivar_spatial(object, ...)
AIC.multivar_sem <- function(object, ...) AIC.multivar_spatial(object, ...)
AIC.multivar_sdem <- function(object, ...) AIC.multivar_spatial(object, ...)


#' BIC extraction
#' 
#' @param object  A 'multivar_spatial' model object
#' @param ...  Additional arguments (passed to methods)
#' @export
BIC.multivar_spatial <- function(object, ...) {
  return(object$fit$BIC)
}

BIC.multivar_sly <- function(object, ...) BIC.multivar_spatial(object, ...)
BIC.multivar_sem <- function(object, ...) BIC.multivar_spatial(object, ...)
BIC.multivar_sdem <- function(object, ...) BIC.multivar_spatial(object, ...)


################################################################################
# 5. Confidence intervals
################################################################################

#' Compute confidence intervals for the parameters
#' 
#' @param object  A 'multivar_spatial' model object
#' @param parm parameter name (optional)
#' @param level confidence level (default: 0.95)
#' @param type target parameter ("beta", "R", "T", "all")
#' @param ...  Additional arguments (passed to methods)
#' @export
confint.multivar_spatial <- function(object, parm = NULL, level = 0.95, 
                                      type = "beta", ...) {
  
  alpha <- 1 - level
  z_crit <- qnorm(1 - alpha/2)
  
  if (type == "beta" || type == "all") {
    if (is.null(object$std_errors$beta)) {
      warning("SE of beta has not been computed")
      return(NULL)
    }
    
    beta <- object$coefficients$beta
    se_beta <- object$std_errors$beta
    
    ci <- data.frame(
      estimate = beta,
      lower = beta - z_crit * se_beta,
      upper = beta + z_crit * se_beta
    )
    
    rownames(ci) <- names(beta)
    colnames(ci) <- c("Estimate", sprintf("%.1f%%", 100*alpha/2), 
                       sprintf("%.1f%%", 100*(1-alpha/2)))
    
    return(ci)
  }
  
  # Confidence intervals for R and T can be implemented similarly
  warning("Currently, only confidence intervals for beta are supported")
  return(NULL)
}

confint.multivar_sly <- function(object, ...) confint.multivar_spatial(object, ...)
confint.multivar_sem <- function(object, ...) confint.multivar_spatial(object, ...)
confint.multivar_sdem <- function(object, ...) confint.multivar_spatial(object, ...)


################################################################################
# 6. Model comparison
################################################################################

#' Build a comparison table of multiple models
#' 
#' @param ...  Additional arguments
#' @param names model names (optional)
#' @return  data.frame with estimation results
compare_models <- function(..., names = NULL) {
  
  models <- list(...)
  n_models <- length(models)
  
  if (is.null(names)) {
    names <- paste0("Model_", 1:n_models)
  }
  
  comparison <- data.frame(
    Model = names,
    Type = character(n_models),
    k = integer(n_models),
    n = integer(n_models),
    Params = integer(n_models),
    LogLik = numeric(n_models),
    AIC = numeric(n_models),
    BIC = numeric(n_models),
    Converged = logical(n_models),
    stringsAsFactors = FALSE
  )
  
  for (i in 1:n_models) {
    m <- models[[i]]
    comparison$Type[i] <- m$model_type
    comparison$k[i] <- m$data_info$k
    comparison$n[i] <- m$data_info$n
    comparison$Params[i] <- m$fit$num_params
    comparison$LogLik[i] <- m$fit$loglik
    comparison$AIC[i] <- m$fit$AIC
    comparison$BIC[i] <- m$fit$BIC
    comparison$Converged[i] <- m$convergence$converged
  }
  
  # Sort by AIC (ascending)
  comparison <- comparison[order(comparison$AIC), ]
  
  # Add ΔAIC and ΔBIC columns
  comparison$Delta_AIC <- comparison$AIC - min(comparison$AIC)
  comparison$Delta_BIC <- comparison$BIC - min(comparison$BIC)
  
  class(comparison) <- c("model_comparison", "data.frame")
  return(comparison)
}


#' Display the model-comparison table
#' 
#' @param x model_comparison object
#' @param ...  Additional arguments (passed to methods)
print.model_comparison <- function(x, ...) {
  
  cat("\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  cat("Model comparison\n")
  cat(paste(rep("=", 80), collapse=""), "\n\n")
  
  # Round numeric values for display
  x$LogLik <- round(x$LogLik, 2)
  x$AIC <- round(x$AIC, 2)
  x$BIC <- round(x$BIC, 2)
  x$Delta_AIC <- round(x$Delta_AIC, 2)
  x$Delta_BIC <- round(x$Delta_BIC, 2)
  
  print(as.data.frame(x), row.names = FALSE)
  
  cat("\nBest model (AIC criterion):", x$Model[1], "\n")
  cat("Best model (BIC criterion):", x$Model[which.min(x$BIC)], "\n")
  cat("\n")
  
  invisible(x)
}


################################################################################
# 7. Post-estimation inference
################################################################################

#' Add a significance test to the estimation result
#' 
#' Call only once after estimation completes; compute vcov, std_errors, inference
#' 
#' @param object  A 'multivar_spatial' model object
#' @param compute_spatial_se whether to compute SE of the spatial parameters (uses numerical Hessian)
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
add_inference <- function(object, compute_spatial_se = TRUE, gamma = 0, verbose = TRUE) {
  
  # Diagonal/VARX/OLS models already have std_errors from individual estimation
  if (object$model_type %in% c("SLY_diagonal", "SEM_diagonal", "SDEM_diagonal",
                                "VARX", "OLS_diagonal")) {
    if (verbose) {
      cat(sprintf("\n%s: Use std_errors from individual estimation (skip add_inference)\n", object$model_type))
    }
    return(object)
  }
  
  if (verbose) {
    cat("\nComputing significance tests...\n")
  }
  
  k <- object$data_info$k
  n <- object$data_info$n
  X <- object$model_data$X
  Sigma <- object$coefficients$Sigma
  
  # ============================================================
  # Variance-covariance matrix of beta (analytic)
  # ============================================================
  
  if (verbose) cat("  Computing the variance-covariance matrix of beta...\n")
  
  # For SEM/SDEM pass Λ and W to account for spatial error filtering (Eq. 16)
  T_for_vcov <- if (object$model_type %in% c("SEM", "SDEM")) object$coefficients$T else NULL
  W_for_vcov <- if (!is.null(T_for_vcov)) object$model_data$W else NULL
  
  Psi <- tryCatch({
    compute_vcov_beta(X, Sigma, k, n, T_mat = T_for_vcov, W = W_for_vcov)
  }, error = function(e) {
    warning("Error computing the variance-covariance matrix of beta: ", e$message)
    NULL
  })
  
  se_beta <- NULL
  if (!is.null(Psi)) {
    se_beta <- sqrt(diag(Psi))
  }
  
  # ============================================================
  # SE of the spatial parameters (numerical Hessian)
  # ============================================================
  
  se_R <- NULL
  se_T <- NULL
  hessian_result <- NULL
  vcov_spatial <- NULL
  
  if (compute_spatial_se) {
    if (verbose) cat("  Computing the SE of the spatial parameters (numerical Hessian)...\n")
    
    # Handling depending on the model type
    if (object$model_type == "SLY" && !is.null(object$coefficients$R)) {
      hessian_result <- compute_hessian_for_R(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          se_R <- matrix(sqrt(diag(vcov_spatial)), nrow = k, ncol = k, byrow = TRUE)
        }
      }
      
    } else if (object$model_type == "SEM" && !is.null(object$coefficients$T)) {
      hessian_result <- compute_hessian_for_T(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          se_T <- matrix(sqrt(diag(vcov_spatial)), nrow = k, ncol = k, byrow = TRUE)
        }
      }
      
    } else if (object$model_type == "SDEM") {
      hessian_result <- compute_hessian_for_RT(object, verbose = FALSE)
      
      if (!is.null(hessian_result)) {
        vcov_spatial <- compute_vcov_from_hessian(hessian_result, gamma = gamma)
        if (!is.null(vcov_spatial)) {
          n_R <- k^2
         se_vec <- sqrt(diag(vcov_spatial))
         se_R <- matrix(se_vec[1:n_R], nrow = k, ncol = k, byrow = TRUE)
         se_T <- matrix(se_vec[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
        }
      }
    }
  }
  
  # ============================================================
  # Build the coefficient table
  # ============================================================
  
  if (verbose) cat("  Building the coefficient table...\n")
  
  coef_table <- create_coefficient_table(object, se_beta, se_R, se_T)
  
  # ============================================================
  # Append to the object
  # ============================================================
  
  object$vcov <- list(
    beta = Psi,
    spatial = vcov_spatial,
    full = NULL  # full version to be implemented later
  )
  
  object$hessian <- list(
    matrix = hessian_result,
    method = "numerical"
  )
  
  object$std_errors <- list(
    beta = se_beta,
    R = se_R,
    T = se_T
  )
  
  object$inference <- list(
    coefficients_table = coef_table
  )
  
  if (verbose) {
    cat("  Done\n\n")
  }
  
  return(object)
}


#' Hessian computation for the R matrix (for SLY)
#' 
#' @param object  A 'multivar_spatial' model object
#' @param eps  Step size for finite-difference approximation
#' @param verbose Logical; print diagnostic messages
#' @return  Numerical Hessian matrix (negative second derivative of log-likelihood)
compute_hessian_for_R <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  R <- object$coefficients$R
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  # Rmatrix: vectorise
  R_vec <- as.vector(t(R))
  n_params <- length(R_vec)
  
  # Negative log-likelihood function (objective for BFGS minimisation)
  neg_loglik <- function(R_vec_temp) {
    R_temp <- matrix(R_vec_temp, nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      compute_log_likelihood(
        R = R_temp, beta = beta, Sigma = Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W,
        k = k, n = n, verbose = 0
      )
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  # numerical Hessian
  H <- compute_hessian_numerical(R_vec, neg_loglik, eps = eps)
  
  return(H)
}


#' Hessian computation for the T matrix (for SEM)
#' 
#' @param object  A 'multivar_spatial' model object
#' @param eps  Step size for finite-difference approximation
#' @param verbose Logical; print diagnostic messages
#' @return  Numerical Hessian matrix (negative second derivative of log-likelihood)
compute_hessian_for_T <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  T_mat <- object$coefficients$T
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  # Tmatrix: vectorise
  T_vec <- as.vector(t(T_mat))
  n_params <- length(T_vec)
  
  # Negative log-likelihood function (objective for BFGS minimisation)
  neg_loglik <- function(T_vec_temp) {
    T_temp <- matrix(T_vec_temp, nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      # compute_log_likelihood_sem is required
      if (exists("compute_log_likelihood_sem")) {
        compute_log_likelihood_sem(
          T_mat = T_temp, beta = beta, Sigma = Sigma,
          y = y, X = X, W = W, eigen_W = eigen_W,
          k = k, n = n, verbose = 0
        )
      } else {
        -Inf
      }
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  H <- compute_hessian_numerical(T_vec, neg_loglik, eps = eps)
  
  return(H)
}


#' Hessian computation for the R and T matrices (for SDEM)
#' 
#' @param object  A 'multivar_spatial' model object
#' @param eps  Step size for finite-difference approximation
#' @param verbose Logical; print diagnostic messages
#' @return  Numerical Hessian matrix (negative second derivative of log-likelihood)
compute_hessian_for_RT <- function(object, eps = 1e-5, verbose = FALSE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  R <- object$coefficients$R
  T_mat <- object$coefficients$T
  beta <- object$coefficients$beta
  Sigma <- object$coefficients$Sigma
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  
  # R, Tmatrix: vectorise
  R_vec <- as.vector(t(R))
  T_vec <- as.vector(t(T_mat))
  param_vec <- c(R_vec, T_vec)
  
  # Negative log-likelihood function (objective for BFGS minimisation)
  neg_loglik <- function(param_temp) {
    n_R <- k^2
    R_temp <- matrix(param_temp[1:n_R], nrow = k, ncol = k, byrow = TRUE)
    T_temp <- matrix(param_temp[(n_R+1):(2*n_R)], nrow = k, ncol = k, byrow = TRUE)
    
    ll <- tryCatch({
      if (exists("compute_log_likelihood_sdem")) {
        compute_log_likelihood_sdem(
          R = R_temp, T_mat = T_temp, beta = beta, Sigma = Sigma,
          y = y, X = X, W = W, eigen_W = eigen_W,
          k = k, n = n, verbose = 0
        )
      } else {
        -Inf
      }
    }, error = function(e) -Inf)
    
    if (!is.finite(ll)) return(1e10)
    return(-ll)
  }
  
  H <- compute_hessian_numerical(param_vec, neg_loglik, eps = eps)
  
  return(H)
}


#' Build the coefficient table
#' 
#' @param object  A 'multivar_spatial' model object
#' @param se_beta SE of beta
#' @param se_R SE matrix of R
#' @param se_T SE matrix of T
#' @return data.frame
create_coefficient_table <- function(object, se_beta, se_R, se_T) {
  
  k <- object$data_info$k
  rows <- list()
  
  # Spatial parameters R (SLY, SDEM)
  if (!is.null(object$coefficients$R)) {
    R <- object$coefficients$R
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("R[%d,%d]", i, j)
        est <- R[i, j]
        se <- if (!is.null(se_R)) se_R[i, j] else NA
        z <- if (!is.na(se) && se > 0) est / se else NA
        p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
        sig <- get_signif_code(p)
        
        rows[[length(rows) + 1]] <- data.frame(
          parameter = param_name,
          estimate = est,
          std_error = se,
          z_value = z,
          p_value = p,
          signif = sig,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Spatial parameters T (SEM, SDEM)
  if (!is.null(object$coefficients$T)) {
    T_mat <- object$coefficients$T
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("T[%d,%d]", i, j)
        est <- T_mat[i, j]
        se <- if (!is.null(se_T)) se_T[i, j] else NA
        z <- if (!is.na(se) && se > 0) est / se else NA
        p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
        sig <- get_signif_code(p)
        
        rows[[length(rows) + 1]] <- data.frame(
          parameter = param_name,
          estimate = est,
          std_error = se,
          z_value = z,
          p_value = p,
          signif = sig,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # regression coefficients beta
  beta <- object$coefficients$beta
  if (!is.null(beta) && !is.null(se_beta)) {
    # coefficient names — retrieve
    beta_names <- names(beta)
    if (is.null(beta_names)) {
      beta_names <- paste0("beta[", 1:length(beta), "]")
    }
    
    for (i in seq_along(beta)) {
      est <- beta[i]
      se <- se_beta[i]
      z <- if (!is.na(se) && se > 0) est / se else NA
      p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
      sig <- get_signif_code(p)
      
      rows[[length(rows) + 1]] <- data.frame(
        parameter = beta_names[i],
        estimate = est,
        std_error = se,
        z_value = z,
        p_value = p,
        signif = sig,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(rows) == 0) {
    return(NULL)
  }
  
  coef_table <- do.call(rbind, rows)
  rownames(coef_table) <- NULL
  
  return(coef_table)
}


################################################################################
# Usage examples
################################################################################

# cat("
# ================================================================================
# Usage examples
# ================================================================================

# [Basic display]
#   print(result)

# [Detailed display]
#   summary(result)

# [Also display the S1 results]
#   summary(result, show_s1 = TRUE)

# [Display without significance test]
#   summary(result, show_inference = FALSE)

# [Additional computation of the significance test]
#   result <- add_inference(result)
#   summary(result)

# [Coefficient extraction]
#   coef(result)                  # All coefficients
#   coef(result, type = 'R')      # Spatial lag matrix only
#   coef(result, type = 'beta0')  # Regression coefficients only

# [variance-covariance matrix]
#   vcov(result)                  # Variance-covariance matrix of beta
#   vcov(result, type = 'spatial')  # Variance-covariance matrix of spatial parameters

# [Goodness-of-fit measures]
#   logLik(result)
#   AIC(result)
#   BIC(result)

# [Residuals]
#   residuals(result)
#   residuals(result, type = 'standardized')

# [fitted values]
#   fitted(result)

# [confidence interval]
# Confidence intervals
# Confidence intervals

# [Model comparison]
#   comp <- compare_models(sly_result, sem_result, sdem_result,
#                          names = c('SLY', 'SEM', 'SDEM'))
#   print(comp)

# ================================================================================
# ")


################################################################################
# START OF FILE: build_output.r
################################################################################

################################################################################
# build_output.r
#
# Helper functions to unify the output format of all fit functions
#
# Purpose:
#   For the object returned by the fit function of each model (SLY/SEM/SDEM/diagonal/VARX/OLS),
#   unify the structure so it is compatible with the S3 methods
#   (add_inference, coef, logLik, residuals, etc.).
#
# Usage:
#
#   Call build_result_object() in the return section of each fit function.
#
# Dependencies:
#   phase2_implementation.r   (compute_information_criteria, compute_residuals)
#   multivar_sem_sdem_v2.r    (compute_residuals_sem, compute_residuals_sdem)
#   spatial_core_functions.r  (compute_fitted_multivar, compute_r_squared_multivar)
#
################################################################################

cat("Loaded build_output.r (unified output-structure helpers)\n")


################################################################################
# 0. helper function: structure extraction from the beta vector
################################################################################

#' Extract the temporal-lag matrix A from the beta vector
#'
#' The last k^2 entries of beta are the elements of the A matrix (row-major)
#' A[i,j] = beta[p0_total + (i-1)*k + j]
#'
#' @param beta  Regression coefficient vector β
#' @param data_list  Output of prepare_data_extended()
#' @param k  Number of response variables K
#' @return  Numeric matrix
extract_alpha <- function(beta, data_list, k) {
  x_vars <- data_list$data_info$x_vars
  include_intercept <- data_list$data_info$include_intercept
  
  # p0_total: number of columns in the explanatory-variable part (intercept + x_vars, for all variables)
  p0_total <- 0
  for (i in 1:k) {
    p0_total <- p0_total + length(x_vars[[i]])
    if (include_intercept) p0_total <- p0_total + 1
  }
  
  # Extract the last k^2 entries of beta
  alpha_idx <- (p0_total + 1):(p0_total + k^2)
  if (max(alpha_idx) > length(beta)) {
    stop("beta vector too short: cannot extract the temporal-lag part")
  }
  
  alpha_vec <- beta[alpha_idx]
  A <- matrix(alpha_vec, nrow = k, ncol = k, byrow = TRUE)
  
  # Set the matrix names
  y_vars <- data_list$data_info$y_vars
  rownames(A) <- y_vars
  colnames(A) <- paste0(y_vars, "_lag")
  
  return(A)
}


#' Structure the beta vector into a per-variable named list
#'
#' @param beta  Regression coefficient vector β
#' @param data_list  Output of prepare_data_extended()
#' @param k  Number of response variables K
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
    
    # Prepend an intercept column (column of ones) if requested
    if (include_intercept) {
      coef_names <- c(coef_names, "(Intercept)")
      coef_vals <- c(coef_vals, beta[idx])
      idx <- idx + 1
    }
    
    # Append the exogenous covariate columns for this response
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

################################################################################
# 1. Main function: build_result_object
################################################################################

#' Generate an object with the unified output structure
#'
#' @param model_type  Model type string (e.g. 'SLY', 'SEM', 'SDEM')
#'   "SLY", "SEM", "SDEM",
#'   "SLY_diagonal", "SEM_diagonal", "SDEM_diagonal",
#'   "VARX", "OLS_diagonal"
#' @param R  K×K spatial lag coefficient matrix
#' @param T_mat  K×K spatial error parameter matrix Λ
#' @param beta  Regression coefficient vector β
#' @param Sigma  K×K error covariance matrix Σ
#' @param loglik unpenalised log-likelihood
#' @param num_params  Number of estimated parameters
#' @param converged convergence flag
#' @param method estimation-method name
#' @param iterations number of iterations (NA allowed)
#' @param data_list  Output of prepare_data_extended()
#' @param gamma  Penalty strength γ ≥ 0
#' @param penalty_value penalty amount
#' @param penalized_loglik penalised log-likelihood
#' @param execution_time execution time (difftime)
#' @param individual_models list of individual-estimation results for the diagonal model (NULL allowed)
#' @param individual_estimates list of individual S1 initial-estimation results (NULL allowed)
#' @param std_errors_R SE matrix of the R matrix (NULL allowed)
#' @param std_errors_T SE matrix of the T matrix (NULL allowed)
#' @param beta0 per-variable beta (list form; if NULL, built automatically from beta)
#' @param alpha temporal-lag matrix A (if NULL, extracted automatically from beta)
#' @param residuals_raw residual vector (computed automatically if NULL)
#' @param residuals_std standardised-residual vector (computed automatically if NULL)
#' @return  S3 model object of class multivar_spatial
#'
build_result_object <- function(
  model_type,
  R = NULL,
  T_mat = NULL,
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
  std_errors_T = NULL,
  beta0 = NULL,
  alpha = NULL,
  residuals_raw = NULL,
  residuals_std = NULL
) {
  
  k <- data_list$k
  n <- data_list$n
  
  include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE
  
  # ==================================================================
  # Information criteria (AIC / BIC using unpenalised log-likelihood)
  # ==================================================================
  ic <- compute_information_criteria(loglik, num_params, k * n)
  
  # ==================================================================
  # Structured β (beta0): extract per-variable coefficient vectors
  # ==================================================================
  if (is.null(beta0) && !is.null(beta) && !is.null(data_list$data_info$x_vars)) {
    beta0 <- tryCatch({
      structure_beta(beta, data_list, k)
    }, error = function(e) NULL)
  }
  
  # ==================================================================
  # Alpha: extract the K×K temporal AR coefficient matrix A
  # ==================================================================
  if (is.null(alpha) && !is.null(beta) && include_time_lag) {
    alpha <- tryCatch({
      extract_alpha(beta, data_list, k)
    }, error = function(e) NULL)
  }
  
  # ==================================================================
  # Residuals: compute raw and standardised residual vectors
  # ==================================================================
  if (is.null(residuals_raw) && !is.null(beta)) {
    residuals_raw <- tryCatch({
      if (model_type %in% c("SLY", "SLY_diagonal")) {
        compute_residuals(R, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else if (model_type %in% c("SEM", "SEM_diagonal")) {
        compute_residuals_sem(T_mat, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else if (model_type %in% c("SDEM", "SDEM_diagonal")) {
        compute_residuals_sdem(R, T_mat, beta, data_list$y, data_list$X, data_list$W, k, n)
      } else {
        # OLS / VARX: y - Xβ
        data_list$y - data_list$X %*% beta
      }
    }, error = function(e) NULL)
  }
  
  # Standardised residuals
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
  
  # ==================================================================
  # Compute averaged pseudo R² (Eq. 26, mstr.pdf):
  # ŷ_k = [(I-R⊗W)^{-1} Xβ]_k,  R²_{k,pseudo} = corr(y_k, ŷ_k)²
  # R̄²_pseudo = (1/K) Σ_k R²_{k,pseudo}
  # ==================================================================
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
        R2_adj <- NA  # no adjusted version for pseudo R^2
      }
    }
  }
  
  # ==================================================================
  # Build penalty section for the result object
  # ==================================================================
  penalty_section <- NULL
  if (!is.null(gamma) && gamma > 0) {
    penalty_section <- list(
      gamma = gamma,
      value = if (!is.null(penalty_value)) penalty_value else 0,
      penalized_loglik = if (!is.null(penalized_loglik)) penalized_loglik else loglik
    )
  }
  
  # ==================================================================
  # Determine the S3 class name based on model type
  # ==================================================================
  class_map <- list(
    "SLY"           = c("multivar_sly", "multivar_spatial"),
    "SEM"           = c("multivar_sem", "multivar_spatial"),
    "SDEM"          = c("multivar_sdem", "multivar_spatial"),
    "SLY_diagonal"  = c("multivar_sly_diagonal", "multivar_sly", "multivar_spatial"),
    "SEM_diagonal"  = c("multivar_sem_diagonal", "multivar_sem", "multivar_spatial"),
    "SDEM_diagonal" = c("multivar_sdem_diagonal", "multivar_sdem", "multivar_spatial"),
    "VARX"          = c("multivar_varx", "multivar_spatial"),
    "OLS_diagonal"  = c("multivar_ols_diagonal", "multivar_ols", "multivar_spatial")
  )
  
  obj_class <- class_map[[model_type]]
  if (is.null(obj_class)) {
    obj_class <- c("multivar_spatial")
    warning(sprintf("Unknown model_type: %s", model_type))
  }
  
  # ==================================================================
  # Build the object
  # ==================================================================
  result <- structure(
    list(
      # Model identification
      model_type = model_type,
      model_description = get_model_description(model_type),
      
      # Estimated coefficients
      coefficients = list(
        R     = R,
        T     = T_mat,
        beta  = beta,
        beta0 = beta0,
        alpha = alpha,
        Sigma = Sigma
      ),
      
      # Goodness of fit
      fit = list(
        loglik     = loglik,
        AIC        = ic$AIC,
        BIC        = ic$BIC,
        num_params = num_params,
        num_obs    = k * n,
        R2         = R2,
        R2_adj     = R2_adj
      ),
      
      # Convergence information
      convergence = list(
        converged  = converged,
        iterations = iterations,
        method     = method,
        message    = ifelse(converged, "Converged", "maximum iterations reached or warning")
      ),
      
      # Penalty information (full model only, gamma>0 case)
      penalty = penalty_section,
      
      # Data information
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
      
      # Residuals
      residuals = list(
        raw          = if (!is.null(residuals_raw)) as.numeric(residuals_raw) else NULL,
        standardized = if (!is.null(residuals_std)) as.numeric(residuals_std) else NULL
      ),
      
      # Model data (used by add_inference etc.)
      model_data = list(
        y       = data_list$y,
        X       = data_list$X,
        W       = data_list$W,
        W_listw = data_list$W_listw,
        y_lag   = data_list$y_lag,
        eigen_W = data_list$eigen_W
      ),
      
      # initial value (S1 estimation result of the full model)
      initial_values = if (!is.null(individual_estimates)) {
        list(individual_estimates = individual_estimates)
      } else {
        NULL
      },
      
      # standard error (already obtained from spatialreg for diagonal models)
      std_errors = list(
        beta = NULL,
        R    = std_errors_R,
        T    = std_errors_T
      ),
      
      # Inference result (added later by add_inference)
      inference = NULL,
      vcov      = NULL,
      hessian   = NULL,
      
      # Individual-model result (diagonal model only)
      individual_models = individual_models,
      
      # Execution metadata
      execution = list(
        time      = execution_time,
        call      = sys.call(-1),
        R_version = R.version.string
      ),
      
      # Original data_list (needed by subsequent processing)
      data_list = data_list
    ),
    class = obj_class
  )
  
  return(result)
}


################################################################################
# 2. Model-name helper
################################################################################

#' Retrieve a description string from the model type
#'
#' @param model_type  Model type string (e.g. 'SLY', 'SEM', 'SDEM')
#' @return description string
get_model_description <- function(model_type) {
  descriptions <- list(
    "SLY"           = "Multivariate Spatial Lag of Y Model (full)",
    "SEM"           = "Multivariate Spatial Error Model (full)",
    "SDEM"          = "Multivariate Spatial Durbin Error Model (full)",
    "SLY_diagonal"  = "Multivariate Spatial Lag of Y Model (diagonal)",
    "SEM_diagonal"  = "Multivariate Spatial Error Model (diagonal)",
    "SDEM_diagonal" = "Multivariate Spatial Durbin Error Model (diagonal)",
    "VARX"          = "Vector Autoregressive with Exogenous Variables Model",
    "OLS_diagonal"  = "OLS Diagonal Model (no spatial dependence)"
  )
  desc <- descriptions[[model_type]]
  if (is.null(desc)) return("Unknown Model")
  return(desc)
}


################################################################################
# 3. for diagonal models: minimal data_list build helper
################################################################################

#' Build a minimal data_list for diagonal models
#'
#' Because diagonal models do not use prepare_data_extended,
#' build a minimal data_list structure to pass to build_result_object.
#'
#' @param y y vector (k*n x 1)
#' @param X X matrix (NULL allowed; may be unnecessary for diagonal models since spatialreg handles it internally)
#' @param W  n×n spatial weight matrix
#' @param W_listw listw object
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_point  Analysis time period (default: max available)
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param include_intercept  Logical; prepend intercept column
#' @param region_var  Name of the region identifier column
#' @param time_var  Name of the time index column
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
  
  # eigenvalues (computed if W is available)
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


################################################################################
# Usage examples
################################################################################
#
# # Inside fit_sly_penalized in penalized_spatial.r:
#
# result <- build_result_object(
#   model_type      = "SLY",
#   R               = R_final,
#   T_mat           = NULL,
#   beta            = beta_final,
#   Sigma           = Sigma_final,
#   loglik          = loglik,
#   num_params      = num_params,
#   converged       = opt_result$converged,
#   method          = "penalized_lbfgsb",
#   iterations      = opt_result$optim_result$counts["function"],
#   data_list       = data_list,
#   gamma           = gamma,
#   penalty_value   = opt_result$penalty,
#   penalized_loglik = opt_result$penalized_loglik,
#   execution_time  = exec_time
# )
#
# # Inside fit_sly_diagonal.r:
#
# result <- build_result_object(
#   model_type        = "SLY_diagonal",
#   R                 = R_diag,
#   T_mat             = NULL,
#   Sigma             = Sigma_diag,
#   loglik            = total_loglik,
#   num_params        = total_params,
#   converged         = TRUE,
#   method            = "lagsarlm_diagonal",
#   data_list         = data_list,
#   individual_models = individual_models,
#   std_errors_R      = R_se,
#   beta0             = beta_list,
#   residuals_raw     = residuals_raw,
#   residuals_std     = residuals_std
# )
#
################################################################################

################################################################################
# START OF FILE: penalized_spatial.r
################################################################################

################################################################################
# penalized_spatial.r
#
# penalised likelihood estimation of the spatial parameters
#
# Purpose:
#   When sample size is small (e.g. n=46 prefectures), spatial parameters (ρ, λ) may
#   Introduce an L2 penalty to mitigate sticking to the boundary values (+/-0.99).
#
# Penalised likelihood:
#   ℓ_pen(θ) = ℓ(θ) - γ/2 * ||θ_spatial||²
#
#   SLY:  ℓ_pen = ℓ - γ/2 * ||R||²_F
#   SEM:  ℓ_pen = ℓ - γ/2 * ||Λ||²_F
#   SDEM: ℓ_pen = ℓ - γ/2 * (||R||²_F + ||Λ||²_F)
#
#   where ||·||_F denotes the Frobenius norm
#
# Note:
#   - Penalty is used only during optimisation; AIC/BIC use the unpenalised likelihood
#   - γ = 0: no penalty (equivalent to standard MLE)
#   - Larger γ → spatial parameters shrink toward 0 (analogous to ridge regression)
#   - γ can be selected by cross-validation or pBIC minimisation
#
# Usage:
#
################################################################################

cat("penalized_spatial.r loaded (penalized estimation of spatial parameters)
")

################################################################################
# 1. penalty function
################################################################################

#' Compute the ridge (L2) penalty term used in the penalised log-likelihood (Eq. 36).
#'
#' The penalised log-likelihood is:
#'   ℓ_p(θ) = ℓ(θ) - (γ/2) θ' D θ
#'
#' where D = diag(0,...,0, 1,...,1) selects only spatial parameters.
#' For each spatial matrix M passed as an argument:
#'   penalty = (γ/2) Σ_i ||M_i||²_F
#'
#' @param ... One or more K×K spatial parameter matrices (R, Λ, or both)
#' @param gamma Non-negative penalty strength γ
#' @return Scalar penalty value (γ/2) Σ ||M_i||²_F
compute_penalty <- function(..., gamma) {
  matrices <- list(...)
  total <- 0
  for (M in matrices) {
    total <- total + sum(M^2)
  }
  return(gamma / 2 * total)
}


################################################################################
# 1b. Hessian / effective-degrees-of-freedom computation for GIC (generalised information criterion)
#
#   GIC(γ) = -2ℓ(θ̂_γ) + 2·df_eff(γ)
#   df_eff(γ) = tr[H(H + 2γI)^{-1}] + k₂
#
#   H: negative Hessian of the unpenalised profile log-likelihood (observed information matrix)
#   k₂: number of non-spatial parameters (β and upper-triangle of Σ)
#
# Depends: numDeriv package
################################################################################

#' Numerical Hessian of the unpenalised profile log-likelihood (SLY model)
#'
#' @param R_hat estimated R matrix
#' @param beta_hat estimated beta vector
#' @param Sigma_hat estimated Sigma matrix
#' @param data_list  Output of prepare_data_extended()
#' @param max_iter_inner maximum number of inner iterations
#' @param tol  Convergence tolerance
#' @return  Scalar log-likelihood or profile log-likelihood
compute_profile_hessian_sly <- function(
  R_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- data_list$W; eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_R <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  neg_profile_loglik <- function(param) {
    R <- param_to_R(param)
    
    # Phase 2: omit check_stationarity for smoother objective function (better Hessian)
    
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


#' Numerical Hessian of the unpenalised profile log-likelihood (SEM model)
compute_profile_hessian_sem <- function(
  T_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- data_list$W; eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_T <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  neg_profile_loglik <- function(param) {
    T_mat <- param_to_T(param)
    
    # Phase 2: omit check_stationarity for smoother objective function (better Hessian)
    
    inner <- tryCatch({
      iterate_beta_sigma_sem(T_mat = T_mat, beta_init = beta_ref, Sigma_init = Sigma_ref,
                              y = y, X = X, W = W, k = k, n = n,
                              max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(inner)) {
      rho_T <- spectral_radius(T_mat, eigen_W)
      return(1e6 * max(rho_T, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_sem(T_mat = T_mat, beta_hat = inner$beta, Sigma_hat = inner$Sigma,
                                      y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
                                      smooth = TRUE)
    }, error = function(e) -Inf)
    if (!is.finite(prof_lik)) {
      rho_T <- spectral_radius(T_mat, eigen_W)
      return(1e6 * max(rho_T, 1))
    }
    
    return(-prof_lik)
  }
  
  theta1_hat <- as.vector(t(T_hat))
  H <- numDeriv::hessian(neg_profile_loglik, theta1_hat)
  return(H)
}


#' Numerical Hessian of the unpenalised profile log-likelihood (SDEM model)
compute_profile_hessian_sdem <- function(
  R_hat, T_hat, beta_hat, Sigma_hat, data_list,
  max_iter_inner = 50, tol_inner = 1e-6
) {
  k <- data_list$k; n <- data_list$n
  y <- data_list$y; X <- data_list$X
  W <- data_list$W; eigen_W <- data_list$eigen_W
  
  beta_ref <- beta_hat; Sigma_ref <- Sigma_hat
  
  param_to_RT <- function(param) {
    R <- matrix(param[1:(k*k)], nrow = k, ncol = k, byrow = TRUE)
    T_mat <- matrix(param[(k*k + 1):(2*k*k)], nrow = k, ncol = k, byrow = TRUE)
    list(R = R, T = T_mat)
  }
  
  neg_profile_loglik <- function(param) {
    RT <- param_to_RT(param)
    
    # Phase 2: improve Hessian accuracy with a smooth objective
    
    inner <- tryCatch({
      iterate_beta_sigma_sdem(R = RT$R, T_mat = RT$T,
                               beta_init = beta_ref, Sigma_init = Sigma_ref,
                               y = y, X = X, W = W, k = k, n = n,
                               max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(inner)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_T <- spectral_radius(RT$T, eigen_W)
      return(1e6 * max(rho_R, rho_T, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_sdem(R = RT$R, T_mat = RT$T,
                                       beta_hat = inner$beta, Sigma_hat = inner$Sigma,
                                       y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
                                       smooth = TRUE)
    }, error = function(e) -Inf)
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_T <- spectral_radius(RT$T, eigen_W)
      return(1e6 * max(rho_R, rho_T, 1))
    }
    
    return(-prof_lik)
  }
  
  theta1_hat <- c(as.vector(t(R_hat)), as.vector(t(T_hat)))
  H <- numDeriv::hessian(neg_profile_loglik, theta1_hat)
  return(H)
}


#' Compute the Generalised Information Criterion (GIC / pAIC / pBIC) (Eq. 25).
#'
#' The penalised information criteria are defined as:
#'   pAIC = -2 ℓ(θ̂) + 2 d_eff(γ)
#'   pBIC = -2 ℓ(θ̂) + log(Kn) · d_eff(γ)
#'
#' where the effective degrees of freedom (Eq. 24) are:
#'   d_eff(γ) = tr[ I_{11}(θ̂) · {I_{11}(θ̂) + γ I_{2K²}}^{-1} ]  +  p + K(K+1)/2
#'
#' The first term is the shrinkage-adjusted dimension of the spatial parameter
#' block (computed from the observed information matrix H via the trace formula),
#' and the second term is the ordinary dimension of the unpenalised (β, Σ) block.
#' When γ = 0, d_eff reduces to the nominal parameter count.
#'
#' @param H      k₁×k₁ negative Hessian of the unpenalised profile log-likelihood
#'               (observed information matrix for spatial parameters θ₁)
#' @param gamma  Penalty strength γ ≥ 0
#' @param loglik Unpenalised maximised log-likelihood ℓ(θ̂)
#' @param k2     Number of unpenalised parameters (β + upper-triangle of Σ)
#' @param n_obs  Effective sample size Kn
#' @return A list with: df_eff, df_spatial, GIC_AIC (= pAIC), GIC_BIC (= pBIC), H
compute_gic <- function(H, gamma, loglik, k2, n_obs) {
  k1 <- nrow(H)
  
  if (gamma == 0) {
    # γ = 0: no shrinkage → df_eff = k1 + k2, coincides with standard AIC/BIC
    df_spatial <- k1
  } else {
    # γ > 0: shrinkage-adjusted effective dimension (Eq. 24)
    # df_spatial = tr[ H · (H + γ I_{k1})^{-1} ]
    H_reg <- H + gamma * diag(k1)
    
    # numerical stability check
    H_reg_cond <- tryCatch(kappa(H_reg), error = function(e) Inf)
    if (H_reg_cond > 1e12) {
      warning(sprintf("GIC: condition number of H + gamma I is large (kappa=%.2e, gamma=%.4g)", H_reg_cond, gamma))
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


################################################################################
# 2. SLY penalisedL-BFGS-B
################################################################################

#' Optimise the R matrix via penalised L-BFGS-B (SLY model)
#'
#' @param gamma  Penalty strength γ ≥ 0
#' @param ...  Additional arguments
optimize_R_lbfgsb_penalized <- function(
  R_init,
  beta_current,
  Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  factr = 1e7,
  pgtol = 1e-5,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- data_list$W
  eigen_W <- data_list$eigen_W
  
  R_to_param <- function(R) as.vector(t(R))
  param_to_R <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  # Penalised objective function (smooth version)
  objective <- function(param) {
    R <- param_to_R(param)
    
    # check_stationarity is not called here; the log-det barrier is sufficient
    
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
    
    # Add the L2 penalty term to the negative profile likelihood (Eq. 36)
    penalty <- compute_penalty(R, gamma = gamma)
    
    return(-prof_lik + penalty)  # Minimise: negative profile likelihood + penalty
  }
  
  param_init <- R_to_param(R_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalised BFGS optimisation start (SLY) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Parameters: %d, gamma = %.4f\n", k*k, gamma))
    init_obj <- objective(param_init)
    cat(sprintf("  Initial penalised objective: %.4f\n", init_obj))
  }
  
  # Step S4: BFGS quasi-Newton update of spatial parameters (Section 3.3)
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  # If BFGS did not converge, fall back to Nelder-Mead for robustness
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS did not converge — retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  R_final <- param_to_R(result$par)
  
  # Compute the final beta, Sigma
  final_inner <- iterate_beta_sigma(
    R = R_final, beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner,
    ridge_eps = 1e-6, verbose = 0)
  
  # Evaluate the unpenalised log-likelihood for AIC/BIC computation
  final_loglik <- compute_log_likelihood(
    R = R_final, beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  # Penalised log-likelihood (stored for reference, not used in IC)
  penalized_loglik <- final_loglik - compute_penalty(R_final, gamma = gamma)
  
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s, iterations: %d\n", ifelse(converged, "yes", "no"), result$counts["function"]))
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalised log-likelihood: %.4f\n", penalized_loglik))
    cat(sprintf("  Penalty value: %.4f\n", compute_penalty(R_final, gamma = gamma)))
    cat("  Final R matrix:\n"); print(round(R_final, 4))
    cat("=== Penalised BFGS complete (SLY) ===\n\n")
  }
  
  return(list(
    R = R_final,
    beta = final_inner$beta,
    Sigma = final_inner$Sigma,
    loglik = final_loglik,           # unpenalised (for AIC/BIC)
    penalized_loglik = penalized_loglik,  # penalised (for reference)
    gamma = gamma,
    penalty = compute_penalty(R_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}


################################################################################
# 3. SEM penalisedL-BFGS-B
################################################################################

optimize_T_lbfgsb_penalized <- function(
  T_init,
  beta_current, Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  factr = 1e7,
  pgtol = 1e-5,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- data_list$W
  eigen_W <- data_list$eigen_W
  
  T_to_param <- function(T_mat) as.vector(t(T_mat))
  param_to_T <- function(param) matrix(param, nrow = k, ncol = k, byrow = TRUE)
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  objective <- function(param) {
    T_mat <- param_to_T(param)
    
    # check_stationarity is not called here; the log-det barrier is sufficient
    
    inner_result <- tryCatch({
      iterate_beta_sigma_sem(
        T_mat = T_mat, beta_init = beta_ref, Sigma_init = Sigma_ref,
        y = y, X = X, W = W, k = k, n = n,
        max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    
    if (is.null(inner_result)) {
      rho_T <- spectral_radius(T_mat, eigen_W)
      return(1e6 * max(rho_T, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_sem(
        T_mat = T_mat, beta_hat = inner_result$beta, Sigma_hat = inner_result$Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
        smooth = TRUE)
    }, error = function(e) -Inf)
    
    if (!is.finite(prof_lik)) {
      rho_T <- spectral_radius(T_mat, eigen_W)
      return(1e6 * max(rho_T, 1))
    }
    
    # ★ penalty term
    penalty <- compute_penalty(T_mat, gamma = gamma)
    
    return(-prof_lik + penalty)
  }
  
  param_init <- T_to_param(T_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalised BFGS optimisation start (SEM) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Parameters: %d, gamma = %.4f\n", k*k, gamma))
  }
  
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS not converged -> retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  T_final <- param_to_T(result$par)
  
  final_inner <- iterate_beta_sigma_sem(
    T_mat = T_final, beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
  
  # unpenalised likelihood
  final_loglik <- compute_log_likelihood_sem(
    T_mat = T_final, beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  penalized_loglik <- final_loglik - compute_penalty(T_final, gamma = gamma)
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s\n", ifelse(converged, "yes", "no")))
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalty value: %.4f\n", compute_penalty(T_final, gamma = gamma)))
    cat("  Final Lambda matrix:\n"); print(round(T_final, 4))
    cat("=== Penalised BFGS complete (SEM) ===\n\n")
  }
  
  return(list(
    T = T_final,
    beta = final_inner$beta,
    Sigma = final_inner$Sigma,
    loglik = final_loglik,
    penalized_loglik = penalized_loglik,
    gamma = gamma,
    penalty = compute_penalty(T_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}


################################################################################
# 4. SDEM penalisedL-BFGS-B
################################################################################

optimize_RT_lbfgsb_penalized <- function(
  R_init, T_init,
  beta_current, Sigma_current,
  data_list,
  gamma = 0,
  max_iter_inner = 50,
  tol_inner = 1e-6,
  max_iter = 100,
  factr = 1e7,
  pgtol = 1e-5,
  verbose = 1
) {
  
  k <- data_list$k
  n <- data_list$n
  y <- data_list$y
  X <- data_list$X
  W <- data_list$W
  eigen_W <- data_list$eigen_W
  
  param_to_RT <- function(param) {
    R <- matrix(param[1:(k*k)], nrow = k, ncol = k, byrow = TRUE)
    T_mat <- matrix(param[(k*k + 1):(2*k*k)], nrow = k, ncol = k, byrow = TRUE)
    list(R = R, T = T_mat)
  }
  RT_to_param <- function(R, T_mat) c(as.vector(t(R)), as.vector(t(T_mat)))
  
  beta_ref <- beta_current
  Sigma_ref <- Sigma_current
  
  objective <- function(param) {
    RT <- param_to_RT(param)
    
    # check_stationarity is not called here; the log-det barrier is sufficient
    
    inner_result <- tryCatch({
      iterate_beta_sigma_sdem(
        R = RT$R, T_mat = RT$T,
        beta_init = beta_ref, Sigma_init = Sigma_ref,
        y = y, X = X, W = W, k = k, n = n,
        max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
    }, error = function(e) NULL)
    
    if (is.null(inner_result)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_T <- spectral_radius(RT$T, eigen_W)
      return(1e6 * max(rho_R, rho_T, 1))
    }
    
    prof_lik <- tryCatch({
      compute_profile_likelihood_sdem(
        R = RT$R, T_mat = RT$T,
        beta_hat = inner_result$beta, Sigma_hat = inner_result$Sigma,
        y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n,
        smooth = TRUE)
    }, error = function(e) -Inf)
    
    if (!is.finite(prof_lik)) {
      rho_R <- spectral_radius(RT$R, eigen_W)
      rho_T <- spectral_radius(RT$T, eigen_W)
      return(1e6 * max(rho_R, rho_T, 1))
    }
    
    # * penalty on both R and Lambda
    penalty <- compute_penalty(RT$R, RT$T, gamma = gamma)
    
    return(-prof_lik + penalty)
  }
  
  param_init <- RT_to_param(R_init, T_init)
  
  if (verbose >= 1) {
    cat("\n", paste(rep("=", 60), collapse=""), "\n")
    cat("=== Penalised BFGS optimisation start (SDEM) ===\n")
    cat(paste(rep("=", 60), collapse=""), "\n")
    cat(sprintf("  Parameters: %d (R: %d, Lambda: %d), gamma = %.4f\n", 2*k*k, k*k, k*k, gamma))
  }
  
  result <- optim(
    par = param_init, fn = objective, method = "BFGS",
    control = list(maxit = max_iter, reltol = 1e-10,
                   trace = ifelse(verbose >= 2, 1, 0)))
  
  if (result$convergence != 0) {
    if (verbose >= 1) cat("  BFGS not converged -> retrying with Nelder-Mead...\n")
    result_nm <- optim(
      par = result$par, fn = objective, method = "Nelder-Mead",
      control = list(maxit = max_iter * 10,
                     trace = ifelse(verbose >= 2, 1, 0)))
    if (result_nm$value < result$value) result <- result_nm
  }
  
  RT_opt <- param_to_RT(result$par)
  R_final <- RT_opt$R; T_final <- RT_opt$T
  
  final_inner <- iterate_beta_sigma_sdem(
    R = R_final, T_mat = T_final,
    beta_init = beta_ref, Sigma_init = Sigma_ref,
    y = y, X = X, W = W, k = k, n = n,
    max_iter = max_iter_inner, tol = tol_inner, verbose = 0)
  
  final_loglik <- compute_log_likelihood_sdem(
    R = R_final, T_mat = T_final,
    beta = final_inner$beta, Sigma = final_inner$Sigma,
    y = y, X = X, W = W, eigen_W = eigen_W, k = k, n = n)
  
  penalized_loglik <- final_loglik - compute_penalty(R_final, T_final, gamma = gamma)
  converged <- (result$convergence == 0)
  
  if (verbose >= 1) {
    cat(sprintf("  Converged: %s\n", ifelse(converged, "yes", "no")))
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", final_loglik))
    cat(sprintf("  Penalty value: %.4f\n", compute_penalty(R_final, T_final, gamma = gamma)))
    cat("  Final R matrix:\n"); print(round(R_final, 4))
    cat("  Final Lambda matrix:\n"); print(round(T_final, 4))
    cat("=== Penalised BFGS complete (SDEM) ===\n\n")
  }
  
  return(list(
    R = R_final, T = T_final,
    beta = final_inner$beta, Sigma = final_inner$Sigma,
    loglik = final_loglik,
    penalized_loglik = penalized_loglik,
    gamma = gamma,
    penalty = compute_penalty(R_final, T_final, gamma = gamma),
    converged = converged,
    optim_result = result
  ))
}


################################################################################
# 5. Combined fit function (penalty-enabled version)
################################################################################

#' Penalised multivariate SLY estimation
#'
#' gamma > 0 gives penalised estimation; gamma = 0 gives ordinary estimation.
#' AIC/BIC are computed from the unpenalised likelihood.
fit_sly_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  R_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalised multivariate SLY estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # S0: Data preparation
  data_list <- prepare_data_extended(
    data_file = data_file, weight_file = weight_file,
    y_vars = y_vars, x_vars = x_vars,
    time_var = time_var, time_point = time_point,
    region_var = region_var, include_intercept = include_intercept,
    include_time_lag = include_time_lag, verbose = FALSE)
  
  k <- data_list$k; n <- data_list$n
  actual_include_time_lag <- data_list$data_info$include_time_lag
  if (is.null(actual_include_time_lag)) actual_include_time_lag <- TRUE
  
  if (verbose) cat(sprintf("  data: k=%d, n=%d\n", k, n))
  
  # ===== S1: Initial estimation =====
  # Estimate individual SAR/SEM models for each response variable y_k
  # to obtain initial diagonal matrices R^(0), Λ^(0), Σ^(0) (Step S1, mstr.pdf)
  if (is.null(R_init)) {
    if (verbose) cat("S1: individualSLY for initial estimation...\n")
    init_result <- initial_estimation_sly_extended(data_list, verbose = ifelse(verbose, 1, 0))
    R <- init_result$R_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    R <- R_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  # ===== S2–S5: Penalised outer loop =====
  # Calls the BFGS-based optimiser which internally runs the inner β/Σ loop
  # (S2) and evaluates the profile likelihood (S3) at each outer iterate (S4)
  opt_result <- optimize_R_lbfgsb_penalized(
    R_init = R, beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  R_final <- opt_result$R
  beta_final <- opt_result$beta
  Sigma_final <- opt_result$Sigma
  loglik <- opt_result$loglik  # unpenalised
  
  # Compute information criteria from the UNPENALISED log-likelihood
  # pAIC and pBIC use d_eff(γ) from compute_gic() (Eqs. 24–25)
  num_params <- ncol(data_list$X) + k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  gamma = %.4f, penalty value = %.4f\n", gamma, opt_result$penalty))
    cat(sprintf("  R diagonal: [%s]\n", paste(sprintf("%.4f", diag(R_final)), collapse=", ")))
    cat(sprintf("  Execution time: %.2f sec\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "SLY",
    R                = R_final,
    beta             = beta_final,
    Sigma            = Sigma_final,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_lbfgsb",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}


#' Penalised multivariate SEM estimation
fit_sem_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  T_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalised multivariate SEM estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  data_list <- prepare_data_extended(
    data_file = data_file, weight_file = weight_file,
    y_vars = y_vars, x_vars = x_vars,
    time_var = time_var, time_point = time_point,
    region_var = region_var, include_intercept = include_intercept,
    include_time_lag = include_time_lag, verbose = FALSE)
  
  k <- data_list$k; n <- data_list$n
  if (verbose) cat(sprintf("  data: k=%d, n=%d\n", k, n))
  
  if (is.null(T_init)) {
    if (verbose) cat("S1: individualSEM for initial estimation...\n")
    init_result <- initial_estimation_sem_extended(data_list, verbose = ifelse(verbose, 1, 0))
    T_mat <- init_result$T_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    T_mat <- T_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  opt_result <- optimize_T_lbfgsb_penalized(
    T_init = T_mat, beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  T_final <- opt_result$T
  loglik <- opt_result$loglik
  
  num_params <- ncol(data_list$X) + k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  gamma = %.4f, penalty value = %.4f\n", gamma, opt_result$penalty))
    cat(sprintf("  Lambda diagonal: [%s]\n", paste(sprintf("%.4f", diag(T_final)), collapse=", ")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "SEM",
    T_mat            = T_final,
    beta             = opt_result$beta,
    Sigma            = opt_result$Sigma,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_lbfgsb",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}


#' Penalised multivariate SDEM estimation
fit_sdem_penalized <- function(
  data_file, weight_file, y_vars, x_vars,
  time_var = "time", time_point = NULL, region_var = "region",
  include_intercept = TRUE, include_time_lag = TRUE,
  R_init = NULL, T_init = NULL, Sigma_init = NULL,
  gamma = 0,
  max_iter_outer = 100, max_iter_inner = 100, tol = 1e-6,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Penalised multivariate SDEM estimation ===\n")
    cat(sprintf("    γ = %.4f\n", gamma))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  data_list <- prepare_data_extended(
    data_file = data_file, weight_file = weight_file,
    y_vars = y_vars, x_vars = x_vars,
    time_var = time_var, time_point = time_point,
    region_var = region_var, include_intercept = include_intercept,
    include_time_lag = include_time_lag, verbose = FALSE)
  
  k <- data_list$k; n <- data_list$n
  if (verbose) cat(sprintf("  data: k=%d, n=%d\n", k, n))
  
  if (is.null(R_init) || is.null(T_init)) {
    if (verbose) cat("S1: individualSDEM for initial estimation...\n")
    init_result <- initial_estimation_sdem_extended(data_list, verbose = ifelse(verbose, 1, 0))
    R <- if(is.null(R_init)) init_result$R_init else R_init
    T_mat <- if(is.null(T_init)) init_result$T_init else T_init
    Sigma <- if(is.null(Sigma_init)) init_result$Sigma_init else Sigma_init
  } else {
    R <- R_init; T_mat <- T_init
    Sigma <- if(is.null(Sigma_init)) diag(1, k) else Sigma_init
  }
  
  beta <- rep(0, ncol(data_list$X))
  
  opt_result <- optimize_RT_lbfgsb_penalized(
    R_init = R, T_init = T_mat,
    beta_current = beta, Sigma_current = Sigma,
    data_list = data_list, gamma = gamma,
    max_iter_inner = max_iter_inner, tol_inner = tol, max_iter = max_iter_outer,
    verbose = ifelse(verbose, 1, 0))
  
  R_final <- opt_result$R; T_final <- opt_result$T
  loglik <- opt_result$loglik
  
  num_params <- ncol(data_list$X) + 2*k^2 + k*(k+1)/2
  ic <- compute_information_criteria(loglik, num_params, k*n)
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n--- Estimation results ---\n")
    cat(sprintf("  Unpenalised log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f, BIC: %.4f\n", ic$AIC, ic$BIC))
    cat(sprintf("  gamma = %.4f, penalty value = %.4f\n", gamma, opt_result$penalty))
    cat("  R:\n"); print(round(R_final, 4))
    cat("  Λ:\n"); print(round(T_final, 4))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  return(build_result_object(
    model_type       = "SDEM",
    R                = R_final,
    T_mat            = T_final,
    beta             = opt_result$beta,
    Sigma            = opt_result$Sigma,
    loglik           = loglik,
    num_params       = num_params,
    converged        = opt_result$converged,
    method           = "penalized_lbfgsb",
    iterations       = opt_result$optim_result$counts["function"],
    data_list        = data_list,
    gamma            = gamma,
    penalty_value    = opt_result$penalty,
    penalized_loglik = opt_result$penalized_loglik,
    execution_time   = exec_time
  ))
}


################################################################################
# 6. Gamma-search utility
################################################################################

#' Run penalised estimation over multiple gamma values and select the optimal gamma by GIC
#'
#' @param fit_func  One of: fit_sly_penalized, fit_sem_penalized, fit_sdem_penalized
#' @param gammas vector of gamma values to try (e.g. c(0, 10^seq(-6, 2, by=0.5)))
#' @param compute_gic_flag whether to compute GIC (numDeriv required when TRUE)
#' @param ...  Additional arguments
#' @return result data frame and the full list of results
compare_gamma <- function(
  fit_func,
  gammas = c(0, 0.5, 1, 2, 5, 10, 20, 50),
  compute_gic_flag = TRUE,
  ...
) {
  
  cat("\n", paste(rep("#", 70), collapse=""), "\n")
  cat("### Gamma search: comparison of penalty strengths (GIC-enabled) ###\n")
  cat(paste(rep("#", 70), collapse=""), "\n\n")
  
  # numDeriv check
  if (compute_gic_flag && !requireNamespace("numDeriv", quietly = TRUE)) {
    warning("numDeriv package not installed. Skipping GIC computation. ")
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
  
  # ==================================================================
  # Phase 1: Estimate over all gamma
  # (each gamma fit is independent -> run them in parallel; the fit
  #  itself — S1 initial estimation, S2–S5 BFGS — is unchanged)
  # ==================================================================
  cat("--- Phase 1: estimation ---\n")
  cat(sprintf("  Fitting %d gamma value(s) in parallel (up to %d cores)...\n",
              length(gammas), MSTR_CORES))
  
  extra_args <- list(...)
  fit_one_gamma <- function(g) {
    tryCatch({
      do.call(fit_func, c(list(gamma = g, verbose = FALSE), extra_args))
    }, error = function(e) {
      structure(list(message = conditionMessage(e)), class = "mstr_fit_error")
    })
  }
  # Give the worker function a minimal, self-contained environment so that the
  # PSOCK backend only ships fit_func and the (already evaluated) extra
  # arguments — not the whole compare_gamma frame.
  environment(fit_one_gamma) <- list2env(
    list(fit_func = fit_func, extra_args = extra_args),
    parent = globalenv()
  )
  
  fit_list <- mstr_parallel_lapply(as.list(gammas), fit_one_gamma)
  
  for (i in seq_along(gammas)) {
    g <- gammas[i]
    res <- fit_list[[i]]
    cat(sprintf("  γ = %.4g (%d/%d) ... ", g, i, length(gammas)))
    
    if (is.null(res) || inherits(res, "try-error")) {
      cat("Error: worker failed\n")
      next
    }
    if (inherits(res, "mstr_fit_error")) {
      cat(sprintf("Error: %s\n", res$message))
      next
    }
    
    results[[as.character(g)]] <- res
    cat(sprintf("loglik = %.4f\n", res$fit$loglik))
  }
  
  # ==================================================================
  # Phase 2: Hessian computation and GIC
  # ==================================================================
  if (compute_gic_flag && length(results) > 0) {
    cat("\n--- Phase 2: Hessian computation and GIC ---\n")
    
    # Retrieve model information from the first result
    first_res <- results[[1]]
    model_type <- first_res$model_type
    k <- first_res$data_info$k
    n_obs <- first_res$fit$num_obs
    
    # k1: number of spatial parameters, k2: number of non-spatial parameters
    if (model_type == "SLY") {
      k1 <- k^2
    } else if (model_type == "SEM") {
      k1 <- k^2
    } else {  # SDEM
      k1 <- 2 * k^2
    }
    k2 <- first_res$fit$num_params - k1
    
    cat(sprintf("  Model: %s, k1=%d (spatial), k2=%d (non-spatial), n_obs=%d\n",
                model_type, k1, k2, n_obs))
  }
  
  # ------------------------------------------------------------------
  # Pre-compute the profile-likelihood Hessian for every gamma IN
  # PARALLEL. The Hessian at each gamma is independent of the others,
  # and compute_profile_hessian_* is exactly the same function as in
  # the sequential version — only the loop is distributed over cores.
  # ------------------------------------------------------------------
  hessian_list <- list()
  if (compute_gic_flag && length(results) > 0) {
    cat(sprintf("  Computing %d Hessian(s) in parallel (up to %d cores)...\n",
                length(results), MSTR_CORES))
    
    compute_one_hessian <- function(g_str) {
      res <- results[[g_str]]
      tryCatch({
        if (res$model_type == "SLY") {
          compute_profile_hessian_sly(
            R_hat = res$coefficients$R,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        } else if (res$model_type == "SEM") {
          compute_profile_hessian_sem(
            T_hat = res$coefficients$T,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        } else {  # SDEM
          compute_profile_hessian_sdem(
            R_hat = res$coefficients$R,
            T_hat = res$coefficients$T,
            beta_hat = res$coefficients$beta,
            Sigma_hat = res$coefficients$Sigma,
            data_list = res$data_list)
        }
      }, error = function(e) {
        structure(list(message = conditionMessage(e)), class = "mstr_hess_error")
      })
    }
    # Minimal, self-contained environment for the PSOCK backend (only the
    # fitted results are shipped, not the whole compare_gamma frame)
    environment(compute_one_hessian) <- list2env(
      list(results = results),
      parent = globalenv()
    )
    
    hessian_list <- mstr_parallel_lapply(as.list(names(results)), compute_one_hessian)
    names(hessian_list) <- names(results)
  }
  
  for (g_str in names(results)) {
    g <- as.numeric(g_str)
    res <- results[[g_str]]
    
    # String representation of the spatial parameters
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
    
    if (res$model_type == "SLY") {
      sp <- format_matrix_params(res$coefficients$R, "R")
    } else if (res$model_type == "SEM") {
      sp <- format_matrix_params(res$coefficients$T, "Λ")
    } else {
      sp <- paste(format_matrix_params(res$coefficients$R, "R"),
                  format_matrix_params(res$coefficients$T, "Λ"))
    }
    
    pen_val <- if (!is.null(res$penalty)) res$penalty$value else 0
    
    # GICcomputation
    df_eff_val <- res$fit$num_params  # default (num_params when GIC not computed)
    gic_aic_val <- res$fit$AIC
    gic_bic_val <- res$fit$BIC
    
    if (compute_gic_flag) {
      cat(sprintf("  γ = %.4g: Computing Hessian...", g))
      
      # Look up the Hessian pre-computed in parallel above
      hessian_result <- hessian_list[[g_str]]
      if (inherits(hessian_result, "mstr_hess_error")) {
        cat(sprintf(" Error: %s", hessian_result$message))
        hessian_result <- NULL
      } else if (inherits(hessian_result, "try-error")) {
        cat(" Error: worker failed")
        hessian_result <- NULL
      }
      
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
          cat(" GICcomputation failed\n")
        }
      } else {
        cat(" Hessiancomputation failed\n")
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
  
  # Sort by gamma
  summary_df <- summary_df[order(summary_df$gamma), ]
  rownames(summary_df) <- NULL
  
  # ==================================================================
  # Output
  # ==================================================================
  cat("\n", paste(rep("=", 90), collapse=""), "\n")
  cat("### Gamma-search result summary ###\n")
  cat(paste(rep("=", 90), collapse=""), "\n\n")
  
  # Numerical summary
  if (compute_gic_flag) {
    summary_numeric <- summary_df[, c("gamma", "loglik", "AIC", "BIC",
                                       "df_eff", "GIC_AIC", "GIC_BIC", "penalty")]
  } else {
    summary_numeric <- summary_df[, c("gamma", "loglik", "AIC", "BIC", "penalty")]
  }
  print(summary_numeric, row.names = FALSE, right = FALSE)
  
  # Display the spatial parameters in matrix form
  cat("\n--- Trajectory of the spatial parameters ---\n")
  for (i in seq_len(nrow(summary_df))) {
    g <- summary_df$gamma[i]
    res <- results[[as.character(g)]]
    if (is.null(res)) next
    cat(sprintf("\nγ = %.4g:\n", g))
    if (res$model_type == "SLY" || res$model_type == "SDEM") {
      if (nrow(res$coefficients$R) == 2) {
        cat(sprintf("  R = [%7.4f %7.4f]\n", res$coefficients$R[1,1], res$coefficients$R[1,2]))
        cat(sprintf("      [%7.4f %7.4f]\n", res$coefficients$R[2,1], res$coefficients$R[2,2]))
      }
    }
    if (res$model_type == "SEM" || res$model_type == "SDEM") {
      if (nrow(res$coefficients$T) == 2) {
        cat(sprintf("  Λ = [%7.4f %7.4f]\n", res$coefficients$T[1,1], res$coefficients$T[1,2]))
        cat(sprintf("      [%7.4f %7.4f]\n", res$coefficients$T[2,1], res$coefficients$T[2,2]))
      }
    }
  }
  
  # Best criterion
  if (nrow(summary_df) > 0) {
    cat("\n--- Selecting the best gamma ---\n")
    
    best_aic_idx <- which.min(summary_df$AIC)
    best_bic_idx <- which.min(summary_df$BIC)
    cat(sprintf("  AIC best:     γ = %.4g (AIC = %.4f)\n",
                summary_df$gamma[best_aic_idx], summary_df$AIC[best_aic_idx]))
    cat(sprintf("  BIC best:     γ = %.4g (BIC = %.4f)\n",
                summary_df$gamma[best_bic_idx], summary_df$BIC[best_bic_idx]))
    
    if (compute_gic_flag) {
      best_gic_aic_idx <- which.min(summary_df$GIC_AIC)
      best_gic_bic_idx <- which.min(summary_df$GIC_BIC)
      cat(sprintf("  GIC_AIC best: γ = %.4g (GIC_AIC = %.4f, df_eff = %.2f)\n",
                  summary_df$gamma[best_gic_aic_idx], summary_df$GIC_AIC[best_gic_aic_idx],
                  summary_df$df_eff[best_gic_aic_idx]))
      cat(sprintf("  GIC_BIC best: γ = %.4g (GIC_BIC = %.4f, df_eff = %.2f)\n",
                  summary_df$gamma[best_gic_bic_idx], summary_df$GIC_BIC[best_gic_bic_idx],
                  summary_df$df_eff[best_gic_bic_idx]))
    }
  }
  
  cat(paste(rep("=", 90), collapse=""), "\n")
  
  return(list(summary = summary_df, results = results))
}


################################################################################
# 7. Two-stage gamma search
################################################################################

#' Two-stage gamma search: coarse search -> fine search around the optimal gamma
#'
#' Stage 1: search coarsely over gammas_coarse and identify the pAIC/pBIC-best gamma
#' Stage 2: refine the search over the intervals adjacent to the best gamma in log10 steps of 0.05
#'
#' Skip conditions (when Stage 2 is not run):
#'   - γ* == 0  (no penalty is best)
#'   - γ* == max(gammas_coarse)  (monotonically decreasing -> boundary solution)
#'
#' @param fit_func  One of: fit_sly_penalized, fit_sem_penalized, fit_sdem_penalized
#' @param gammas_coarse coarse gamma grid for Stage 1 (e.g. c(0, 10^seq(-2, 4, by=0.5)))
#' @param fine_log10_step log10 step width for Stage 2 (default: 0.05)
#' @param compute_gic_flag whether to compute GIC
#' @param ...  Additional arguments
#' @return same format as compare_gamma: list(summary, results)
compare_gamma_twostage <- function(
  fit_func,
  gammas_coarse = c(0, 10^seq(-2, 4, by = 0.5)),
  fine_log10_step = 0.05,
  compute_gic_flag = TRUE,
  ...
) {
  
  cat("\n", paste(rep("#", 70), collapse=""), "\n")
  cat("### Two-stage gamma search ###\n")
  cat(paste(rep("#", 70), collapse=""), "\n\n")
  
  # ================================================================
  # Stage 1: coarse search
  # ================================================================
  cat("========== Stage 1: coarse search ==========\n")
  cat(sprintf("  Grid: %d points (%.4g to %.4g)\n",
              length(gammas_coarse), min(gammas_coarse), max(gammas_coarse)))
  
  stage1 <- compare_gamma(
    fit_func = fit_func,
    gammas = gammas_coarse,
    compute_gic_flag = compute_gic_flag,
    ...
  )
  
  s1 <- stage1$summary
  
  if (nrow(s1) == 0) {
    warning("No valid results were obtained in Stage 1")
    return(stage1)
  }
  
  # ================================================================
  # Identify the pAIC/pBIC-best gamma from Stage 1
  # ================================================================
  gamma_max <- max(gammas_coarse)
  
  # pAIC best gamma
  if ("GIC_AIC" %in% colnames(s1) && any(!is.na(s1$GIC_AIC))) {
    idx_aic <- which.min(s1$GIC_AIC)
    g_star_aic <- s1$gamma[idx_aic]
    criterion_aic <- "GIC_AIC"
  } else {
    idx_aic <- which.min(s1$AIC)
    g_star_aic <- s1$gamma[idx_aic]
    criterion_aic <- "AIC"
  }
  
  # pBIC best gamma
  if ("GIC_BIC" %in% colnames(s1) && any(!is.na(s1$GIC_BIC))) {
    idx_bic <- which.min(s1$GIC_BIC)
    g_star_bic <- s1$gamma[idx_bic]
    criterion_bic <- "GIC_BIC"
  } else {
    idx_bic <- which.min(s1$BIC)
    g_star_bic <- s1$gamma[idx_bic]
    criterion_bic <- "BIC"
  }
  
  cat(sprintf("\n  Stage 1 Result: %s best gamma* = %.4g, %s best gamma* = %.4g\n",
              criterion_aic, g_star_aic, criterion_bic, g_star_bic))
  
  # ================================================================
  # Stage 2: refined-grid construction
  # ================================================================
  # Extract only the positive values from the coarse grid and sort (for log10 transform)
  gammas_positive <- sort(gammas_coarse[gammas_coarse > 0])
  
  build_fine_grid <- function(g_star, criterion_name) {
    # Skip condition 1: γ* == 0
    if (g_star == 0) {
      cat(sprintf("  %s: γ* = 0 → skip Stage 2 (no penalty is best)\n",
                  criterion_name))
      return(numeric(0))
    }
    # Skip condition 2: γ* == max(gammas_coarse)
    if (g_star >= gamma_max) {
      cat(sprintf("  %s: γ* = %.4g = max(grid) -> skip Stage 2 (boundary solution)\n",
                  criterion_name, g_star))
      return(numeric(0))
    }
    
    # Locate g_star within gammas_positive and retrieve its neighbours
    pos <- which(abs(gammas_positive - g_star) < 1e-12)
    if (length(pos) == 0) {
      # If g_star is not found in gammas_positive (g_star==0 is handled above)
      cat(sprintf("  %s: γ* = %.4g is not in the grid -> skip Stage 2\n",
                  criterion_name, g_star))
      return(numeric(0))
    }
    pos <- pos[1]
    
    # Lower bound: previous positive gamma, or g_star/10 if none
    if (pos > 1) {
      g_low <- gammas_positive[pos - 1]
    } else {
      g_low <- g_star / 10
    }
    # Upper bound: next gamma, or gamma_max if none
    if (pos < length(gammas_positive)) {
      g_high <- gammas_positive[pos + 1]
    } else {
      g_high <- gamma_max
    }
    
    # Refine in log10 steps
    log10_low  <- log10(g_low)
    log10_high <- log10(g_high)
    fine_grid <- 10^seq(log10_low, log10_high, by = fine_log10_step)
    
    cat(sprintf("  %s: gamma* = %.4g -> refine [%.4g, %.4g] (%d points)\n",
                criterion_name, g_star, g_low, g_high, length(fine_grid)))
    
    return(fine_grid)
  }
  
  cat("\n========== Stage 2: preparing the refined search ==========\n")
  fine_aic <- build_fine_grid(g_star_aic, criterion_aic)
  fine_bic <- build_fine_grid(g_star_bic, criterion_bic)
  
  # Merge and de-duplicate
  fine_all <- sort(unique(c(fine_aic, fine_bic)))
  
  # Exclude gammas already computed in Stage 1 (accounting for floating-point rounding error)
  is_already_computed <- sapply(fine_all, function(g) {
    any(abs(s1$gamma - g) / max(g, 1e-10) < 1e-6)
  })
  gammas_stage2 <- fine_all[!is_already_computed]
  
  if (length(gammas_stage2) == 0) {
    cat("  no additional search needed -> return the Stage 1 results as is\n")
    return(stage1)
  }
  
  cat(sprintf("  Stage 2 additional gamma: %d points\n", length(gammas_stage2)))
  
  # ================================================================
  # Stage 2: running the refined search
  # ================================================================
  cat("\n========== Stage 2: running the refined search ==========\n")
  
  stage2 <- compare_gamma(
    fit_func = fit_func,
    gammas = gammas_stage2,
    compute_gic_flag = compute_gic_flag,
    ...
  )
  
  # ================================================================
  # Combined Stage 1 + Stage 2
  # ================================================================
  cat("\n========== Result merging ==========\n")
  
  # Merge the summaries (sorted by gamma)
  combined_summary <- rbind(stage1$summary, stage2$summary)
  combined_summary <- combined_summary[order(combined_summary$gamma), ]
  rownames(combined_summary) <- NULL
  
  # Merge the results
  combined_results <- c(stage1$results, stage2$results)
  
  # ================================================================
  # Display the combined (joint) estimation results
  # ================================================================
  cat("\n", paste(rep("=", 90), collapse=""), "\n")
  cat("### Two-stage gamma search: final result summary ###\n")
  cat(sprintf("  Stage 1: %d points, Stage 2: %d points, total: %d points\n",
              nrow(stage1$summary), nrow(stage2$summary), nrow(combined_summary)))
  cat(paste(rep("=", 90), collapse=""), "\n\n")
  
  # Numerical summary
  if (compute_gic_flag && "GIC_AIC" %in% colnames(combined_summary)) {
    summary_numeric <- combined_summary[, c("gamma", "loglik", "AIC", "BIC",
                                             "df_eff", "GIC_AIC", "GIC_BIC", "penalty")]
  } else {
    summary_numeric <- combined_summary[, c("gamma", "loglik", "AIC", "BIC", "penalty")]
  }
  print(summary_numeric, row.names = FALSE, right = FALSE)
  
  # Display the best gamma
  if (nrow(combined_summary) > 0) {
    cat("\n--- Selecting the best gamma (after two-stage search) ---\n")
    
    best_aic_idx <- which.min(combined_summary$AIC)
    best_bic_idx <- which.min(combined_summary$BIC)
    cat(sprintf("  AIC best:     γ = %.4g (AIC = %.4f)\n",
                combined_summary$gamma[best_aic_idx], combined_summary$AIC[best_aic_idx]))
    cat(sprintf("  BIC best:     γ = %.4g (BIC = %.4f)\n",
                combined_summary$gamma[best_bic_idx], combined_summary$BIC[best_bic_idx]))
    
    if (compute_gic_flag && "GIC_AIC" %in% colnames(combined_summary)) {
      best_gic_aic_idx <- which.min(combined_summary$GIC_AIC)
      best_gic_bic_idx <- which.min(combined_summary$GIC_BIC)
      cat(sprintf("  GIC_AIC best: γ = %.4g (GIC_AIC = %.4f, df_eff = %.2f)\n",
                  combined_summary$gamma[best_gic_aic_idx],
                  combined_summary$GIC_AIC[best_gic_aic_idx],
                  combined_summary$df_eff[best_gic_aic_idx]))
      cat(sprintf("  GIC_BIC best: γ = %.4g (GIC_BIC = %.4f, df_eff = %.2f)\n",
                  combined_summary$gamma[best_gic_bic_idx],
                  combined_summary$GIC_BIC[best_gic_bic_idx],
                  combined_summary$df_eff[best_gic_bic_idx]))
    }
  }
  
  cat(paste(rep("=", 90), collapse=""), "\n")
  
  return(list(summary = combined_summary, results = combined_results))
}


################################################################################
# Usage examples
################################################################################
# 
# # Usage example with real data:
#
# # source the required files
#
# # --- Gamma search with the SLY model ---
# gamma_search_sly <- compare_gamma(
#   fit_func = fit_sly_penalized,
#   gammas = c(0, 1, 5, 10, 20, 50),
#   data_file = "simulated_data_1111_n100_T5.csv",
#   weight_file = "spatial_weights_n100.csv",
#   time_point = 2,
#   y_vars = c("y1", "y2"),
#    x_vars =  list(
#    y1 = c("x_common1", "x_common2", "x_specific1_1"),
#    y2 = c("x_common1", "x_common2", "x_specific2_1")
#    ),
#   include_time_lag = TRUE
# )
#
# # --- Single run at a specific gamma ---
# result_sly <- fit_sly_penalized(
#   data_file = "simulated_data_1111_n100_T5.csv",
#   weight_file = "spatial_weights_n100.csv",
#   time_point = 2,
#   y_vars = c("y1", "y2"),
#    x_vars =  list(
#    y1 = c("x_common1", "x_common2", "x_specific1_1"),
#    y2 = c("x_common1", "x_common2", "x_specific2_1")
#    ),
#   gamma = 5,
#   verbose = TRUE
# )
#
# # --- Gamma search with the SDEM model ---
# gamma_search_sdem <- compare_gamma(
#   fit_func = fit_sdem_penalized,
#   gammas = c(0, 1, 5, 10, 20, 50),
#   data_file = "simulated_data_1111_n100_T5.csv",
#   weight_file = "spatial_weights_n100.csv",
#   time_point = 2,
#   y_vars = c("y1", "y2"),
#    x_vars =  list(
#    y1 = c("x_common1", "x_common2", "x_specific1_1"),
#    y2 = c("x_common1",  "x_common2", "x_specific2_1")
#    ),
#   include_time_lag = TRUE
# )
# # --- Gamma search with the SEM model ---
# gamma_search_sem <- compare_gamma(
#   fit_func = fit_sem_penalized,
#   gammas = c(0, 1, 5, 10, 20, 50),
#   data_file = "simulated_data_1111_n100_T5.csv",
#   weight_file = "spatial_weights_n100.csv",
#   time_point = 2,
#   y_vars = c("y1", "y2"),
#    x_vars =  list(
#    y1 = c("x_common1", "x_common2",  "x_specific1_1"),
#    y2 = c("x_common1", "x_common2", "x_specific2_1")
#    ),
#   include_time_lag = TRUE
# )
################################################################################

################################################################################
# START OF FILE: fit_sly_diagonal.r
################################################################################

################################################################################
# fit_sly_diagonal.r
# 
# Diagonal-constrained multivariate SLY estimation (d0dd model)
# 
# Estimate y1 and y2 separately with lagsarlm()
# - R matrix: diagonal (ρ₁₁ and ρ₂₂ only, no cross-variable spatial effects)
# - Lambda matrix: 0 (no spatial error dependence)
# - A matrix: diagonal (AR(1) coefficients are independent, no cross-equation lags)
# - Sigma matrix: diagonal (errors are independent across responses)
#
# Usage:
#
#   result <- fit_sly_diagonal(
#     data_file = "simulated_data_d0dd_n400_T5.csv",
#     weight_file = "spatial_weights_n400.csv",
#     y_vars = c("y1", "y2"),
#     x_vars = list(
#       y1 = c("x_common1", "x_common2", "x_specific1_1"),
#       y2 = c("x_common1", "x_common2", "x_specific2_1")
#     ),
#     include_time_lag = TRUE,
#     verbose = TRUE
#   )
#
################################################################################

cat("Loaded fit_sly_diagonal.r (diagonal-constrained SLY estimation, d0dd model)\n")

################################################################################
# Main function: fit_sly_diagonal
################################################################################

#' Diagonal-constrained multivariate SLY estimation (d0dd model)
#' 
#' Estimate each variable individually with lagsarlm() and combine the results
#' 
#' @param data_file  Path to the panel data CSV file
#' @param weight_file  Path to the spatial weight matrix CSV file
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
#' 
fit_sly_diagonal <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained SLY estimation (d0dd model) ===\n")
    cat(sprintf("    temporal lag: %s(no cross terms)\n", ifelse(include_time_lag, "with", "without")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Package availability check
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("spdep package is required: install.packages('spdep')")
  }
  
  k <- length(y_vars)
  
  #############################################################################
  # S0: Data preparation
  #############################################################################
  
  if (verbose) cat("\nS0: Data preparation...\n")
  
  # Load data file from CSV
  data <- read.csv(data_file)
  W_matrix <- as.matrix(read.csv(weight_file, header = TRUE))
  
  # Determine the analysis time point
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  # Extract data for the current period t
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  # Extract data for the previous period t-1 (AR(1) lag variables)
  if (include_time_lag) {
    data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
    if (nrow(data_t_lag) != n) {
      stop("temporal lagdata's number of regions do not match")
    }
  }
  
  # Convert spatial weight matrix to listw format
  W_listw <- spdep::mat2listw(W_matrix, style = "W")
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Analysis time point: t = %d\n", time_point))
    cat(sprintf("  Temporal lag: %s\n", ifelse(include_time_lag, "with (no cross terms)", "without")))
  }
  
  #############################################################################
  # S1-S2: Individual SAR estimation (lagsarlm) for each response variable
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Individual SLY estimation (lagsarlm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Initialise storage for per-variable estimation results
  R_diag <- matrix(0, k, k)      # diagonal R matrix
  Sigma_diag <- matrix(0, k, k)  # diagonal Sigma matrix
  
  individual_models <- list()     # per-variable model results
  beta_list <- list()             # regression-coefficient list
  residuals_list <- list()        # residual list
  fitted_list <- list()           # fitted-value list
  
  total_loglik <- 0               # total log-likelihood
  total_params <- 0               # total number of parameters
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    if (verbose) {
      cat(sprintf("\n--- Estimating %s ---\n", var_name))
    }
    
    # Construct data frame for spatialreg::lagsarlm / errorsarlm
    df_i <- data.frame(
      y = data_t[[var_name]]
    )
    
    # Append the exogenous covariate columns for this response
    for (xvar in xi_vars) {
      df_i[[xvar]] <- data_t[[xvar]]
    }
    
    # Add own-lag AR(1) variable only (no cross-variable lag terms)
    if (include_time_lag) {
      lag_var_name <- paste0(var_name, "_lag")
      df_i[[lag_var_name]] <- data_t_lag[[var_name]]
    }
    
    # Build regression formula string
    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)
    
    # Estimation via lagsarlm() (spatial-lag model)
    tryCatch({
      model_i <- spatialreg::lagsarlm(
        formula = formula_i,
        data = df_i,
        listw = W_listw,
        method = "eigen",
        quiet = !verbose
      )
      
      # Extract estimation results from the model object
      R_diag[i, i] <- model_i$rho       # spatial-lag parameter
      Sigma_diag[i, i] <- model_i$s2    # error variance
      
      # Regression coefficients
      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef
      
      # Extract residuals and fitted values
      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))
      
      # Extract model summary for SE and significance codes
      s <- summary(model_i)
      
      # Log-likelihood and parameter count
      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i
      
      # Spatial parameter significance test
      rho <- model_i$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      
      # Likelihood-ratio test (H_0: spatial parameter = 0)
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # Build regression coefficient table
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Compute R² statistics
      y_i <- df_i$y
      fitted_i <- fitted(model_i)
      residuals_i <- residuals(model_i)
      
      r2_results <- compute_r_squared(
        y = y_i,
        fitted = as.numeric(fitted_i),
        residuals = as.numeric(residuals_i),
        loglik = loglik_i,
        num_params = num_params_i
      )
      
      # Save individual model estimation results
      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = get_signif_code(rho_p),
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model_i$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params_i,
          num_obs = length(y_i)
        )
      )
      
      # Verbose output
      if (verbose) {
        cat("  [Spatial parameters]\n")
        cat(sprintf("    ρ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    rho, rho_se, rho_z, rho_p, get_signif_code(rho_p)))
        
        # Likelihood-ratio test (H_0: spatial parameter = 0)
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "Reject rho=0", "Fail to reject rho=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f\n", model_i$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²     = %.4f\n", r2_results$R2_adj))
        cat(sprintf("    log-likelihood   = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: lagsarlm() failed when estimating variable %s\n", var_name))
      cat("Error message:", e$message, "\n")
      stop("Estimation failed")
    })
  }
  
  #############################################################################
  # S3: Aggregate individual results into a joint model object
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Merge the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Compute information criteria
  ic <- compute_information_criteria(total_loglik, total_params, k * n)
  
  # Vectorise residuals and fitted values: [y1 block; y2 block; ...]
  residuals_raw <- unlist(residuals_list)
  fitted_vals <- unlist(fitted_list)
  
  # Standardised residuals
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_list[[y_vars[i]]] / sqrt(Sigma_diag[i, i])
  }
  
  # Assemble the stacked response vector y
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- data_t[[y_vars[i]]]
  }
  
  # Compute mean R² and mean adjusted R² across responses
  r2_list <- sapply(individual_models, function(m) m$fit$R2)
  r2_adj_list <- sapply(individual_models, function(m) m$fit$R2_adj)
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  # Structure estimated coefficients (T = NULL for SLY)
  coefficients <- list(
    R = R_diag,
    T = NULL,  # no T matrix for SLY
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  # Structure standard errors
  R_se <- matrix(0, k, k)
  for (i in 1:k) {
    var_name <- y_vars[i]
    R_se[i, i] <- individual_models[[var_name]]$spatial_params$rho_se
  }
  
  std_errors <- list(
    R = R_se,
    T = NULL  # no T matrix for SLY
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  #############################################################################
  # Print estimation results
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[diagonal R matrix (spatial lag)]\n")
    print(round(R_diag, 4))
    
    cat("\n[diagonal Sigma matrix (error variance)]\n")
    print(round(Sigma_diag, 4))
    
    cat("\n[goodness-of-fit (combined)]\n")
    cat(sprintf("  total log-likelihood: %.4f\n", total_loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  meanR²: %.4f\n", r2_mean))
    cat(sprintf("  meanAdj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  #############################################################################
  # Build unified result object (standardised output format)
  #############################################################################
  
  # minimal data_list for build_result_object
  dl <- build_data_list_from_parts(
    y = y_vec, W = W_matrix, W_listw = W_listw,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = include_time_lag,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "SLY_diagonal",
    R                 = R_diag,
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "lagsarlm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    std_errors_R      = R_se,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}


################################################################################
# Standard S3 methods
################################################################################

#' Print method for the model object
print.multivar_sly_diagonal <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained SLY (d0dd model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  cat("\n[diagonal R matrix (spatial lag ρ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$R)), collapse=", ")))
  
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  meanR²: %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for detailed output.\n\n")
  
  invisible(x)
}


#' Summary method with detailed inference tables
summary.multivar_sly_diagonal <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained SLY - detailed results (d0dd model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Model specification section
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  # Per-variable estimation results
  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. estimation result of %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    # Spatial parameters
    cat("\nspatial parameters:\n")
    sp <- model_i$spatial_params
    cat(sprintf("  ρ (spatial lag)   = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$rho, sp$rho_se, sp$rho_z, sp$rho_p, sp$rho_signif))
    
    # Likelihood-ratio test (H_0: spatial parameter = 0)
    if (!is.na(sp$LR_statistic)) {
      LR_result <- ifelse(sp$LR_p_value < 0.05, "reject rho=0", "cannot reject rho=0")
      cat(sprintf("\nlikelihood-ratio test:\n"))
      cat(sprintf("  LR statistic = %.2f, p-value = %.4f (%s)\n", 
                  sp$LR_statistic, sp$LR_p_value, LR_result))
    }
    
    # Error variance
    cat(sprintf("\nerror variance: σ² = %.6f\n", model_i$sigma2))
    
    # Regression coefficients
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    # Goodness of fit
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  log-likelihood     = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))
  }
  
  # Combined (joint) estimation results
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. combinedgoodness-of-fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ndiagonal R matrix:\n")
  print(round(x$coefficients$R, digits))
  
  cat("\ndiagonal Sigma matrix:\n")
  print(round(x$coefficients$Sigma, digits))
  
  cat("\ncombinedgoodness-of-fit:\n")
  cat(sprintf("  total log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  meanR²       = %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}


#' Extract estimated coefficients
coef.multivar_sly_diagonal <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(object$coefficients$R)
  } else if (type == "T") {
    return(NULL)  # no T matrix for SLY
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract log-likelihood
logLik.multivar_sly_diagonal <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' Extract AIC
AIC.multivar_sly_diagonal <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' Extract BIC
BIC.multivar_sly_diagonal <- function(object, ...) {
  return(object$fit$BIC)
}


#' Extract residuals
residuals.multivar_sly_diagonal <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract fitted values
fitted.multivar_sly_diagonal <- function(object, ...) {
  
  fitted_vals <- numeric(object$data_info$k * object$data_info$n)
  
  for (i in 1:object$data_info$k) {
    var_name <- object$data_info$y_vars[i]
    idx <- ((i - 1) * object$data_info$n + 1):(i * object$data_info$n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}


################################################################################
# Helper functions (shared across model fitting files) (shared with the other *_diagonal.r)
################################################################################

#' Return significance code string for a p-value
if (!exists("get_signif_code")) {
  get_signif_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  }
}


#' Compute R² statistics
if (!exists("compute_r_squared")) {
  compute_r_squared <- function(y, fitted, residuals, loglik = NULL, num_params = NULL) {
    
    n <- length(y)
    
    # Standard SS-based R²
    ss_res <- sum(residuals^2)
    ss_tot <- sum((y - mean(y))^2)
    R2 <- 1 - ss_res / ss_tot
    
    # Adjusted R²
    if (!is.null(num_params)) {
      R2_adj <- 1 - (ss_res / (n - num_params)) / (ss_tot / (n - 1))
    } else {
      R2_adj <- NA
    }
    
    # Correlation-based R² = corr(y, ŷ)²
    R2_cor <- cor(y, fitted)^2
    
    # Pseudo R² (likelihood-ratio based, McFadden-type)
    if (!is.null(loglik)) {
      ll_null <- sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE))
      R2_pseudo <- 1 - exp(-2/n * (loglik - ll_null))
    } else {
      R2_pseudo <- NA
    }
    
    return(list(
      R2 = R2,
      R2_adj = R2_adj,
      R2_cor = R2_cor,
      R2_pseudo = R2_pseudo
    ))
  }
}


#' Compute AIC and BIC from log-likelihood and parameter count
if (!exists("compute_information_criteria")) {
  compute_information_criteria <- function(loglik, num_params, num_obs) {
    
    AIC <- -2 * loglik + 2 * num_params
    BIC <- -2 * loglik + log(num_obs) * num_params
    
    return(list(
      AIC = AIC,
      BIC = BIC
    ))
  }
}


################################################################################
# Usage examples
################################################################################

cat("
================================================================================
fit_sly_diagonal.r - Usage example
================================================================================

[Basic usage]
#   # source('multivar_sly_phase1_1.R')  # commented out in integrated version
#   # source('phase2_implementation.r')  # commented out in integrated version
#   # source('spatial_core_functions.r')  # commented out in integrated version
#   # source('spatial_output_functions.r')  # commented out in integrated version
#   # source('fit_sly_diagonal.r')  # commented out in integrated version

  result <- fit_sly_diagonal(
    data_file = 'simulated_data_d0dd_n400_T5.csv',
    weight_file = 'spatial_weights_n400.csv',
    y_vars = c('y1', 'y2'),
    x_vars = list(
      y1 = c('x_common1', 'x_common2', 'x_specific1_1'),
      y2 = c('x_common1', 'x_common2', 'x_specific2_1')
    ),
    include_time_lag = TRUE,
    verbose = TRUE
  )

[Checking the result]
  print(result)
  summary(result)

[Extracting coefficients]
  coef(result)
  coef(result, type = 'R')
  coef(result, type = 'beta0')

[goodness-of-fit measures]
  logLik(result)
  AIC(result)
  BIC(result)

[residuals]
  residuals(result)
  residuals(result, type = 'standardized')

[Accessing individual models]
  result$individual_models$y1$spatial_params
  result$individual_models$y2$coefficients

================================================================================
")

################################################################################
# START OF FILE: fit_sem_diagonal.r
################################################################################

################################################################################
# fit_sem_diagonal.r
# 
# Diagonal-constrained multivariate SEM estimation (0ddd model)
# 
# Estimate y1 and y2 separately with errorsarlm()
# - R matrix: 0 (no spatial lag)
# - Lambda matrix: diagonal (λ₁₁ and λ₂₂ only)  
# - A matrix: diagonal (AR(1) coefficients are independent, no cross-equation lags)
# - Sigma matrix: diagonal (errors are independent across responses)
#
# Usage:
#
#   result <- fit_sem_diagonal(
#     data_file = "simulated_data_0ddd_n400_T5.csv",
#     weight_file = "spatial_weights_n400.csv",
#     y_vars = c("y1", "y2"),
#     x_vars = list(
#       y1 = c("x_common1", "x_common2", "x_specific1_1"),
#       y2 = c("x_common1", "x_common2", "x_specific2_1")
#     ),
#     include_time_lag = TRUE,
#     verbose = TRUE
#   )
#
################################################################################

cat("Loaded fit_sem_diagonal.r (diagonal-constrained SEM estimation, 0ddd model)\n")

################################################################################
# Main function: fit_sem_diagonal
################################################################################

#' Diagonal-constrained multivariate SEM estimation (0ddd model)
#' 
#' Estimate each variable individually with errorsarlm() and combine the results
#' 
#' @param data_file  Path to the panel data CSV file
#' @param weight_file  Path to the spatial weight matrix CSV file
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
#' 
fit_sem_diagonal <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained SEM estimation (0ddd model) ===\n")
    cat(sprintf("    temporal lag: %s(no cross terms)\n", ifelse(include_time_lag, "with", "without")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Package availability check
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("spdep package is required: install.packages('spdep')")
  }
  
  k <- length(y_vars)
  
  #############################################################################
  # S0: Data preparation
  #############################################################################
  
  if (verbose) cat("\nS0: Data preparation...\n")
  
  # Load data file from CSV
  data <- read.csv(data_file)
  W_matrix <- as.matrix(read.csv(weight_file, header = TRUE))
  
  # Determine the analysis time point
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  # Extract data for the current period t
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  # Extract data for the previous period t-1 (AR(1) lag variables)
  if (include_time_lag) {
    data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
    if (nrow(data_t_lag) != n) {
      stop("temporal lagdata's number of regions do not match")
    }
  }
  
  # Convert spatial weight matrix to listw format
  W_listw <- spdep::mat2listw(W_matrix, style = "W")
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Analysis time point: t = %d\n", time_point))
    cat(sprintf("  Temporal lag: %s\n", ifelse(include_time_lag, "with (no cross terms)", "without")))
  }
  
  #############################################################################
  # S1-S2: Individual SEM estimation (errorsarlm) for each response variable
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Individual SEM estimation (errorsarlm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Initialise storage for per-variable estimation results
  T_diag <- matrix(0, k, k)      # diagonal T matrix
  Sigma_diag <- matrix(0, k, k)  # diagonal Sigma matrix
  
  individual_models <- list()     # per-variable model results
  beta_list <- list()             # regression-coefficient list
  residuals_list <- list()        # residual list
  fitted_list <- list()           # fitted-value list
  
  total_loglik <- 0               # total log-likelihood
  total_params <- 0               # total number of parameters
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    if (verbose) {
      cat(sprintf("\n--- Estimating %s ---\n", var_name))
    }
    
    # Construct data frame for spatialreg::lagsarlm / errorsarlm
    df_i <- data.frame(
      y = data_t[[var_name]]
    )
    
    # Append the exogenous covariate columns for this response
    for (xvar in xi_vars) {
      df_i[[xvar]] <- data_t[[xvar]]
    }
    
    # Add own-lag AR(1) variable only (no cross-variable lag terms)
    if (include_time_lag) {
      lag_var_name <- paste0(var_name, "_lag")
      df_i[[lag_var_name]] <- data_t_lag[[var_name]]
    }
    
    # Build regression formula string
    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)
    
    # Estimation via errorsarlm() (spatial-error model)
    tryCatch({
      model_i <- spatialreg::errorsarlm(
        formula = formula_i,
        data = df_i,
        listw = W_listw,
        method = "eigen",
        quiet = !verbose
      )
      
      # Extract estimation results from the model object
      T_diag[i, i] <- model_i$lambda    # spatial-error parameter
      Sigma_diag[i, i] <- model_i$s2    # error variance
      
      # Regression coefficients
      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef
      
      # Extract residuals and fitted values
      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))
      
      # Extract model summary for SE and significance codes
      s <- summary(model_i)
      
      # Log-likelihood and parameter count
      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i
      
      # Spatial parameter significance test
      lambda <- model_i$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      
      # Likelihood-ratio test (H_0: spatial parameter = 0)
      LR_stat <- if (!is.null(s$LR1)) s$LR1$statistic else NA
      LR_p <- if (!is.null(s$LR1)) s$LR1$p.value else NA
      
      # Build regression coefficient table
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Compute R² statistics
      y_i <- df_i$y
      fitted_i <- fitted(model_i)
      residuals_i <- residuals(model_i)
      
      r2_results <- compute_r_squared(
        y = y_i,
        fitted = as.numeric(fitted_i),
        residuals = as.numeric(residuals_i),
        loglik = loglik_i,
        num_params = num_params_i
      )
      
      # Save individual model estimation results
      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        spatial_params = list(
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = get_signif_code(lambda_p),
          LR_statistic = LR_stat,
          LR_p_value = LR_p
        ),
        
        coefficients = coef_table,
        sigma2 = model_i$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params_i,
          num_obs = length(y_i)
        )
      )
      
      # Verbose output
      if (verbose) {
        cat("  [Spatial parameters]\n")
        cat(sprintf("    λ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    lambda, lambda_se, lambda_z, lambda_p, get_signif_code(lambda_p)))
        
        # Likelihood-ratio test (H_0: spatial parameter = 0)
        if (!is.na(LR_stat)) {
          LR_result <- ifelse(LR_p < 0.05, "Reject lambda=0", "Fail to reject lambda=0")
          cat("  [Likelihood-ratio test]\n")
          cat(sprintf("    LR statistic = %.2f, p-value = %.4f (%s)\n", LR_stat, LR_p, LR_result))
        }
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f\n", model_i$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²     = %.4f\n", r2_results$R2_adj))
        cat(sprintf("    log-likelihood   = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: errorsarlm() failed when estimating variable %s\n", var_name))
      cat("Error message:", e$message, "\n")
      stop("Estimation failed")
    })
  }
  
  #############################################################################
  # S3: Aggregate individual results into a joint model object
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Merge the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Compute information criteria
  ic <- compute_information_criteria(total_loglik, total_params, k * n)
  
  # Vectorise residuals and fitted values: [y1 block; y2 block; ...]
  residuals_raw <- unlist(residuals_list)
  fitted_vals <- unlist(fitted_list)
  
  # Standardised residuals
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_list[[y_vars[i]]] / sqrt(Sigma_diag[i, i])
  }
  
  # Assemble the stacked response vector y
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- data_t[[y_vars[i]]]
  }
  
  # Compute mean R² and mean adjusted R² across responses
  r2_list <- sapply(individual_models, function(m) m$fit$R2)
  r2_adj_list <- sapply(individual_models, function(m) m$fit$R2_adj)
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  # Structure the coefficients (R=0 for SEM)
  coefficients <- list(
    R = NULL,  # no R matrix for SEM
    T = T_diag,
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  # Structure standard errors
  T_se <- matrix(0, k, k)
  for (i in 1:k) {
    var_name <- y_vars[i]
    T_se[i, i] <- individual_models[[var_name]]$spatial_params$lambda_se
  }
  
  std_errors <- list(
    R = NULL,  # no R matrix for SEM
    T = T_se
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  #############################################################################
  # Print estimation results
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[diagonal T matrix (spatial error)]\n")
    print(round(T_diag, 4))
    
    cat("\n[diagonal Sigma matrix (error variance)]\n")
    print(round(Sigma_diag, 4))
    
    cat("\n[goodness-of-fit (combined)]\n")
    cat(sprintf("  total log-likelihood: %.4f\n", total_loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  meanR²: %.4f\n", r2_mean))
    cat(sprintf("  meanAdj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  #############################################################################
  # Build unified result object (standardised output format)
  #############################################################################
  
  dl <- build_data_list_from_parts(
    y = y_vec, W = W_matrix, W_listw = W_listw,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = include_time_lag,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "SEM_diagonal",
    T_mat             = T_diag,
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "errorsarlm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    std_errors_T      = T_se,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}


################################################################################
# Standard S3 methods
################################################################################

#' Print method for the model object
print.multivar_sem_diagonal <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained SEM (0ddd model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  cat("\n[diagonal T matrix (spatial error λ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$T)), collapse=", ")))
  
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  meanR²: %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for detailed output.\n\n")
  
  invisible(x)
}


#' Summary method with detailed inference tables
summary.multivar_sem_diagonal <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained SEM - detailed results (0ddd model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Model specification section
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  # Per-variable estimation results
  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. estimation result of %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    # Spatial parameters
    cat("\nspatial parameters:\n")
    sp <- model_i$spatial_params
    cat(sprintf("  λ (spatial error)   = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$lambda, sp$lambda_se, sp$lambda_z, sp$lambda_p, sp$lambda_signif))
    
    # Likelihood-ratio test (H_0: spatial parameter = 0)
    if (!is.na(sp$LR_statistic)) {
      LR_result <- ifelse(sp$LR_p_value < 0.05, "reject lambda=0", "cannot reject lambda=0")
      cat(sprintf("\nlikelihood-ratio test:\n"))
      cat(sprintf("  LR statistic = %.2f, p-value = %.4f (%s)\n", 
                  sp$LR_statistic, sp$LR_p_value, LR_result))
    }
    
    # Error variance
    cat(sprintf("\nerror variance: σ² = %.6f\n", model_i$sigma2))
    
    # Regression coefficients
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    # Goodness of fit
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  log-likelihood     = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))
  }
  
  # Combined (joint) estimation results
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. combinedgoodness-of-fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ndiagonal T matrix:\n")
  print(round(x$coefficients$T, digits))
  
  cat("\ndiagonal Sigma matrix:\n")
  print(round(x$coefficients$Sigma, digits))
  
  cat("\ncombinedgoodness-of-fit:\n")
  cat(sprintf("  total log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  meanR²       = %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}


#' Extract estimated coefficients
coef.multivar_sem_diagonal <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(NULL)  # no R matrix for SEM
  } else if (type == "T") {
    return(object$coefficients$T)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract log-likelihood
logLik.multivar_sem_diagonal <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' Extract AIC
AIC.multivar_sem_diagonal <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' Extract BIC
BIC.multivar_sem_diagonal <- function(object, ...) {
  return(object$fit$BIC)
}


#' Extract residuals
residuals.multivar_sem_diagonal <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract fitted values
fitted.multivar_sem_diagonal <- function(object, ...) {
  
  fitted_vals <- numeric(object$data_info$k * object$data_info$n)
  
  for (i in 1:object$data_info$k) {
    var_name <- object$data_info$y_vars[i]
    idx <- ((i - 1) * object$data_info$n + 1):(i * object$data_info$n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}


################################################################################
# Helper functions (shared across model fitting files) (shared with fit_sdem_diagonal.r)
################################################################################

#' Return significance code string for a p-value
if (!exists("get_signif_code")) {
  get_signif_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  }
}


#' Compute R² statistics
if (!exists("compute_r_squared")) {
  compute_r_squared <- function(y, fitted, residuals, loglik = NULL, num_params = NULL) {
    
    n <- length(y)
    
    # Standard SS-based R²
    ss_res <- sum(residuals^2)
    ss_tot <- sum((y - mean(y))^2)
    R2 <- 1 - ss_res / ss_tot
    
    # Adjusted R²
    if (!is.null(num_params)) {
      R2_adj <- 1 - (ss_res / (n - num_params)) / (ss_tot / (n - 1))
    } else {
      R2_adj <- NA
    }
    
    # Correlation-based R² = corr(y, ŷ)²
    R2_cor <- cor(y, fitted)^2
    
    # Pseudo R² (likelihood-ratio based, McFadden-type)
    if (!is.null(loglik)) {
      ll_null <- sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE))
      R2_pseudo <- 1 - exp(-2/n * (loglik - ll_null))
    } else {
      R2_pseudo <- NA
    }
    
    return(list(
      R2 = R2,
      R2_adj = R2_adj,
      R2_cor = R2_cor,
      R2_pseudo = R2_pseudo
    ))
  }
}


#' Compute AIC and BIC from log-likelihood and parameter count
if (!exists("compute_information_criteria")) {
  compute_information_criteria <- function(loglik, num_params, num_obs) {
    
    AIC <- -2 * loglik + 2 * num_params
    BIC <- -2 * loglik + log(num_obs) * num_params
    
    return(list(
      AIC = AIC,
      BIC = BIC
    ))
  }
}


################################################################################
# Usage examples
################################################################################

cat("
================================================================================
fit_sem_diagonal.r - Usage example
================================================================================

[Basic usage]
#   # source('multivar_sly_phase1_1.R')  # commented out in integrated version
#   # source('phase2_implementation.r')  # commented out in integrated version
#   # source('spatial_core_functions.r')  # commented out in integrated version
#   # source('spatial_output_functions.r')  # commented out in integrated version
#   # source('fit_sem_diagonal.r')  # commented out in integrated version

  result <- fit_sem_diagonal(
    data_file = 'simulated_data_0ddd_n400_T5.csv',
    weight_file = 'spatial_weights_n400.csv',
    y_vars = c('y1', 'y2'),
    x_vars = list(
      y1 = c('x_common1', 'x_common2', 'x_specific1_1'),
      y2 = c('x_common1', 'x_common2', 'x_specific2_1')
    ),
    include_time_lag = TRUE,
    verbose = TRUE
  )

[Checking the result]
  print(result)
  summary(result)

[Extracting coefficients]
  coef(result)
  coef(result, type = 'T')
  coef(result, type = 'beta0')

[goodness-of-fit measures]
  logLik(result)
  AIC(result)
  BIC(result)

[residuals]
  residuals(result)
  residuals(result, type = 'standardized')

[Accessing individual models]
  result$individual_models$y1$spatial_params
  result$individual_models$y2$coefficients

================================================================================
")

################################################################################
# START OF FILE: fit_sdem_diagonal.r
################################################################################

################################################################################
# fit_sdem_diagonal.r
# 
# Diagonal-constrained multivariate SDEM estimation (dddd model)
# 
# Estimate y1 and y2 separately with sacsarlm()
# - R matrix: diagonal (ρ₁₁ and ρ₂₂ only, no cross-variable spatial effects)
# - Lambda matrix: diagonal (λ₁₁ and λ₂₂ only)  
# - A matrix: diagonal (AR(1) coefficients are independent, no cross-equation lags)
# - Sigma matrix: diagonal (errors are independent across responses)
#
# Usage:
#
#   result <- fit_sdem_diagonal(
#     data_file = "simulated_data_1111_n400_T5.csv",
#     weight_file = "spatial_weights_n400.csv",
#     y_vars = c("y1", "y2"),
#     x_vars = list(
#       y1 = c("x_common1", "x_common2", "x_specific1_1"),
#       y2 = c("x_common1", "x_common2", "x_specific2_1")
#     ),
#     include_time_lag = TRUE,
#     verbose = TRUE
#   )
#
################################################################################

cat("Loaded fit_sdem_diagonal.r (diagonal-constrained SDEM estimation)\n")

################################################################################
# Main function: fit_sdem_diagonal
################################################################################

#' Diagonal-constrained multivariate SDEM estimation (dddd model)
#' 
#' Estimate each variable individually with sacsarlm() and combine the results
#' 
#' @param data_file  Path to the panel data CSV file
#' @param weight_file  Path to the spatial weight matrix CSV file
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param include_time_lag  Logical; include AR(1) lag columns
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
#' 
fit_sdem_diagonal <- function(
  data_file,
  weight_file,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained SDEM estimation (dddd model) ===\n")
    cat(sprintf("    temporal lag: %s(no cross terms)\n", ifelse(include_time_lag, "with", "without")))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Package availability check
  if (!requireNamespace("spatialreg", quietly = TRUE)) {
    stop("spatialreg package is required: install.packages('spatialreg')")
  }
  if (!requireNamespace("spdep", quietly = TRUE)) {
    stop("spdep package is required: install.packages('spdep')")
  }
  
  k <- length(y_vars)
  
  #############################################################################
  # S0: Data preparation
  #############################################################################
  
  if (verbose) cat("\nS0: Data preparation...\n")
  
  # Load data file from CSV
  data <- read.csv(data_file)
  W_matrix <- as.matrix(read.csv(weight_file, header = TRUE))
  
  # Determine the analysis time point
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  # Extract data for the current period t
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  # Extract data for the previous period t-1 (AR(1) lag variables)
  if (include_time_lag) {
    data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
    if (nrow(data_t_lag) != n) {
      stop("temporal lagdata's number of regions do not match")
    }
  }
  
  # Convert spatial weight matrix to listw format
  W_listw <- spdep::mat2listw(W_matrix, style = "W")
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Analysis time point: t = %d\n", time_point))
    cat(sprintf("  Temporal lag: %s\n", ifelse(include_time_lag, "with (no cross terms)", "without")))
  }
  
  #############################################################################
  # S1-S2: estimate each variable individually with sacsarlm()
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Individual SDEM estimation (sacsarlm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Initialise storage for per-variable estimation results
  R_diag <- matrix(0, k, k)      # diagonal R matrix
  T_diag <- matrix(0, k, k)      # diagonal T matrix
  Sigma_diag <- matrix(0, k, k)  # diagonal Sigma matrix
  
  individual_models <- list()     # per-variable model results
  beta_list <- list()             # regression-coefficient list
  residuals_list <- list()        # residual list
  fitted_list <- list()           # fitted-value list
  
  total_loglik <- 0               # total log-likelihood
  total_params <- 0               # total number of parameters
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    if (verbose) {
      cat(sprintf("\n--- Estimating %s ---\n", var_name))
    }
    
    # Construct data frame for spatialreg::lagsarlm / errorsarlm
    df_i <- data.frame(
      y = data_t[[var_name]]
    )
    
    # Append the exogenous covariate columns for this response
    for (xvar in xi_vars) {
      df_i[[xvar]] <- data_t[[xvar]]
    }
    
    # Add own-lag AR(1) variable only (no cross-variable lag terms)
    if (include_time_lag) {
      lag_var_name <- paste0(var_name, "_lag")
      df_i[[lag_var_name]] <- data_t_lag[[var_name]]
    }
    
    # Build regression formula string
    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)
    
    # Estimation via sacsarlm()
    tryCatch({
      model_i <- spatialreg::sacsarlm(
        formula = formula_i,
        data = df_i,
        listw = W_listw,
        method = "eigen",
        quiet = !verbose
      )
      
      # Extract estimation results from the model object
      R_diag[i, i] <- model_i$rho       # spatial-lag parameter
      T_diag[i, i] <- model_i$lambda    # spatial-error parameter
      Sigma_diag[i, i] <- model_i$s2    # error variance
      
      # Regression coefficients
      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef
      
      # Extract residuals and fitted values
      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))
      
      # Extract model summary for SE and significance codes
      s <- summary(model_i)
      
      # Log-likelihood and parameter count
      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  # add only sigma^2 (spatial params are included in coef)
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i
      
      # Spatial parameter significance test
      rho <- model_i$rho
      rho_se <- s$rho.se
      rho_z <- rho / rho_se
      rho_p <- 2 * pnorm(-abs(rho_z))
      
      lambda <- model_i$lambda
      lambda_se <- s$lambda.se
      lambda_z <- lambda / lambda_se
      lambda_p <- 2 * pnorm(-abs(lambda_z))
      
      # Build regression coefficient table
      coef_mat <- s$Coef
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, 1],
        std_error = coef_mat[, 2],
        z_value = coef_mat[, 3],
        p_value = coef_mat[, 4],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Compute R² statistics
      y_i <- df_i$y
      fitted_i <- fitted(model_i)
      residuals_i <- residuals(model_i)
      
      r2_results <- compute_r_squared(
        y = y_i,
        fitted = as.numeric(fitted_i),
        residuals = as.numeric(residuals_i),
        loglik = loglik_i,
        num_params = num_params_i
      )
      
      # Save individual model estimation results
      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        spatial_params = list(
          rho = rho,
          rho_se = rho_se,
          rho_z = rho_z,
          rho_p = rho_p,
          rho_signif = get_signif_code(rho_p),
          lambda = lambda,
          lambda_se = lambda_se,
          lambda_z = lambda_z,
          lambda_p = lambda_p,
          lambda_signif = get_signif_code(lambda_p)
        ),
        
        coefficients = coef_table,
        sigma2 = model_i$s2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = r2_results$R2,
          R2_adj = r2_results$R2_adj,
          R2_cor = r2_results$R2_cor,
          R2_pseudo = r2_results$R2_pseudo,
          num_params = num_params_i,
          num_obs = length(y_i)
        )
      )
      
      # Verbose output
      if (verbose) {
        cat("  [Spatial parameters]\n")
        cat(sprintf("    ρ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    rho, rho_se, rho_z, rho_p, get_signif_code(rho_p)))
        cat(sprintf("    λ = %.4f  (SE: %.4f, z: %.2f, p: %.4f %s)\n",
                    lambda, lambda_se, lambda_z, lambda_p, get_signif_code(lambda_p)))
        
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f\n", model_i$s2))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$z_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", r2_results$R2))
        cat(sprintf("    Adj.R²     = %.4f\n", r2_results$R2_adj))
        cat(sprintf("    log-likelihood   = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: sacsarlm() failed when estimating variable %s\n", var_name))
      cat("Error message:", e$message, "\n")
      stop("Estimation failed")
    })
  }
  
  #############################################################################
  # S3: Aggregate individual results into a joint model object
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Merge the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Compute information criteria
  ic <- compute_information_criteria(total_loglik, total_params, k * n)
  
  # Vectorise residuals and fitted values: [y1 block; y2 block; ...]
  residuals_raw <- unlist(residuals_list)
  fitted_vals <- unlist(fitted_list)
  
  # Standardised residuals
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_list[[y_vars[i]]] / sqrt(Sigma_diag[i, i])
  }
  
  # Assemble the stacked response vector y
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- data_t[[y_vars[i]]]
  }
  
  # Compute mean R² and mean adjusted R² across responses
  r2_list <- sapply(individual_models, function(m) m$fit$R2)
  r2_adj_list <- sapply(individual_models, function(m) m$fit$R2_adj)
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  # Structure the coefficients
  coefficients <- list(
    R = R_diag,
    T = T_diag,
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  # Structure standard errors
  R_se <- matrix(0, k, k)
  T_se <- matrix(0, k, k)
  for (i in 1:k) {
    var_name <- y_vars[i]
    R_se[i, i] <- individual_models[[var_name]]$spatial_params$rho_se
    T_se[i, i] <- individual_models[[var_name]]$spatial_params$lambda_se
  }
  
  std_errors <- list(
    R = R_se,
    T = T_se
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  #############################################################################
  # Print estimation results
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[diagonal R matrix (spatial lag)]\n")
    print(round(R_diag, 4))
    
    cat("\n[diagonal T matrix (spatial error)]\n")
    print(round(T_diag, 4))
    
    cat("\n[diagonal Sigma matrix (error variance)]\n")
    print(round(Sigma_diag, 4))
    
    cat("\n[goodness-of-fit (combined)]\n")
    cat(sprintf("  total log-likelihood: %.4f\n", total_loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  meanR²: %.4f\n", r2_mean))
    cat(sprintf("  meanAdj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  #############################################################################
  # Build unified result object (standardised output format)
  #############################################################################
  
  dl <- build_data_list_from_parts(
    y = y_vec, W = W_matrix, W_listw = W_listw,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = include_time_lag,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "SDEM_diagonal",
    R                 = R_diag,
    T_mat             = T_diag,
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "sacsarlm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    std_errors_R      = R_se,
    std_errors_T      = T_se,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}


################################################################################
# Standard S3 methods
################################################################################

#' Print method for the model object
print.multivar_sdem_diagonal <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained SDEM (dddd model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  cat("\n[diagonal R matrix (spatial lag ρ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$R)), collapse=", ")))
  
  cat("\n[diagonal T matrix (spatial error λ)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$T)), collapse=", ")))
  
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  meanR²: %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for detailed output.\n\n")
  
  invisible(x)
}


#' Summary method with detailed inference tables
summary.multivar_sdem_diagonal <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained SDEM - detailed results (dddd model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Model specification section
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  cat(sprintf("  temporal lag: %s\n", ifelse(x$data_info$include_time_lag, "with (no cross terms)", "without")))
  
  # Per-variable estimation results
  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. estimation result of %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    # Spatial parameters
    cat("\nspatial parameters:\n")
    sp <- model_i$spatial_params
    cat(sprintf("  ρ (spatial lag)   = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$rho, sp$rho_se, sp$rho_z, sp$rho_p, sp$rho_signif))
    cat(sprintf("  λ (spatial error)   = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                sp$lambda, sp$lambda_se, sp$lambda_z, sp$lambda_p, sp$lambda_signif))
    
    # Error variance
    cat(sprintf("\nerror variance: σ² = %.6f\n", model_i$sigma2))
    
    # Regression coefficients
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    # Goodness of fit
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  log-likelihood     = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))
  }
  
  # Combined (joint) estimation results
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. combinedgoodness-of-fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ndiagonal R matrix:\n")
  print(round(x$coefficients$R, digits))
  
  cat("\ndiagonal T matrix:\n")
  print(round(x$coefficients$T, digits))
  
  cat("\ndiagonal Sigma matrix:\n")
  print(round(x$coefficients$Sigma, digits))
  
  cat("\ncombinedgoodness-of-fit:\n")
  cat(sprintf("  total log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  meanR²       = %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}


#' Extract estimated coefficients
coef.multivar_sdem_diagonal <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(object$coefficients$R)
  } else if (type == "T") {
    return(object$coefficients$T)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract log-likelihood
logLik.multivar_sdem_diagonal <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' Extract AIC
AIC.multivar_sdem_diagonal <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' Extract BIC
BIC.multivar_sdem_diagonal <- function(object, ...) {
  return(object$fit$BIC)
}


#' Extract residuals
residuals.multivar_sdem_diagonal <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract fitted values
fitted.multivar_sdem_diagonal <- function(object, ...) {
  
  fitted_vals <- numeric(object$data_info$k * object$data_info$n)
  
  for (i in 1:object$data_info$k) {
    var_name <- object$data_info$y_vars[i]
    idx <- ((i - 1) * object$data_info$n + 1):(i * object$data_info$n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}


################################################################################
# Helper functions (shared across model fitting files)
################################################################################

#' Return significance code string for a p-value
get_signif_code <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  if (p < 0.1) return(".")
  return("")
}


#' Compute R² statistics
compute_r_squared <- function(y, fitted, residuals, loglik = NULL, num_params = NULL) {
  
  n <- length(y)
  
  # Standard SS-based R²
  ss_res <- sum(residuals^2)
  ss_tot <- sum((y - mean(y))^2)
  R2 <- 1 - ss_res / ss_tot
  
  # Adjusted R²
  if (!is.null(num_params)) {
    R2_adj <- 1 - (ss_res / (n - num_params)) / (ss_tot / (n - 1))
  } else {
    R2_adj <- NA
  }
  
  # Correlation-based R² = corr(y, ŷ)²
  R2_cor <- cor(y, fitted)^2
  
  # Pseudo R² (likelihood-ratio based, McFadden-type)
  if (!is.null(loglik)) {
    ll_null <- sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE))
    R2_pseudo <- 1 - exp(-2/n * (loglik - ll_null))
  } else {
    R2_pseudo <- NA
  }
  
  return(list(
    R2 = R2,
    R2_adj = R2_adj,
    R2_cor = R2_cor,
    R2_pseudo = R2_pseudo
  ))
}


#' Compute AIC and BIC from log-likelihood and parameter count
compute_information_criteria <- function(loglik, num_params, num_obs) {
  
  AIC <- -2 * loglik + 2 * num_params
  BIC <- -2 * loglik + log(num_obs) * num_params
  
  return(list(
    AIC = AIC,
    BIC = BIC
  ))
}


################################################################################
# Usage examples
################################################################################

# cat("
# ================================================================================
# Usage examples
# ================================================================================

# [Basic usage]
#   source('multivar_sly_phase1_1.R')
#   source('phase2_implementation.r')
#   source('spatial_core_functions.r')
#   source('spatial_output_functions.r')
#   source('fit_sdem_diagonal.r')

#   result <- fit_sdem_diagonal(
#     data_file = 'simulated_data_dddd_n400_T5.csv',
#     weight_file = 'spatial_weights_n400.csv',
#     y_vars = c('y1', 'y2'),
#     x_vars = list(
#       y1 = c('x_common1', 'x_common2', 'x_specific1_1'),
#       y2 = c('x_common1', 'x_common2', 'x_specific2_1')
#     ),
#     include_time_lag = TRUE,
#     verbose = TRUE
#   )

# [Check the result]
#   print(result)
#   summary(result)

# [Coefficient extraction]
#   coef(result)
#   coef(result, type = 'R')
#   coef(result, type = 'T')
#   coef(result, type = 'beta0')

# [Goodness-of-fit measures]
#   logLik(result)
#   AIC(result)
#   BIC(result)

# [Residuals]
#   residuals(result)
#   residuals(result, type = 'standardized')

# [Accessing individual models]
#   result$individual_models$y1$spatial_params
#   result$individual_models$y2$coefficients

# ================================================================================
# ")

################################################################################
# START OF FILE: fit_varx.r
################################################################################

################################################################################
# fit_varx.r
# 
# VARX model estimation (0011 model)
# 
# no spatial correlation; temporal-lag cross terms present; error correlation present
# Uses SUR estimation from the systemfit package
#
# Model equation:
#   y1_t = α11 * y1_{t-1} + α12 * y2_{t-1} + X1 * β1 + ε1
#   y2_t = α21 * y1_{t-1} + α22 * y2_{t-1} + X2 * β2 + ε2
#   
#   ε ~ N(0, Σ)  where Σ is non-diagonal (correlated errors)
#
# Parameters:
#   - R matrix: 0 (no spatial lag)
#   - Lambda matrix: 0 (no spatial error)
#   - A matrix: off-diagonal (α₁₁, α₁₂, α₂₁, α₂₂ — full K²)
#   - Sigma matrix: full off-diagonal (cross-equation error correlations)
#
# Usage:
#
#   result <- fit_varx(
#     data_file = "simulated_data_0011_n400_T5.csv",
#     weight_file = "spatial_weights_n400.csv",
#     y_vars = c("y1", "y2"),
#     x_vars = list(
#       y1 = c("x_common1", "x_common2", "x_specific1_1"),
#       y2 = c("x_common1", "x_common2", "x_specific2_1")
#     ),
#     verbose = TRUE
#   )
#
################################################################################

cat("Loaded fit_varx.r (VARX model estimation, 0011 model)\n")

################################################################################
# Main function: fit_varx
################################################################################

#' VARX model estimation (0011 model)
#' 
#' Uses SUR estimation from the systemfit package
#' 
#' @param data_file  Path to the panel data CSV file
#' @param weight_file  Path to the spatial weight matrix CSV file
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param method estimation method ("SUR", "OLS", "WLS", "3SLS")
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
#' 
fit_varx <- function(
  data_file,
  weight_file = NULL,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  method = "SUR",
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== VARX model estimation (0011 model) ===\n")
    cat(sprintf("    estimation method: %s\n", method))
    cat(sprintf("    Temporal lag: with cross terms\n"))
    cat(sprintf("    Error correlation: yes\n"))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  # Package availability check
  if (!requireNamespace("systemfit", quietly = TRUE)) {
    stop("The systemfit package is required: install.packages('systemfit')")
  }
  
  k <- length(y_vars)
  
  #############################################################################
  # S0: Data preparation
  #############################################################################
  
  if (verbose) cat("\nS0: Data preparation...\n")
  
  # Load data file from CSV
  data <- read.csv(data_file)
  
  # Determine the analysis time point
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  # Extract data for the current period t
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  # Extract data for the previous period t-1 (AR(1) lag variables)
  data_t_lag <- data[data[[time_var]] == (time_point - 1), ]
  if (nrow(data_t_lag) != n) {
    stop("temporal lagdata's number of regions do not match")
  }
  
  # Sort rows by region ID to ensure consistent ordering
  if (region_var %in% names(data_t)) {
    data_t <- data_t[order(data_t[[region_var]]), ]
    data_t_lag <- data_t_lag[order(data_t_lag[[region_var]]), ]
  }
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Analysis time point: t = %d\n", time_point))
  }
  
  #############################################################################
  # Construct data frame for spatialreg estimation
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Create the data frame (for SUR)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Construct data frame for spatialreg estimation
  df_sur <- data.frame(row.names = 1:n)
  
  # Add the response variable
  for (i in 1:k) {
    var_name <- y_vars[i]
    df_sur[[var_name]] <- data_t[[var_name]]
  }
  
  # Add temporal-lag variables (lags of all variables)
  for (i in 1:k) {
    var_name <- y_vars[i]
    lag_name <- paste0(var_name, "_lag")
    df_sur[[lag_name]] <- data_t_lag[[var_name]]
  }
  
  # Add exogenous regressors (variables specific to each equation)
  all_x_vars <- unique(unlist(x_vars))
  for (xvar in all_x_vars) {
    df_sur[[xvar]] <- data_t[[xvar]]
  }
  
  if (verbose) {
    cat(sprintf("  Data frame: %d rows x %d columns\n", nrow(df_sur), ncol(df_sur)))
    cat(sprintf("  Variables: %s\n", paste(names(df_sur), collapse=", ")))
  }
  
  #############################################################################
  # S2: Building the formula
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Build the formula\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  formula_list <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    # Temporal-lag variables (including lags of all y variables)
    lag_vars <- paste0(y_vars, "_lag")
    
    # Right-hand-side variables
    rhs_vars <- c(lag_vars, xi_vars)
    
    # Build the formula
    if (include_intercept) {
      formula_str <- paste(var_name, "~", paste(rhs_vars, collapse = " + "))
    } else {
      formula_str <- paste(var_name, "~ -1 +", paste(rhs_vars, collapse = " + "))
    }
    
    formula_list[[var_name]] <- as.formula(formula_str)
    
    if (verbose) {
      cat(sprintf("  %s: %s\n", var_name, formula_str))
    }
  }
  
  #############################################################################
  # S3: SUR estimation
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("SUR estimation (%s)\n", method))
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Estimation via systemfit
  tryCatch({
    sur_result <- systemfit::systemfit(
      formula = formula_list,
      method = method,
      data = df_sur
    )
    
    if (verbose) {
      cat("  ✓ estimation complete\n")
      cat("  Class of the result object:", class(sur_result), "\n")
      cat("  Elements of the result object:", paste(names(sur_result), collapse=", "), "\n")
    }
    
  }, error = function(e) {
    cat(sprintf("\nError: systemfit() failed\n"))
    cat("Error message:", e$message, "\n")
    stop("Estimation failed")
  })
  
  #############################################################################
  # S4: result extraction
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Extracting results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # residuals() and fitted() return a data.frame
  res_df <- residuals(sur_result)
  fit_df <- fitted(sur_result)
  
  # Convert the data.frame to a matrix
  residuals_mat <- as.matrix(res_df)
  fitted_mat <- as.matrix(fit_df)
  
  # Retrieve column names (check they are in the same order as y_vars)
  res_colnames <- colnames(residuals_mat)
  
  if (verbose) {
    cat(sprintf("  residual matrix: %d x %d\n", nrow(residuals_mat), ncol(residuals_mat)))
    cat(sprintf("  Column names: %s\n", paste(res_colnames, collapse=", ")))
  }
  
  # Reorder if the column names do not match y_vars
  if (!is.null(res_colnames) && all(y_vars %in% res_colnames)) {
    residuals_mat <- residuals_mat[, y_vars, drop = FALSE]
    fitted_mat <- fitted_mat[, y_vars, drop = FALSE]
  }
  
  # Estimate the error covariance matrix Sigma
  Sigma <- (t(residuals_mat) %*% residuals_mat) / n
  rownames(Sigma) <- y_vars
  colnames(Sigma) <- y_vars
  
  # log-likelihood
  log_det_Sigma <- determinant(Sigma, logarithm = TRUE)$modulus[1]
  loglik <- -0.5 * k * n * log(2 * pi) - 0.5 * n * log_det_Sigma - 0.5 * n * k
  
  # number of parameters computation
  num_params_per_eq <- sapply(1:k, function(i) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    # intercept (if requested) + k temporal-lag terms + explanatory variables
    n_params <- ifelse(include_intercept, 1, 0) + k + length(xi_vars)
    return(n_params)
  })
  num_beta <- sum(num_params_per_eq)
  num_sigma <- k * (k + 1) / 2  # number of independent parameters of Sigma
  total_params <- num_beta + num_sigma
  
  # information criterion
  ic <- compute_information_criteria(loglik, total_params, k * n)
  
  if (verbose) cat("  information criterion computed\n")
  
  # A matrix (temporal-lag coefficients) extraction
  # Retrieve coefficients from summary
  sur_summary <- summary(sur_result)
  
  if (verbose) cat("  summary retrieved\n")
  
  A_mat <- matrix(0, k, k)
  A_se <- matrix(0, k, k)
  A_z <- matrix(0, k, k)
  A_p <- matrix(0, k, k)
  rownames(A_mat) <- y_vars
  colnames(A_mat) <- y_vars
  rownames(A_se) <- y_vars
  colnames(A_se) <- y_vars
  
  # Retrieve all coefficients
  # coef(sur_summary) may return a list of summary.systemfit.equation objects
  # Use sur_result$coefficients directly
  all_coef_vec <- coef(sur_result)
  all_se_vec <- sqrt(diag(vcov(sur_result)))
  
  if (verbose) {
    cat("  Length of the coefficient vector:", length(all_coef_vec), "\n")
    cat("  coefficient names:\n")
    print(names(all_coef_vec))
  }
  
  # Extract the A matrix from the coefficient names
  coef_names <- names(all_coef_vec)
  
  for (i in 1:k) {
    eq_prefix <- y_vars[i]
    
    for (j in 1:k) {
      lag_name <- paste0(y_vars[j], "_lag")
      # systemfit coefficient names have the form "eqname_varname"
      coef_name <- paste0(eq_prefix, "_", lag_name)
      
      if (coef_name %in% coef_names) {
        A_mat[i, j] <- all_coef_vec[coef_name]
        A_se[i, j] <- all_se_vec[coef_name]
        A_z[i, j] <- all_coef_vec[coef_name] / all_se_vec[coef_name]
        A_p[i, j] <- 2 * pnorm(-abs(A_z[i, j]))
      }
    }
  }
  
  if (verbose) cat("  A-matrix extraction complete\n")
  
  # regression coefficients (beta) extraction
  beta_list <- list()
  coef_tables <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    eq_prefix <- var_name
    
    # Extract this equation's coefficients
    eq_coef_idx <- grep(paste0("^", eq_prefix, "_"), coef_names)
    eq_coef_names <- coef_names[eq_coef_idx]
    eq_coef_vals <- all_coef_vec[eq_coef_idx]
    eq_se_vals <- all_se_vec[eq_coef_idx]
    
    # Remove the equation prefix from the coefficient names
    short_names <- sub(paste0("^", eq_prefix, "_"), "", eq_coef_names)
    
    # Extract the non-temporal-lag coefficients
    lag_names <- paste0(y_vars, "_lag")
    non_lag_mask <- !short_names %in% lag_names
    
    beta_coef <- eq_coef_vals[non_lag_mask]
    names(beta_coef) <- short_names[non_lag_mask]
    beta_list[[var_name]] <- beta_coef
    
    # Compute z-values and p-values
    z_vals <- eq_coef_vals / eq_se_vals
    p_vals <- 2 * pnorm(-abs(z_vals))
    
    # Coefficient table (overall)
    coef_table <- data.frame(
      parameter = short_names,
      estimate = eq_coef_vals,
      std_error = eq_se_vals,
      z_value = z_vals,
      p_value = p_vals,
      stringsAsFactors = FALSE
    )
    coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
    rownames(coef_table) <- NULL
    coef_tables[[var_name]] <- coef_table
  }
  
  if (verbose) cat("  regression-coefficient extraction complete\n")
  
  # Residuals and fitted values (residuals_mat and fitted_mat already created)
  residuals_raw <- as.vector(residuals_mat)
  fitted_vals <- as.vector(fitted_mat)
  
  # Standardised residuals
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_mat[, i] / sqrt(Sigma[i, i])
  }
  
  # y vector
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- df_sur[[y_vars[i]]]
  }
  
  # Compute R² statistics
  r2_list <- numeric(k)
  r2_adj_list <- numeric(k)
  names(r2_list) <- y_vars
  names(r2_adj_list) <- y_vars
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    y_i <- df_sur[[var_name]]
    fitted_i <- fitted_mat[, i]
    residuals_i <- residuals_mat[, i]
    
    r2_results <- compute_r_squared(
      y = y_i,
      fitted = fitted_i,
      residuals = residuals_i,
      loglik = NULL,
      num_params = num_params_per_eq[i]
    )
    
    r2_list[i] <- r2_results$R2
    r2_adj_list[i] <- r2_results$R2_adj
  }
  
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  #############################################################################
  # S5: construction of individual-model results
  #############################################################################
  
  individual_models <- list()
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    
    individual_models[[var_name]] <- list(
      # Temporal AR(1) coefficient matrix A
      time_lag_params = list(
        alpha = A_mat[i, ],
        alpha_se = A_se[i, ],
        alpha_z = A_z[i, ],
        alpha_p = A_p[i, ],
        alpha_signif = sapply(A_p[i, ], get_signif_code)
      ),
      
      # Build regression coefficient table
      coefficients = coef_tables[[var_name]],
      
      # Error variance
      sigma2 = Sigma[i, i],
      
      # Extract residuals and fitted values
      residuals = residuals_mat[, i],
      fitted = fitted_mat[, i],
      
      # Goodness of fit
      fit = list(
        R2 = r2_list[i],
        R2_adj = r2_adj_list[i],
        num_params = num_params_per_eq[i],
        num_obs = n
      )
    )
  }
  
  #############################################################################
  # Print estimation results
  #############################################################################
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[Amatrix (temporal-lag coefficient)]\n")
    print(round(A_mat, 4))
    
    cat("\n[Amatrix's standard error]\n")
    print(round(A_se, 4))
    
    cat("\n[p-values of the A matrix]\n")
    print(round(A_p, 4))
    
    cat("\n[Σmatrix (error covariance)]\n")
    print(round(Sigma, 6))
    
    # error correlation matrix
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    rownames(corr_mat) <- y_vars
    colnames(corr_mat) <- y_vars
    cat("\n[error correlation matrix]\n")
    print(round(corr_mat, 4))
    
    cat("\n[goodness-of-fit (combined)]\n")
    cat(sprintf("  log-likelihood: %.4f\n", loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  meanR²: %.4f\n", r2_mean))
    cat(sprintf("  meanAdj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  #############################################################################
  # Build unified result object (standardised output format)
  #############################################################################
  
  dl <- build_data_list_from_parts(
    y = y_vec, W = NULL, W_listw = NULL,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = TRUE,
    include_intercept = include_intercept,
    region_var = region_var, time_var = time_var
  )
  
  result <- build_result_object(
    model_type        = "VARX",
    Sigma             = Sigma,
    loglik            = loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = method,
    data_list         = dl,
    beta0             = beta_list,
    alpha             = A_mat,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  # Add VARX-specific fields
  result$coefficients$A <- A_mat
  result$std_errors$A <- A_se
  result$inference <- list(A_z = A_z, A_p = A_p)
  result$sur_result <- sur_result
  result$data_info$time_lag_cross <- TRUE
  
  return(result)
}


################################################################################
# Standard S3 methods
################################################################################

#' Print method for the model object
print.multivar_varx <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("VARX model (0011 model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Temporal lag: with cross terms\n"))
  cat(sprintf("  Error correlation: yes\n"))
  
  cat("\n[Amatrix (temporal-lag coefficient)]\n")
  print(round(x$coefficients$A, digits))
  
  cat("\n[Σmatrix (error covariance)diagonal elements]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$Sigma)), collapse=", ")))
  
  # error correlation
  Sigma <- x$coefficients$Sigma
  if (nrow(Sigma) > 1) {
    D <- diag(1/sqrt(diag(Sigma)))
    corr_mat <- D %*% Sigma %*% D
    cat(sprintf("  error correlation(1,2): %.4f\n", corr_mat[1,2]))
  }
  
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  meanR²: %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for detailed output.\n\n")
  
  invisible(x)
}


#' Summary method with detailed inference tables
summary.multivar_varx <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("VARX model - detailed results (0011 model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Model specification section
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  cat(sprintf("  estimation method: %s\n", x$convergence$method))
  
  # Temporal AR(1) coefficient matrix A (Amatrix)
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[2. temporal-lag coefficient (Amatrix)]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ncoefficient matrix:\n")
  print(round(x$coefficients$A, digits))
  
  cat("\nstandard error:\n")
  print(round(x$std_errors$A, digits))
  
  cat("\np-value:\n")
  print(round(x$inference$A_p, digits))
  
  cat("\ndetails (rows: response variables, columns: lag variables):\n")
  k <- x$data_info$k
  y_vars <- x$data_info$y_vars
  
  for (i in 1:k) {
    for (j in 1:k) {
      cat(sprintf("  α[%s,%s_lag] = %7.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  y_vars[i], y_vars[j],
                  x$coefficients$A[i, j],
                  x$std_errors$A[i, j],
                  x$inference$A_z[i, j],
                  x$inference$A_p[i, j],
                  get_signif_code(x$inference$A_p[i, j])))
    }
  }
  
  # Details of each equation
  for (i in 1:k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. details of the %s equation]\n", i + 2, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, z: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$z_value, row$p_value, row$signif))
    }
    
    cat(sprintf("\nerror variance: σ² = %.6f\n", model_i$sigma2))
    
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²       = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²   = %.4f\n", model_i$fit$R2_adj))
  }
  
  # Error covariance matrix Σ
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. error covariance matrix Σ]\n", k + 3))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ncovariance matrix:\n")
  print(round(x$coefficients$Sigma, digits + 2))
  
  # correlation matrix
  Sigma <- x$coefficients$Sigma
  D <- diag(1/sqrt(diag(Sigma)))
  corr_mat <- D %*% Sigma %*% D
  rownames(corr_mat) <- y_vars
  colnames(corr_mat) <- y_vars
  cat("\ncorrelation matrix:\n")
  print(round(corr_mat, digits))
  
  # Goodness of fit
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. combinedgoodness-of-fit]\n", k + 4))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("\n  log-likelihood     = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  meanR²       = %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}


#' Extract estimated coefficients
coef.multivar_varx <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(NULL)
  } else if (type == "T") {
    return(NULL)
  } else if (type == "A") {
    return(object$coefficients$A)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract log-likelihood
logLik.multivar_varx <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' Extract AIC
AIC.multivar_varx <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' Extract BIC
BIC.multivar_varx <- function(object, ...) {
  return(object$fit$BIC)
}


#' Extract residuals
residuals.multivar_varx <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract fitted values
fitted.multivar_varx <- function(object, ...) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars
  
  fitted_vals <- numeric(k * n)
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    idx <- ((i - 1) * n + 1):(i * n)
    fitted_vals[idx] <- object$individual_models[[var_name]]$fitted
  }
  
  return(fitted_vals)
}


################################################################################
# Helper functions (shared across model fitting files)
################################################################################

#' Return significance code string for a p-value
if (!exists("get_signif_code")) {
  get_signif_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  }
}


#' Compute R² statistics
if (!exists("compute_r_squared")) {
  compute_r_squared <- function(y, fitted, residuals, loglik = NULL, num_params = NULL) {
    
    n <- length(y)
    
    ss_res <- sum(residuals^2)
    ss_tot <- sum((y - mean(y))^2)
    R2 <- 1 - ss_res / ss_tot
    
    if (!is.null(num_params)) {
      R2_adj <- 1 - (ss_res / (n - num_params)) / (ss_tot / (n - 1))
    } else {
      R2_adj <- NA
    }
    
    R2_cor <- cor(y, fitted)^2
    
    if (!is.null(loglik)) {
      ll_null <- sum(dnorm(y, mean = mean(y), sd = sd(y), log = TRUE))
      R2_pseudo <- 1 - exp(-2/n * (loglik - ll_null))
    } else {
      R2_pseudo <- NA
    }
    
    return(list(
      R2 = R2,
      R2_adj = R2_adj,
      R2_cor = R2_cor,
      R2_pseudo = R2_pseudo
    ))
  }
}


#' Compute AIC and BIC from log-likelihood and parameter count
if (!exists("compute_information_criteria")) {
  compute_information_criteria <- function(loglik, num_params, num_obs) {
    
    AIC <- -2 * loglik + 2 * num_params
    BIC <- -2 * loglik + log(num_obs) * num_params
    
    return(list(
      AIC = AIC,
      BIC = BIC
    ))
  }
}


################################################################################
# Usage examples
################################################################################

cat("
================================================================================
fit_varx.r - Usage example
================================================================================

[Basic usage]
#   # source('fit_varx.r')  # commented out in integrated version

  result <- fit_varx(
    data_file = 'simulated_data_0011_n400_T5.csv',
    weight_file = 'spatial_weights_n400.csv',
    y_vars = c('y1', 'y2'),
    x_vars = list(
      y1 = c('x_common1', 'x_common2', 'x_specific1_1'),
      y2 = c('x_common1', 'x_common2', 'x_specific2_1')
    ),
    verbose = TRUE
  )

[Checking the result]
  print(result)
  summary(result)

[Extracting coefficients]
  coef(result)
  coef(result, type = 'A')       # Temporal AR(1) coefficient matrix Amatrix
  coef(result, type = 'Sigma')   # Error covariance matrix Σ
  coef(result, type = 'beta0')   # regression coefficients

[goodness-of-fit measures]
  logLik(result)
  AIC(result)
  BIC(result)

[residuals]
  residuals(result)
  residuals(result, type = 'standardized')

[Accessing individual equations]
  result$individual_models$y1$coefficients
  result$individual_models$y2$time_lag_params

[Accessing the systemfit result]
  summary(result$sur_result)

================================================================================
")

################################################################################
# START OF FILE: fit_ols_diagonal.r
################################################################################

################################################################################
# fit_ols_diagonal.r
# 
# Diagonal-constrained OLS estimation (000d model)
# 
# Estimate y1 and y2 separately with lm()
# - R matrix: 0 (no spatial lag)
# - Lambda matrix: 0 (no spatial error dependence)
# - A matrix: 0 (no temporal AR(1) component)
# - Sigma matrix: diagonal (errors are independent across responses)
#
# Usage:
#
#   result <- fit_ols_diagonal(
#     data_file = "simulated_data_000d_n400_T5.csv",
#     y_vars = c("y1", "y2"),
#     x_vars = list(
#       y1 = c("x_common1", "x_common2", "x_specific1_1"),
#       y2 = c("x_common1", "x_common2", "x_specific2_1")
#     ),
#     verbose = TRUE
#   )
#
################################################################################

cat("Loaded fit_ols_diagonal.r (diagonal-constrained OLS estimation, 000d model)\n")

################################################################################
# Main function: fit_ols_diagonal
################################################################################

#' Diagonal-constrained OLS estimation (000d model)
#' 
#' Estimate each variable individually with lm() and combine the results
#' 
#' @param data_file  Path to the panel data CSV file
#' @param weight_file  Path to the spatial weight matrix CSV file
#' @param y_vars  Character vector of K response variable names
#' @param x_vars  Named list of exogenous covariate names per response
#' @param time_var  Name of the time index column
#' @param time_point  Time period to use (default: max available)
#' @param region_var  Name of the region identifier column
#' @param include_intercept  Logical; include intercept column
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
#' 
fit_ols_diagonal <- function(
  data_file,
  weight_file = NULL,
  y_vars,
  x_vars,
  time_var = "time",
  time_point = NULL,
  region_var = "region",
  include_intercept = TRUE,
  verbose = TRUE
) {
  
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n", paste(rep("=", 70), collapse=""), "\n")
    cat("=== Diagonal-constrained OLS estimation (000d model) ===\n")
    cat(sprintf("    temporal lag: none\n"))
    cat(sprintf("    spatial correlation: none\n"))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  k <- length(y_vars)
  
  #############################################################################
  # S0: Data preparation
  #############################################################################
  
  if (verbose) cat("\nS0: Data preparation...\n")
  
  # Load data file from CSV
  data <- read.csv(data_file)
  
  # Determine the analysis time point
  if (is.null(time_point)) {
    time_point <- max(data[[time_var]])
  }
  
  # Extract data for the current period t
  data_t <- data[data[[time_var]] == time_point, ]
  n <- nrow(data_t)
  
  # Sort rows by region ID to ensure consistent ordering
  if (region_var %in% names(data_t)) {
    data_t <- data_t[order(data_t[[region_var]]), ]
  }
  
  if (verbose) {
    cat(sprintf("  Data: k=%d variables, n=%d regions\n", k, n))
    cat(sprintf("  Analysis time point: t = %d\n", time_point))
  }
  
  #############################################################################
  # S1-S2: Individual OLS estimation (lm) for each response variable
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Individual OLS estimation (lm)\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Initialise storage for per-variable estimation results
  Sigma_diag <- matrix(0, k, k)  # diagonal Sigma matrix
  
  individual_models <- list()     # per-variable model results
  beta_list <- list()             # regression-coefficient list
  residuals_list <- list()        # residual list
  fitted_list <- list()           # fitted-value list
  
  total_loglik <- 0               # total log-likelihood
  total_params <- 0               # total number of parameters
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    xi_vars <- x_vars[[var_name]]
    if (is.null(xi_vars)) xi_vars <- x_vars[[i]]
    
    if (verbose) {
      cat(sprintf("\n--- Estimating %s ---\n", var_name))
    }
    
    # Construct data frame for spatialreg::lagsarlm / errorsarlm
    df_i <- data.frame(
      y = data_t[[var_name]]
    )
    
    # Append the exogenous covariate columns for this response
    for (xvar in xi_vars) {
      df_i[[xvar]] <- data_t[[xvar]]
    }
    
    # Build regression formula string
    if (include_intercept) {
      formula_str <- "y ~ ."
    } else {
      formula_str <- "y ~ . - 1"
    }
    formula_i <- as.formula(formula_str)
    
    # Estimation via lm()
    tryCatch({
      model_i <- lm(
        formula = formula_i,
        data = df_i
      )
      
      # Extract estimation results from the model object
      s <- summary(model_i)
      sigma2 <- s$sigma^2
      Sigma_diag[i, i] <- sigma2
      
      # Regression coefficients
      beta_coef <- coef(model_i)
      beta_list[[var_name]] <- beta_coef
      
      # Extract residuals and fitted values
      residuals_list[[var_name]] <- as.numeric(residuals(model_i))
      fitted_list[[var_name]] <- as.numeric(fitted(model_i))
      
      # Log-likelihood and parameter count
      loglik_i <- as.numeric(logLik(model_i))
      num_params_i <- length(coef(model_i)) + 1  # β + σ²
      
      total_loglik <- total_loglik + loglik_i
      total_params <- total_params + num_params_i
      
      # Build regression coefficient table
      coef_mat <- s$coefficients
      coef_table <- data.frame(
        parameter = rownames(coef_mat),
        estimate = coef_mat[, "Estimate"],
        std_error = coef_mat[, "Std. Error"],
        t_value = coef_mat[, "t value"],
        p_value = coef_mat[, "Pr(>|t|)"],
        stringsAsFactors = FALSE
      )
      coef_table$signif <- sapply(coef_table$p_value, get_signif_code)
      rownames(coef_table) <- NULL
      
      # Compute R² statistics
      R2 <- s$r.squared
      R2_adj <- s$adj.r.squared
      
      # F-test for overall significance of the regression
      f_stat <- s$fstatistic
      if (!is.null(f_stat)) {
        f_value <- f_stat[1]
        f_df1 <- f_stat[2]
        f_df2 <- f_stat[3]
        f_p <- pf(f_value, f_df1, f_df2, lower.tail = FALSE)
      } else {
        f_value <- NA
        f_df1 <- NA
        f_df2 <- NA
        f_p <- NA
      }
      
      # Save individual model estimation results
      individual_models[[var_name]] <- list(
        model = model_i,
        summary = s,
        
        coefficients = coef_table,
        sigma2 = sigma2,
        
        fit = list(
          loglik = loglik_i,
          AIC = AIC(model_i),
          BIC = BIC(model_i),
          R2 = R2,
          R2_adj = R2_adj,
          F_statistic = f_value,
          F_df1 = f_df1,
          F_df2 = f_df2,
          F_p_value = f_p,
          num_params = num_params_i,
          num_obs = length(df_i$y),
          df_residual = s$df[2]
        )
      )
      
      # Verbose output
      if (verbose) {
        cat("  [Error variance]\n")
        cat(sprintf("    σ² = %.6f (σ = %.4f)\n", sigma2, sqrt(sigma2)))
        
        cat("  [Regression coefficients]\n")
        for (j in 1:nrow(coef_table)) {
          row <- coef_table[j, ]
          cat(sprintf("    %-20s = %8.4f  (SE: %.4f, t: %6.2f, p: %.4f %s)\n",
                      row$parameter, row$estimate, row$std_error,
                      row$t_value, row$p_value, row$signif))
        }
        
        cat("  [Goodness of fit]\n")
        cat(sprintf("    R²         = %.4f\n", R2))
        cat(sprintf("    Adj.R²     = %.4f\n", R2_adj))
        cat(sprintf("    log-likelihood   = %.4f\n", loglik_i))
        cat(sprintf("    AIC        = %.4f\n", AIC(model_i)))
        cat(sprintf("    BIC        = %.4f\n", BIC(model_i)))
        
        if (!is.na(f_value)) {
          cat(sprintf("    Fstatistic    = %.2f (df1=%d, df2=%d, p=%.4f)\n",
                      f_value, f_df1, f_df2, f_p))
        }
      }
      
    }, error = function(e) {
      cat(sprintf("\nError: lm() failed when estimating variable %s\n", var_name))
      cat("Error message:", e$message, "\n")
      stop("Estimation failed")
    })
  }
  
  #############################################################################
  # S3: Aggregate individual results into a joint model object
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat("Merge the results\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
  }
  
  # Compute information criteria
  ic <- compute_information_criteria(total_loglik, total_params, k * n)
  
  # Vectorise residuals and fitted values: [y1 block; y2 block; ...]
  residuals_raw <- unlist(residuals_list)
  fitted_vals <- unlist(fitted_list)
  
  # Standardised residuals
  residuals_std <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    residuals_std[idx] <- residuals_list[[y_vars[i]]] / sqrt(Sigma_diag[i, i])
  }
  
  # Assemble the stacked response vector y
  y_vec <- numeric(k * n)
  for (i in 1:k) {
    idx <- ((i - 1) * n + 1):(i * n)
    y_vec[idx] <- data_t[[y_vars[i]]]
  }
  
  # Compute mean R² and mean adjusted R² across responses
  r2_list <- sapply(individual_models, function(m) m$fit$R2)
  r2_adj_list <- sapply(individual_models, function(m) m$fit$R2_adj)
  r2_mean <- mean(r2_list)
  r2_adj_mean <- mean(r2_adj_list)
  
  # Structure the coefficients (for OLS, R, T, A are all absent)
  coefficients <- list(
    R = NULL,      # no spatial lag
    T = NULL,      # no spatial error
    A = NULL,      # no temporal lag
    Sigma = Sigma_diag,
    beta0 = beta_list
  )
  
  end_time <- Sys.time()
  exec_time <- difftime(end_time, start_time, units = "secs")
  
  #############################################################################
  # Print estimation results
  #############################################################################
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    cat("=== Estimation result summary ===\n")
    cat(paste(rep("=", 70), collapse=""), "\n")
    
    cat("\n[diagonal Sigma matrix (error variance)]\n")
    print(round(Sigma_diag, 6))
    
    cat("\n[goodness-of-fit (combined)]\n")
    cat(sprintf("  total log-likelihood: %.4f\n", total_loglik))
    cat(sprintf("  AIC: %.4f\n", ic$AIC))
    cat(sprintf("  BIC: %.4f\n", ic$BIC))
    cat(sprintf("  meanR²: %.4f\n", r2_mean))
    cat(sprintf("  meanAdj.R²: %.4f\n", r2_adj_mean))
    cat(sprintf("  Number of parameters: %d\n", total_params))
    
    cat(sprintf("\nExecution time: %.2f s\n", as.numeric(exec_time)))
    cat(paste(rep("=", 70), collapse=""), "\n")
  }
  
  #############################################################################
  # Build unified result object (standardised output format)
  #############################################################################
  
  dl <- build_data_list_from_parts(
    y = y_vec,
    k = k, n = n, y_vars = y_vars, x_vars = x_vars,
    time_point = time_point,
    include_time_lag = FALSE,
    include_intercept = include_intercept
  )
  
  result <- build_result_object(
    model_type        = "OLS_diagonal",
    Sigma             = Sigma_diag,
    loglik            = total_loglik,
    num_params        = total_params,
    converged         = TRUE,
    method            = "lm_diagonal",
    data_list         = dl,
    beta0             = beta_list,
    residuals_raw     = residuals_raw,
    residuals_std     = residuals_std,
    individual_models = individual_models,
    execution_time    = exec_time
  )
  
  return(result)
}


################################################################################
# Standard S3 methods
################################################################################

#' Print method for the model object
print.multivar_ols_diagonal <- function(x, digits = 4, ...) {
  
  cat("\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Diagonal-constrained OLS (000d model)\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  cat("\n[Model specification]\n")
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  temporal lag: none\n"))
  cat(sprintf("  spatial correlation: none\n"))
  
  cat("\n[diagonal Sigma matrix (error variance)]\n")
  cat(sprintf("  Diagonal elements: %s\n", 
              paste(sprintf("%.4f", diag(x$coefficients$Sigma)), collapse=", ")))
  
  cat("\n[Goodness-of-fit]\n")
  cat(sprintf("  log-likelihood: %.2f\n", x$fit$loglik))
  cat(sprintf("  AIC: %.2f\n", x$fit$AIC))
  cat(sprintf("  BIC: %.2f\n", x$fit$BIC))
  cat(sprintf("  meanR²: %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²: %.4f\n", x$fit$R2_adj))
  
  cat("\nUse summary() for detailed output.\n\n")
  
  invisible(x)
}


#' Summary method with detailed inference tables
summary.multivar_ols_diagonal <- function(object, digits = 4, ...) {
  
  x <- object
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  cat("Diagonal-constrained OLS - detailed results (000d model)\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  # Model specification section
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat("[1. Model information]\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat(sprintf("  Response variables: %s\n", paste(x$data_info$y_vars, collapse=", ")))
  cat(sprintf("  number of variables (k): %d\n", x$data_info$k))
  cat(sprintf("  number of regions (n): %d\n", x$data_info$n))
  cat(sprintf("  Number of observations: %d\n", x$fit$num_obs))
  cat(sprintf("  Number of parameters: %d\n", x$fit$num_params))
  cat(sprintf("  Time point used: t = %d\n", x$data_info$time_point_used))
  
  # Per-variable estimation results
  for (i in 1:x$data_info$k) {
    var_name <- x$data_info$y_vars[i]
    model_i <- x$individual_models[[var_name]]
    
    cat("\n")
    cat(paste(rep("-", 70), collapse=""), "\n")
    cat(sprintf("[%d. estimation result of %s]\n", i + 1, var_name))
    cat(paste(rep("-", 70), collapse=""), "\n")
    
    # Error variance
    cat(sprintf("\nerror variance: σ² = %.6f (σ = %.4f)\n", 
                model_i$sigma2, sqrt(model_i$sigma2)))
    
    # Regression coefficients
    cat("\nRegression coefficients:\n")
    coef_table <- model_i$coefficients
    for (j in 1:nrow(coef_table)) {
      row <- coef_table[j, ]
      cat(sprintf("  %-20s = %8.4f  (SE: %.4f, t: %6.2f, p: %.4f %s)\n",
                  row$parameter, row$estimate, row$std_error,
                  row$t_value, row$p_value, row$signif))
    }
    
    # Goodness of fit
    cat("\nGoodness of fit:\n")
    cat(sprintf("  R²           = %.4f\n", model_i$fit$R2))
    cat(sprintf("  Adj.R²       = %.4f\n", model_i$fit$R2_adj))
    cat(sprintf("  log-likelihood     = %.4f\n", model_i$fit$loglik))
    cat(sprintf("  AIC          = %.4f\n", model_i$fit$AIC))
    cat(sprintf("  BIC          = %.4f\n", model_i$fit$BIC))
    
    # F-test for overall significance of the regression
    if (!is.na(model_i$fit$F_statistic)) {
      cat(sprintf("  Fstatistic      = %.2f (df1=%d, df2=%d, p=%.4f)\n",
                  model_i$fit$F_statistic, 
                  model_i$fit$F_df1, 
                  model_i$fit$F_df2, 
                  model_i$fit$F_p_value))
    }
  }
  
  # Combined (joint) estimation results
  cat("\n")
  cat(paste(rep("-", 70), collapse=""), "\n")
  cat(sprintf("[%d. combinedgoodness-of-fit]\n", x$data_info$k + 2))
  cat(paste(rep("-", 70), collapse=""), "\n")
  
  cat("\ndiagonal Sigma matrix:\n")
  print(round(x$coefficients$Sigma, digits + 2))
  
  cat("\ncombinedgoodness-of-fit:\n")
  cat(sprintf("  total log-likelihood = %.4f\n", x$fit$loglik))
  cat(sprintf("  AIC          = %.4f\n", x$fit$AIC))
  cat(sprintf("  BIC          = %.4f\n", x$fit$BIC))
  cat(sprintf("  meanR²       = %.4f\n", x$fit$R2))
  cat(sprintf("  meanAdj.R²   = %.4f\n", x$fit$R2_adj))
  
  cat("\n")
  cat(paste(rep("=", 70), collapse=""), "\n")
  
  invisible(x)
}


#' Extract estimated coefficients
coef.multivar_ols_diagonal <- function(object, type = "all", ...) {
  
  if (type == "all") {
    return(object$coefficients)
  } else if (type == "R") {
    return(NULL)
  } else if (type == "T") {
    return(NULL)
  } else if (type == "A") {
    return(NULL)
  } else if (type == "Sigma") {
    return(object$coefficients$Sigma)
  } else if (type == "beta0") {
    return(object$coefficients$beta0)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract log-likelihood
logLik.multivar_ols_diagonal <- function(object, ...) {
  
  val <- object$fit$loglik
  attr(val, "df") <- object$fit$num_params
  attr(val, "nobs") <- object$fit$num_obs
  class(val) <- "logLik"
  
  return(val)
}


#' Extract AIC
AIC.multivar_ols_diagonal <- function(object, ..., k = 2) {
  return(object$fit$AIC)
}


#' Extract BIC
BIC.multivar_ols_diagonal <- function(object, ...) {
  return(object$fit$BIC)
}


#' Extract residuals
residuals.multivar_ols_diagonal <- function(object, type = "raw", ...) {
  
  if (type == "raw") {
    return(object$residuals$raw)
  } else if (type == "standardized") {
    return(object$residuals$standardized)
  } else {
    stop("Unknown type: ", type)
  }
}


#' Extract fitted values
fitted.multivar_ols_diagonal <- function(object, ...) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars
  
  fitted_vals <- numeric(k * n)
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    idx <- ((i - 1) * n + 1):(i * n)
    fitted_vals[idx] <- fitted(object$individual_models[[var_name]]$model)
  }
  
  return(fitted_vals)
}


################################################################################
# Helper functions (shared across model fitting files)
################################################################################

#' Return significance code string for a p-value
if (!exists("get_signif_code")) {
  get_signif_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01) return("**")
    if (p < 0.05) return("*")
    if (p < 0.1) return(".")
    return("")
  }
}


#' Compute AIC and BIC from log-likelihood and parameter count
if (!exists("compute_information_criteria")) {
  compute_information_criteria <- function(loglik, num_params, num_obs) {
    
    AIC <- -2 * loglik + 2 * num_params
    BIC <- -2 * loglik + log(num_obs) * num_params
    
    return(list(
      AIC = AIC,
      BIC = BIC
    ))
  }
}


################################################################################
# Usage examples
################################################################################

cat("
================================================================================
fit_ols_diagonal.r - Usage example
================================================================================

[Basic usage]
#   # source('fit_ols_diagonal.r')  # commented out in integrated version

  result <- fit_ols_diagonal(
    data_file = 'simulated_data_000d_n400_T5.csv',
    y_vars = c('y1', 'y2'),
    x_vars = list(
      y1 = c('x_common1', 'x_common2', 'x_specific1_1'),
      y2 = c('x_common1', 'x_common2', 'x_specific2_1')
    ),
    verbose = TRUE
  )

[Checking the result]
  print(result)
  summary(result)

[Extracting coefficients]
  coef(result)
  coef(result, type = 'Sigma')
  coef(result, type = 'beta0')

[goodness-of-fit measures]
  logLik(result)
  AIC(result)
  BIC(result)

[residuals]
  residuals(result)
  residuals(result, type = 'standardized')

[Accessing individual models]
  result$individual_models$y1$coefficients
  result$individual_models$y1$fit$R2

================================================================================
")

################################################################################
# START OF FILE: full_hessian_inference.r
################################################################################

################################################################################
# full_hessian_inference.r
#
# Significance test for beta based on the joint full-parameter Hessian (see Appendix A, mstr.pdf)
#
# Overview:
#   compute_vcov_beta() treats (R, Λ) as known (fixed at estimates) and
#   computes the conditional variance Ψ of β (Eq. 16, mstr.pdf).
#
#   This file computes the joint Hessian w.r.t. θ = (β', vec(R)', vec(Λ)')'
#   numerically, and derives a sandwich variance estimator (Appendix A) that
#   propagates estimation uncertainty from (R, Λ) into the variance of β̂.
#
# Theory:
#   Penalised log-likelihood: ℓ_p(θ) = ℓ(θ) - (γ/2)θ'Dθ  (Eq. 36)
#   
#   H(θ̂)  = -∂²ℓ/∂θ∂θ'|_{θ=θ̂}   (unpenalised Hessian)
#   H_p(θ̂) = H(θ̂) + γD             (penalised Hessian)
#   
#   Avar(θ̂) ≈ H_p^{-1} H H_p^{-1}      (sandwich estimator; see Appendix A)
#   
#   D = diag(0,...,0, 1,...,1)  (0 for β block, 1 for spatial block)
#
# Usage:
#   result <- add_full_inference(result, gamma = 5)
#
# Dependencies:
#   multivar_sly_phase1_1.R, phase2_implementation.r,
#   multivar_sem_sdem_v2.r, spatial_core_functions.r,
#   penalized_spatial.r, numDeriv package
#
################################################################################

cat("Loaded full_hessian_inference.r (joint full-parameter Hessian)\n")


################################################################################
# 1. log-likelihood with Sigma profiled out (function of theta = (beta, spatial))
################################################################################

#' SLY: Sigma-profiled log-likelihood -- function of theta = (beta, vec(R))
#'
#' Solve analytically for hat{Sigma} = (1/n)[z'_k z_l] from z = (I-R(x)W)y - X*beta,
#' and return ell(theta) = -Kn/2 log(2*pi*e) + log|I-R(x)W| - n/2 log|hat{Sigma}|
#'
#' @param theta theta = c(beta, vec(R)') -- vector of length p + K^2
#' @param y  Kn×1 stacked response vector
#' @param X  Kn×p composite design matrix
#' @param W    n×n row-normalised spatial weight matrix
#' @param eigen_W  Pre-computed eigenvalues of W (vector)
#' @param k  Number of response variables K
#' @param n  Number of regions
#' @param p dimension of beta
#' @return  Scalar log-likelihood or profile log-likelihood
full_loglik_sly <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  R_vec <- theta[(p + 1):(p + k^2)]
  R <- matrix(R_vec, nrow = k, ncol = k, byrow = TRUE)
  
  # log|I - R⊗W|
  log_det_R <- tryCatch(
    log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_R)) return(-1e10)
  
  # residuals z = (I-R⊗W)y - Xβ
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  z <- y - RWy - X %*% beta
  
  # Σ̂ = (1/n) [z'_k z_l]
  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }
  
  # Check positive-definiteness of Sigma
  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  # ℓ = -Kn/2 log(2πe) + log|I-R⊗W| - n/2 log|Σ̂|
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_R - (n / 2) * log_det_Sigma
  return(loglik)
}


#' SEM: Sigma-profiled log-likelihood -- function of theta = (beta, vec(Lambda))
full_loglik_sem <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  T_vec <- theta[(p + 1):(p + k^2)]
  T_mat <- matrix(T_vec, nrow = k, ncol = k, byrow = TRUE)
  
  # log|I - Λ⊗W|
  log_det_T <- tryCatch(
    log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_T)) return(-1e10)
  
  # residuals z = (I-Λ⊗W)(y - Xβ)
  residual_raw <- y - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Σ̂
  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }
  
  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_T - (n / 2) * log_det_Sigma
  return(loglik)
}


#' SDEM: Sigma-profiled log-likelihood -- function of theta = (beta, vec(R), vec(Lambda))
full_loglik_sdem <- function(theta, y, X, W, eigen_W, k, n, p) {
  
  beta <- theta[1:p]
  R_vec <- theta[(p + 1):(p + k^2)]
  T_vec <- theta[(p + k^2 + 1):(p + 2 * k^2)]
  R <- matrix(R_vec, nrow = k, ncol = k, byrow = TRUE)
  T_mat <- matrix(T_vec, nrow = k, ncol = k, byrow = TRUE)
  
  # log|I - R⊗W| + log|I - Λ⊗W|
  log_det_R <- tryCatch(
    log_det_spatial(R, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  log_det_T <- tryCatch(
    log_det_spatial(T_mat, eigen_W, k, n, verbose = FALSE, smooth = TRUE),
    error = function(e) -Inf
  )
  if (!is.finite(log_det_R) || !is.finite(log_det_T)) return(-1e10)
  
  # residuals z = (I-Λ⊗W){(I-R⊗W)y - Xβ}
  RWy <- compute_RW_times_y(R, W, y, k, n, verbose = FALSE)
  residual_raw <- y - RWy - X %*% beta
  z <- compute_I_minus_TW_times_v(T_mat, W, residual_raw, k, n)
  
  # Σ̂
  Sigma_hat <- matrix(0, k, k)
  for (i in 1:k) {
    i_idx <- ((i - 1) * n + 1):(i * n)
    for (j in i:k) {
      j_idx <- ((j - 1) * n + 1):(j * n)
      Sigma_hat[i, j] <- sum(z[i_idx] * z[j_idx]) / n
      if (i != j) Sigma_hat[j, i] <- Sigma_hat[i, j]
    }
  }
  
  eig_S <- eigen(Sigma_hat, only.values = TRUE)$values
  if (min(Re(eig_S)) <= 0) {
    Sigma_hat <- Sigma_hat + diag(1e-8, k)
  }
  
  log_det_Sigma <- determinant(Sigma_hat, logarithm = TRUE)$modulus[1]
  if (!is.finite(log_det_Sigma)) return(-1e10)
  
  loglik <- -(k * n / 2) * log(2 * pi * exp(1)) + log_det_R + log_det_T - (n / 2) * log_det_Sigma
  return(loglik)
}


################################################################################
# 2. Joint full-parameter Hessian computation
################################################################################

#' Numerical Hessian for all parameters theta = (beta, spatial)
#'
#' Compute the second derivative of -ell(theta) using numDeriv::hessian
#'
#' @param object  A 'multivar_spatial' model object
#' @param verbose Logical; print diagnostic messages
#' @return list(H, theta_hat, p_beta, p_spatial, dim_theta) or NULL
compute_full_hessian <- function(object, verbose = TRUE) {
  
  k <- object$data_info$k
  n <- object$data_info$n
  y <- object$model_data$y
  X <- object$model_data$X
  W <- object$model_data$W
  eigen_W <- object$model_data$eigen_W
  beta <- object$coefficients$beta
  model_type <- object$model_type
  
  p <- length(beta)
  
  # Construct theta_hat and select the model-specific likelihood function
  if (model_type == "SLY") {
    R_vec <- as.vector(t(object$coefficients$R))
    theta_hat <- c(beta, R_vec)
    p_spatial <- k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_sly(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else if (model_type == "SEM") {
    T_vec <- as.vector(t(object$coefficients$T))
    theta_hat <- c(beta, T_vec)
    p_spatial <- k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_sem(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else if (model_type == "SDEM") {
    R_vec <- as.vector(t(object$coefficients$R))
    T_vec <- as.vector(t(object$coefficients$T))
    theta_hat <- c(beta, R_vec, T_vec)
    p_spatial <- 2 * k^2
    
    neg_loglik <- function(theta) {
      -full_loglik_sdem(theta, y, X, W, eigen_W, k, n, p)
    }
    
  } else {
    warning(sprintf("Model type '%s' is not supported for the full-parameter Hessian", model_type))
    return(NULL)
  }
  
  dim_theta <- length(theta_hat)
  
  if (verbose) {
    cat(sprintf("  Dimension of theta: %d (beta: %d, spatial: %d)\n", dim_theta, p, p_spatial))
    cat("  Computing numDeriv::hessian...\n")
  }
  
  # Numerical Hessian computation
  H <- tryCatch({
    numDeriv::hessian(neg_loglik, theta_hat)
  }, error = function(e) {
    warning("Error computing the full-parameter Hessian: ", e$message)
    NULL
  })
  
  if (is.null(H)) return(NULL)
  
  # Symmetrise (absorbs numerical error)
  H <- (H + t(H)) / 2
  
  # Check positive-definiteness of Sigma
  eig_H <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  min_eig <- min(eig_H)
  
  if (verbose) {
    cat(sprintf("  Hessianeigenvalue: min = %.4e, max = %.4e\n",
                min_eig, max(eig_H)))
  }
  
  if (min_eig <= 0) {
    warning(sprintf("Hessian is not positive definite (minimum eigenvalue = %.2e). Applying ridge correction", min_eig))
    ridge <- abs(min_eig) + 1e-6
    H <- H + diag(ridge, dim_theta)
  }
  
  return(list(
    H         = H,
    theta_hat = theta_hat,
    p_beta    = p,
    p_spatial = p_spatial,
    dim_theta = dim_theta
  ))
}


################################################################################
# 3. Sandwich variance estimation (see Appendix A) -- with the D matrix
################################################################################

#' Sandwich variance of all parameters
#'
#' D = diag(0,...,0, 1,...,1): no penalty on beta, penalty on the spatial parameters
#' Avar(θ̂) = (H + γD)⁻¹ H (H + γD)⁻¹
#'
#' @param H unpenalisedHessian (dim_theta × dim_theta)
#' @param gamma  Penalty strength γ ≥ 0
#' @param p_beta dimension of beta
#' @param p_spatial dimension of the spatial parameters
#' @return  Numeric matrix
compute_full_sandwich_vcov <- function(H, gamma, p_beta, p_spatial) {
  
  dim_theta <- nrow(H)
  
  if (gamma == 0) {
    # No penalty: ordinary inverse information matrix
    vcov <- tryCatch(solve(H), error = function(e) {
      warning("Error inverting the full Hessian. Applying regularisation: ", e$message)
      eig <- eigen(H, symmetric = TRUE)
      ridge <- max(abs(min(eig$values)), 1e-6)
      solve(H + diag(ridge, dim_theta))
    })
  } else {
    # D matrix: 0 for beta, 1 for the spatial parameters (the D matrix in Eq. 36, mstr.pdf)
    D_diag <- c(rep(0, p_beta), rep(1, p_spatial))
    D <- diag(D_diag)
    
    H_pen <- H + gamma * D
    
    H_pen_inv <- tryCatch(solve(H_pen), error = function(e) {
      warning("Error inverting (H + gamma D). Applying regularisation: ", e$message)
      solve(H_pen + diag(1e-6, dim_theta))
    })
    
    # sandwich: H_p⁻¹ H H_p⁻¹
    vcov <- H_pen_inv %*% H %*% H_pen_inv
  }
  
  return(vcov)
}


################################################################################
# 4. Construction of the Z statistic (see Appendix A)
################################################################################

#' Compute the Z statistic (see Appendix A) for all parameters
#'
#' Z_j = θ̂_j / sqrt([H_p⁻¹ H H_p⁻¹]_jj)
#'
#' @param theta_hat theta_hat vector
#' @param vcov_full sandwich variance of all parameters
#' @return data.frame (index, estimate, se, z, p_value)
compute_full_z_statistics <- function(theta_hat, vcov_full) {
  
  se <- sqrt(pmax(diag(vcov_full), 0))
  z <- ifelse(se > 0, theta_hat / se, NA_real_)
  p_value <- ifelse(!is.na(z), 2 * pnorm(-abs(z)), NA_real_)
  
  data.frame(
    index    = seq_along(theta_hat),
    estimate = theta_hat,
    se       = se,
    z_value  = z,
    p_value  = p_value,
    stringsAsFactors = FALSE
  )
}


################################################################################
# 5. Main function: add_full_inference
################################################################################

#' Add a full-parameter-Hessian-based significance test to the estimation result
#'
#' In addition to the Psi-based SE computed by the existing add_inference(),
#' this adds SE from the joint full-parameter Hessian -> sandwich estimator
#' as a "full" slot.
#'
#' @param object  A 'multivar_spatial' model object
#' @param gamma  Penalty strength γ ≥ 0
#' @param verbose Logical; print diagnostic messages
#' @return  S3 model object of class multivar_spatial
add_full_inference <- function(object, gamma = NULL, verbose = TRUE) {
  
  model_type <- object$model_type
  
  # Supported-model check
  if (!(model_type %in% c("SLY", "SEM", "SDEM"))) {
    if (verbose) cat(sprintf("%s: the full-parameter Hessian is supported only for SLY/SEM/SDEM\n", model_type))
    return(object)
  }
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("Joint full-parameter-Hessian significance test (see Appendix A)\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
  }
  
  # Retrieve the penalty strength γ
  if (is.null(gamma)) {
    if (!is.null(object$penalty$gamma)) {
      gamma <- object$penalty$gamma
    } else {
      gamma <- 0
    }
  }
  if (verbose) cat(sprintf("  γ = %g\n", gamma))
  
  k <- object$data_info$k
  n <- object$data_info$n
  y_vars <- object$data_info$y_vars
  
  # ------------------------------------------------------------------
  # Step 1: Full-parameter Hessian computation
  # ------------------------------------------------------------------
  if (verbose) cat("\n[Step 1] Computing the full-parameter Hessian H(theta_hat)\n")
  
  hess_result <- compute_full_hessian(object, verbose = verbose)
  
  if (is.null(hess_result)) {
    warning("Failed to compute the full-parameter Hessian")
    return(object)
  }
  
  H          <- hess_result$H
  theta_hat  <- hess_result$theta_hat
  p_beta     <- hess_result$p_beta
  p_spatial  <- hess_result$p_spatial
  dim_theta  <- hess_result$dim_theta
  
  # ------------------------------------------------------------------
  # Step 2: sandwich variance (see Appendix A)
  # ------------------------------------------------------------------
  if (verbose) cat("\n[Step 2] Sandwich variance Avar(theta_hat) = H_p^{-1} H H_p^{-1}\n")
  
  vcov_full <- tryCatch({
    compute_full_sandwich_vcov(H, gamma, p_beta, p_spatial)
  }, error = function(e) {
    warning("Error computing the sandwich variance: ", e$message)
    NULL
  })
  
  if (is.null(vcov_full)) return(object)
  
  # ------------------------------------------------------------------
  # Step 3: Z statistic (see Appendix A)
  # ------------------------------------------------------------------
  if (verbose) cat("\n[Step 3] Construction of the Z statistic\n")
  
  z_stats <- compute_full_z_statistics(theta_hat, vcov_full)
  
  # ------------------------------------------------------------------
  # Step 4: beta-block extraction
  # ------------------------------------------------------------------
  se_beta_full <- z_stats$se[1:p_beta]
  z_beta_full  <- z_stats$z_value[1:p_beta]
  p_beta_full  <- z_stats$p_value[1:p_beta]
  
  # spatial-parameter-block extraction
  spatial_idx <- (p_beta + 1):dim_theta
  se_spatial_full <- z_stats$se[spatial_idx]
  z_spatial_full  <- z_stats$z_value[spatial_idx]
  p_spatial_full  <- z_stats$p_value[spatial_idx]
  
  # ------------------------------------------------------------------
  # Step 5: Build a comparison table against the existing Psi-based SE
  # ------------------------------------------------------------------
  if (verbose) {
    cat("\n[Step 4] Comparison of the results\n")
    cat(paste(rep("-", 70), collapse = ""), "\n")
    
    # Existing Psi-based SE
    T_for_vcov <- if (model_type %in% c("SEM", "SDEM")) object$coefficients$T else NULL
    W_for_vcov <- if (!is.null(T_for_vcov)) object$model_data$W else NULL
    
    Psi <- tryCatch(
      compute_vcov_beta(object$model_data$X, object$coefficients$Sigma, k, n,
                        T_mat = T_for_vcov, W = W_for_vcov),
      error = function(e) NULL
    )
    se_beta_psi <- if (!is.null(Psi)) sqrt(pmax(diag(Psi), 0)) else rep(NA, p_beta)
    
    # beta comparison table
    beta_names <- names(object$coefficients$beta)
    if (is.null(beta_names)) beta_names <- paste0("β[", 1:p_beta, "]")
    
    cat("\n  Standard error comparison for beta:\n")
    cat(sprintf("  %-30s %12s %12s %8s\n", "Parameter", "SE(Psi)", "SE(Full)", "Ratio"))
    cat(paste0("  ", paste(rep("-", 62), collapse = ""), "\n"))
    
    for (i in 1:p_beta) {
      ratio <- if (!is.na(se_beta_psi[i]) && se_beta_psi[i] > 0) {
        se_beta_full[i] / se_beta_psi[i]
      } else {
        NA
      }
      sig_full <- get_signif_code(p_beta_full[i])
      cat(sprintf("  %-30s %12.6f %12.6f %7.3f %s\n",
                  beta_names[i], se_beta_psi[i], se_beta_full[i],
                  ifelse(is.na(ratio), NA, ratio), sig_full))
    }
    
    # spatial-parameter comparison table
    cat("\n  Spatial parameters:\n")
    cat(sprintf("  %-20s %12s %12s %12s %8s\n",
                "Parameter", "Estimate", "SE(Full)", "z-value", "Signif"))
    cat(paste0("  ", paste(rep("-", 64), collapse = ""), "\n"))
    
    sp_names <- character(0)
    if (model_type %in% c("SLY", "SDEM")) {
      for (i in 1:k) for (j in 1:k) {
        sp_names <- c(sp_names, sprintf("R[%d,%d]", i, j))
      }
    }
    if (model_type %in% c("SEM", "SDEM")) {
      for (i in 1:k) for (j in 1:k) {
        sp_names <- c(sp_names, sprintf("Λ[%d,%d]", i, j))
      }
    }
    
    for (i in seq_along(spatial_idx)) {
      sig <- get_signif_code(p_spatial_full[i])
      cat(sprintf("  %-20s %12.6f %12.6f %12.4f %s\n",
                  sp_names[i], theta_hat[spatial_idx[i]],
                  se_spatial_full[i], z_spatial_full[i], sig))
    }
  }
  
  # ------------------------------------------------------------------
  # Step 6: Store into the object
  # ------------------------------------------------------------------
  
  # Reshape the spatial-parameter SE into a matrix
  se_R_full <- NULL
  se_T_full <- NULL
  
  if (model_type == "SLY") {
    se_R_full <- matrix(se_spatial_full, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "SEM") {
    se_T_full <- matrix(se_spatial_full, nrow = k, ncol = k, byrow = TRUE)
  } else if (model_type == "SDEM") {
    se_R_full <- matrix(se_spatial_full[1:(k^2)], nrow = k, ncol = k, byrow = TRUE)
    se_T_full <- matrix(se_spatial_full[(k^2 + 1):(2 * k^2)], nrow = k, ncol = k, byrow = TRUE)
  }
  
  object$vcov$full <- vcov_full
  
  object$hessian$full <- list(
    matrix    = H,
    theta_hat = theta_hat,
    p_beta    = p_beta,
    p_spatial = p_spatial,
    method    = "numerical_full"
  )
  
  object$std_errors$beta_full <- se_beta_full
  object$std_errors$R_full    <- se_R_full
  object$std_errors$T_full    <- se_T_full
  
  object$inference$full <- list(
    z_statistics = z_stats,
    gamma        = gamma,
    method       = "full_hessian_sandwich"
  )
  
  if (verbose) {
    cat("\n")
    cat(paste(rep("=", 70), collapse = ""), "\n")
    cat("Added a full-parameter-Hessian-based significance test\n")
    cat("  Access via:\n")
    cat("    object$std_errors$beta_full  — full-Hessian SE of beta\n")
    cat("    object$std_errors$R_full     — SE matrix of R\n")
    cat("    object$std_errors$T_full     — SE matrix of Lambda\n")
    cat("    object$vcov$full             — variance-covariance matrix of all parameters\n")
    cat("    object$inference$full        — Z-statistic table\n")
    cat(paste(rep("=", 70), collapse = ""), "\n\n")
  }
  
  return(object)
}


################################################################################
# 6. CSV-output helper: CSV including full-Hessian results
################################################################################

#' Generate a coefficient table including full-Hessian SE
#'
#' @param object object after add_full_inference
#' @return data.frame (parameter, estimate, se_psi, se_full, z_psi, z_full, ...)
export_full_inference_table <- function(object) {
  
  if (is.null(object$inference$full)) {
    warning("No full-Hessian inference result. Run add_full_inference() first")
    return(NULL)
  }
  
  k <- object$data_info$k
  n <- object$data_info$n
  model_type <- object$model_type
  beta <- object$coefficients$beta
  p <- length(beta)
  
  # Psi-based SE
  T_for_vcov <- if (model_type %in% c("SEM", "SDEM")) object$coefficients$T else NULL
  W_for_vcov <- if (!is.null(T_for_vcov)) object$model_data$W else NULL
  
  Psi <- tryCatch(
    compute_vcov_beta(object$model_data$X, object$coefficients$Sigma, k, n,
                      T_mat = T_for_vcov, W = W_for_vcov),
    error = function(e) NULL
  )
  se_psi <- if (!is.null(Psi)) sqrt(pmax(diag(Psi), 0)) else rep(NA, p)
  
  # Full SE
  se_full <- object$std_errors$beta_full
  
  # Parameter names
  beta_names <- names(beta)
  if (is.null(beta_names)) beta_names <- paste0("beta[", 1:p, "]")
  
  # beta table
  z_psi <- ifelse(se_psi > 0, beta / se_psi, NA)
  z_full <- ifelse(se_full > 0, beta / se_full, NA)
  p_psi <- ifelse(!is.na(z_psi), 2 * pnorm(-abs(z_psi)), NA)
  p_full <- ifelse(!is.na(z_full), 2 * pnorm(-abs(z_full)), NA)
  
  tbl <- data.frame(
    parameter  = beta_names,
    estimate   = beta,
    se_psi     = se_psi,
    se_full    = se_full,
    z_psi      = z_psi,
    z_full     = z_full,
    p_psi      = p_psi,
    p_full     = p_full,
    signif_psi = sapply(p_psi, get_signif_code),
    signif_full = sapply(p_full, get_signif_code),
    stringsAsFactors = FALSE
  )
  
  # Add spatial-parameter rows
  z_stats <- object$inference$full$z_statistics
  spatial_start <- p + 1
  spatial_end <- nrow(z_stats)
  
  sp_names <- character(0)
  if (model_type %in% c("SLY", "SDEM")) {
    for (i in 1:k) for (j in 1:k) sp_names <- c(sp_names, sprintf("R[%d,%d]", i, j))
  }
  if (model_type %in% c("SEM", "SDEM")) {
    for (i in 1:k) for (j in 1:k) sp_names <- c(sp_names, sprintf("Lambda[%d,%d]", i, j))
  }
  
  sp_rows <- z_stats[spatial_start:spatial_end, ]
  
  sp_tbl <- data.frame(
    parameter   = sp_names,
    estimate    = sp_rows$estimate,
    se_psi      = NA_real_,
    se_full     = sp_rows$se,
    z_psi       = NA_real_,
    z_full      = sp_rows$z_value,
    p_psi       = NA_real_,
    p_full      = sp_rows$p_value,
    signif_psi  = "",
    signif_full = sapply(sp_rows$p_value, get_signif_code),
    stringsAsFactors = FALSE
  )
  
  result_tbl <- rbind(tbl, sp_tbl)
  rownames(result_tbl) <- NULL
  
  return(result_tbl)
}


################################################################################
# Usage examples
################################################################################
#
# # --- Prerequisite: an estimated object result_sdem exists ---
#
#
# # Run estimation (e.g. SDEM)
# result <- fit_sdem_penalized(
#   data_file = "simulated_data_1111_n100_T5.csv",
#   weight_file = "spatial_weights_n100.csv",
#   time_point = 2,
#   y_vars = c("y1", "y2"),
#   x_vars = list(
#     y1 = c("x_common1", "x_common2", "x_specific1_1"),
#     y2 = c("x_common1", "x_common2", "x_specific2_1")
#   ),
#   include_time_lag = TRUE,
#   gamma = 5
# )
#
# # Add full-parameter-Hessian inference
# result <- add_full_inference(result, gamma = 5)
#
# # Check the result
# result$std_errors$beta_full   # full-Hessian SE
# result$std_errors$beta        # Psi-based SE (for comparison)
#
# # CSV output of the comparison table
# tbl <- export_full_inference_table(result)
# write.csv(tbl, "full_vs_psi_comparison.csv", row.names = FALSE)
#
################################################################################

################################################################################
# START OF FILE: export_results_csv.r
################################################################################

################################################################################
# export_results_csv.r
#
# Function to write the estimation result to CSV
#
# Output format:
#   parameter, estimate, std_error, z_value, p_value, signif
#   as a long-format table; goodness-of-fit rows are appended at the bottom.
#
# Dependencies:
#   build_output.r, spatial_output_functions.r, spatial_core_functions.r,
#   penalized_spatial.r (compute_profile_hessian_*, compute_gic)
#
################################################################################

cat("Loaded export_results_csv.r\n")

################################################################################
# 1. Sandwich variance computation (Appendix A)
################################################################################

#' Compute the sandwich variance of the spatial parameters
#'
#' Penalised-estimation case with gamma > 0:
#'   Var(θ̂_γ) = (H + γI)⁻¹ H (H + γI)⁻¹
#' γ = 0  case:
#'   Var(theta_hat) = H^{-1} (ordinary inverse information matrix)
#'
#' @param H negative Hessian of the unpenalised profile log-likelihood (k1 x k1)
#' @param gamma  Penalty strength γ ≥ 0
#' @return  Numeric matrix
compute_sandwich_vcov <- function(H, gamma) {
  k1 <- nrow(H)
  
  if (gamma == 0) {
    # ordinary inverse information matrix
    vcov <- tryCatch({
      solve(H)
    }, error = function(e) {
      warning("Error inverting the Hessian. Applying regularisation: ", e$message)
      eigen_H <- eigen(H, symmetric = TRUE)
      ridge <- max(abs(min(eigen_H$values)), 1e-6)
      solve(H + diag(ridge, k1))
    })
  } else {
    # Sandwich variance: (H + gamma I)^{-1} H (H + gamma I)^{-1}
    H_pen <- H + gamma * diag(k1)
    
    H_pen_inv <- tryCatch({
      solve(H_pen)
    }, error = function(e) {
      warning("Error inverting (H + gamma I). Applying regularisation: ", e$message)
      solve(H_pen + diag(1e-6, k1))
    })
    
    vcov <- H_pen_inv %*% H %*% H_pen_inv
  }
  
  # Check that the diagonal elements are not negative
  diag_vcov <- diag(vcov)
  if (any(diag_vcov < 0)) {
    warning("The variance-covariance matrix has negative diagonal elements")
  }
  
  return(vcov)
}


################################################################################
# 2. Standard error of the error covariance matrix Σ
################################################################################

#' Compute the asymptotic SE of the error covariance matrix Sigma
#'
#' Normal-distribution case: Var(hat{sigma}_ij) = (sigma_ii sigma_jj + sigma_ij^2) / n
#'
#' @param Sigma  K×K error covariance matrix Σ
#' @param n  Number of regions
#' @return  Numeric matrix
compute_se_sigma <- function(Sigma, n) {
  k <- nrow(Sigma)
  se_Sigma <- matrix(NA, k, k)
  
  for (i in 1:k) {
    for (j in 1:k) {
      var_sigma_ij <- (Sigma[i, i] * Sigma[j, j] + Sigma[i, j]^2) / n
      se_Sigma[i, j] <- sqrt(max(var_sigma_ij, 0))
    }
  }
  
  return(se_Sigma)
}


################################################################################
# 3. Main function: export_results_csv
################################################################################

#' Write the estimation result to a CSV file
#'
#' @param result object returned by build_result_object()
#' @param output_file  Path for the output CSV file
#' @param gamma  Penalty strength γ ≥ 0
#' @param compute_gic_flag whether to compute GIC
#' @param verbose Logical; print diagnostic messages
#' @return  data.frame with estimation results
export_results_csv <- function(
  result,
  output_file = "estimation_results.csv",
  gamma = NULL,
  compute_gic_flag = TRUE,
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== CSV output start ===\n")
  
  # ================================================================
  # Get basic information
  # ================================================================
  k <- result$data_info$k
  n <- result$data_info$n
  model_type <- result$model_type
  
  # Retrieve the penalty strength γ
  if (is.null(gamma)) {
    if (!is.null(result$penalty$gamma)) {
      gamma <- result$penalty$gamma
    } else {
      gamma <- 0
    }
  }
  
  include_time_lag <- result$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE
  
  y_vars <- result$data_info$y_vars
  x_vars <- result$data_info$x_vars
  
  # ================================================================
  # SE of beta (analytic -- exact regardless of gamma)
  # ================================================================
  if (verbose) cat("  Computing the SE of beta...\n")
  
  se_beta <- NULL
  # For SEM/SDEM pass Λ and W to account for spatial error filtering (Eq. 16)
  T_for_vcov <- if (model_type %in% c("SEM", "SDEM")) result$coefficients$T else NULL
  W_for_vcov <- if (!is.null(T_for_vcov)) result$model_data$W else NULL
  
  vcov_beta <- tryCatch({
    compute_vcov_beta(result$model_data$X, result$coefficients$Sigma, k, n,
                      T_mat = T_for_vcov, W = W_for_vcov)
  }, error = function(e) {
    warning("Error computing the variance-covariance matrix of beta: ", e$message)
    NULL
  })
  
  if (!is.null(vcov_beta)) {
    se_beta <- sqrt(pmax(diag(vcov_beta), 0))
  }
  
  # ================================================================
  # SE of the spatial parameters (sandwich variance)
  # ================================================================
  se_R <- NULL
  se_T <- NULL
  
  if (model_type %in% c("SLY", "SEM", "SDEM")) {
    if (verbose) cat("  Computing the SE of the spatial parameters (sandwich variance)...\n")
    
    H_spatial <- tryCatch({
      if (model_type == "SLY") {
        compute_profile_hessian_sly(
          R_hat = result$coefficients$R,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      } else if (model_type == "SEM") {
        compute_profile_hessian_sem(
          T_hat = result$coefficients$T,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      } else {  # SDEM
        compute_profile_hessian_sdem(
          R_hat = result$coefficients$R,
          T_hat = result$coefficients$T,
          beta_hat = result$coefficients$beta,
          Sigma_hat = result$coefficients$Sigma,
          data_list = result$data_list)
      }
    }, error = function(e) {
      warning("Error in Hessian computation: ", e$message)
      NULL
    })
    
    if (!is.null(H_spatial)) {
      vcov_spatial <- tryCatch({
        compute_sandwich_vcov(H_spatial, gamma)
      }, error = function(e) {
        warning("Error computing the sandwich variance: ", e$message)
        NULL
      })
      
      if (!is.null(vcov_spatial)) {
        se_spatial <- sqrt(pmax(diag(vcov_spatial), 0))
        
        if (model_type == "SLY") {
          se_R <- matrix(se_spatial, nrow = k, ncol = k, byrow = TRUE)
        } else if (model_type == "SEM") {
          se_T <- matrix(se_spatial, nrow = k, ncol = k, byrow = TRUE)
        } else {  # SDEM
          n_R <- k^2
          se_R <- matrix(se_spatial[1:n_R], nrow = k, ncol = k, byrow = TRUE)
          se_T <- matrix(se_spatial[(n_R + 1):(2 * n_R)], nrow = k, ncol = k, byrow = TRUE)
        }
      }
    }
  }
  
  # ================================================================
  # SE of Sigma
  # ================================================================
  if (verbose) cat("  Computing the SE of Sigma...\n")
  se_Sigma <- compute_se_sigma(result$coefficients$Sigma, n)
  
  # ================================================================
  # Assemble output rows
  # ================================================================
  if (verbose) cat("  Assembling the CSV table...\n")
  
  rows <- list()
  
  add_row <- function(param_name, estimate, se = NA) {
    z <- if (!is.na(se) && se > 0) estimate / se else NA
    p <- if (!is.na(z)) 2 * pnorm(-abs(z)) else NA
    sig <- if (!is.na(p)) get_signif_code(p) else ""
    
    # estimated value + significance markers (e.g. "0.234567***")
    if (!is.na(estimate)) {
      est_sig <- paste0(sprintf("%.6f", estimate), sig)
    } else {
      est_sig <- ""
    }
    
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = param_name,
      estimate_signif = est_sig,
      estimate = estimate,
      std_error = se,
      z_value = z,
      p_value = p,
      signif = sig,
      stringsAsFactors = FALSE
    )
  }
  
  # For goodness-of-fit measures (no SE/z/p)
  add_fit_row <- function(param_name, value) {
    # Goodness-of-fit measures have no significance markers
    est_sig <- if (!is.na(value)) sprintf("%.6f", value) else ""
    
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = param_name,
      estimate_signif = est_sig,
      estimate = value,
      std_error = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      signif = "",
      stringsAsFactors = FALSE
    )
  }
  
  # ----------------------------------------------------------------
  # (a) R matrix -- SLY, SDEM only
  # ----------------------------------------------------------------
  if (model_type %in% c("SLY", "SDEM") && !is.null(result$coefficients$R)) {
    R_mat <- result$coefficients$R
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("R[%d,%d]", i, j)
        se_val <- if (!is.null(se_R)) se_R[i, j] else NA
        add_row(param_name, R_mat[i, j], se_val)
      }
    }
  } else {
    # no R matrix -> fill with NA
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("R[%d,%d]", i, j), NA)
      }
    }
  }
  
  # ----------------------------------------------------------------
  # (b) Lambda matrix -- SEM, SDEM only
  # ----------------------------------------------------------------
  if (model_type %in% c("SEM", "SDEM") && !is.null(result$coefficients$T)) {
    T_mat <- result$coefficients$T
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("Lambda[%d,%d]", i, j)
        se_val <- if (!is.null(se_T)) se_T[i, j] else NA
        add_row(param_name, T_mat[i, j], se_val)
      }
    }
  } else {
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("Lambda[%d,%d]", i, j), NA)
      }
    }
  }
  
  # ----------------------------------------------------------------
  # (c) Amatrix (temporal lag)
  # ----------------------------------------------------------------
  if (include_time_lag && !is.null(result$coefficients$alpha)) {
    alpha <- result$coefficients$alpha
    
    # The SE of alpha is taken from the corresponding part of the beta vector
    se_alpha <- NULL
    if (!is.null(se_beta) && !is.null(result$coefficients$beta)) {
      # Compute the total number of beta0 and identify the index of the alpha part
      n_beta0 <- 0
      for (i in 1:k) {
        n_beta0 <- n_beta0 + length(x_vars[[i]])
        if (result$data_info$include_intercept) n_beta0 <- n_beta0 + 1
      }
      alpha_idx <- (n_beta0 + 1):(n_beta0 + k^2)
      if (max(alpha_idx) <= length(se_beta)) {
        se_alpha_vec <- se_beta[alpha_idx]
        se_alpha <- matrix(se_alpha_vec, nrow = k, ncol = k, byrow = TRUE)
      }
    }
    
    for (i in 1:k) {
      for (j in 1:k) {
        param_name <- sprintf("A[%d,%d]", i, j)
        se_val <- if (!is.null(se_alpha)) se_alpha[i, j] else NA
        add_row(param_name, alpha[i, j], se_val)
      }
    }
  } else {
    for (i in 1:k) {
      for (j in 1:k) {
        add_fit_row(sprintf("A[%d,%d]", i, j), NA)
      }
    }
  }
  
  # ----------------------------------------------------------------
  # (d) Sigmamatrix
  # ----------------------------------------------------------------
  Sigma_mat <- result$coefficients$Sigma
  for (i in 1:k) {
    for (j in 1:k) {
      param_name <- sprintf("Sigma[%d,%d]", i, j)
      se_val <- se_Sigma[i, j]
      add_row(param_name, Sigma_mat[i, j], se_val)
    }
  }
  
  # ----------------------------------------------------------------
  # (e) β regression coefficients (per-variable)
  # ----------------------------------------------------------------
  beta_vec <- result$coefficients$beta
  beta0 <- result$coefficients$beta0
  include_intercept <- result$data_info$include_intercept
  
  # Track the beta indices
  beta_idx <- 1
  
  for (i in 1:k) {
    var_name <- y_vars[i]
    
    # Prepend an intercept column (column of ones) if requested
    if (include_intercept) {
      param_name <- sprintf("beta_intercept_%s", var_name)
      est <- beta_vec[beta_idx]
      se_val <- if (!is.null(se_beta) && beta_idx <= length(se_beta)) se_beta[beta_idx] else NA
      add_row(param_name, est, se_val)
      beta_idx <- beta_idx + 1
    }
    
    # Append the exogenous covariate columns for this response
    for (x_name in x_vars[[i]]) {
      param_name <- sprintf("beta_%s_%s", x_name, var_name)
      est <- beta_vec[beta_idx]
      se_val <- if (!is.null(se_beta) && beta_idx <= length(se_beta)) se_beta[beta_idx] else NA
      add_row(param_name, est, se_val)
      beta_idx <- beta_idx + 1
    }
  }
  
  # Temporal AR(1) coefficient matrix A (the trailing part of the beta vector has already been output as the A matrix)
  # → Skip the beta indices
  if (include_time_lag) {
    beta_idx <- beta_idx + k^2
  }
  
  # ----------------------------------------------------------------
  # (f) goodness-of-fitmeasures
  # ----------------------------------------------------------------
  
  # AIC, BIC
  add_fit_row("AIC", result$fit$AIC)
  add_fit_row("BIC", result$fit$BIC)
  
  # GIC (full model and compute_gic_flag case)
  if (compute_gic_flag && model_type %in% c("SLY", "SEM", "SDEM")) {
    if (!is.null(H_spatial)) {
      # number of non-spatial parameters k2
      n_beta0 <- 0
      for (i in 1:k) {
        n_beta0 <- n_beta0 + length(x_vars[[i]])
        if (include_intercept) n_beta0 <- n_beta0 + 1
      }
      if (include_time_lag) n_beta0 <- n_beta0 + k^2
      k2 <- n_beta0 + k * (k + 1) / 2  # β + vech(Σ)
      n_obs <- k * n
      
      gic_result <- tryCatch({
        compute_gic(H = H_spatial, gamma = gamma,
                    loglik = result$fit$loglik, k2 = k2, n_obs = n_obs)
      }, error = function(e) {
        warning("Error in GIC computation: ", e$message)
        NULL
      })
      
      if (!is.null(gic_result)) {
        add_fit_row("GIC_AIC", gic_result$GIC_AIC)
        add_fit_row("GIC_BIC", gic_result$GIC_BIC)
        add_fit_row("df_eff", gic_result$df_eff)
      } else {
        add_fit_row("GIC_AIC", NA)
        add_fit_row("GIC_BIC", NA)
        add_fit_row("df_eff", NA)
      }
    } else {
      add_fit_row("GIC_AIC", NA)
      add_fit_row("GIC_BIC", NA)
      add_fit_row("df_eff", NA)
    }
  }
  
  # Pseudo R² (per-variable + mean) (Eq. 26, mstr.pdf)
  if (!is.null(result$coefficients$beta)) {
    fitted_vals <- tryCatch({
      compute_fitted_multivar(
        model_type = model_type,
        R = result$coefficients$R,
        beta = result$coefficients$beta,
        X = result$model_data$X,
        W = result$model_data$W,
        k = k, n = n
      )
    }, error = function(e) {
      warning("Error computing y_hat: ", e$message)
      NULL
    })
    
    r2_result <- NULL
    if (!is.null(fitted_vals)) {
      r2_result <- tryCatch({
        compute_r_squared_multivar(
          y = result$model_data$y,
          fitted = fitted_vals,
          k = k, n = n
        )
      }, error = function(e) {
        warning("Error computing R^2: ", e$message)
        NULL
      })
    }
    
    if (!is.null(r2_result)) {
      for (i in 1:k) {
        add_fit_row(sprintf("R2_%s", y_vars[i]), r2_result$R2_pseudo_individual[i])
      }
      add_fit_row("R2", r2_result$R2_pseudo_mean)
    } else {
      for (i in 1:k) add_fit_row(sprintf("R2_%s", y_vars[i]), NA)
      add_fit_row("R2", NA)
    }
  } else {
    for (i in 1:k) add_fit_row(sprintf("R2_%s", y_vars[i]), NA)
    add_fit_row("R2", result$fit$R2)
  }
  
  # Others
  add_fit_row("n_param", result$fit$num_params)
  add_fit_row("loglik", result$fit$loglik)
  add_fit_row("n_obs", result$fit$num_obs)
  add_fit_row("n_region", n)
  add_fit_row("gamma", gamma)
  
  # ================================================================
  # Build data.frame and write CSV
  # ================================================================
  
  output_df <- do.call(rbind, rows)
  rownames(output_df) <- NULL
  
  # Write out CSV
  write.csv(output_df, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output destination: %s\n", output_file))
    cat(sprintf("  Number of rows: %d\n", nrow(output_df)))
    cat("=== CSV output complete ===\n\n")
  }
  
  invisible(output_df)
}


################################################################################
# 4. Batch CSV export for multiple models
################################################################################

#' Write the estimation results of multiple models to a single CSV
#'
#' Arrange each model's results column-wise
#'
#' @param results_list named list (model name -> result object)
#' @param output_file  Path for the output CSV file
#' @param gammas gamma for each model (taken from result if NULL)
#' @param verbose Logical; print diagnostic messages
#' @return  data.frame with estimation results
export_multiple_models_csv <- function(
  results_list,
  output_file = "model_comparison.csv",
  gammas = NULL,
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== Batch CSV output for multiple models: start ===\n")
  
  model_names <- names(results_list)
  n_models <- length(results_list)
  
  # Get each model's results in CSV form
  all_dfs <- list()
  
  for (i in seq_along(results_list)) {
    model_name <- model_names[i]
    result <- results_list[[i]]
    g <- if (!is.null(gammas)) gammas[i] else NULL
    
    if (verbose) cat(sprintf("  Processing %s...\n", model_name))
    
    df <- export_results_csv(
      result,
      output_file = tempfile(),  # do not write individual files
      gamma = g,
      verbose = FALSE
    )
    
    # Qualify column names with the model name (except the parameter column)
    colnames(df)[2:7] <- paste0(colnames(df)[2:7], ".", model_name)
    
    all_dfs[[model_name]] <- df
  }
  
  # Merge tables on the 'parameter' column
  merged <- all_dfs[[1]]
  if (n_models > 1) {
    for (i in 2:n_models) {
      merged <- merge(merged, all_dfs[[i]], by = "parameter",
                      all = TRUE, sort = FALSE)
    }
  }
  
  # Preserve the original order
  param_order <- all_dfs[[1]]$parameter
  merged <- merged[match(param_order, merged$parameter), ]
  rownames(merged) <- NULL
  
  write.csv(merged, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output destination: %s\n", output_file))
    cat(sprintf("  Number of models: %d, number of parameter rows: %d\n", n_models, nrow(merged)))
    cat("=== Batch CSV output for multiple models: complete ===\n\n")
  }
  
  invisible(merged)
}


################################################################################
# Usage examples
################################################################################
#
# # Single-model output
# result_sdem <- fit_sdem_penalized(...)
# export_results_csv(result_sdem, "sdem_results.csv", verbose = TRUE)
#
# # Multi-model comparison output
# export_multiple_models_csv(
#   results_list = list(
#     SDEM = all_results$SDEM,
#     SEM  = all_results$SEM,
#     SLY  = all_results$SLY
#   ),
#   output_file = "model_comparison.csv"
# )
#
################################################################################

################################################################################
# START OF FILE: experiment_output_functions.r
################################################################################

################################################################################
# experiment_output_functions.r
#
# Table-generation functions for Phase 4 of the numerical experiment
#
# Functions provided:
#   - extract_params_uniform()      : Extract parameters in unified format across model types
#   - build_comparison_table()      : Parameter comparison table across all models (Outputs 2-3)
#   - build_model_selection_table() : Model selection summary table (Output 4)
#   - build_bias_table()            : Bias table vs. true parameters (Output 5)
#   - build_se_comparison_table()   : Psi vs Hessian SE comparison table (Output 6)
#
# Dependencies:
#   spatial_core_functions.r, full_hessian_inference.r
#
################################################################################

cat("Loaded experiment_output_functions.r\n")


################################################################################
# 0. Helper: significance code function
################################################################################

if (!exists("get_signif_code")) {
  get_signif_code <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    if (p < 0.1)   return(".")
    return("")
  }
}


################################################################################
# 1. Unified parameter extraction across all model types
################################################################################

#' Extract all parameters from an estimation-result object in a uniform format
#'
#' For all model types (full/diagonal/VARX/OLS).
#' Parameter names are output in a fixed order so rows align across models.
#'
#' @param result multivar_spatial object
#' @param model_id model-ID string (e.g. "1111", "d0dd")
#' @param gamma  Penalty strength γ ≥ 0
#' @param pAIC pAIC value (equals AIC if NULL)
#' @param pBIC pBIC value (equals BIC if NULL)
#' @param d_eff effective degrees of freedom (num_params if NULL)
#' @return data.frame (parameter, estimate, se_psi, se_hessian, signif_psi, signif_hessian)
extract_params_uniform <- function(result, model_id, gamma = 0,
                                    pAIC = NULL, pBIC = NULL, d_eff = NULL) {
  
  k <- result$data_info$k
  n <- result$data_info$n
  model_type <- result$model_type
  y_vars <- result$data_info$y_vars
  x_vars <- result$data_info$x_vars
  include_time_lag <- result$data_info$include_time_lag
  if (is.null(include_time_lag)) include_time_lag <- FALSE
  include_intercept <- result$data_info$include_intercept
  if (is.null(include_intercept)) include_intercept <- TRUE
  
  # ================================================================
  # Safe row-addition function
  # ================================================================
  rows <- list()
  
  safe_scalar <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    x <- x[1]
    if (is.na(x)) return(NA_real_)
    return(as.numeric(x))
  }
  
  safe_signif <- function(est, se) {
    est <- safe_scalar(est)
    se  <- safe_scalar(se)
    if (is.na(est) || is.na(se) || se <= 0) return("")
    p <- 2 * pnorm(-abs(est / se))
    return(get_signif_code(p))
  }
  
  add_row <- function(param, est, se_p = NA_real_, se_h = NA_real_) {
    est  <- safe_scalar(est)
    se_p <- safe_scalar(se_p)
    se_h <- safe_scalar(se_h)
    rows[[length(rows) + 1]] <<- data.frame(
      parameter      = as.character(param),
      estimate       = est,
      se_psi         = se_p,
      se_hessian     = se_h,
      signif_psi     = safe_signif(est, se_p),
      signif_hessian = safe_signif(est, se_h),
      stringsAsFactors = FALSE
    )
  }
  
  add_fit <- function(param, val) {
    val <- safe_scalar(val)
    rows[[length(rows) + 1]] <<- data.frame(
      parameter = as.character(param), estimate = val,
      se_psi = NA_real_, se_hessian = NA_real_,
      signif_psi = "", signif_hessian = "",
      stringsAsFactors = FALSE
    )
  }
  
  # ================================================================
  # Pre-collect SE information
  # ================================================================
  is_full_model <- model_type %in% c("SLY", "SEM", "SDEM")
  is_diag_model <- model_type %in% c("SLY_diagonal", "SEM_diagonal", "SDEM_diagonal")
  
  # Full model: Psi-based SE (for beta; set by add_inference)
  se_beta_psi  <- result$std_errors$beta       # add_inference
  se_beta_hess <- result$std_errors$beta_full   # add_full_inference
  
  # Full model: spatial parameters Hessian SE
  se_R_hess <- result$std_errors$R_full
  se_T_hess <- result$std_errors$T_full
  
  # Diagonal/Full model: profile-Hessian SE (set by add_inference)
  se_R_prof <- result$std_errors$R
  se_T_prof <- result$std_errors$T
  
  # Full model: Recompute Psi-based SE (if not set)
  if (is_full_model && is.null(se_beta_psi) && !is.null(result$model_data$X)) {
    T_for_vcov <- if (model_type %in% c("SEM", "SDEM")) result$coefficients$T else NULL
    W_for_vcov <- if (!is.null(T_for_vcov) && !is.null(result$model_data)) result$model_data$W else NULL
    Psi <- tryCatch(
      compute_vcov_beta(result$model_data$X, result$coefficients$Sigma, k, n,
                        T_mat = T_for_vcov, W = W_for_vcov),
      error = function(e) NULL
    )
    if (!is.null(Psi)) se_beta_psi <- sqrt(pmax(diag(Psi), 0))
  }
  
  # ================================================================
  # (a) beta regression coefficients (intercept + explanatory variables)
  # ================================================================
  beta_vec <- result$coefficients$beta  # full models (SLY/SEM/SDEM) only
  beta0    <- result$coefficients$beta0 # diagonal model / VARX / OLS
  
  if (!is.null(beta_vec)) {
    # --- Full model / VARX: extract from the unified beta vector ---
    idx <- 1
    for (i in 1:k) {
      vn <- y_vars[i]
      if (include_intercept) {
        pname <- sprintf("beta_intercept_%s", vn)
        se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
        se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
        add_row(pname, beta_vec[idx], se_p, se_h)
        idx <- idx + 1
      }
      for (xn in x_vars[[i]]) {
        pname <- sprintf("beta_%s_%s", xn, vn)
        se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
        se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
        add_row(pname, beta_vec[idx], se_p, se_h)
        idx <- idx + 1
      }
    }
    
    # A matrix (temporal lag) -- trailing part of the beta vector
    if (include_time_lag && idx <= length(beta_vec)) {
      for (i in 1:k) {
        for (j in 1:k) {
          pname <- sprintf("A[%d,%d]", i, j)
          se_p <- if (!is.null(se_beta_psi) && idx <= length(se_beta_psi)) se_beta_psi[idx] else NA
          se_h <- if (!is.null(se_beta_hess) && idx <= length(se_beta_hess)) se_beta_hess[idx] else NA
          add_row(pname, beta_vec[idx], se_p, se_h)
          idx <- idx + 1
        }
      }
    } else {
      # no temporal lag, or VARX (A stored separately)
      # VARX: obtained from coefficients$A
      A_mat <- result$coefficients$A
      alpha_mat <- result$coefficients$alpha
      A_use <- if (!is.null(A_mat)) A_mat else alpha_mat
      for (i in 1:k) {
        for (j in 1:k) {
          pname <- sprintf("A[%d,%d]", i, j)
          est <- if (!is.null(A_use)) A_use[i, j] else NA
          add_row(pname, est, NA, NA)
        }
      }
    }
    
  } else if (!is.null(beta0)) {
    # --- Diagonal / OLS model: extract from the beta0 list ---
    # beta0 is a list of named numeric vectors; it may include temporal-lag terms (e.g. "y1_lag").
    # Temporal-lag terms are output separately as the A matrix, so they are excluded here.
    
    lag_pattern <- "^y[0-9]+_lag$"  # "y1_lag", "y2_lag", etc.
    
    for (i in 1:k) {
      vn <- y_vars[i]
      coef_i <- beta0[[vn]]
      if (is.null(coef_i)) coef_i <- beta0[[i]]
      if (is.null(coef_i)) {
        # Empty case: output NA for intercept + explanatory variables
        if (include_intercept) add_row(sprintf("beta_intercept_%s", vn), NA, NA, NA)
        for (xn in x_vars[[i]]) add_row(sprintf("beta_%s_%s", xn, vn), NA, NA, NA)
        next
      }
      
      cnames <- names(coef_i)
      im <- result$individual_models[[vn]]
      ct <- if (!is.null(im)) im$coefficients else NULL  # spatialreg coefficient table
      
      # Output non-lag terms
      for (j in seq_along(coef_i)) {
        cn <- cnames[j]
        if (grepl(lag_pattern, cn)) next  # skip lag terms
        
        # Parameter-name normalisation
        if (cn == "(Intercept)") {
          pname <- sprintf("beta_intercept_%s", vn)
        } else {
          pname <- sprintf("beta_%s_%s", cn, vn)
        }
        
        # SE: Obtained from the individual_models coefficient table
        se_p <- NA
        if (!is.null(ct)) {
          match_row <- ct[ct$parameter == cn, ]
          if (nrow(match_row) > 0) se_p <- match_row$std_error[1]
        }
        
        add_row(pname, coef_i[j], se_p, NA)
      }
    }
    
    # Get the A matrix: VARX -> coefficients$A/$alpha; diagonal model -> lag terms in beta0
    A_diag <- matrix(NA, k, k)
    A_se   <- matrix(NA, k, k)
    
    # (1) If coefficients$A or $alpha exists, prefer it (for VARX)
    A_from_result <- result$coefficients$A
    if (is.null(A_from_result)) A_from_result <- result$coefficients$alpha
    
    if (!is.null(A_from_result)) {
      A_diag <- A_from_result
      if (!is.null(result$std_errors$A)) {
        A_se <- result$std_errors$A
      }
    } else {
      # (2) diagonal model: extract diagonal elements from the lag terms in beta0
      # no temporal lag for OLS
      for (i in 1:k) {
        vn <- y_vars[i]
        coef_i <- beta0[[vn]]
        if (is.null(coef_i)) coef_i <- beta0[[i]]
        if (is.null(coef_i)) next
        
        im <- result$individual_models[[vn]]
        ct <- if (!is.null(im)) im$coefficients else NULL
        
        for (j in 1:k) {
          lag_name <- sprintf("y%d_lag", j)
          if (lag_name %in% names(coef_i)) {
            A_diag[i, j] <- coef_i[[lag_name]]
            if (!is.null(ct)) {
              match_row <- ct[ct$parameter == lag_name, ]
              if (nrow(match_row) > 0) A_se[i, j] <- match_row$std_error[1]
            }
          }
        }
      }
    }
    
    for (i in 1:k) for (j in 1:k) {
      add_row(sprintf("A[%d,%d]", i, j), A_diag[i, j], A_se[i, j], NA)
    }
    
  } else {
    # both beta and beta0 are NULL: all NA
    for (i in 1:k) {
      vn <- y_vars[i]
      if (include_intercept) add_row(sprintf("beta_intercept_%s", vn), NA, NA, NA)
      for (xn in x_vars[[i]]) add_row(sprintf("beta_%s_%s", xn, vn), NA, NA, NA)
    }
    for (i in 1:k) for (j in 1:k) add_row(sprintf("A[%d,%d]", i, j), NA, NA, NA)
  }
  
  # ================================================================
  # (b) Rmatrix
  # ================================================================
  R_mat <- result$coefficients$R
  for (i in 1:k) for (j in 1:k) {
    est  <- if (!is.null(R_mat)) R_mat[i, j] else NA
    se_p <- if (!is.null(se_R_prof)) se_R_prof[i, j] else NA
    se_h <- if (!is.null(se_R_hess)) se_R_hess[i, j] else NA
    add_row(sprintf("R[%d,%d]", i, j), est, se_p, se_h)
  }
  
  # ================================================================
  # (c) Λmatrix
  # ================================================================
  T_mat <- result$coefficients$T
  for (i in 1:k) for (j in 1:k) {
    est  <- if (!is.null(T_mat)) T_mat[i, j] else NA
    se_p <- if (!is.null(se_T_prof)) se_T_prof[i, j] else NA
    se_h <- if (!is.null(se_T_hess)) se_T_hess[i, j] else NA
    add_row(sprintf("Lambda[%d,%d]", i, j), est, se_p, se_h)
  }
  
  # ================================================================
  # (d) Sigma matrix (upper triangle only)
  # ================================================================
  Sigma_mat <- result$coefficients$Sigma
  for (i in 1:k) for (j in i:k) {
    est <- if (!is.null(Sigma_mat)) Sigma_mat[i, j] else NA
    add_row(sprintf("Sigma[%d,%d]", i, j), est, NA, NA)
  }
  
  # ================================================================
  # (e) goodness-of-fitmeasures
  # ================================================================
  add_fit("loglik",     result$fit$loglik)
  add_fit("AIC",        result$fit$AIC)
  add_fit("BIC",        result$fit$BIC)
  add_fit("pAIC",       if (!is.null(pAIC)) pAIC else result$fit$AIC)
  add_fit("pBIC",       if (!is.null(pBIC)) pBIC else result$fit$BIC)
  add_fit("d_eff",      if (!is.null(d_eff)) d_eff else result$fit$num_params)
  add_fit("gamma",      gamma)
  add_fit("n_params",   result$fit$num_params)
  
  # pseudo_R2 / adj_R2: if result$fit$R2 is NULL, compute from individual_models
  r2_mean <- result$fit$R2
  r2_adj_mean <- result$fit$R2_adj
  
  if ((is.null(r2_mean) || is.na(r2_mean)) && !is.null(result$individual_models)) {
    r2_vals <- c()
    r2_adj_vals <- c()
    for (i in 1:k) {
      vn <- y_vars[i]
      im <- result$individual_models[[vn]]
      if (!is.null(im) && !is.null(im$fit)) {
        r2_i <- if (!is.null(im$fit$R2_pseudo)) im$fit$R2_pseudo else im$fit$R2
        r2_vals <- c(r2_vals, r2_i)
        if (!is.null(im$fit$R2_adj) && !is.na(im$fit$R2_adj)) {
          r2_adj_vals <- c(r2_adj_vals, im$fit$R2_adj)
        }
      }
    }
    if (length(r2_vals) > 0) r2_mean <- mean(r2_vals, na.rm = TRUE)
    if (length(r2_adj_vals) > 0) r2_adj_mean <- mean(r2_adj_vals, na.rm = TRUE)
  }
  
  add_fit("pseudo_R2",  r2_mean)
  add_fit("adj_R2",     r2_adj_mean)
  
  # per-variable R^2
  if (!is.null(result$individual_models)) {
    for (i in 1:k) {
      vn <- y_vars[i]
      im <- result$individual_models[[vn]]
      if (!is.null(im) && !is.null(im$fit)) {
        r2_val <- if (!is.null(im$fit$R2_pseudo)) im$fit$R2_pseudo else im$fit$R2
        add_fit(sprintf("pseudo_R2_%s", vn), r2_val)
      }
    }
  }
  
  # ================================================================
  # Combine and return
  # ================================================================
  df <- do.call(rbind, rows)
  rownames(df) <- NULL
  return(df)
}


################################################################################
# 2. Cross-model parameter comparison table
################################################################################

#' Collect all models' estimation results into a wide comparison table
#'
#' @param results_list  Named list mapping model ID to result object
#' @param gamma_info named list (model ID -> list(gamma, pAIC, pBIC, d_eff))
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return data.frame (invisible)
build_comparison_table <- function(
  results_list,
  gamma_info = NULL,
  output_file = "comparison_table.csv",
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== Building the cross-model parameter comparison table ===\n")
  
  model_ids <- names(results_list)
  n_models <- length(results_list)
  
  # Extract each model's parameters in a uniform format
  all_dfs <- list()
  for (mid in model_ids) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    if (verbose) cat(sprintf("  Processing %s...\n", mid))
    
    df <- extract_params_uniform(
      result, mid,
      gamma = gi$gamma,
      pAIC  = gi$pAIC,
      pBIC  = gi$pBIC,
      d_eff = gi$d_eff
    )
    
    # Qualify column names with the model ID
    colnames(df)[2:6] <- paste0(colnames(df)[2:6], ".", mid)
    all_dfs[[mid]] <- df
  }
  
  # Merge tables on the 'parameter' column
  merged <- all_dfs[[1]]
  if (n_models > 1) {
    for (i in 2:n_models) {
      merged <- merge(merged, all_dfs[[i]], by = "parameter",
                      all = TRUE, sort = FALSE)
    }
  }
  
  # Preserve the original parameter order
  param_order <- all_dfs[[1]]$parameter
  merged <- merged[match(param_order, merged$parameter), ]
  rownames(merged) <- NULL
  
  write.csv(merged, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d rows x %d models)\n", output_file, nrow(merged), n_models))
  }
  
  invisible(merged)
}


################################################################################
# 3. Model-selection summary table
################################################################################

#' Build a summary table condensing model-selection measures to one row per model
#'
#' @param results_list  Named list mapping model ID to result object
#' @param gamma_info_pAIC gamma info under the pAIC criterion
#' @param gamma_info_pBIC gamma info under the pBIC criterion
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return data.frame (invisible)
build_model_selection_table <- function(
  results_pAIC,
  results_pBIC = NULL,
  gamma_info_pAIC = NULL,
  gamma_info_pBIC = NULL,
  output_file = "model_selection_summary.csv",
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== Building the model-selection summary table ===\n")
  
  model_ids <- names(results_pAIC)
  
  rows <- list()
  for (mid in model_ids) {
    r_aic <- results_pAIC[[mid]]
    gi_aic <- if (!is.null(gamma_info_pAIC[[mid]])) gamma_info_pAIC[[mid]] else list(gamma = 0)
    gi_bic <- if (!is.null(gamma_info_pBIC) && !is.null(gamma_info_pBIC[[mid]])) gamma_info_pBIC[[mid]] else gi_aic
    
    # If a pBIC-criterion result exists, retrieve the loglik for BIC
    r_bic <- if (!is.null(results_pBIC[[mid]])) results_pBIC[[mid]] else r_aic
    
    # R2: if result$fit$R2 is NULL, compute from individual_models
    r2_val <- r_aic$fit$R2
    r2_adj_val <- r_aic$fit$R2_adj
    if ((is.null(r2_val) || is.na(r2_val)) && !is.null(r_aic$individual_models)) {
      r2_tmp <- c(); r2_adj_tmp <- c()
      for (vn in r_aic$data_info$y_vars) {
        im <- r_aic$individual_models[[vn]]
        if (!is.null(im) && !is.null(im$fit)) {
          r2_i <- if (!is.null(im$fit$R2_pseudo)) im$fit$R2_pseudo else im$fit$R2
          if (!is.null(r2_i) && !is.na(r2_i)) r2_tmp <- c(r2_tmp, r2_i)
          if (!is.null(im$fit$R2_adj) && !is.na(im$fit$R2_adj)) r2_adj_tmp <- c(r2_adj_tmp, im$fit$R2_adj)
        }
      }
      if (length(r2_tmp) > 0) r2_val <- mean(r2_tmp, na.rm = TRUE)
      if (length(r2_adj_tmp) > 0) r2_adj_val <- mean(r2_adj_tmp, na.rm = TRUE)
    }
    
    rows[[mid]] <- data.frame(
      model_id       = mid,
      n_params       = r_aic$fit$num_params,
      d_eff_pAIC     = if (!is.null(gi_aic$d_eff)) gi_aic$d_eff else r_aic$fit$num_params,
      d_eff_pBIC     = if (!is.null(gi_bic$d_eff)) gi_bic$d_eff else r_aic$fit$num_params,
      loglik_pAIC    = r_aic$fit$loglik,
      loglik_pBIC    = r_bic$fit$loglik,
      AIC            = r_aic$fit$AIC,
      pAIC           = if (!is.null(gi_aic$pAIC)) gi_aic$pAIC else r_aic$fit$AIC,
      BIC            = r_aic$fit$BIC,
      pBIC           = if (!is.null(gi_bic$pBIC)) gi_bic$pBIC else r_aic$fit$BIC,
      gamma_pAIC     = gi_aic$gamma,
      gamma_pBIC     = gi_bic$gamma,
      R2_mean        = if (!is.null(r2_val)) r2_val else NA,
      adj_R2_mean    = if (!is.null(r2_adj_val)) r2_adj_val else NA,
      stringsAsFactors = FALSE
    )
  }
  
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL
  
  # ΔAIC, ΔpAIC, ΔBIC, ΔpBIC  add
  tbl$delta_AIC  <- tbl$AIC  - min(tbl$AIC,  na.rm = TRUE)
  tbl$delta_pAIC <- tbl$pAIC - min(tbl$pAIC, na.rm = TRUE)
  tbl$delta_BIC  <- tbl$BIC  - min(tbl$BIC,  na.rm = TRUE)
  tbl$delta_pBIC <- tbl$pBIC - min(tbl$pBIC, na.rm = TRUE)
  
  # Sort by AIC
  tbl <- tbl[order(tbl$AIC), ]
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d models)\n", output_file, nrow(tbl)))
    cat(sprintf("  AIC best:  %s (%.2f)\n", tbl$model_id[1], tbl$AIC[1]))
    cat(sprintf("  pAIC best: %s (%.2f)\n",
                tbl$model_id[which.min(tbl$pAIC)], min(tbl$pAIC)))
    cat(sprintf("  BIC best:  %s (%.2f)\n",
                tbl$model_id[which.min(tbl$BIC)], min(tbl$BIC)))
    cat(sprintf("  pBIC best: %s (%.2f)\n",
                tbl$model_id[which.min(tbl$pBIC)], min(tbl$pBIC)))
  }
  
  invisible(tbl)
}


################################################################################
# 4. Comparison table against true values (bias assessment)
################################################################################

#' Build a table assessing the bias of estimates against true values
#'
#' @param results_list  Named list mapping model ID to result object
#' @param true_params TRUE_PARAMS list
#' @param gamma_info gamma info (passed to extract_params_uniform)
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return data.frame (invisible)
build_bias_table <- function(
  results_list,
  true_params,
  gamma_info = NULL,
  output_file = "bias_comparison.csv",
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== Building the true-value comparison table (bias assessment) ===\n")
  
  # Vectorise the true values (parameter name -> value mapping)
  tv <- list()
  # β
  tv[["beta_intercept_y1"]] <- true_params$beta_intercept_y1
  tv[["beta_intercept_y2"]] <- true_params$beta_intercept_y2
  tv[["beta_x_common1_y1"]] <- true_params$beta_common1_y1
  tv[["beta_x_common1_y2"]] <- true_params$beta_common1_y2
  tv[["beta_x_common2_y1"]] <- true_params$beta_common2_y1
  tv[["beta_x_common2_y2"]] <- true_params$beta_common2_y2
  tv[["beta_x_specific1_1_y1"]] <- true_params$beta_specific1_1
  tv[["beta_x_specific2_1_y2"]] <- true_params$beta_specific2_1
  # R
  for (i in 1:2) for (j in 1:2) tv[[sprintf("R[%d,%d]", i, j)]] <- true_params$R[i, j]
  # Lambda
  for (i in 1:2) for (j in 1:2) tv[[sprintf("Lambda[%d,%d]", i, j)]] <- true_params$Lambda[i, j]
  # A
  for (i in 1:2) for (j in 1:2) tv[[sprintf("A[%d,%d]", i, j)]] <- true_params$A[i, j]
  # Sigma
  tv[["Sigma[1,1]"]] <- true_params$Sigma[1, 1]
  tv[["Sigma[1,2]"]] <- true_params$Sigma[1, 2]
  tv[["Sigma[2,2]"]] <- true_params$Sigma[2, 2]
  
  param_names <- names(tv)
  
  # Extract the estimated values for each model
  model_ids <- names(results_list)
  
  # True column
  tbl <- data.frame(parameter = param_names, true = unlist(tv),
                    stringsAsFactors = FALSE)
  rownames(tbl) <- NULL
  
  for (mid in model_ids) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    df <- extract_params_uniform(result, mid, gamma = gi$gamma,
                                  pAIC = gi$pAIC, pBIC = gi$pBIC, d_eff = gi$d_eff)
    
    # Parameter-name matching
    est_col <- rep(NA_real_, length(param_names))
    bias_col <- rep(NA_real_, length(param_names))
    
    for (p_idx in seq_along(param_names)) {
      pn <- param_names[p_idx]
      # Match against the parameter names output by extract_params_uniform
      match_idx <- which(df$parameter == pn)
      if (length(match_idx) == 0) {
        # Handle name variants (e.g. beta_x_common1_y1 vs beta_common1_y1)
        # Try pattern matching
        pn_alt <- gsub("^beta_x_", "beta_", pn)
        match_idx <- which(df$parameter == pn_alt)
      }
      if (length(match_idx) > 0) {
        val <- df$estimate[match_idx[1]]
        if (!is.na(val)) {
          est_col[p_idx] <- val
          bias_col[p_idx] <- val - tv[[pn]]
        }
      }
    }
    
    tbl[[paste0("est_", mid)]] <- est_col
    tbl[[paste0("bias_", mid)]] <- bias_col
  }
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d parameters x %d models)\n",
                output_file, length(param_names), length(model_ids)))
  }
  
  invisible(tbl)
}


################################################################################
# 5. Ψ vs Hessian SE comparison table
################################################################################

#' Build a table comparing Psi-based and Hessian-based SE for the full model
#'
#' @param results_list list of full-model estimation results (after add_full_inference)
#' @param gamma_info gamma info
#' @param criterion "pAIC" or "pBIC"
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return data.frame (invisible)
build_se_comparison_table <- function(
  results_list,
  gamma_info = NULL,
  criterion = "pAIC",
  output_file = "se_comparison.csv",
  verbose = TRUE
) {
  
  if (verbose) cat(sprintf("\n=== Building the Psi vs Hessian SE comparison table (%s criterion) ===\n", criterion))
  
  all_rows <- list()
  
  for (mid in names(results_list)) {
    result <- results_list[[mid]]
    gi <- if (!is.null(gamma_info[[mid]])) gamma_info[[mid]] else list(gamma = 0)
    
    df <- extract_params_uniform(result, mid, gamma = gi$gamma,
                                  pAIC = gi$pAIC, pBIC = gi$pBIC, d_eff = gi$d_eff)
    
    # Exclude goodness-of-fit-measure rows (out of scope for SE comparison)
    fit_params <- c("loglik", "AIC", "BIC", "pAIC", "pBIC", "d_eff",
                    "gamma", "n_params", "pseudo_R2", "adj_R2",
                    "pseudo_R2_y1", "pseudo_R2_y2")
    df <- df[!df$parameter %in% fit_params, ]
    
    # SE_ratio computation (beta coefficients only)
    df$SE_ratio <- ifelse(
      !is.na(df$se_psi) & !is.na(df$se_hessian) & df$se_psi > 0,
      df$se_hessian / df$se_psi,
      NA
    )
    
    df$model_id <- mid
    df$gamma <- gi$gamma
    df$criterion <- criterion
    
    all_rows[[mid]] <- df
  }
  
  tbl <- do.call(rbind, all_rows)
  rownames(tbl) <- NULL
  
  # Tidy the column order
  tbl <- tbl[, c("parameter", "model_id", "gamma", "criterion", "estimate",
                  "se_psi", "se_hessian", "signif_psi", "signif_hessian", "SE_ratio")]
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Output: %s (%d rows)\n", output_file, nrow(tbl)))
  }
  
  invisible(tbl)
}


################################################################################
# 6. CSV output of the gamma-search trajectory table
################################################################################

#' Write the result of compare_gamma() to CSV
#'
#' @param gamma_result return value of compare_gamma()
#' @param model_id model-ID string
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return summary data.frame (invisible)
export_gamma_search_csv <- function(gamma_result, model_id,
                                     output_file, verbose = TRUE) {
  
  summary_df <- gamma_result$summary
  results <- gamma_result$results
  
  # Add detail columns for the spatial parameters
  k <- NULL
  extra_cols <- list()
  
  for (i in seq_len(nrow(summary_df))) {
    g_str <- as.character(summary_df$gamma[i])
    res <- results[[g_str]]
    
    if (is.null(res)) {
      extra_cols[[i]] <- data.frame(
        R11 = NA, R12 = NA, R21 = NA, R22 = NA,
        L11 = NA, L12 = NA, L21 = NA, L22 = NA,
        stringsAsFactors = FALSE
      )
      next
    }
    
    if (is.null(k)) k <- res$data_info$k
    
    R_mat <- res$coefficients$R
    T_mat <- res$coefficients$T
    
    row_data <- data.frame(
      R11 = if (!is.null(R_mat)) R_mat[1, 1] else NA,
      R12 = if (!is.null(R_mat)) R_mat[1, 2] else NA,
      R21 = if (!is.null(R_mat)) R_mat[2, 1] else NA,
      R22 = if (!is.null(R_mat)) R_mat[2, 2] else NA,
      L11 = if (!is.null(T_mat)) T_mat[1, 1] else NA,
      L12 = if (!is.null(T_mat)) T_mat[1, 2] else NA,
      L21 = if (!is.null(T_mat)) T_mat[2, 1] else NA,
      L22 = if (!is.null(T_mat)) T_mat[2, 2] else NA,
      stringsAsFactors = FALSE
    )
    extra_cols[[i]] <- row_data
  }
  
  extra_df <- do.call(rbind, extra_cols)
  out_df <- cbind(summary_df[, c("gamma", "loglik", "penalty",
                                   "AIC", "BIC", "df_eff",
                                   "GIC_AIC", "GIC_BIC")],
                  extra_df)
  
  write.csv(out_df, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) {
    cat(sprintf("  Gamma-search trajectory: %s -> %s (%d rows)\n", model_id, output_file, nrow(out_df)))
  }
  
  invisible(out_df)
}


################################################################################
# 7. Optimal-gamma selection summary
################################################################################

#' Build a summary of the optimal gamma from all gamma-search results
#'
#' @param gamma_searches named list (model ID -> compare_gamma result)
#' @param output_file  Path for the output CSV file
#' @param verbose Logical; print diagnostic messages
#' @return data.frame (model_id, gamma_pAIC, pAIC, d_eff_pAIC, gamma_pBIC, pBIC, d_eff_pBIC)
build_gamma_optimal_summary <- function(
  gamma_searches,
  output_file = "gamma_optimal_summary.csv",
  verbose = TRUE
) {
  
  if (verbose) cat("\n=== Optimal-gamma selection summary ===\n")
  
  rows <- list()
  for (mid in names(gamma_searches)) {
    gs <- gamma_searches[[mid]]
    s <- gs$summary
    
    # Use GIC_AIC if computed, otherwise AIC
    if ("GIC_AIC" %in% colnames(s) && any(!is.na(s$GIC_AIC))) {
      idx_aic <- which.min(s$GIC_AIC)
      paic_val <- s$GIC_AIC[idx_aic]
    } else {
      idx_aic <- which.min(s$AIC)
      paic_val <- s$AIC[idx_aic]
    }
    
    if ("GIC_BIC" %in% colnames(s) && any(!is.na(s$GIC_BIC))) {
      idx_bic <- which.min(s$GIC_BIC)
      pbic_val <- s$GIC_BIC[idx_bic]
    } else {
      idx_bic <- which.min(s$BIC)
      pbic_val <- s$BIC[idx_bic]
    }
    
    rows[[mid]] <- data.frame(
      model_id    = mid,
      gamma_pAIC  = s$gamma[idx_aic],
      pAIC        = paic_val,
      d_eff_pAIC  = if ("df_eff" %in% colnames(s)) s$df_eff[idx_aic] else NA,
      loglik_pAIC = s$loglik[idx_aic],
      gamma_pBIC  = s$gamma[idx_bic],
      pBIC        = pbic_val,
      d_eff_pBIC  = if ("df_eff" %in% colnames(s)) s$df_eff[idx_bic] else NA,
      loglik_pBIC = s$loglik[idx_bic],
      stringsAsFactors = FALSE
    )
    
    if (verbose) {
      cat(sprintf("  %s: γ*_pAIC=%.4g (pAIC=%.2f), γ*_pBIC=%.4g (pBIC=%.2f)\n",
                  mid, s$gamma[idx_aic], paic_val, s$gamma[idx_bic], pbic_val))
    }
  }
  
  tbl <- do.call(rbind, rows)
  rownames(tbl) <- NULL
  
  write.csv(tbl, file = output_file, row.names = FALSE, na = "")
  
  if (verbose) cat(sprintf("  Output: %s\n", output_file))
  
  invisible(tbl)
}

cat("experiment_output_functions.r: all functions loaded\n")

################################################################################
# START OF FILE: run_experiment_n400.r
################################################################################

################################################################################
# run_experiment_n400.r
#
# Numerical-experiment run script (n=400)
#
# For synthetic data generated by DGP = MGNS (1111),
# estimate 11 models and output gamma search, significance tests, and comparison tables.
#
# Usage:
#
################################################################################

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat("Numerical experiment: n = 400\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
total_start_time <- Sys.time()

################################################################################
# Phase 0: configuration
################################################################################

cat("\n### Phase 0: configuration ###\n")

N <- 400
DATA_FILE   <- "../simulated_data_1111_n400_T5.csv"
WEIGHT_FILE <- "../spatial_weights_n400.csv"
OUTPUT_DIR  <- "output_n400"

# Gamma-search grid (two-stage: coarse -> fine)
GAMMAS_COARSE <- c(0, 10^seq(-2, 4, by = 0.5))

# Common arguments
Y_VARS <- c("y1", "y2")
X_VARS <- list(
  y1 = c("x_common1", "x_common2", "x_specific1_1"),
  y2 = c("x_common1", "x_common2", "x_specific2_1")
)
TIME_VAR    <- "time"
TIME_POINT  <- 2
REGION_VAR  <- "region"

# True values
TRUE_PARAMS <- list(
  beta_intercept_y1 =  0.20,  beta_intercept_y2 = -0.10,
  beta_common1_y1   =  1.20,  beta_common1_y2   =  0.80,
  beta_common2_y1   = -0.60,  beta_common2_y2   = -0.30,
  beta_specific1_1  =  0.80,
  beta_specific2_1  =  0.50,
  R = matrix(c(0.45, 0.15, 0.10, 0.35), nrow = 2, byrow = TRUE),
  Lambda = matrix(c(0.22, 0.08, 0.05, 0.17), nrow = 2, byrow = TRUE),
  A = matrix(c(0.55, 0.05, 0.10, 0.65), nrow = 2, byrow = TRUE),
  Sigma = matrix(c(0.10, 0.03, 0.03, 0.10), nrow = 2, byrow = TRUE)
)

# Create the output directory
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Load source files
cat("  Loading source files...\n")
required_packages <- c("spdep", "spatialreg", "Matrix", "numDeriv")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("Package %s is required: install.packages('%s')", pkg, pkg))
  }
}


# Check file existence
stopifnot(file.exists(DATA_FILE))
stopifnot(file.exists(WEIGHT_FILE))

cat(sprintf("  Data: %s\n", DATA_FILE))
cat(sprintf("  Weight matrix: %s\n", WEIGHT_FILE))
cat(sprintf("  Output destination: %s/\n", OUTPUT_DIR))
cat(sprintf("  Gamma grid (coarse): %d points (%.4g to %.4g)\n",
            length(GAMMAS_COARSE), min(GAMMAS_COARSE), max(GAMMAS_COARSE)))

################################################################################
# Phase 1: gamma search (6 full models)
################################################################################

cat("\n")
cat(paste(rep("#", 80), collapse = ""), "\n")
cat("### Phase 1: Tuning parameter search (grid search for optimal γ) ###\n")
cat(paste(rep("#", 80), collapse = ""), "\n")

phase1_start <- Sys.time()

# Gamma-search configuration: model_id -> (fit_func, include_time_lag)
gamma_search_configs <- list(
  "1111" = list(func = fit_sdem_penalized, time_lag = TRUE,  label = "SDEM(1111)"),
  "0111" = list(func = fit_sem_penalized,  time_lag = TRUE,  label = "SEM(0111)"),
  "1011" = list(func = fit_sly_penalized,  time_lag = TRUE,  label = "SLY(1011)"),
  "1101" = list(func = fit_sdem_penalized, time_lag = FALSE, label = "SDEM(1101)"),
  "0101" = list(func = fit_sem_penalized,  time_lag = FALSE, label = "SEM(0101)"),
  "1001" = list(func = fit_sly_penalized,  time_lag = FALSE, label = "SLY(1001)")
)

gamma_searches <- list()

for (mid in names(gamma_search_configs)) {
  cfg <- gamma_search_configs[[mid]]
  
  cat(sprintf("\n--- Gamma search for %s ---\n", cfg$label))
  
  gs <- tryCatch({
    compare_gamma_twostage(
      fit_func = cfg$func,
      gammas_coarse = GAMMAS_COARSE,
      compute_gic_flag = TRUE,
      data_file = DATA_FILE,
      weight_file = WEIGHT_FILE,
      y_vars = Y_VARS,
      x_vars = X_VARS,
      time_var = TIME_VAR,
      time_point = TIME_POINT,
      region_var = REGION_VAR,
      include_time_lag = cfg$time_lag
    )
  }, error = function(e) {
    cat(sprintf("  FAIL - error: %s\n", e$message))
    NULL
  })
  
  if (!is.null(gs)) {
    gamma_searches[[mid]] <- gs
    
    # CSV output (Output 1)
    export_gamma_search_csv(
      gs, mid,
      output_file = file.path(OUTPUT_DIR, sprintf("gamma_search_%s.csv", mid)),
      verbose = TRUE
    )
  }
}

# Optimal-gamma summary (summary of Output 1)
gamma_summary <- build_gamma_optimal_summary(
  gamma_searches,
  output_file = file.path(OUTPUT_DIR, "gamma_optimal_summary.csv"),
  verbose = TRUE
)

phase1_time <- difftime(Sys.time(), phase1_start, units = "mins")
cat(sprintf("\nPhase 1 complete: %.1f min\n", as.numeric(phase1_time)))

################################################################################
# Phase 2: estimation of all 11 models
################################################################################

cat("\n")
cat(paste(rep("#", 80), collapse = ""), "\n")
cat("### Phase 2: estimation of all 11 models ###\n")
cat(paste(rep("#", 80), collapse = ""), "\n")

phase2_start <- Sys.time()

# --- 2a: full models -- take the optimal-gamma result from the gamma-search results ---

results_pAIC <- list()   # pAIC criterion
results_pBIC <- list()   # pBIC criterion
gamma_info_pAIC <- list() # gamma info (pAIC criterion)
gamma_info_pBIC <- list() # gamma info (pBIC criterion)

for (mid in names(gamma_searches)) {
  gs <- gamma_searches[[mid]]
  s <- gs$summary
  
  # Identify the pAIC-optimal gamma
  if ("GIC_AIC" %in% colnames(s) && any(!is.na(s$GIC_AIC))) {
    idx_aic <- which.min(s$GIC_AIC)
    g_aic <- s$gamma[idx_aic]
    paic_val <- s$GIC_AIC[idx_aic]
    deff_aic <- if ("df_eff" %in% colnames(s)) s$df_eff[idx_aic] else s$AIC[idx_aic]
  } else {
    idx_aic <- which.min(s$AIC)
    g_aic <- s$gamma[idx_aic]
    paic_val <- s$AIC[idx_aic]
    deff_aic <- NA
  }
  
  # Identify the pBIC-optimal gamma
  if ("GIC_BIC" %in% colnames(s) && any(!is.na(s$GIC_BIC))) {
    idx_bic <- which.min(s$GIC_BIC)
    g_bic <- s$gamma[idx_bic]
    pbic_val <- s$GIC_BIC[idx_bic]
    deff_bic <- if ("df_eff" %in% colnames(s)) s$df_eff[idx_bic] else s$BIC[idx_bic]
  } else {
    idx_bic <- which.min(s$BIC)
    g_bic <- s$gamma[idx_bic]
    pbic_val <- s$BIC[idx_bic]
    deff_bic <- NA
  }
  
  # Get the result object
  res_aic <- gs$results[[as.character(g_aic)]]
  res_bic <- gs$results[[as.character(g_bic)]]
  
  results_pAIC[[mid]] <- res_aic
  results_pBIC[[mid]] <- res_bic
  
  gamma_info_pAIC[[mid]] <- list(
    gamma = g_aic, pAIC = paic_val,
    pBIC = if ("GIC_BIC" %in% colnames(s)) s$GIC_BIC[idx_aic] else s$BIC[idx_aic],
    d_eff = deff_aic
  )
  gamma_info_pBIC[[mid]] <- list(
    gamma = g_bic,
    pAIC = if ("GIC_AIC" %in% colnames(s)) s$GIC_AIC[idx_bic] else s$AIC[idx_bic],
    pBIC = pbic_val, d_eff = deff_bic
  )
  
  cat(sprintf("  %s: γ*_pAIC=%.4g, γ*_pBIC=%.4g\n", mid, g_aic, g_bic))
}

# --- 2b: diagonal model ---
cat("\n--- Diagonal-model estimation ---\n")

common_diag_args <- list(
  data_file = DATA_FILE, weight_file = WEIGHT_FILE,
  y_vars = Y_VARS, x_vars = X_VARS,
  time_var = TIME_VAR, time_point = TIME_POINT,
  region_var = REGION_VAR, include_time_lag = TRUE,
  verbose = TRUE
)

diag_models <- list(
  "d0dd" = list(func = fit_sly_diagonal,  label = "SLY diagonal"),
  "0ddd" = list(func = fit_sem_diagonal,  label = "SEM diagonal"),
  "dddd" = list(func = fit_sdem_diagonal, label = "SDEM diagonal")
)

for (mid in names(diag_models)) {
  dm <- diag_models[[mid]]
  cat(sprintf("\n--- %s (%s) ---\n", mid, dm$label))
  
  result <- tryCatch({
    do.call(dm$func, common_diag_args)
  }, error = function(e) {
    cat(sprintf("  FAIL - error: %s\n", e$message))
    NULL
  })
  
  if (!is.null(result)) {
    results_pAIC[[mid]] <- result
    results_pBIC[[mid]] <- result  # diagonal models are gamma-independent
    gamma_info_pAIC[[mid]] <- list(gamma = 0, pAIC = result$fit$AIC,
                                    pBIC = result$fit$BIC, d_eff = result$fit$num_params)
    gamma_info_pBIC[[mid]] <- gamma_info_pAIC[[mid]]
    cat(sprintf("  ✓ %s done (AIC=%.2f)\n", mid, result$fit$AIC))
  }
}

# --- Model 0011: VARX (R = Λ = 0, full A and Σ) ---
cat("\n--- VARX (0011) ---\n")
result_varx <- tryCatch({
  fit_varx(
    data_file = DATA_FILE, weight_file = WEIGHT_FILE,
    y_vars = Y_VARS, x_vars = X_VARS,
    time_var = TIME_VAR, time_point = TIME_POINT,
    region_var = REGION_VAR, verbose = TRUE
  )
}, error = function(e) {
  cat(sprintf("  FAIL - error: %s\n", e$message))
  NULL
})

if (!is.null(result_varx)) {
  results_pAIC[["0011"]] <- result_varx
  results_pBIC[["0011"]] <- result_varx
  gamma_info_pAIC[["0011"]] <- list(gamma = 0, pAIC = result_varx$fit$AIC,
                                     pBIC = result_varx$fit$BIC, d_eff = result_varx$fit$num_params)
  gamma_info_pBIC[["0011"]] <- gamma_info_pAIC[["0011"]]
  cat(sprintf("  ✓ 0011 done (AIC=%.2f)\n", result_varx$fit$AIC))
}

# --- 2d: OLSdiagonal ---
cat("\n--- OLS diagonal (000d) ---\n")
result_ols <- tryCatch({
  fit_ols_diagonal(
    data_file = DATA_FILE,
    y_vars = Y_VARS, x_vars = X_VARS,
    time_var = TIME_VAR, time_point = TIME_POINT,
    region_var = REGION_VAR, verbose = TRUE
  )
}, error = function(e) {
  cat(sprintf("  FAIL - error: %s\n", e$message))
  NULL
})

if (!is.null(result_ols)) {
  results_pAIC[["000d"]] <- result_ols
  results_pBIC[["000d"]] <- result_ols
  gamma_info_pAIC[["000d"]] <- list(gamma = 0, pAIC = result_ols$fit$AIC,
                                     pBIC = result_ols$fit$BIC, d_eff = result_ols$fit$num_params)
  gamma_info_pBIC[["000d"]] <- gamma_info_pAIC[["000d"]]
  cat(sprintf("  ✓ 000d done (AIC=%.2f)\n", result_ols$fit$AIC))
}

phase2_time <- difftime(Sys.time(), phase2_start, units = "mins")
cat(sprintf("\nPhase 2 complete: %.1f min\n", as.numeric(phase2_time)))
cat(sprintf("  Number of estimated models: pAIC=%d, pBIC=%d\n",
            length(results_pAIC), length(results_pBIC)))

################################################################################
# Phase 3: Significance testing — compute Ψ-based and full-Hessian standard
# errors for β, R, Λ, and add inference tables to each estimated model object.
################################################################################

cat("\n")
cat(paste(rep("#", 80), collapse = ""), "\n")
cat("### Phase 3: significance test ###\n")
cat(paste(rep("#", 80), collapse = ""), "\n")

phase3_start <- Sys.time()

# Full spatial models from Table 1 (all have non-trivial R or Λ or both)
full_model_ids <- c("1111", "0111", "1011", "1101", "0101", "1001")

# --- pAIC-criterion version ---
cat("\n--- Significance testing (pAIC-optimal models) ---\n")

# A helper that runs both inference steps for one model. The computations
# (Psi-based SE, profile-Hessian SE, full-Hessian sandwich SE) are identical
# to the sequential version; models are simply processed on separate cores.
mstr_run_inference <- function(obj, g) {
  err_inf  <- NULL
  err_full <- NULL
  
  # Psi-based SE + profile-Hessian SE
  obj <- tryCatch({
    add_inference(obj, compute_spatial_se = TRUE, gamma = g, verbose = FALSE)
  }, error = function(e) {
    err_inf <<- conditionMessage(e)
    obj
  })
  
  # Full Hessian sandwich SE (see Appendix A)
  obj <- tryCatch({
    add_full_inference(obj, gamma = g, verbose = FALSE)
  }, error = function(e) {
    err_full <<- conditionMessage(e)
    obj
  })
  
  list(object = obj, err_inference = err_inf, err_full_inference = err_full)
}

paic_ids <- full_model_ids[sapply(full_model_ids,
                                  function(m) !is.null(results_pAIC[[m]]))]
cat(sprintf("  Running inference for %d model(s) in parallel (up to %d cores)...\n",
            length(paic_ids), MSTR_CORES))

# Pass each model object through the job list so the task is self-contained
# (required for the PSOCK backend on Windows; harmless on the fork backend)
paic_jobs <- lapply(paic_ids, function(mid) {
  list(mid = mid,
       obj = results_pAIC[[mid]],
       g   = gamma_info_pAIC[[mid]]$gamma)
})

paic_inference <- mstr_parallel_lapply(paic_jobs, function(job) {
  mstr_run_inference(job$obj, job$g)
})
names(paic_inference) <- paic_ids

for (mid in paic_ids) {
  g <- gamma_info_pAIC[[mid]]$gamma
  cat(sprintf("\n  %s (γ=%.4g): ", mid, g))
  out <- paic_inference[[mid]]
  if (is.null(out) || inherits(out, "try-error")) {
    cat("worker failed — keeping pre-inference result ")
  } else {
    if (!is.null(out$err_inference)) {
      cat(sprintf("add_inference error: %s ", out$err_inference))
    }
    if (!is.null(out$err_full_inference)) {
      cat(sprintf("add_full_inference error: %s ", out$err_full_inference))
    }
    results_pAIC[[mid]] <- out$object
  }
  cat("Done\n")
}
rm(paic_inference)

# --- pBIC-criterion version ---
cat("\n--- Significance testing (pBIC-optimal models) ---\n")

# First decide which models can reuse the pAIC result (same gamma), then run
# the remaining models in parallel.
pbic_pending <- character(0)
for (mid in full_model_ids) {
  if (is.null(results_pBIC[[mid]])) next
  g <- gamma_info_pBIC[[mid]]$gamma
  
  # Skip if gamma equals the pAIC one (same object)
  if (!is.null(gamma_info_pAIC[[mid]]) &&
      g == gamma_info_pAIC[[mid]]$gamma) {
    results_pBIC[[mid]] <- results_pAIC[[mid]]
    cat(sprintf("  %s: gamma*_pBIC == gamma*_pAIC — reusing pAIC result\n", mid))
    next
  }
  
  pbic_pending <- c(pbic_pending, mid)
}

if (length(pbic_pending) > 0) {
  cat(sprintf("  Running inference for %d model(s) in parallel (up to %d cores)...\n",
              length(pbic_pending), MSTR_CORES))
  
  pbic_jobs <- lapply(pbic_pending, function(mid) {
    list(mid = mid,
         obj = results_pBIC[[mid]],
         g   = gamma_info_pBIC[[mid]]$gamma)
  })
  
  pbic_inference <- mstr_parallel_lapply(pbic_jobs, function(job) {
    mstr_run_inference(job$obj, job$g)
  })
  names(pbic_inference) <- pbic_pending
  
  for (mid in pbic_pending) {
    g <- gamma_info_pBIC[[mid]]$gamma
    cat(sprintf("\n  %s (γ=%.4g): ", mid, g))
    out <- pbic_inference[[mid]]
    if (is.null(out) || inherits(out, "try-error")) {
      cat("worker failed — keeping pre-inference result ")
    } else {
      if (!is.null(out$err_inference)) {
        cat(sprintf("add_inference error: %s ", out$err_inference))
      }
      if (!is.null(out$err_full_inference)) {
        cat(sprintf("add_full_inference error: %s ", out$err_full_inference))
      }
      results_pBIC[[mid]] <- out$object
    }
    cat("Done\n")
  }
  rm(pbic_inference)
}

# --- add_inference for diagonal models and VARX ---
cat("\n--- Significance testing (non-full spatial models) ---\n")
non_full_ids <- c("d0dd", "0ddd", "dddd", "0011", "000d")
for (mid in non_full_ids) {
  if (is.null(results_pAIC[[mid]])) next
  cat(sprintf("  %s: ", mid))
  results_pAIC[[mid]] <- tryCatch({
    add_inference(results_pAIC[[mid]], compute_spatial_se = TRUE, gamma = 0, verbose = FALSE)
  }, error = function(e) {
    cat(sprintf("Error: %s ", e$message))
    results_pAIC[[mid]]
  })
  results_pBIC[[mid]] <- results_pAIC[[mid]]
  cat("Done\n")
}

phase3_time <- difftime(Sys.time(), phase3_start, units = "mins")
cat(sprintf("\nPhase 3 complete: %.1f min\n", as.numeric(phase3_time)))

################################################################################
# Phase 4: Output table generation
# Produce CSV comparison tables (pAIC- and pBIC-optimal estimates),
# model selection summary, bias table vs. true parameters, and SE comparison.
################################################################################

cat("\n")
cat(paste(rep("#", 80), collapse = ""), "\n")
cat("### Phase 4: output-table generation ###\n")
cat(paste(rep("#", 80), collapse = ""), "\n")

phase4_start <- Sys.time()

# --- Output 2: cross-model parameter comparison table (pAIC criterion) ---
cat("\n--- Output 2: Parameter comparison table (pAIC-optimal) ---\n")
build_comparison_table(
  results_list = results_pAIC,
  gamma_info = gamma_info_pAIC,
  output_file = file.path(OUTPUT_DIR, "comparison_table_pAIC.csv"),
  verbose = TRUE
)

# --- Output 3: cross-model parameter comparison table (pBIC criterion) ---
cat("\n--- Output 3: Parameter comparison table (pBIC-optimal) ---\n")
build_comparison_table(
  results_list = results_pBIC,
  gamma_info = gamma_info_pBIC,
  output_file = file.path(OUTPUT_DIR, "comparison_table_pBIC.csv"),
  verbose = TRUE
)

# --- Output 4: model-selection summary table ---
cat("\n--- Output 4: Model selection summary table ---\n")
build_model_selection_table(
  results_pAIC = results_pAIC,
  results_pBIC = results_pBIC,
  gamma_info_pAIC = gamma_info_pAIC,
  gamma_info_pBIC = gamma_info_pBIC,
  output_file = file.path(OUTPUT_DIR, "model_selection_summary.csv"),
  verbose = TRUE
)

# --- Output 5: true-value comparison table ---
cat("\n--- Output 5: Comparison with true parameters (bias assessment) ---\n")
build_bias_table(
  results_list = results_pAIC,
  true_params = TRUE_PARAMS,
  gamma_info = gamma_info_pAIC,
  output_file = file.path(OUTPUT_DIR, "bias_comparison.csv"),
  verbose = TRUE
)

# --- Output 6: Psi vs Hessian SE comparison table ---
cat("\n--- Output 6: Psi vs Hessian SE comparison table (pAIC) ---\n")
full_results_pAIC <- results_pAIC[intersect(full_model_ids, names(results_pAIC))]
build_se_comparison_table(
  results_list = full_results_pAIC,
  gamma_info = gamma_info_pAIC,
  criterion = "pAIC",
  output_file = file.path(OUTPUT_DIR, "se_comparison_full_models_pAIC.csv"),
  verbose = TRUE
)

cat("\n--- Output 6: Psi vs Hessian SE comparison table (pBIC) ---\n")
full_results_pBIC <- results_pBIC[intersect(full_model_ids, names(results_pBIC))]
build_se_comparison_table(
  results_list = full_results_pBIC,
  gamma_info = gamma_info_pBIC,
  criterion = "pBIC",
  output_file = file.path(OUTPUT_DIR, "se_comparison_full_models_pBIC.csv"),
  verbose = TRUE
)

phase4_time <- difftime(Sys.time(), phase4_start, units = "mins")
cat(sprintf("\nPhase 4 complete: %.1f min\n", as.numeric(phase4_time)))

################################################################################
# Phase 5: Save results and report completion
################################################################################

cat("\n")
cat(paste(rep("#", 80), collapse = ""), "\n")
cat("### Phase 5: save and completion report ###\n")
cat(paste(rep("#", 80), collapse = ""), "\n")

# RDS save
saveRDS(results_pAIC,    file.path(OUTPUT_DIR, "results_pAIC.rds"))
saveRDS(results_pBIC,    file.path(OUTPUT_DIR, "results_pBIC.rds"))
saveRDS(gamma_searches,  file.path(OUTPUT_DIR, "gamma_searches.rds"))

# Shut down the parallel workers (no-op on fork backend / single core)
mstr_stop_cluster()

cat("\nRDS save complete:\n")
cat(sprintf("  %s/results_pAIC.rds (%d model)\n", OUTPUT_DIR, length(results_pAIC)))
cat(sprintf("  %s/results_pBIC.rds (%d model)\n", OUTPUT_DIR, length(results_pBIC)))
cat(sprintf("  %s/gamma_searches.rds (%d model)\n", OUTPUT_DIR, length(gamma_searches)))

# List of output files
cat("\nList of output files:\n")
files <- list.files(OUTPUT_DIR, full.names = TRUE)
for (f in files) {
  finfo <- file.info(f)
  cat(sprintf("  %s (%.1f KB)\n", basename(f), finfo$size / 1024))
}

# Execution-time summary
total_time <- difftime(Sys.time(), total_start_time, units = "mins")
cat(sprintf("\nExecution-time summary:\n"))
  cat(sprintf("  Phase 1 (gamma search):   %.1f min\n", as.numeric(phase1_time)))
  cat(sprintf("  Phase 2 (model fitting):  %.1f min\n", as.numeric(phase2_time)))
  cat(sprintf("  Phase 3 (inference):      %.1f min\n", as.numeric(phase3_time)))
  cat(sprintf("  Phase 4 (output tables):  %.1f min\n", as.numeric(phase4_time)))
  cat(sprintf("  Total:                    %.1f min\n", as.numeric(total_time)))

cat("\n")
cat(paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("Numerical experiment n=%d complete\n", N))
cat(paste(rep("=", 80), collapse = ""), "\n")

################################################################################
# Program-wide elapsed-time measurement (END marker + save to file)
# ---------------------------------------------------------------------------
# Reached only when the whole program runs to completion. We compare the
# current time against the `.program_start_time` recorded at the very top to
# obtain the total wall-clock elapsed time, and use the matching proc.time()
# snapshot for cumulative CPU usage. The result is written to a small text
# file so it is preserved independently of the console log.
################################################################################
.program_end_time <- Sys.time()
.program_end_proc <- proc.time()

# Wall-clock elapsed time (the headline "elapsed" figure), in seconds.
.elapsed_secs <- as.numeric(difftime(.program_end_time,
                                      .program_start_time, units = "secs"))

# CPU time consumed since the start marker: user + system, and the
# proc.time() "elapsed" component (also wall-clock, as a cross-check).
.proc_delta    <- .program_end_proc - .program_start_proc
.cpu_user_secs <- as.numeric(.proc_delta["user.self"])
.cpu_sys_secs  <- as.numeric(.proc_delta["sys.self"])

# Format the wall-clock elapsed time as HH:MM:SS for easy reading.
.hh <-  .elapsed_secs %/% 3600
.mm <- (.elapsed_secs %%  3600) %/% 60
.ss <-  .elapsed_secs %%  60
.elapsed_hms <- sprintf("%02d:%02d:%05.2f", as.integer(.hh), as.integer(.mm), .ss)

# Choose an output location: prefer OUTPUT_DIR if it is defined, otherwise the
# current working directory. The start timestamp is embedded in the file name
# so repeated runs accumulate rather than overwrite one another.
.timing_dir <- if (exists("OUTPUT_DIR") && dir.exists(OUTPUT_DIR)) {
  OUTPUT_DIR
} else {
  "."
}
.timing_file <- file.path(
  .timing_dir,
  sprintf("elapsed_time_%s.txt", format(.program_start_time, "%Y%m%d_%H%M%S"))
)

# Assemble the report and write it to the file (wrapped in tryCatch so a write
# failure cannot mask the just-completed run).
.timing_lines <- c(
  "Program elapsed-time report",
  strrep("-", 40),
  sprintf("Start time      : %s", format(.program_start_time, "%Y-%m-%d %H:%M:%S")),
  sprintf("End time        : %s", format(.program_end_time,   "%Y-%m-%d %H:%M:%S")),
  sprintf("Elapsed (wall)  : %s (%.2f sec = %.2f min)",
          .elapsed_hms, .elapsed_secs, .elapsed_secs / 60),
  sprintf("CPU user time   : %.2f sec", .cpu_user_secs),
  sprintf("CPU system time : %.2f sec", .cpu_sys_secs)
)

tryCatch({
  writeLines(.timing_lines, .timing_file)
  cat(sprintf("\n[timing] Total elapsed (wall-clock): %s (%.2f min)\n",
              .elapsed_hms, .elapsed_secs / 60))
  cat(sprintf("[timing] Elapsed-time report saved to: %s\n", .timing_file))
}, error = function(e) {
  # If the file cannot be written, still report the timing to the console.
  cat(sprintf("\n[timing] Total elapsed (wall-clock): %s (%.2f min)\n",
              .elapsed_hms, .elapsed_secs / 60))
  cat(sprintf("[timing] WARNING: could not write timing file (%s)\n", e$message))
})
