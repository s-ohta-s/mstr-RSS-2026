# ============================================================
# implement-B.r
#
# Complete R implementation of the bivariate Multivariate
# Spatio-Temporal Regression (MSTR) / Multivariate General
# Nesting Spatial (MGNS) model described in mstr.pdf.
#
# Statistical map to the current equation numbering in mstr.pdf:
#   * Equation (4) defines the stacked response y_t and the K x K matrices
#     R, A, and Lambda.  This program uses K = 2.
#   * Equations (5) and (6) define the MGNS structural equation and the
#     spatial-error process, respectively.
#   * Equation (7) absorbs the time-lag term (A kron I_n)y_{t-1} into the
#     composite design matrix X_t and coefficient vector beta.
#   * Equation (8) is the standard MGNS form implemented by this program.
#   * Equations (9) and (10) give the log-likelihood, quadratic form Q, and
#     transformed residual vector z.
#   * Equations (11)-(15) give the likelihood equations and the updates for
#     beta and Sigma; Equation (16) gives Var(beta_hat).
#   * Step S1 in Subsection 3.3 uses the univariate GNS model in Equation (17).
#     Steps S2-S5 alternate the beta/Sigma updates and BFGS optimization of
#     R and Lambda subject to |rho_kl| < 1 and |lambda_kl| < 1.
#   * Equations (18)-(23) describe determinant evaluation from eigenvalues and
#     trace-based elementary symmetric polynomials.
#   * Table 1 in Subsection 3.5 defines the 11 model IDs in the order
#     R / Lambda / A / Sigma, with 1 = full, d = diagonal, and 0 = zero.
#   * Equation (24) defines the total effective degrees of freedom, and
#     Equation (25) defines pAIC and pBIC.  Equation (26) defines the averaged
#     pseudo coefficient of determination.
#   * Appendix A, Equations (35)-(44), gives the penalized profile-likelihood,
#     profile-Hessian, influence-matrix, and effective-df justification used
#     for the spatial-parameter block.
#
# Expected input files in DATA_DIR:
#   simulated_data_<DATA_ID>_n<N_VALUE>_T<T_VALUE>.csv
#   spatial_weights_n<N_VALUE>.csv
#   Optional: gamma_grid.csv
#   Optional: insert-template_new.xlsx
#
# Main output files in OUT_DIR:
#   weight_matrix_diagnostics.csv
#   gamma_grid_used.csv
#   gamma_path_<MODEL_ID>.csv
#   run_summary_n<N_VALUE>.csv
#   estimated_parameters_<MODEL_ID>_pAIC.xlsx
#   estimated_parameters_<MODEL_ID>_pBIC.xlsx
#   true_value_evaluation_1111.csv / .xlsx, when model 1111 is available
#   elapsed_time_n<N_VALUE>.csv
# ============================================================

rm(list = ls(all.names = TRUE))

# ------------------------------------------------------------
# Whole-program elapsed-time measurement
# ------------------------------------------------------------
# Record both the wall-clock start time and R's elapsed process time immediately
# after clearing the workspace.  The process-time counter is used for the actual
# elapsed duration because it is monotonic during the R session and is therefore
# not affected by a system-clock adjustment while this long script is running.
PROGRAM_START_TIME <- Sys.time()
PROGRAM_START_ELAPSED <- unname(proc.time()[["elapsed"]])

suppressPackageStartupMessages({
  library(Matrix)
  library(openxlsx)
  library(parallel)
})
if (!requireNamespace("numDeriv", quietly = TRUE)) {
  warning("Package 'numDeriv' is not installed; optimHess will be used as a fallback for profile theta Hessians.")
}

# ------------------------------------------------------------
# User settings and reproducibility controls
# ------------------------------------------------------------
# Edit this block when changing the simulated data size, target time, input
# directory, output behavior, gamma grid source, optimization limits, or
# parallel execution.  The statistical functions below are written so that the
# same code covers n = 100, 400, and 900, and all 11 model restrictions.

DATA_ID       <- "1111"
N_VALUE       <- 400L     # change to 100L / 400L / 900L
T_VALUE       <- 5L
TARGET_TIME   <- 2L
DATA_DIR      <- "."

# Output root selection.  Interactive runs ask whether outputs may be written
# under the current working directory.  Non-interactive runs should set an
# environment variable when a specific output root is required.
#
# Non-interactive runs cannot answer readline() prompts.  In that case, set one
# of the following environment variables before source()/Rscript execution:
#   Sys.setenv(MSTR_USE_CURRENT_DIR = "yes")  # use getwd()
#   Sys.setenv(MSTR_OUT_ROOT = "C:/mstr_out") # force another output root
#' Choose the root directory for all generated outputs.
#'
#' Interactive runs ask whether the current working directory should be used.
#' Non-interactive runs can be controlled through MSTR_USE_CURRENT_DIR or
#' MSTR_OUT_ROOT.  This keeps the script portable across Windows, macOS, and
#' Linux paths without hard-coding a machine-specific output directory.
choose_output_root <- function() {
  env_root <- Sys.getenv("MSTR_OUT_ROOT", unset = "")
  if (nzchar(env_root)) return(path.expand(env_root))

  current_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  short_root <- if (.Platform$OS.type == "windows") {
    user_profile <- Sys.getenv("USERPROFILE", unset = "")
    if (nzchar(user_profile)) file.path(user_profile, "mstr_out") else file.path(tempdir(), "mstr_out")
  } else {
    file.path(getwd(), "mstr_out")
  }

  use_current_env <- tolower(Sys.getenv("MSTR_USE_CURRENT_DIR", unset = ""))
  if (use_current_env %in% c("yes", "y", "true", "1")) return(current_root)
  if (use_current_env %in% c("no", "n", "false", "0")) return(short_root)

  if (interactive()) {
    message("Current working directory:")
    message("  ", current_root)
    message("Proposed output directory:")
    message("  ", file.path(current_root, sprintf("results_n%d", N_VALUE)))
    ans <- readline("カレントディレクトリ配下に出力してよいですか？ [Y/n]: ")
    ans <- trimws(tolower(ans))
    if (ans %in% c("", "y", "yes", "はい", "ha", "h")) {
      return(current_root)
    }

    alt <- readline(sprintf("別の出力ルートを入力してください。空欄なら短い既定パス [%s] を使います: ", short_root))
    alt <- trimws(alt)
    if (nzchar(alt)) return(path.expand(alt))
    return(short_root)
  }

  message("Non-interactive run: no output prompt was possible.")
  message("Using short output root: ", short_root)
  message("To use the current directory in non-interactive runs, set MSTR_USE_CURRENT_DIR=yes.")
  short_root
}

OUT_ROOT      <- choose_output_root()

# The Excel template is optional.  When present it can provide a preferred
# workbook layout; otherwise the result workbooks are generated from scratch.
TEMPLATE_FILE <- file.path(DATA_DIR, "insert-template_new.xlsx")
OUT_DIR       <- file.path(OUT_ROOT, sprintf("results_n%d", N_VALUE))

# Standard input file locations.  The execution block also checks the parent
# directory so that the script can be run either from DATA_DIR or from a
# subdirectory below it.
DATA_FILE     <- file.path(DATA_DIR, sprintf("../simulated_data_%s_n%d_T%d.csv", DATA_ID, N_VALUE, T_VALUE))
WEIGHTS_FILE  <- file.path(DATA_DIR, sprintf("../spatial_weights_n%d.csv", N_VALUE))

# Gamma grid source.  Use GAMMA_SOURCE = "code" to force the vector below,
# "csv" to force GAMMA_CSV_FILE, "xlsx" to force TEMPLATE_FILE/gamma_grid,
# or "auto" to use CSV if present, otherwise code, otherwise xlsx.
#
# Default gamma grid for penalized spatial-parameter selection.  The grid starts
# with gamma = 0 and then follows a log10 sequence from 0.01 through 10000.
# Keeping the default source as "code" makes the run reproducible unless the
# user explicitly switches to a CSV or workbook-provided grid.
GAMMA_SOURCE   <- "code"
GAMMA_CSV_FILE <- file.path(DATA_DIR, "gamma_grid.csv")
GAMMA_VALUES   <- sort(unique(c(0, 10^seq(-2, 4, by = 0.5))))

# Row-standardization check for W.  Rows whose sums differ from 1 by more than
# W_ROW_TOL are normalized.  Zero rows are left unchanged and reported.
W_ROW_STANDARDIZE <- TRUE
W_ROW_TOL         <- 1e-8

MODEL_IDS     <- c("0011","1011","0111","1111","1001","0101","1101",
                   "000d","d0dd","0ddd","dddd")

# Gamma-penalty policy for information criteria.  In this implementation, the
# gamma-dependent effective degrees of freedom in Equation (24) are applied to
# models with full spatial-parameter blocks.  Models without full spatial
# matrices, including diagonal systems, are evaluated without the gamma
# penalty; their reported pAIC/pBIC therefore equal ordinary AIC/BIC.
GAMMA_PENALTY_MODEL_IDS <- c("1111", "0111", "1011", "1101", "0101", "1001")

#' Return whether a model ID uses gamma penalization for spatial parameters.
#'
#' The Equation (24) effective-df shrinkage and Equation (25) pAIC/pBIC
#' formulas are applied only to models with full spatial-parameter blocks in
#' this implementation.  Diagonal or non-spatial models are compared by
#' ordinary AIC/BIC.
uses_gamma_penalty_model <- function(model_id) {
  as.character(model_id) %in% GAMMA_PENALTY_MODEL_IDS
}

#' Build the gamma grid used for a given model.
#'
#' Full spatial models receive the user-specified gamma grid.  Other models are
#' evaluated once at gamma = 0 because no gamma shrinkage is applied to them.
gamma_grid_for_model <- function(model_id, gamma_values) {
  vals <- sort(unique(as.numeric(gamma_values)))
  vals <- vals[is.finite(vals) & vals >= 0]
  if (length(vals) == 0L) stop("gamma grid is empty for model ", model_id)
  if (uses_gamma_penalty_model(model_id)) return(vals)
  0
}

# parallel settings
USE_PARALLEL                <- TRUE
N_CORES                     <- NULL    # NULL = auto
RETRY_FAILED_SEQUENTIALLY   <- TRUE
RESERVE_CORES               <- 1L

# optimization settings
MAXIT_MAIN                  <- 220L
MAXIT_RETRY                 <- 320L
REL_TOL                     <- 1e-8


# -------------------------------
# Profile-likelihood effective-df controls
# -------------------------------
# DEFF_EVALUATION_MODE:
#   "candidate" = fast default.  Compute exact profile theta Hessians only for
#                 gamma values that can still be optimal by a rigorous lower-bound
#                 screen.  The selected pAIC/pBIC values are exact; non-candidate
#                 gamma rows keep pAIC/pBIC as NA and show bounds.
#   "all"       = compute exact profile-Hessian deff/pAIC/pBIC for every gamma grid row.
DEFF_EVALUATION_MODE        <- "candidate"
DEFF_SCREEN_INITIAL_TOP     <- 3L
DEFF_SCREEN_BATCH_SIZE      <- 4L
DEFF_SCREEN_TOL             <- 1e-10

# Numerical Hessian evaluation controls.  The Hessian objective refits beta and Sigma with
# max_iter_inner around 50; this keeps the same idea while being faster than a
# full optimization loop inside every Hessian evaluation.
HESSIAN_MAX_GLS_ITER        <- 50L
HESSIAN_GLS_TOL             <- 1e-6
HESSIAN_THETA_BOUND         <- 0.995
HESSIAN_RIDGE               <- 1e-8
USE_NUMDERIV_HESSIAN        <- TRUE   # TRUE uses numDeriv::hessian when available; otherwise optimHess is used.

# Output speed controls.  The CSV summary/gamma paths are always fast.  Workbook
# output is done only after all models finish, not inside parallel workers.
WRITE_GAMMA_PATH_CSV        <- TRUE
WRITE_XLSX                  <- TRUE
USE_TEMPLATE_FOR_XLSX       <- FALSE  # TRUE preserves template formatting but is slower.
WRITE_TRUE_VALUE_XLSX       <- TRUE



# ------------------------------------------------------------
# Safe output and rerun utilities
# ------------------------------------------------------------
# These helpers make the long estimation run robust to missing directories,
# locked CSV/XLSX files, and reruns from different working directories.

#' Return a compact timestamp used in safe fallback file names.
now_stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

#' Insert a suffix before a file extension.
#'
#' Used when an output file is locked or already open, so the script can save a
#' timestamped fallback file instead of failing after a long model run.
insert_file_suffix <- function(path, suffix) {
  ext <- tools::file_ext(path)
  stem <- if (nzchar(ext)) sub(paste0("\\.", ext, "$"), "", path) else path
  if (nzchar(ext)) paste0(stem, suffix, ".", ext) else paste0(stem, suffix)
}

#' Ensure that an output directory exists and is writable.
#'
#' The function creates the directory recursively if necessary and performs a
#' write/delete test before the main outputs are generated.
ensure_writable_dir <- function(out_dir, label = "output directory") {
  if (is.null(out_dir) || length(out_dir) != 1L || is.na(out_dir) || !nzchar(out_dir)) {
    stop(label, " is empty or invalid.")
  }
  out_dir <- path.expand(out_dir)

  if (!dir.exists(out_dir)) {
    ok <- dir.create(out_dir, recursive = TRUE, showWarnings = TRUE)
    if (!ok && !dir.exists(out_dir)) {
      stop("Cannot create ", label, ": ", normalizePath(out_dir, mustWork = FALSE))
    }
  }

  test_file <- file.path(out_dir, paste0("_write_test_", now_stamp(), ".tmp"))
  ok <- tryCatch({
    writeLines("ok", test_file, useBytes = TRUE)
    unlink(test_file, force = TRUE)
    TRUE
  }, error = function(e) {
    message("Write test failed for ", label, ": ", conditionMessage(e))
    FALSE
  })
  if (!ok) {
    stop("Cannot write to ", label, ": ", normalizePath(out_dir, mustWork = FALSE))
  }

  normalizePath(out_dir, winslash = "/", mustWork = TRUE)
}

#' Write a CSV file through a temporary file and fallback path.
#'
#' This protects results from partial writes and handles common cases where a
#' previous CSV is open in Excel or another application.
safe_write_csv <- function(x, file, row.names = FALSE, ..., fallback_on_lock = TRUE) {
  parent <- dirname(file)
  parent <- ensure_writable_dir(parent, label = "parent directory for CSV")
  file <- file.path(parent, basename(file))

  tmp <- tempfile(pattern = paste0(".", tools::file_path_sans_ext(basename(file)), "_"),
                  tmpdir = parent, fileext = ".csv")
  on.exit(if (file.exists(tmp)) unlink(tmp, force = TRUE), add = TRUE)

  utils::write.csv(x, file = tmp, row.names = row.names, ...)

  copy_to <- function(target) {
    if (file.exists(target)) unlink(target, force = TRUE)
    isTRUE(file.copy(tmp, target, overwrite = TRUE)) && file.exists(target)
  }

  ok <- tryCatch(copy_to(file), error = function(e) {
    message("CSV write failed for ", normalizePath(file, mustWork = FALSE), ": ", conditionMessage(e))
    FALSE
  })

  if (!ok && isTRUE(fallback_on_lock)) {
    alt <- insert_file_suffix(file, paste0("_", now_stamp()))
    ok <- tryCatch(copy_to(alt), error = function(e) {
      message("Fallback CSV write failed for ", normalizePath(alt, mustWork = FALSE), ": ", conditionMessage(e))
      FALSE
    })
    if (ok) {
      warning("Could not overwrite the requested CSV path. Wrote fallback file instead: ",
              normalizePath(alt, mustWork = FALSE))
      return(invisible(alt))
    }
  }

  if (!ok) stop("Could not write CSV file: ", normalizePath(file, mustWork = FALSE))
  invisible(file)
}

#' Save an openxlsx workbook with a fallback file name when needed.
#'
#' Workbooks are often locked by spreadsheet software.  This function preserves
#' the model run by saving a timestamped fallback workbook if overwrite fails.
safe_save_workbook <- function(wb, file, overwrite = TRUE, fallback_on_lock = TRUE) {
  parent <- dirname(file)
  parent <- ensure_writable_dir(parent, label = "parent directory for workbook")
  file <- file.path(parent, basename(file))

  ok <- tryCatch({
    openxlsx::saveWorkbook(wb, file, overwrite = overwrite)
    file.exists(file)
  }, error = function(e) {
    message("Workbook write failed for ", normalizePath(file, mustWork = FALSE), ": ", conditionMessage(e))
    FALSE
  })

  if (!ok && isTRUE(fallback_on_lock)) {
    alt <- insert_file_suffix(file, paste0("_", now_stamp()))
    ok <- tryCatch({
      openxlsx::saveWorkbook(wb, alt, overwrite = TRUE)
      file.exists(alt)
    }, error = function(e) {
      message("Fallback workbook write failed for ", normalizePath(alt, mustWork = FALSE), ": ", conditionMessage(e))
      FALSE
    })
    if (ok) {
      warning("Could not overwrite the requested workbook path. Wrote fallback file instead: ",
              normalizePath(alt, mustWork = FALSE))
      return(invisible(alt))
    }
  }

  if (!ok) stop("Could not write workbook file: ", normalizePath(file, mustWork = FALSE))
  invisible(file)
}

# ------------------------------------------------------------
# Numerical and reporting utilities
# ------------------------------------------------------------
# The likelihood calculations repeatedly solve GLS normal equations and Hessian
# systems.  The small helpers here centralize safe linear algebra and reporting
# conventions.

#' Convert a p-value into conventional significance symbols.
#'
#' The output is used only for reporting parameter tables; it does not affect
#' optimization, likelihood values, or information criteria.
stars_from_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.10)  return(".")
  ""
}

#' Solve a symmetric linear system with a ridge fallback.
#'
#' GLS normal equations and Hessian inversions can be nearly singular.  The
#' matrix is symmetrized and, if direct solving fails, a small ridge is added.
safe_solve <- function(A, b = NULL, ridge = 1e-8) {
  A2 <- as.matrix(A)
  A2 <- (A2 + t(A2)) / 2
  d <- nrow(A2)
  if (is.null(b)) {
    out <- tryCatch(solve(A2), error = function(e) NULL)
    if (!is.null(out)) return(out)
    return(solve(A2 + diag(ridge, d)))
  }
  out <- tryCatch(solve(A2, b), error = function(e) NULL)
  if (!is.null(out)) return(out)
  solve(A2 + diag(ridge, d), b)
}

#' Solve a general linear system with a ridge fallback.
#'
#' This variant does not symmetrize the matrix, so it is appropriate for systems
#' such as (I - R kron W) y = X beta when producing fitted values.
safe_solve_general <- function(A, b = NULL, ridge = 1e-8) {
  A2 <- as.matrix(A)
  d <- nrow(A2)
  if (is.null(b)) {
    out <- tryCatch(solve(A2), error = function(e) NULL)
    if (!is.null(out)) return(out)
    return(solve(A2 + diag(ridge, d)))
  }
  out <- tryCatch(solve(A2, b), error = function(e) NULL)
  if (!is.null(out)) return(out)
  solve(A2 + diag(ridge, d), b)
}

#' Select a conservative number of parallel workers.
#'
#' Larger spatial grids require more memory per worker, so the function reserves
#' cores and scales down parallelism as N_VALUE increases.
choose_auto_cores <- function(n_value, n_models, requested_cores = NULL, reserve_cores = 1L) {
  detected <- tryCatch(parallel::detectCores(logical = TRUE), error = function(e) 1L)
  if (is.na(detected) || detected < 1L) detected <- 1L

  if (!is.null(requested_cores)) {
    return(max(1L, min(as.integer(requested_cores), as.integer(n_models), detected)))
  }

  usable <- max(1L, detected - as.integer(reserve_cores))

  if (n_value <= 100L) {
    chosen <- usable
  } else if (n_value <= 400L) {
    chosen <- max(1L, floor(usable * 0.75))
  } else {
    chosen <- max(1L, floor(usable * 0.50))
  }

  max(1L, min(chosen, as.integer(n_models)))
}

# ------------------------------------------------------------
# Model-restriction helpers for the 11 MSTR specifications
# ------------------------------------------------------------
# Table 1 in Subsection 3.5 assigns each four-character model ID to restrictions
# on R, Lambda, A, and Sigma.  These masks determine which parameters are
# optimized and which are fixed at zero.

#' Convert one model-ID character into a logical 2 x 2 free-parameter mask.
#'
#' Table 1 of mstr.pdf uses 1 for a full K x K matrix, d for a diagonal
#' K x K matrix, and 0 for the zero matrix.
make_mask <- function(code, K = 2L) {
  if (code == "1") return(matrix(TRUE, K, K))
  if (code == "d") return(diag(TRUE, K))
  if (code == "0") return(matrix(FALSE, K, K))
  stop("Unknown code: ", code)
}

#' Decode a four-character model ID into R, Lambda, A, and Sigma restrictions.
#'
#' The character order follows Table 1: R (spatial autoregression), Lambda
#' (spatial error), A (temporal autoregression), and Sigma (cross-equation
#' innovation covariance).
decode_model_id <- function(model_id) {
  chars <- strsplit(model_id, "")[[1]]
  stopifnot(length(chars) == 4)
  list(
    R_code = chars[1],
    L_code = chars[2],
    A_code = chars[3],
    S_code = chars[4],
    R_mask = make_mask(chars[1], 2L),
    L_mask = make_mask(chars[2], 2L)
  )
}

#' Return the regression-coefficient names implied by the A restriction.
#'
#' In the composite design representation of Equation (7), A = 1 includes all
#' lagged responses in each equation, A = d includes only own-lag terms, and
#' A = 0 excludes time-AR terms from X_t.
beta_names_by_model <- function(A_code) {
  if (A_code == "1") {
    return(c("beta_intercept_y1", "beta_common1_y1", "beta_common2_y1", "beta_specific1_1",
             "alpha[1,1]", "alpha[1,2]",
             "beta_intercept_y2", "beta_common1_y2", "beta_common2_y2", "beta_specific2_1",
             "alpha[2,1]", "alpha[2,2]"))
  }
  if (A_code == "d") {
    return(c("beta_intercept_y1", "beta_common1_y1", "beta_common2_y1", "beta_specific1_1",
             "alpha[1,1]",
             "beta_intercept_y2", "beta_common1_y2", "beta_common2_y2", "beta_specific2_1",
             "alpha[2,2]"))
  }
  if (A_code == "0") {
    return(c("beta_intercept_y1", "beta_common1_y1", "beta_common2_y1", "beta_specific1_1",
             "beta_intercept_y2", "beta_common1_y2", "beta_common2_y2", "beta_specific2_1"))
  }
  stop("Unknown A_code")
}

#' Define the unified output order for all reported parameters.
#'
#' Parameters that are fixed at zero by a model restriction are still included
#' in the output table so that all model workbooks share a comparable layout.
full_parameter_order <- function() {
  c("gamma",
    "rho[1,1]", "rho[1,2]", "rho[2,1]", "rho[2,2]",
    "lambda[1,1]", "lambda[1,2]", "lambda[2,1]", "lambda[2,2]",
    "alpha[1,1]", "alpha[1,2]", "alpha[2,1]", "alpha[2,2]",
    "sigma[1,1]", "sigma[1,2]", "sigma[2,1]", "sigma[2,2]",
    "beta_intercept_y1", "beta_intercept_y2",
    "beta_common1_y1", "beta_common1_y2",
    "beta_common2_y1", "beta_common2_y2",
    "beta_specific1_1", "beta_specific2_1")
}

#' Pack free R and Lambda entries into a vector theta.
#'
#' This vector is theta_1 in Subsection 3.6 and Appendix A: the free entries
#' of R and Lambda.  It is used in BFGS optimization, the profile Hessian in
#' Equation (38), and the effective-df calculations in Equations (24) and
#' (41)-(44).
pack_spatial_theta <- function(R, L, spec) {
  out <- c()
  for (mat_name in c("R", "L")) {
    M <- if (mat_name == "R") R else L
    mask <- if (mat_name == "R") spec$R_mask else spec$L_mask
    for (i in 1:2) for (j in 1:2) if (mask[i, j]) out <- c(out, M[i, j])
  }
  out
}

#' Convert unconstrained optimization variables into R and Lambda matrices.
#'
#' The tanh transformation enforces the constraints |rho_kl| < 1 and
#' |lambda_kl| < 1 stated in Step S4 of Subsection 3.3 while allowing
#' unconstrained BFGS updates.
unpack_spatial_u <- function(u, spec) {
  R <- matrix(0, 2, 2)
  L <- matrix(0, 2, 2)
  idx <- 1L
  for (mat_name in c("R", "L")) {
    mask <- if (mat_name == "R") spec$R_mask else spec$L_mask
    M <- matrix(0, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        if (mask[i, j]) {
          M[i, j] <- tanh(u[idx])
          idx <- idx + 1L
        }
      }
    }
    if (mat_name == "R") R <- M else L <- M
  }
  list(R = R, L = L)
}

#' Return human-readable names for the free spatial parameters.
#'
#' The names are used to align estimates, standard errors, p-values, and output
#' columns across the 11 restricted MSTR specifications.
spatial_param_names <- function(spec) {
  out <- c()
  for (pair in list(list(prefix = "rho", mask = spec$R_mask),
                    list(prefix = "lambda", mask = spec$L_mask))) {
    for (i in 1:2) for (j in 1:2) {
      if (pair$mask[i, j]) out <- c(out, sprintf("%s[%d,%d]", pair$prefix, i, j))
    }
  }
  out
}

#' Construct multi-start initial values for spatial-parameter optimization.
#'
#' The profile likelihood can have local optima.  This function combines
#' univariate initial values, zero starts, warm starts from adjacent gamma
#' values, sign perturbations, and optional aggressive starts for retries.
build_multistart_theta <- function(spec, init_R, init_L,
                                   warm_theta = NULL,
                                   aggressive = FALSE) {
  base <- pack_spatial_theta(init_R, init_L, spec)
  m <- length(base)
  if (m == 0L) return(list(numeric(0)))

  add_start <- function(store, th) {
    th2 <- pmax(pmin(as.numeric(th), 0.90), -0.90)
    key <- paste(sprintf("%.6f", th2), collapse = "|")
    if (!key %in% names(store)) store[[key]] <- th2
    store
  }

  starts <- list()
  if (!is.null(warm_theta) && length(warm_theta) == m) {
    starts <- add_start(starts, warm_theta)
    starts <- add_start(starts, warm_theta + 0.02)
    starts <- add_start(starts, warm_theta - 0.02)
  }

  starts <- add_start(starts, base)
  starts <- add_start(starts, rep(0, m))
  starts <- add_start(starts, base + 0.05)
  starts <- add_start(starts, base - 0.05)
  starts <- add_start(starts, 0.5 * base)
  starts <- add_start(starts, -base)

  spatial_names <- spatial_param_names(spec)
  diag_like <- base
  if (length(spatial_names) > 0) {
    off_idx <- grepl("[1,2]", spatial_names, fixed = TRUE) | grepl("[2,1]", spatial_names, fixed = TRUE)
    if (any(off_idx)) diag_like[off_idx] <- 0
    starts <- add_start(starts, diag_like)
  }

  if (spec$R_code == "1" && spec$L_code == "1") {
    signs <- expand.grid(rep(list(c(-1, 1)), m), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    max_extra <- min(8L, nrow(signs))
    for (ii in seq_len(max_extra)) {
      starts <- add_start(starts, 0.15 * as.numeric(signs[ii, ]))
    }
  } else {
    starts <- add_start(starts, rep(0.15, m))
    starts <- add_start(starts, rep(-0.15, m))
    starts <- add_start(starts, rep(c(0.15, -0.15), length.out = m))
  }

  if (isTRUE(aggressive)) {
    coarse_grid <- c(-0.60, -0.35, -0.20, 0.20, 0.35, 0.60)
    for (val in coarse_grid) starts <- add_start(starts, rep(val, m))
    starts <- add_start(starts, rep(c(-0.60, 0.60), length.out = m))
    starts <- add_start(starts, rep(c(0.60, -0.60), length.out = m))
    starts <- add_start(starts, rep(c(-0.35, 0.35), length.out = m))
    starts <- add_start(starts, rep(c(0.35, -0.35), length.out = m))
    starts <- add_start(starts, pmax(pmin(base + 0.15, 0.90), -0.90))
    starts <- add_start(starts, pmax(pmin(base - 0.15, 0.90), -0.90))
    starts <- add_start(starts, pmax(pmin(1.5 * base, 0.90), -0.90))
    starts <- add_start(starts, pmax(pmin(-1.5 * base, 0.90), -0.90))
  }

  unname(starts)
}

#' Compute log |I - theta W| from the eigenvalues of W.
#'
#' This is used for the univariate GNS initialization in Equation (17).  It is
#' the scalar counterpart of the eigenvalue determinant calculation described
#' in Subsection 3.4.
logdet_univ <- function(theta, eigW, eps = 1e-10) {
  vals <- 1 - eigW * theta
  if (any(abs(vals) < eps)) return(-Inf)
  sum(log(abs(vals)))
}

#' Compute log |I - M kron W| for a 2 x 2 spatial matrix M.
#'
#' For K = 2 and each eigenvalue omega_i of W, the Kronecker determinant
#' reduces to |I_2 - omega_i M|.  The product over i is the K = 2 form of
#' Equation (18); Equations (19)-(23) give the corresponding trace expansion
#' for general K.
logdet_kron2 <- function(M, eigW, eps = 1e-10) {
  out <- 0
  I2 <- diag(2)
  for (w in eigW) {
    d <- determinant(I2 - w * M, logarithm = TRUE)
    if (abs(d$modulus) == Inf || d$sign == 0) return(-Inf)
    out <- out + as.numeric(d$modulus)
  }
  out
}

# ------------------------------------------------------------
# Data preparation and cached objects
# ------------------------------------------------------------
# This block reads panel data and W, validates row normalization, constructs
# equation-specific design matrices, and caches eigenvalues/initial values used
# by all model fits.

#' Check and, when necessary, row-standardize the spatial weight matrix W.
#'
#' Section 1 and Equation (5) use W as a row-normalized adjacency matrix.
#' This function records diagnostics before and after normalization and leaves
#' zero-neighbor rows unchanged.
standardize_W_rows <- function(W, tol = W_ROW_TOL, standardize = W_ROW_STANDARDIZE) {
  W <- as.matrix(W)
  storage.mode(W) <- "numeric"
  if (any(!is.finite(W))) stop("W contains NA/NaN/Inf values after numeric conversion.")

  row_sums_before <- rowSums(W)
  zero_rows <- which(abs(row_sums_before) <= tol)
  nonzero_rows <- setdiff(seq_len(nrow(W)), zero_rows)
  max_dev_before <- if (length(nonzero_rows) > 0L) {
    max(abs(row_sums_before[nonzero_rows] - 1))
  } else {
    NA_real_
  }

  standardized <- FALSE
  if (isTRUE(standardize) && length(nonzero_rows) > 0L &&
      is.finite(max_dev_before) && max_dev_before > tol) {
    W[nonzero_rows, ] <- sweep(W[nonzero_rows, , drop = FALSE], 1,
                               row_sums_before[nonzero_rows], `/`)
    standardized <- TRUE
  }

  row_sums_after <- rowSums(W)
  max_dev_after <- if (length(nonzero_rows) > 0L) {
    max(abs(row_sums_after[nonzero_rows] - 1))
  } else {
    NA_real_
  }

  attr(W, "row_standardization") <- data.frame(
    n = nrow(W),
    zero_row_count = length(zero_rows),
    min_row_sum_before = if (length(row_sums_before) > 0L) min(row_sums_before) else NA_real_,
    max_row_sum_before = if (length(row_sums_before) > 0L) max(row_sums_before) else NA_real_,
    max_abs_deviation_before = max_dev_before,
    row_standardized = standardized,
    min_row_sum_after = if (length(row_sums_after) > 0L) min(row_sums_after) else NA_real_,
    max_row_sum_after = if (length(row_sums_after) > 0L) max(row_sums_after) else NA_real_,
    max_abs_deviation_after = max_dev_after,
    tolerance = tol,
    stringsAsFactors = FALSE
  )
  W
}

#' Read a numeric square matrix from a possibly labeled CSV file.
#'
#' Header rows, row-name columns, and nonnumeric label artifacts are removed so
#' that spatial-weight files exported from spreadsheet software can be used.
safe_numeric_matrix_from_csv <- function(path) {
  raw <- read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  M <- suppressWarnings(apply(raw, 2, as.numeric))
  M <- as.matrix(M)

  # Drop header rows read as data, for example "V1,V2,...".
  all_na_rows <- which(rowSums(is.na(M)) == ncol(M))
  if (length(all_na_rows) > 0L) M <- M[-all_na_rows, , drop = FALSE]

  # Drop row-name columns, if present.
  all_na_cols <- which(colSums(is.na(M)) == nrow(M))
  if (length(all_na_cols) > 0L) M <- M[, -all_na_cols, drop = FALSE]

  # Handle one extra row or column caused by labels.
  if (nrow(M) == ncol(M) + 1L) M <- M[-1, , drop = FALSE]
  if (ncol(M) == nrow(M) + 1L) M <- M[, -1, drop = FALSE]

  if (nrow(M) != ncol(M)) {
    stop(sprintf("W must be square after reading. Got %d x %d from %s",
                 nrow(M), ncol(M), path))
  }
  if (any(!is.finite(M))) stop("W contains non-numeric entries after header/label removal: ", path)
  M
}

#' Read the simulated panel data, spatial weights, W eigenvalues, and template path.
#'
#' The eigenvalues of W are cached because they are repeatedly needed for the
#' spatial determinant terms in the likelihood.
read_inputs <- function(data_file, weights_file, template_file = TEMPLATE_FILE) {
  dat <- read.csv(data_file, stringsAsFactors = FALSE)
  W <- safe_numeric_matrix_from_csv(weights_file)
  W <- standardize_W_rows(W, tol = W_ROW_TOL, standardize = W_ROW_STANDARDIZE)
  W_stats <- attr(W, "row_standardization")

  template_path <- if (!is.null(template_file) && length(template_file) == 1L &&
                       !is.na(template_file) && file.exists(template_file)) {
    normalizePath(template_file, mustWork = TRUE)
  } else {
    NA_character_
  }

  list(
    dat = dat,
    W = W,
    W_stats = W_stats,
    eigW = Re(eigen(W, symmetric = FALSE, only.values = TRUE)$values),
    template_path = template_path
  )
}

#' Read and validate the gamma grid from code, CSV, workbook, or auto mode.
#'
#' Gamma values are filtered to finite nonnegative values and sorted uniquely
#' before model-specific grids are constructed.
read_gamma_values <- function(gamma_values = GAMMA_VALUES,
                              gamma_csv_file = GAMMA_CSV_FILE,
                              template_file = TEMPLATE_FILE,
                              source = GAMMA_SOURCE) {
  source <- match.arg(source, c("auto", "code", "csv", "xlsx"))

  read_from_csv <- function(path) {
    if (!file.exists(path)) stop("Gamma CSV file not found: ", normalizePath(path, mustWork = FALSE))
    gdat_h <- read.csv(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    col <- which(tolower(names(gdat_h)) == "gamma")
    if (ncol(gdat_h) >= 1L && length(col) > 0L) {
      return(as.numeric(gdat_h[[col[1]]]))
    }
    # Headerless CSV is also accepted; use the first column and do not drop row 1.
    gdat <- read.csv(path, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(gdat) < 1L) stop("Gamma CSV has no columns: ", path)
    as.numeric(gdat[[1]])
  }

  read_from_xlsx <- function(path) {
    if (!file.exists(path)) stop("Gamma xlsx/template file not found: ", normalizePath(path, mustWork = FALSE))
    gamma_values <- openxlsx::readWorkbook(path, sheet = "gamma_grid",
                                           cols = 1, rows = 2:10000,
                                           colNames = FALSE)[[1]]
    as.numeric(gamma_values)
  }

  vals <- NULL
  gamma_source_used <- NULL
  if (source == "csv") {
    vals <- read_from_csv(gamma_csv_file)
    gamma_source_used <- "csv"
  } else if (source == "xlsx") {
    vals <- read_from_xlsx(template_file)
    gamma_source_used <- "xlsx"
  } else if (source == "code") {
    vals <- gamma_values
    gamma_source_used <- "code"
  } else {
    if (!is.null(gamma_csv_file) && file.exists(gamma_csv_file)) {
      vals <- read_from_csv(gamma_csv_file)
      gamma_source_used <- "csv"
    } else if (!is.null(gamma_values) && length(gamma_values) > 0L) {
      vals <- gamma_values
      gamma_source_used <- "code"
    } else if (!is.null(template_file) && file.exists(template_file)) {
      vals <- read_from_xlsx(template_file)
      gamma_source_used <- "xlsx"
    } else {
      stop("No gamma grid is available. Set GAMMA_VALUES or provide gamma_grid.csv.")
    }
  }

  vals <- suppressWarnings(as.numeric(vals))
  vals <- sort(unique(vals[is.finite(vals) & vals >= 0]))
  if (length(vals) == 0L) stop("Gamma grid is empty after filtering finite gamma >= 0.")
  attr(vals, "source_used") <- gamma_source_used
  vals
}

#' Build the response vectors and equation-specific design matrices for one time point.
#'
#' The lagged responses at TARGET_TIME - 1 are appended according to the A-code,
#' implementing the composite X_t beta representation in Equation (7), which
#' leads to the standard MGNS form in Equation (8).
build_period_data <- function(dat, A_code, target_time = TARGET_TIME) {
  if (!(target_time %in% dat$time)) stop("target_time is not found in data: ", target_time)
  if (!((target_time - 1L) %in% dat$time)) stop("lagged time is not found in data: ", target_time - 1L)

  cur <- dat[dat$time == target_time, ]
  lag <- dat[dat$time == (target_time - 1L), ]
  cur <- cur[order(cur$region), ]
  lag <- lag[order(lag$region), ]

  X1_base <- cbind(1, cur$x_common1, cur$x_common2, cur$x_specific1_1)
  X2_base <- cbind(1, cur$x_common1, cur$x_common2, cur$x_specific2_1)

  if (A_code == "1") {
    X1 <- cbind(X1_base, lag$y1, lag$y2)
    X2 <- cbind(X2_base, lag$y1, lag$y2)
  } else if (A_code == "d") {
    X1 <- cbind(X1_base, lag$y1)
    X2 <- cbind(X2_base, lag$y2)
  } else if (A_code == "0") {
    X1 <- X1_base
    X2 <- X2_base
  } else {
    stop("Unknown A_code")
  }

  list(list(
    y1 = as.numeric(cur$y1),
    y2 = as.numeric(cur$y2),
    X1 = as.matrix(X1),
    X2 = as.matrix(X2)
  ))
}

# ------------------------------------------------------------
# Step S1 initialization: univariate GNS fits (Equation (17))
# ------------------------------------------------------------
# Step S1 in Subsection 3.3 starts from separate univariate GNS fits defined by
# Equation (17).  Their diagonal rho/lambda/sigma estimates initialize R,
# Lambda, and Sigma for the multivariate optimization.

#' Estimate a univariate GNS model for initialization of one response equation.
#'
#' This implements Step S1(a)-(b) of Subsection 3.3 using Equation (17): it
#' obtains diagonal starting values for rho, lambda, and sigma before fitting
#' the restricted multivariate model.
fit_univariate_gns <- function(eq_index, periods, W, eigW, control = list()) {
  n <- nrow(W)
  TT <- length(periods)
  N <- n * TT
  spec_col <- if (eq_index == 1) "X1" else "X2"

  y_list <- lapply(periods, function(z) if (eq_index == 1) z$y1 else z$y2)
  X_list <- lapply(periods, function(z) z[[spec_col]])

  profile_ll <- function(rho, lambda) {
    if (max(abs(rho), abs(lambda)) >= 0.999999) return(list(ok = FALSE))
    H <- diag(n) - rho * W
    L <- diag(n) - lambda * W
    ld <- TT * (logdet_univ(rho, eigW) + logdet_univ(lambda, eigW))
    if (!is.finite(ld)) return(list(ok = FALSE))

    p <- ncol(X_list[[1]])
    M <- matrix(0, p, p)
    b <- rep(0, p)
    for (tt in seq_len(TT)) {
      LX <- L %*% X_list[[tt]]
      Ly <- L %*% (H %*% y_list[[tt]])
      M <- M + crossprod(LX)
      b <- b + crossprod(LX, Ly)
    }
    # Conditional GLS estimate for the univariate GNS model in Equation (17).
    beta <- as.numeric(safe_solve(M, b))
    rss <- 0
    for (tt in seq_len(TT)) {
      z <- L %*% (H %*% y_list[[tt]] - X_list[[tt]] %*% beta)
      rss <- rss + sum(z^2)
    }
    sigma2 <- rss / N
    if (!is.finite(sigma2) || sigma2 <= 0) return(list(ok = FALSE))
    ll <- -N / 2 * (log(2 * pi) + 1 + log(sigma2)) + ld
    list(ok = TRUE, logLik = ll, beta = beta, sigma2 = sigma2)
  }

  obj <- function(u) {
    rho <- tanh(u[1])
    lambda <- tanh(u[2])
    ans <- profile_ll(rho, lambda)
    if (!ans$ok) return(1e100)
    -ans$logLik
  }

  starts <- list(c(0, 0), c(0.2, 0.2), c(-0.2, 0.2), c(0.3, -0.1))
  best <- NULL
  for (st in starts) {
    fit <- optim(st, obj, method = "BFGS",
                 control = modifyList(list(maxit = 200, reltol = 1e-8), control))
    if (is.null(best) || fit$value < best$value) best <- fit
  }

  rho_hat <- tanh(best$par[1])
  lambda_hat <- tanh(best$par[2])
  prof <- profile_ll(rho_hat, lambda_hat)
  list(rho = rho_hat, lambda = lambda_hat, beta = prof$beta,
       sigma2 = prof$sigma2, logLik = prof$logLik)
}

#' Assemble data objects shared across all model IDs.
#'
#' Data, W, eigenvalues, gamma grid, period-specific design matrices, and
#' univariate initial values are cached once to keep the 11-model run efficient.
build_common_inputs <- function(data_file = DATA_FILE,
                                weights_file = WEIGHTS_FILE,
                                template_file = TEMPLATE_FILE,
                                target_time = TARGET_TIME) {
  inp <- read_inputs(data_file, weights_file, template_file)
  gamma_values <- read_gamma_values(gamma_values = GAMMA_VALUES,
                                    gamma_csv_file = GAMMA_CSV_FILE,
                                    template_file = template_file,
                                    source = GAMMA_SOURCE)

  periods_cache <- list(
    "0" = build_period_data(inp$dat, "0", target_time = target_time),
    "d" = build_period_data(inp$dat, "d", target_time = target_time),
    "1" = build_period_data(inp$dat, "1", target_time = target_time)
  )

  init_cache <- list()
  for (A_code in c("0", "d", "1")) {
    periods <- periods_cache[[A_code]]
    uni1 <- fit_univariate_gns(1, periods, inp$W, inp$eigW)
    uni2 <- fit_univariate_gns(2, periods, inp$W, inp$eigW)

    R0 <- matrix(0, 2, 2)
    L0 <- matrix(0, 2, 2)
    S0 <- diag(2)
    R0[1, 1] <- uni1$rho
    R0[2, 2] <- uni2$rho
    L0[1, 1] <- uni1$lambda
    L0[2, 2] <- uni2$lambda
    S0[1, 1] <- uni1$sigma2
    S0[2, 2] <- uni2$sigma2

    init_cache[[A_code]] <- list(R0 = R0, L0 = L0, S0 = S0)
  }

  list(
    dat = inp$dat,
    W = inp$W,
    W_stats = inp$W_stats,
    eigW = inp$eigW,
    template_path = inp$template_path,
    gamma_values = gamma_values,
    gamma_source = attr(gamma_values, "source_used"),
    periods_cache = periods_cache,
    init_cache = init_cache,
    target_time = target_time
  )
}

# ------------------------------------------------------------
# Profile likelihood for fixed R and Lambda: Equations (9)-(16), Steps S2-S3
# ------------------------------------------------------------
# For fixed R and Lambda, the MGNS model in Equation (8) becomes a multivariate
# GLS problem with correlated innovations.  Step S2 alternates the beta update
# in Equations (13)-(14) and the Sigma update in Equation (15); Step S3 evaluates
# the profiled log-likelihood derived from Equations (9)-(10).

#' Profile beta and Sigma for fixed R and Lambda by iterative GLS.
#'
#' This function implements the likelihood equations (11)-(12), the beta
#' update (13)-(14), the Sigma update (15), and Var(beta_hat) in Equation (16).
#' It evaluates the log-likelihood and residual quadratic form in Equations
#' (9)-(10), and obtains fitted means from Equation (8) for the pseudo R-squared
#' calculation in Equation (26).
profile_beta_sigma <- function(R, L, periods, W, eigW, spec, Sigma_init = NULL,
                               max_gls_iter = 150, tol = 1e-8) {
  n <- nrow(W)
  TT <- length(periods)
  N <- n * TT
  p1 <- ncol(periods[[1]]$X1)
  p2 <- ncol(periods[[1]]$X2)
  p <- p1 + p2

  # The two Jacobian log-determinants are the second and third terms of
  # Equation (9); logdet_kron2 evaluates them using Equation (18).
  ldH <- TT * logdet_kron2(R, eigW)
  ldL <- TT * logdet_kron2(L, eigW)
  if (!is.finite(ldH) || !is.finite(ldL)) return(list(ok = FALSE))

  Sigma <- if (is.null(Sigma_init)) diag(2) else Sigma_init
  if (spec$S_code == "d") Sigma <- diag(diag(Sigma))

  stabilize_sigma <- function(S) {
    S2 <- (as.matrix(S) + t(as.matrix(S))) / 2
    if (spec$S_code == "d") S2 <- diag(diag(S2))
    ev <- tryCatch(eigen(S2, symmetric = TRUE, only.values = TRUE)$values,
                   error = function(e) c(-Inf, -Inf))
    min_ev <- min(ev)
    if (!is.finite(min_ev)) min_ev <- -Inf
    if (min_ev <= 1e-8) {
      bump <- max(1e-8, 1e-6 - min_ev)
      S2 <- S2 + diag(bump, nrow(S2))
      if (spec$S_code == "d") S2 <- diag(diag(S2))
    }
    S2
  }

  Sigma <- stabilize_sigma(Sigma)
  last_beta <- NULL
  last_Sigma <- NULL

  make_transformed <- function(X1, X2, q1, q2, L) {
    WX1 <- W %*% X1
    WX2 <- W %*% X2
    Wq1 <- W %*% q1
    Wq2 <- W %*% q2

    U11 <- X1 - L[1, 1] * WX1
    U12 <- -L[1, 2] * WX2
    U21 <- -L[2, 1] * WX1
    U22 <- X2 - L[2, 2] * WX2

    U1 <- cbind(U11, U12)
    U2 <- cbind(U21, U22)
    v1 <- q1 - L[1, 1] * Wq1 - L[1, 2] * Wq2
    v2 <- q2 - L[2, 1] * Wq1 - L[2, 2] * Wq2
    list(U1 = U1, U2 = U2, v1 = v1, v2 = v2)
  }

  # Step S2(c)-(d): alternate the beta update in Equations (13)-(14) and
  # the Sigma update in Equation (15) until convergence.
  for (iter in seq_len(max_gls_iter)) {
    Sigma <- stabilize_sigma(Sigma)
    Sinv <- tryCatch(safe_solve(Sigma), error = function(e) NULL)
    if (is.null(Sinv)) return(list(ok = FALSE))

    M <- matrix(0, p, p)
    b <- rep(0, p)

    for (tt in seq_len(TT)) {
      y1 <- periods[[tt]]$y1; y2 <- periods[[tt]]$y2
      X1 <- periods[[tt]]$X1; X2 <- periods[[tt]]$X2
      Wy1 <- W %*% y1; Wy2 <- W %*% y2
      q1 <- y1 - R[1, 1] * Wy1 - R[1, 2] * Wy2
      q2 <- y2 - R[2, 1] * Wy1 - R[2, 2] * Wy2

      trn <- make_transformed(X1, X2, q1, q2, L)
      O1 <- Sinv[1, 1] * trn$U1 + Sinv[1, 2] * trn$U2
      O2 <- Sinv[2, 1] * trn$U1 + Sinv[2, 2] * trn$U2
      M <- M + crossprod(trn$U1, O1) + crossprod(trn$U2, O2)
      b <- b + crossprod(trn$U1, Sinv[1, 1] * trn$v1 + Sinv[1, 2] * trn$v2) +
        crossprod(trn$U2, Sinv[2, 1] * trn$v1 + Sinv[2, 2] * trn$v2)
    }

    # GLS solution beta_hat = B y from Equations (13)-(14).
    beta <- as.numeric(safe_solve(M, b))
    beta1 <- beta[seq_len(p1)]
    beta2 <- beta[p1 + seq_len(p2)]

    S <- matrix(0, 2, 2)
    for (tt in seq_len(TT)) {
      y1 <- periods[[tt]]$y1; y2 <- periods[[tt]]$y2
      X1 <- periods[[tt]]$X1; X2 <- periods[[tt]]$X2
      Wy1 <- W %*% y1; Wy2 <- W %*% y2
      q1 <- y1 - R[1, 1] * Wy1 - R[1, 2] * Wy2
      q2 <- y2 - R[2, 1] * Wy1 - R[2, 2] * Wy2
      r1 <- q1 - X1 %*% beta1
      r2 <- q2 - X2 %*% beta2
      Wr1 <- W %*% r1; Wr2 <- W %*% r2
      z1 <- r1 - L[1, 1] * Wr1 - L[1, 2] * Wr2
      z2 <- r2 - L[2, 1] * Wr1 - L[2, 2] * Wr2
      S <- S + matrix(c(sum(z1 * z1), sum(z1 * z2), sum(z2 * z1), sum(z2 * z2)), 2, 2, byrow = TRUE)
    }

    # Equation (15), pooled over TT periods; when TT = 1, N = n exactly.
    Sigma_new <- S / N
    Sigma_new <- (Sigma_new + t(Sigma_new)) / 2
    if (spec$S_code == "d") Sigma_new <- diag(diag(Sigma_new))
    Sigma_new <- stabilize_sigma(Sigma_new)

    ev <- eigen(Sigma_new, symmetric = TRUE, only.values = TRUE)$values
    if (any(!is.finite(ev)) || min(ev) <= 0) return(list(ok = FALSE))

    if (!is.null(last_beta) && max(abs(beta - last_beta), abs(Sigma_new - last_Sigma)) < tol) {
      Sigma <- Sigma_new
      break
    }
    last_beta <- beta
    last_Sigma <- Sigma_new
    Sigma <- Sigma_new
  }

  Sinv <- safe_solve(Sigma)
  # Q and z are evaluated according to Equation (10).
  Q <- 0
  M_beta <- matrix(0, p, p)
  fitted_all_y1 <- c(); fitted_all_y2 <- c()
  actual_all_y1 <- c(); actual_all_y2 <- c()

  for (tt in seq_len(TT)) {
    y1 <- periods[[tt]]$y1; y2 <- periods[[tt]]$y2
    X1 <- periods[[tt]]$X1; X2 <- periods[[tt]]$X2
    beta1 <- beta[seq_len(p1)]
    beta2 <- beta[p1 + seq_len(p2)]
    Wy1 <- W %*% y1; Wy2 <- W %*% y2
    q1 <- y1 - R[1, 1] * Wy1 - R[1, 2] * Wy2
    q2 <- y2 - R[2, 1] * Wy1 - R[2, 2] * Wy2

    trn <- make_transformed(X1, X2, q1, q2, L)
    O1 <- Sinv[1, 1] * trn$U1 + Sinv[1, 2] * trn$U2
    O2 <- Sinv[2, 1] * trn$U1 + Sinv[2, 2] * trn$U2
    M_beta <- M_beta + crossprod(trn$U1, O1) + crossprod(trn$U2, O2)

    r1 <- q1 - X1 %*% beta1
    r2 <- q2 - X2 %*% beta2
    Wr1 <- W %*% r1; Wr2 <- W %*% r2
    z1 <- r1 - L[1, 1] * Wr1 - L[1, 2] * Wr2
    z2 <- r2 - L[2, 1] * Wr1 - L[2, 2] * Wr2
    Q <- Q + sum(cbind(z1, z2) * (cbind(z1, z2) %*% Sinv))

    rhs <- c(X1 %*% beta1, X2 %*% beta2)
    H_big <- diag(2 * n) - kronecker(R, W)
    yhat <- as.numeric(safe_solve_general(H_big, rhs))
    fitted_all_y1 <- c(fitted_all_y1, yhat[1:n])
    fitted_all_y2 <- c(fitted_all_y2, yhat[(n + 1):(2 * n)])
    actual_all_y1 <- c(actual_all_y1, y1)
    actual_all_y2 <- c(actual_all_y2, y2)
  }

  ldS <- determinant(Sigma, logarithm = TRUE)
  if (ldS$sign <= 0) return(list(ok = FALSE))
  # Equation (9), with K = 2 and N = n * TT.
  logLik <- -(N * 2 / 2) * log(2 * pi) + ldH + ldL - (N / 2) * as.numeric(ldS$modulus) - 0.5 * Q
  # Equation (16): covariance matrix of beta_hat.
  Psi_beta <- safe_solve(M_beta)

  list(ok = TRUE, R = R, L = L, beta = beta, Sigma = Sigma,
       logLik = logLik, Psi_beta = Psi_beta,
       fitted_y1 = fitted_all_y1, fitted_y2 = fitted_all_y2,
       actual_y1 = actual_all_y1, actual_y2 = actual_all_y2,
       p_beta = p, n_eff = N)
}

# ------------------------------------------------------------
# One-model / one-gamma estimation: Appendix A and Steps S4-S5
# ------------------------------------------------------------
# This block maximizes the penalized profile log-likelihood in Equation (36)
# over the free spatial parameters.  BFGS implements Steps S4-S5 of Subsection
# 3.3, with beta and Sigma profiled as in Equation (35).

#' Fit one restricted MSTR/MGNS model at one gamma value.
#'
#' The function decodes the Table 1 restrictions, profiles beta and Sigma as
#' in Equation (35), maximizes the penalized objective in Equation (36) by
#' BFGS, and returns the quantities required for Equations (24)-(26).
fit_mstr_model <- function(model_id, common_inputs,
                           gamma = 0,
                           warm_theta = NULL,
                           control_optim = list(maxit = MAXIT_MAIN, reltol = REL_TOL),
                           aggressive_starts = FALSE) {
  spec <- decode_model_id(model_id)
  use_gamma_penalty <- uses_gamma_penalty_model(model_id)
  gamma_fit <- if (isTRUE(use_gamma_penalty)) gamma else 0
  A_code <- spec$A_code
  periods <- common_inputs$periods_cache[[A_code]]
  W <- common_inputs$W
  eigW <- common_inputs$eigW
  init_vals <- common_inputs$init_cache[[A_code]]

  p_beta <- length(beta_names_by_model(spec$A_code))
  sigma_free <- if (spec$S_code == "1") 3L else 2L
  beta_sigma_count <- p_beta + sigma_free
  m_spatial <- sum(spec$R_mask) + sum(spec$L_mask)
  spatial_names <- spatial_param_names(spec)

  if (m_spatial == 0L) {
    prof <- profile_beta_sigma(matrix(0, 2, 2), matrix(0, 2, 2), periods, W, eigW, spec, Sigma_init = init_vals$S0)
    if (!prof$ok) stop("Profile evaluation failed for model ", model_id)
    return(list(model_id = model_id, gamma = gamma_fit, spec = spec,
                uses_gamma_penalty = use_gamma_penalty,
                spatial_names = character(0), spatial_theta = numeric(0),
                spatial_vcov = matrix(0, 0, 0), I11 = matrix(0, 0, 0),
                u_hat = numeric(0), fit = prof,
                beta_sigma_count = beta_sigma_count,
                nominal_param_count = beta_sigma_count,
                free_spatial_param_count = 0L,
                optim = list(convergence = NA_integer_, converged = NA,
                             n_starts = 0L, n_successful_starts = 0L,
                             n_converged_starts = 0L, value = NA_real_,
                             counts = NA, message = NA_character_, retry_used = FALSE)))
  }

  starts <- build_multistart_theta(spec, init_vals$R0, init_vals$L0,
                                   warm_theta = warm_theta,
                                   aggressive = aggressive_starts)

  # Negative unpenalized profile log-likelihood, Equation (35).
  obj_unpen <- function(u) {
    mats <- unpack_spatial_u(u, spec)
    prof <- profile_beta_sigma(mats$R, mats$L, periods, W, eigW, spec, Sigma_init = init_vals$S0)
    if (!prof$ok || !is.finite(prof$logLik)) return(1e100)
    -prof$logLik
  }

  # Negative penalized profile log-likelihood, Equation (36):
  #   -{ell_c(theta_1) - gamma * theta_1' theta_1 / 2}.
  obj_pen <- function(u) {
    mats <- unpack_spatial_u(u, spec)
    prof <- profile_beta_sigma(mats$R, mats$L, periods, W, eigW, spec, Sigma_init = init_vals$S0)
    if (!prof$ok || !is.finite(prof$logLik)) return(1e100)
    theta <- pack_spatial_theta(mats$R, mats$L, spec)
    penalty <- 0.5 * gamma_fit * sum(theta^2)
    -(prof$logLik - penalty)
  }

  best <- NULL
  n_successful_starts <- 0L
  n_converged_starts <- 0L
  for (th in starts) {
    u0 <- atanh(pmax(pmin(th, 0.98), -0.98))
    fit <- tryCatch(optim(u0, obj_pen, method = "BFGS", control = control_optim),
                    error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value)) {
      n_successful_starts <- n_successful_starts + 1L
      if (identical(as.integer(fit$convergence), 0L)) n_converged_starts <- n_converged_starts + 1L
      if (is.null(best) || fit$value < best$value) best <- fit
    }
  }
  if (is.null(best)) stop("All multi-start optimizations failed for model ", model_id, " at gamma = ", gamma_fit)

  u_hat <- best$par
  mats <- unpack_spatial_u(u_hat, spec)
  prof <- profile_beta_sigma(mats$R, mats$L, periods, W, eigW, spec, Sigma_init = init_vals$S0)
  if (!prof$ok) stop("Final profile evaluation failed for model ", model_id, " at gamma = ", gamma_fit)

  # Do not compute the expensive profile Hessian in Equation (38) here.
  # Equation (24) and Appendix A, Equations (41)-(43), use this Hessian for
  # pAIC/pBIC, so it is evaluated later
  # only for the gamma values needed for exact selection (or for all values
  # when DEFF_EVALUATION_MODE = "all").
  theta_hat <- pack_spatial_theta(mats$R, mats$L, spec)
  vcov_theta <- matrix(NA_real_, length(theta_hat), length(theta_hat))
  I11 <- matrix(NA_real_, length(theta_hat), length(theta_hat))

  list(model_id = model_id, gamma = gamma_fit, spec = spec,
       uses_gamma_penalty = use_gamma_penalty,
       spatial_names = spatial_names, spatial_theta = theta_hat,
       spatial_vcov = vcov_theta, I11 = I11, u_hat = u_hat,
       fit = prof, beta_sigma_count = beta_sigma_count,
       nominal_param_count = beta_sigma_count + m_spatial,
       free_spatial_param_count = m_spatial,
       optim = list(convergence = as.integer(best$convergence),
                    converged = identical(as.integer(best$convergence), 0L),
                    n_starts = length(starts),
                    n_successful_starts = n_successful_starts,
                    n_converged_starts = n_converged_starts,
                    value = best$value,
                    counts = best$counts,
                    message = if (is.null(best$message)) NA_character_ else as.character(best$message),
                    retry_used = FALSE))
}

# ------------------------------------------------------------
# Inference and model criteria: Equations (16), (24)-(26), and (38)-(44)
# ------------------------------------------------------------
# This block constructs standard errors, the profile-Hessian effective degrees
# of freedom, pAIC/pBIC, ordinary AIC/BIC, and pseudo R-squared summaries using
# the current equation numbering in mstr.pdf.

#' Approximate standard errors for the free entries of Sigma.
#'
#' The table respects whether Sigma is full or diagonal and is used only for
#' reporting covariance-parameter p-values.
sigma_se_table <- function(Sigma, n_eff, S_code) {
  out <- list(
    `sigma[1,1]` = NA_real_,
    `sigma[1,2]` = NA_real_,
    `sigma[2,1]` = NA_real_,
    `sigma[2,2]` = NA_real_
  )
  out[["sigma[1,1]"]] <- sqrt(2 * Sigma[1, 1]^2 / n_eff)
  out[["sigma[2,2]"]] <- sqrt(2 * Sigma[2, 2]^2 / n_eff)
  if (S_code == "1") {
    se12 <- sqrt((Sigma[1, 2]^2 + Sigma[1, 1] * Sigma[2, 2]) / n_eff)
    out[["sigma[1,2]"]] <- se12
    out[["sigma[2,1]"]] <- se12
  }
  out
}

#' Convert a constrained theta vector into R and Lambda matrices.
#'
#' Unlike unpack_spatial_u(), this function expects theta already on the
#' coefficient scale and is used by profile-Hessian calculations.
unpack_spatial_theta <- function(theta, spec) {
  R <- matrix(0, 2, 2)
  L <- matrix(0, 2, 2)
  idx <- 1L
  for (mat_name in c("R", "L")) {
    mask <- if (mat_name == "R") spec$R_mask else spec$L_mask
    M <- matrix(0, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        if (mask[i, j]) {
          M[i, j] <- theta[idx]
          idx <- idx + 1L
        }
      }
    }
    if (mat_name == "R") R <- M else L <- M
  }
  list(R = R, L = L)
}

#' Return the index of the smallest finite value in a numeric vector.
#'
#' NA is returned when no finite value is available, preventing accidental
#' selection from failed gamma fits.
which_min_finite <- function(x) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  if (!any(ok)) return(NA_integer_)
  which(ok)[which.min(x[ok])]
}

#' Compute the profile-likelihood Hessian for free spatial parameters.
#'
#' At each theta_1 value, beta and Sigma are re-estimated as in Equation (35).
#' The returned matrix is the negative Hessian I_c(theta_1) in Equation (38),
#' denoted I11 in the main-text effective-df formula (24).
theta_hessian_profile <- function(fit_obj, common_inputs,
                                  max_gls_iter = HESSIAN_MAX_GLS_ITER,
                                  tol = HESSIAN_GLS_TOL,
                                  theta_bound = HESSIAN_THETA_BOUND) {
  m <- fit_obj$free_spatial_param_count
  if (m == 0L) return(matrix(0, 0, 0))

  spec <- fit_obj$spec
  periods <- common_inputs$periods_cache[[spec$A_code]]
  W <- common_inputs$W
  eigW <- common_inputs$eigW
  theta_hat <- as.numeric(fit_obj$spatial_theta)
  Sigma_start <- fit_obj$fit$Sigma

  obj_theta_unpen <- function(theta) {
    theta <- as.numeric(theta)
    if (length(theta) != m || any(!is.finite(theta)) || any(abs(theta) >= theta_bound)) {
      return(1e100)
    }
    mats <- unpack_spatial_theta(theta, spec)
    prof <- tryCatch(
      profile_beta_sigma(mats$R, mats$L, periods, W, eigW, spec,
                         Sigma_init = Sigma_start,
                         max_gls_iter = max_gls_iter, tol = tol),
      error = function(e) NULL
    )
    if (is.null(prof) || !isTRUE(prof$ok) || !is.finite(prof$logLik)) return(1e100)
    -prof$logLik
  }

  H <- tryCatch({
    if (isTRUE(USE_NUMDERIV_HESSIAN) && requireNamespace("numDeriv", quietly = TRUE)) {
      numDeriv::hessian(obj_theta_unpen, theta_hat)
    } else {
      optimHess(theta_hat, obj_theta_unpen)
    }
  }, error = function(e) NULL)
  if (is.null(H) || any(!is.finite(H))) {
    stop("profile theta-space Hessian failed for model ", fit_obj$model_id,
         ", gamma = ", fit_obj$gamma)
  }
  H <- (H + t(H)) / 2
  H
}

#' Calculate the effective degrees of freedom for pAIC and pBIC.
#'
#' The spatial contribution is the equivalent trace form in Equation (43),
#' trace{I_c (I_c + gamma I)^(-1)}.  Equation (24) adds the full contribution
#' of the unpenalized regression and covariance parameters.
calc_deff <- function(I11, beta_sigma_count, gamma) {
  # Spatial effective df: Equation (43), equivalent to Equation (42).
  # For a symmetric profile Hessian, the trace equals sum(ev / (ev + gamma)).
  # At gamma = 0, Equation (44) gives the number of free spatial parameters.
  m_spatial <- nrow(I11)
  if (m_spatial == 0L) return(beta_sigma_count)
  H <- (I11 + t(I11)) / 2
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values

  # Numerical Hessians can have very small negative eigenvalues.  Treat them as
  # zero; large negative eigenvalues are floored as well, with a warning.
  if (any(ev < -1e-6, na.rm = TRUE)) {
    warning(sprintf("Theta Hessian has negative eigenvalues; min = %.3e. Floored at 0 for deff.", min(ev, na.rm = TRUE)))
  }
  ev <- pmax(ev, 0)

  trace_part <- if (gamma == 0) {
    m_spatial
  } else {
    sum(ev / (ev + gamma))
  }
  beta_sigma_count + trace_part
}

#' Attach Hessian-based effective degrees of freedom to one fitted model.
#'
#' For penalized spatial models, this computes the Equation (38) profile
#' Hessian and the total effective df in Equation (24).  For models to which
#' this implementation does not apply the gamma penalty, deff is the nominal
#' parameter count and no shrinkage is used.
attach_profile_deff <- function(fit_i, common_inputs, gamma) {
  m <- fit_i$free_spatial_param_count

  # For diagonal / non-full models, do not apply the gamma effective-df
  # shrinkage.  Their pAIC/pBIC are ordinary AIC/BIC, so d_eff is the
  # nominal parameter count including any diagonal spatial parameters.
  if (!uses_gamma_penalty_model(fit_i$model_id)) {
    if (m == 0L) {
      fit_i$I11 <- matrix(0, 0, 0)
      fit_i$spatial_vcov <- matrix(0, 0, 0)
    }
    return(list(fit = fit_i,
                deff = fit_i$nominal_param_count,
                hessian_ok = TRUE,
                hessian_error = "gamma penalty not used; pAIC/pBIC equal AIC/BIC"))
  }

  if (m == 0L) {
    fit_i$I11 <- matrix(0, 0, 0)
    fit_i$spatial_vcov <- matrix(0, 0, 0)
    deff_i <- calc_deff(fit_i$I11, fit_i$beta_sigma_count, gamma)
    return(list(fit = fit_i, deff = deff_i, hessian_ok = TRUE, hessian_error = ""))
  }

  out <- tryCatch({
    H <- theta_hessian_profile(fit_i, common_inputs)
    fit_i$I11 <- H

    H_for_vcov <- H
    ev <- eigen((H_for_vcov + t(H_for_vcov)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (!all(is.finite(ev)) || min(ev) <= HESSIAN_RIDGE) {
      H_for_vcov <- H_for_vcov + diag(max(HESSIAN_RIDGE, HESSIAN_RIDGE - min(ev, na.rm = TRUE)), nrow(H_for_vcov))
    }
    fit_i$spatial_vcov <- safe_solve(H_for_vcov)
    fit_i$spatial_vcov <- (fit_i$spatial_vcov + t(fit_i$spatial_vcov)) / 2

    deff_i <- calc_deff(H, fit_i$beta_sigma_count, gamma)
    list(fit = fit_i, deff = deff_i, hessian_ok = TRUE, hessian_error = "")
  }, error = function(e) {
    list(fit = fit_i, deff = NA_real_, hessian_ok = FALSE,
         hessian_error = conditionMessage(e))
  })
  out
}


#' Assemble estimates, p-values, and significance marks for workbook output.
#'
#' The output uses a fixed parameter order, filling restricted parameters with
#' zero estimates and blank significance marks for easy cross-model comparison.
assemble_parameter_table <- function(fit_obj) {
  spec <- fit_obj$spec
  fit <- fit_obj$fit
  n_eff <- fit$n_eff

  vals <- setNames(rep(0, length(full_parameter_order()) - 1L), full_parameter_order()[-1])
  pval <- setNames(rep(NA_real_, length(full_parameter_order()) - 1L), full_parameter_order()[-1])
  sig  <- setNames(rep("", length(full_parameter_order()) - 1L), full_parameter_order()[-1])

  if (length(fit_obj$spatial_names) > 0) {
    for (i in seq_along(fit_obj$spatial_names)) {
      nm <- fit_obj$spatial_names[i]
      vals[nm] <- fit_obj$spatial_theta[i]
      se <- sqrt(max(fit_obj$spatial_vcov[i, i], 0))
      if (is.finite(se) && se > 0) {
        z <- vals[nm] / se
        pval[nm] <- 2 * pnorm(abs(z), lower.tail = FALSE)
        sig[nm] <- stars_from_p(pval[nm])
      }
    }
  }

  beta_names <- beta_names_by_model(spec$A_code)
  beta_hat <- fit$beta
  beta_se <- sqrt(pmax(diag(fit$Psi_beta), 0))
  for (i in seq_along(beta_names)) {
    nm <- beta_names[i]
    vals[nm] <- beta_hat[i]
    if (is.finite(beta_se[i]) && beta_se[i] > 0) {
      z <- beta_hat[i] / beta_se[i]
      pval[nm] <- 2 * pnorm(abs(z), lower.tail = FALSE)
      sig[nm] <- stars_from_p(pval[nm])
    }
  }

  Sigma <- fit$Sigma
  sigma_vals <- c(`sigma[1,1]` = Sigma[1, 1], `sigma[1,2]` = Sigma[1, 2],
                  `sigma[2,1]` = Sigma[2, 1], `sigma[2,2]` = Sigma[2, 2])
  sigma_se <- sigma_se_table(Sigma, n_eff, spec$S_code)
  for (nm in names(sigma_vals)) {
    vals[nm] <- sigma_vals[nm]
    se <- sigma_se[[nm]]
    if (is.finite(se) && se > 0) {
      z <- sigma_vals[nm] / se
      pval[nm] <- 2 * pnorm(abs(z), lower.tail = FALSE)
      sig[nm] <- stars_from_p(pval[nm])
    }
  }

  list(values = vals, pvalues = pval, sig = sig)
}

#' Compute likelihood-based model-selection metrics for a fitted model.
#'
#' The function evaluates pAIC and pBIC from Equation (25) and the response-wise
#' and averaged pseudo R-squared values from Equation (26).  It also returns
#' ordinary AIC/BIC and the parameter-count components used in the run summary.
selected_fit_metrics <- function(fit_obj, gamma, deff) {
  fit <- fit_obj$fit
  ll <- fit$logLik
  # Equation (25) uses log(Kn); here K = 2 and fit$n_eff = n * TT.
  n_eff_total <- 2 * fit$n_eff
  nominal_k <- fit_obj$nominal_param_count
  free_spatial_param_count <- fit_obj$free_spatial_param_count
  beta_sigma_count <- fit_obj$beta_sigma_count
  pseudo_R2_y1 <- cor(fit$actual_y1, fit$fitted_y1)^2
  pseudo_R2_y2 <- cor(fit$actual_y2, fit$fitted_y2)^2
  # Equation (26), specialized to K = 2.
  pseudo_R2_avg <- mean(c(pseudo_R2_y1, pseudo_R2_y2))

  list(
    gamma = gamma,
    logLik = ll,
    deff = deff,
    free_spatial_param_count = free_spatial_param_count,
    beta_sigma_count = beta_sigma_count,
    nominal_param_count = nominal_k,
    pAIC = -2 * ll + 2 * deff,
    pBIC = -2 * ll + log(n_eff_total) * deff,
    AIC = -2 * ll + 2 * nominal_k,
    BIC = -2 * ll + log(n_eff_total) * nominal_k,
    pseudo_R2_y1 = pseudo_R2_y1,
    pseudo_R2_y2 = pseudo_R2_y2,
    pseudo_R2_avg = pseudo_R2_avg
  )
}

# ------------------------------------------------------------
# Gamma-path fitting and exact-selection screening
# ------------------------------------------------------------
# Full spatial models are fitted over the gamma grid.  Warm starts reduce
# runtime, and lower/upper bounds on pAIC/pBIC avoid unnecessary Hessian work
# when DEFF_EVALUATION_MODE = "candidate".

#' Fit one model ID over its gamma path using warm starts.
#'
#' Warm starts propagate the previous gamma estimate to the next grid point.
#' Exact deff calculations can be evaluated for every gamma or only for
#' candidates that can still be optimal under rigorous lower-bound screening.
fit_over_gamma_grid <- function(model_id, common_inputs) {
  gamma_values <- gamma_grid_for_model(model_id, common_inputs$gamma_values)
  fits <- vector("list", length(gamma_values))
  names(fits) <- as.character(gamma_values)

  warm_theta <- NULL
  for (i in seq_along(gamma_values)) {
    g <- gamma_values[i]

    retry_used <- FALSE
    fit_i <- tryCatch(
      fit_mstr_model(model_id, common_inputs,
                     gamma = g,
                     warm_theta = warm_theta,
                     control_optim = list(maxit = MAXIT_MAIN, reltol = REL_TOL),
                     aggressive_starts = FALSE),
      error = function(e) NULL
    )

    if (is.null(fit_i)) {
      retry_used <- TRUE
      fit_i <- fit_mstr_model(model_id, common_inputs,
                              gamma = g,
                              warm_theta = warm_theta,
                              control_optim = list(maxit = MAXIT_RETRY, reltol = REL_TOL),
                              aggressive_starts = TRUE)
    }
    if (!is.null(fit_i$optim)) fit_i$optim$retry_used <- retry_used

    # Metrics will be completed after profile-Hessian deff is computed.
    met_i <- selected_fit_metrics(fit_i, g, NA_real_)
    fits[[i]] <- list(gamma = g, fit = fit_i, metrics = met_i,
                      hessian_ok = NA, hessian_error = "",
                      deff_exact_computed = FALSE)

    if (length(fit_i$spatial_theta) > 0) warm_theta <- fit_i$spatial_theta
  }

  nfit <- length(fits)
  ll_vec <- vapply(fits, function(z) z$metrics$logLik, numeric(1))
  beta_sigma_vec <- vapply(fits, function(z) z$metrics$beta_sigma_count, numeric(1))
  m_spatial_vec <- vapply(fits, function(z) z$metrics$free_spatial_param_count, numeric(1))
  n_eff_total_vec <- vapply(fits, function(z) 2 * z$fit$fit$n_eff, numeric(1))
  bic_pen_vec <- log(n_eff_total_vec)

  use_gamma_penalty <- uses_gamma_penalty_model(model_id)

  # For gamma-penalized full spatial models, beta_sigma_count <= deff <=
  # beta_sigma_count + m_spatial gives rigorous screening bounds.  For
  # diagonal / non-full models, pAIC/pBIC are ordinary AIC/BIC, so the lower
  # and upper bounds are both the nominal-count AIC/BIC.
  if (isTRUE(use_gamma_penalty)) {
    pAIC_lower <- -2 * ll_vec + 2 * beta_sigma_vec
    pAIC_upper <- -2 * ll_vec + 2 * (beta_sigma_vec + m_spatial_vec)
    pBIC_lower <- -2 * ll_vec + bic_pen_vec * beta_sigma_vec
    pBIC_upper <- -2 * ll_vec + bic_pen_vec * (beta_sigma_vec + m_spatial_vec)
  } else {
    nominal_vec <- beta_sigma_vec + m_spatial_vec
    pAIC_lower <- pAIC_upper <- -2 * ll_vec + 2 * nominal_vec
    pBIC_lower <- pBIC_upper <- -2 * ll_vec + bic_pen_vec * nominal_vec
  }

  computed <- rep(FALSE, nfit)

  compute_exact_for_index <- function(i) {
    if (isTRUE(computed[i])) return(invisible(NULL))
    g <- fits[[i]]$gamma
    ans <- attach_profile_deff(fits[[i]]$fit, common_inputs, g)
    fits[[i]]$fit <<- ans$fit
    fits[[i]]$metrics <<- selected_fit_metrics(ans$fit, g, ans$deff)
    fits[[i]]$hessian_ok <<- ans$hessian_ok
    fits[[i]]$hessian_error <<- ans$hessian_error
    fits[[i]]$deff_exact_computed <<- isTRUE(ans$hessian_ok) || ans$fit$free_spatial_param_count == 0L
    computed[i] <<- TRUE
    invisible(NULL)
  }

  mode <- match.arg(DEFF_EVALUATION_MODE, choices = c("candidate", "all"))

  if (mode == "all") {
    for (i in seq_len(nfit)) {
      message(sprintf("  [%s] profile deff/Hessian for gamma %d/%d: %g",
                      model_id, i, nfit, fits[[i]]$gamma))
      compute_exact_for_index(i)
    }
  } else {
    # No-spatial models are exact without Hessian.
    initial_idx <- which(m_spatial_vec == 0L)

    # Start with the best lower-bound candidates for both criteria.
    top_n <- max(1L, min(as.integer(DEFF_SCREEN_INITIAL_TOP), nfit))
    initial_idx <- unique(c(
      initial_idx,
      head(order(pAIC_lower), top_n),
      head(order(pBIC_lower), top_n)
    ))

    for (i in initial_idx) {
      message(sprintf("  [%s] profile deff/Hessian initial candidate: gamma=%g",
                      model_id, fits[[i]]$gamma))
      compute_exact_for_index(i)
    }

    repeat {
      exact_pAIC <- vapply(fits, function(z) z$metrics$pAIC, numeric(1))
      exact_pBIC <- vapply(fits, function(z) z$metrics$pBIC, numeric(1))
      best_aic <- suppressWarnings(min(exact_pAIC, na.rm = TRUE))
      best_bic <- suppressWarnings(min(exact_pBIC, na.rm = TRUE))
      if (!is.finite(best_aic)) best_aic <- Inf
      if (!is.finite(best_bic)) best_bic <- Inf

      need <- which(!computed & (
        pAIC_lower <= best_aic + DEFF_SCREEN_TOL |
          pBIC_lower <= best_bic + DEFF_SCREEN_TOL
      ))

      if (length(need) == 0L) break

      # Evaluate only a small batch, then re-screen.
      rank_score <- pmin(pAIC_lower[need] - best_aic,
                         pBIC_lower[need] - best_bic)
      need <- need[order(rank_score, pAIC_lower[need], pBIC_lower[need])]
      batch_n <- max(1L, as.integer(DEFF_SCREEN_BATCH_SIZE))
      batch <- head(need, batch_n)

      for (i in batch) {
        message(sprintf("  [%s] profile deff/Hessian screened candidate: gamma=%g",
                        model_id, fits[[i]]$gamma))
        compute_exact_for_index(i)
      }
    }
  }

  gamma_tbl <- data.frame(
    gamma = vapply(fits, function(z) z$gamma, numeric(1)),
    logLik = vapply(fits, function(z) z$metrics$logLik, numeric(1)),
    deff = vapply(fits, function(z) z$metrics$deff, numeric(1)),
    deff_exact_computed = vapply(fits, function(z) isTRUE(z$deff_exact_computed), logical(1)),
    hessian_ok = vapply(fits, function(z) if (is.na(z$hessian_ok)) FALSE else isTRUE(z$hessian_ok), logical(1)),
    hessian_error = vapply(fits, function(z) if (is.null(z$hessian_error)) "" else as.character(z$hessian_error), character(1)),
    free_spatial_param_count = vapply(fits, function(z) z$metrics$free_spatial_param_count, numeric(1)),
    beta_sigma_count = vapply(fits, function(z) z$metrics$beta_sigma_count, numeric(1)),
    nominal_param_count = vapply(fits, function(z) z$metrics$nominal_param_count, numeric(1)),
    pAIC = vapply(fits, function(z) z$metrics$pAIC, numeric(1)),
    pAIC_lower_bound = pAIC_lower,
    pAIC_upper_bound = pAIC_upper,
    pBIC = vapply(fits, function(z) z$metrics$pBIC, numeric(1)),
    pBIC_lower_bound = pBIC_lower,
    pBIC_upper_bound = pBIC_upper,
    AIC = vapply(fits, function(z) z$metrics$AIC, numeric(1)),
    BIC = vapply(fits, function(z) z$metrics$BIC, numeric(1)),
    optim_convergence = vapply(fits, function(z) if (is.null(z$fit$optim$convergence)) NA_integer_ else as.integer(z$fit$optim$convergence), integer(1)),
    best_optim_converged = vapply(fits, function(z) if (is.null(z$fit$optim$converged)) NA else isTRUE(z$fit$optim$converged), logical(1)),
    n_starts = vapply(fits, function(z) if (is.null(z$fit$optim$n_starts)) NA_integer_ else as.integer(z$fit$optim$n_starts), integer(1)),
    n_successful_starts = vapply(fits, function(z) if (is.null(z$fit$optim$n_successful_starts)) NA_integer_ else as.integer(z$fit$optim$n_successful_starts), integer(1)),
    n_converged_starts = vapply(fits, function(z) if (is.null(z$fit$optim$n_converged_starts)) NA_integer_ else as.integer(z$fit$optim$n_converged_starts), integer(1)),
    retry_used = vapply(fits, function(z) if (is.null(z$fit$optim$retry_used)) FALSE else isTRUE(z$fit$optim$retry_used), logical(1)),
    stringsAsFactors = FALSE
  )

  list(model_id = model_id,
       template_path = common_inputs$template_path,
       fits = fits,
       gamma_tbl = gamma_tbl,
       deff_evaluation_mode = mode)
}

# -------------------------------
# Workbook writing
# -------------------------------

#' Create a result workbook, optionally starting from a template.
#'
#' If no template is supplied, a clean workbook is created and populated with
#' standardized result sheets.
create_output_workbook <- function(template_path = NA_character_) {
  if (!is.null(template_path) && length(template_path) == 1L &&
      !is.na(template_path) && file.exists(template_path)) {
    return(openxlsx::loadWorkbook(template_path))
  }
  openxlsx::createWorkbook()
}

#' Replace a workbook sheet with a data frame.
#'
#' This keeps repeated runs deterministic by removing any existing sheet before
#' writing the new data.
write_replace_sheet <- function(wb, sheet, x) {
  if (sheet %in% openxlsx::sheets(wb)) openxlsx::removeWorksheet(wb, sheet)
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeData(wb, sheet = sheet, x = x,
                      startRow = 1, startCol = 1,
                      colNames = TRUE, rowNames = FALSE)
  invisible(wb)
}

#' Populate and save the workbook for the pAIC- or pBIC-selected model.
#'
#' The workbook contains parameter estimates, model-fit metrics, and the full
#' gamma grid diagnostics for the selected model ID.
fill_workbook <- function(template_path, gamma_fit_obj,
                          selected_by = c("pAIC", "pBIC"), out_file) {
  selected_by <- match.arg(selected_by)
  template_to_use <- if (isTRUE(USE_TEMPLATE_FOR_XLSX)) template_path else NA_character_
  wb <- create_output_workbook(template_to_use)

  tbl <- gamma_fit_obj$gamma_tbl
  pick_idx <- if (selected_by == "pAIC") which_min_finite(tbl$pAIC) else which_min_finite(tbl$pBIC)
  if (!is.finite(pick_idx)) stop("No finite ", selected_by, " value is available for workbook output.")
  best_fit <- gamma_fit_obj$fits[[pick_idx]]$fit
  best_metrics <- gamma_fit_obj$fits[[pick_idx]]$metrics

  param_tab <- assemble_parameter_table(best_fit)
  param_order <- full_parameter_order()
  param_no_gamma <- param_order[-1]

  est_value_vec <- c(gamma = best_metrics$gamma,
                     unname(as.numeric(param_tab$values[param_no_gamma])))
  pvalue_vec <- c(NA_real_, unname(as.numeric(param_tab$pvalues[param_no_gamma])))
  sig_vec <- c("", unname(as.character(param_tab$sig[param_no_gamma])))

  est_df <- data.frame(
    Parameter = param_order,
    `Estimated value` = est_value_vec,
    `p-value` = pvalue_vec,
    `significance level` = sig_vec,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_replace_sheet(wb, "est_parameters", est_df)

  opt <- best_fit$optim
  model_fit_df <- data.frame(
    Metric = c("selected_by", "gamma", "logLik", "deff",
               "free_spatial_param_count", "beta_sigma_count", "nominal_param_count",
               "pAIC", "pBIC",
               "pseudo_R2_y1", "pseudo_R2_y2", "pseudo_R2_avg", "AIC", "BIC",
               "optim_convergence", "optim_converged", "n_starts",
               "n_successful_starts", "n_converged_starts", "retry_used"),
    Value = c(selected_by,
              best_metrics$gamma,
              best_metrics$logLik,
              best_metrics$deff,
              best_metrics$free_spatial_param_count,
              best_metrics$beta_sigma_count,
              best_metrics$nominal_param_count,
              best_metrics$pAIC,
              best_metrics$pBIC,
              best_metrics$pseudo_R2_y1,
              best_metrics$pseudo_R2_y2,
              best_metrics$pseudo_R2_avg,
              best_metrics$AIC,
              best_metrics$BIC,
              if (is.null(opt$convergence)) NA else opt$convergence,
              if (is.null(opt$converged)) NA else opt$converged,
              if (is.null(opt$n_starts)) NA else opt$n_starts,
              if (is.null(opt$n_successful_starts)) NA else opt$n_successful_starts,
              if (is.null(opt$n_converged_starts)) NA else opt$n_converged_starts,
              if (is.null(opt$retry_used)) FALSE else opt$retry_used),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_replace_sheet(wb, "model_fit", model_fit_df)

  write_replace_sheet(wb, "gamma_grid", tbl)
  safe_save_workbook(wb, out_file, overwrite = TRUE)
}

# ------------------------------------------------------------
# True-value evaluation for the 1111 data-generating process
# ------------------------------------------------------------
# When the full MGNS model is fitted, the selected estimates can be compared
# with known simulation truth to report differences, RMSE, and convergence.

TRUE_PARAMS_1111 <- list(
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

#' Return the true parameter values used for the 1111 data-generating model.
#'
#' These values support simulation diagnostics by comparing selected estimates
#' with the known data-generating parameters.
true_parameter_values_1111 <- function(true_params = TRUE_PARAMS_1111) {
  c(
    `rho[1,1]` = true_params$R[1, 1],
    `rho[1,2]` = true_params$R[1, 2],
    `rho[2,1]` = true_params$R[2, 1],
    `rho[2,2]` = true_params$R[2, 2],
    `lambda[1,1]` = true_params$Lambda[1, 1],
    `lambda[1,2]` = true_params$Lambda[1, 2],
    `lambda[2,1]` = true_params$Lambda[2, 1],
    `lambda[2,2]` = true_params$Lambda[2, 2],
    `alpha[1,1]` = true_params$A[1, 1],
    `alpha[1,2]` = true_params$A[1, 2],
    `alpha[2,1]` = true_params$A[2, 1],
    `alpha[2,2]` = true_params$A[2, 2],
    `sigma[1,1]` = true_params$Sigma[1, 1],
    `sigma[1,2]` = true_params$Sigma[1, 2],
    `sigma[2,1]` = true_params$Sigma[2, 1],
    `sigma[2,2]` = true_params$Sigma[2, 2],
    beta_intercept_y1 = true_params$beta_intercept_y1,
    beta_intercept_y2 = true_params$beta_intercept_y2,
    beta_common1_y1 = true_params$beta_common1_y1,
    beta_common1_y2 = true_params$beta_common1_y2,
    beta_common2_y1 = true_params$beta_common2_y1,
    beta_common2_y2 = true_params$beta_common2_y2,
    beta_specific1_1 = true_params$beta_specific1_1,
    beta_specific2_1 = true_params$beta_specific2_1
  )
}

#' Extract estimates from the pAIC- or pBIC-selected gamma row.
#'
#' The output is aligned to the unified parameter order so that true-value
#' errors and RMSE summaries can be computed directly.
extract_selected_estimates <- function(gamma_fit_obj, selected_by = c("pAIC", "pBIC")) {
  selected_by <- match.arg(selected_by)
  tbl <- gamma_fit_obj$gamma_tbl
  pick_idx <- if (selected_by == "pAIC") which_min_finite(tbl$pAIC) else which_min_finite(tbl$pBIC)
  best_fit <- gamma_fit_obj$fits[[pick_idx]]$fit
  param_tab <- assemble_parameter_table(best_fit)
  param_order <- full_parameter_order()[-1]
  setNames(as.numeric(param_tab$values[param_order]), param_order)
}

#' Create the true-value evaluation table for the full MGNS model 1111.
#'
#' The table reports estimates selected by pAIC and pBIC, differences from the
#' known data-generating values, RMSEs, and convergence diagnostics.
build_true_evaluation_1111 <- function(gamma_fit_obj, out_dir = OUT_DIR,
                                       true_params = TRUE_PARAMS_1111) {
  out_dir <- ensure_writable_dir(out_dir)
  true_vals <- true_parameter_values_1111(true_params)
  bias_tables <- list()
  rmse_rows <- list()

  for (crit in c("pAIC", "pBIC")) {
    tbl <- gamma_fit_obj$gamma_tbl
    pick_idx <- if (crit == "pAIC") which_min_finite(tbl$pAIC) else which_min_finite(tbl$pBIC)
    est_vals <- extract_selected_estimates(gamma_fit_obj, crit)
    params <- intersect(names(true_vals), names(est_vals))

    df <- data.frame(
      selected_by = crit,
      parameter = params,
      true_value = as.numeric(true_vals[params]),
      estimated_value = as.numeric(est_vals[params]),
      difference = as.numeric(est_vals[params] - true_vals[params]),
      abs_difference = abs(as.numeric(est_vals[params] - true_vals[params])),
      squared_error = as.numeric((est_vals[params] - true_vals[params])^2),
      stringsAsFactors = FALSE
    )
    bias_tables[[crit]] <- df

    # Count convergence in two ways: across selected best optimizers for each gamma,
    # and across all multi-start attempts used along the gamma path.
    best_conv_count <- sum(tbl$best_optim_converged %in% TRUE, na.rm = TRUE)
    rmse_rows[[crit]] <- data.frame(
      selected_by = crit,
      selected_gamma = tbl$gamma[pick_idx],
      selected_gamma_index = pick_idx,
      gamma_grid_points = nrow(tbl),
      RMSE = sqrt(mean(df$squared_error, na.rm = TRUE)),
      n_parameters_compared = nrow(df),
      best_optimizer_convergence_zero_count = best_conv_count,
      selected_optim_converged = tbl$best_optim_converged[pick_idx],
      selected_optim_convergence_code = tbl$optim_convergence[pick_idx],
      total_multistart_attempts = sum(tbl$n_starts, na.rm = TRUE),
      total_successful_starts = sum(tbl$n_successful_starts, na.rm = TRUE),
      total_converged_starts = sum(tbl$n_converged_starts, na.rm = TRUE),
      retry_used_count = sum(tbl$retry_used %in% TRUE, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  bias_tbl <- do.call(rbind, bias_tables)
  rmse_tbl <- do.call(rbind, rmse_rows)

  bias_csv <- file.path(out_dir, "true_value_bias_1111.csv")
  rmse_csv <- file.path(out_dir, "true_value_rmse_convergence_1111.csv")
  safe_write_csv(bias_tbl, bias_csv, row.names = FALSE)
  safe_write_csv(rmse_tbl, rmse_csv, row.names = FALSE)

  xlsx_file <- NA_character_
  if (isTRUE(WRITE_TRUE_VALUE_XLSX)) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "parameter_bias")
    openxlsx::writeData(wb, "parameter_bias", bias_tbl)
    openxlsx::addWorksheet(wb, "rmse_convergence")
    openxlsx::writeData(wb, "rmse_convergence", rmse_tbl)
    xlsx_file <- file.path(out_dir, "true_value_evaluation_1111.xlsx")
    safe_save_workbook(wb, xlsx_file, overwrite = TRUE)
    message("Saved true-value evaluation workbook: ", normalizePath(xlsx_file, mustWork = FALSE))
  }

  invisible(list(parameter_bias = bias_tbl, rmse_convergence = rmse_tbl,
                 files = c(bias_csv = bias_csv, rmse_csv = rmse_csv, xlsx = xlsx_file)))
}

# -------------------------------
# Runners
# -------------------------------

#' Run the complete gamma-path estimation for one model ID.
#'
#' Errors are caught and returned in a structured object so that failure of one
#' model does not prevent the remaining specifications from being estimated.
run_one_model <- function(model_id, common_inputs, out_dir = OUT_DIR) {
  # Compute only.  Workbook/CSV writing is done after all models finish, so
  # parallel workers do not compete for slow Excel/template I/O.
  gamma_fit <- fit_over_gamma_grid(model_id, common_inputs)
  invisible(gamma_fit)
}

#' Write gamma paths and selected-parameter workbooks for all successful models.
#'
#' CSV gamma paths are written first; Excel workbooks are then produced for the
#' pAIC- and pBIC-selected fits when workbook output is enabled.
write_model_outputs <- function(results_all, common_inputs, out_dir = OUT_DIR) {
  out_dir <- ensure_writable_dir(out_dir)

  for (id in names(results_all)) {
    res <- results_all[[id]]
    if (!isTRUE(res$ok) || is.null(res$fit)) next
    gamma_fit <- res$fit

    if (isTRUE(WRITE_GAMMA_PATH_CSV)) {
      gamma_csv <- file.path(out_dir, sprintf("gamma_path_%s.csv", id))
      safe_write_csv(gamma_fit$gamma_tbl, gamma_csv, row.names = FALSE)
      message("Saved gamma path: ", normalizePath(gamma_csv, mustWork = FALSE))
    }

    if (isTRUE(WRITE_XLSX)) {
      out_aic <- file.path(out_dir, sprintf("estimated_parameters_%s_pAIC.xlsx", id))
      out_bic <- file.path(out_dir, sprintf("estimated_parameters_%s_pBIC.xlsx", id))
      fill_workbook(common_inputs$template_path, gamma_fit, selected_by = "pAIC", out_file = out_aic)
      fill_workbook(common_inputs$template_path, gamma_fit, selected_by = "pBIC", out_file = out_bic)
      message("Saved: ", normalizePath(out_aic, mustWork = FALSE))
      message("Saved: ", normalizePath(out_bic, mustWork = FALSE))
    }
  }
  invisible(TRUE)
}


#' Estimate all requested model IDs, optionally in parallel.
#'
#' The function exports shared data and helper functions to PSOCK workers, then
#' optionally retries failed models sequentially with the same structured output.
run_all_models <- function(common_inputs,
                           model_ids = MODEL_IDS,
                           out_dir = OUT_DIR,
                           parallel = USE_PARALLEL,
                           n_cores = N_CORES,
                           retry_failed_sequentially = RETRY_FAILED_SEQUENTIALLY) {
  model_ids <- unique(as.character(model_ids))
  if (length(model_ids) == 0L) stop("model_ids is empty")

  out_dir <- ensure_writable_dir(out_dir)

  if (is.null(n_cores)) {
    n_cores <- choose_auto_cores(N_VALUE, length(model_ids), NULL, reserve_cores = RESERVE_CORES)
  }
  n_cores <- max(1L, min(as.integer(n_cores), length(model_ids)))

  worker_fun <- function(id) {
    message("--- Running model ID = ", id, " ---")
    tryCatch(
      {
        fit <- run_one_model(id, common_inputs, out_dir = out_dir)
        list(ok = TRUE, model_id = id, fit = fit, error = NULL)
      },
      error = function(e) {
        msg <- conditionMessage(e)
        message("*** Model ", id, " failed: ", msg)
        list(ok = FALSE, model_id = id, fit = NULL, error = msg)
      }
    )
  }

  use_parallel <- isTRUE(parallel) && length(model_ids) >= 2L && n_cores >= 2L

  if (!use_parallel) {
    results <- lapply(model_ids, worker_fun)
    names(results) <- model_ids
    failures <- vapply(results, function(x) !isTRUE(x$ok), logical(1))
    if (any(failures)) warning("Some models failed: ", paste(names(results)[failures], collapse = ", "))
    return(invisible(results))
  }

  os_type <- .Platform$OS.type

  if (os_type != "windows") {
    results <- parallel::mclapply(model_ids, worker_fun, mc.cores = n_cores)
    names(results) <- model_ids
    failures <- vapply(results, function(x) !isTRUE(x$ok), logical(1))
    if (any(failures) && isTRUE(retry_failed_sequentially)) {
      failed_ids <- names(results)[failures]
      message("Retrying failed models sequentially: ", paste(failed_ids, collapse = ", "))
      retry_res <- lapply(failed_ids, worker_fun)
      names(retry_res) <- failed_ids
      for (id in failed_ids) if (isTRUE(retry_res[[id]]$ok)) results[[id]] <- retry_res[[id]]
      failures <- vapply(results, function(x) !isTRUE(x$ok), logical(1))
    }
    if (any(failures)) warning("Some models failed: ", paste(names(results)[failures], collapse = ", "))
    return(invisible(results))
  }

  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages({
      library(Matrix)
      library(openxlsx)
      library(parallel)
    })
    NULL
  })
  parallel::clusterExport(
    cl,
    varlist = c(
      "OUT_DIR", "MAXIT_MAIN", "MAXIT_RETRY", "REL_TOL", "DEFF_EVALUATION_MODE", "DEFF_SCREEN_INITIAL_TOP", "DEFF_SCREEN_BATCH_SIZE", "DEFF_SCREEN_TOL", "HESSIAN_MAX_GLS_ITER", "HESSIAN_GLS_TOL", "HESSIAN_THETA_BOUND", "HESSIAN_RIDGE", "USE_NUMDERIV_HESSIAN",
      "stars_from_p", "safe_solve", "safe_solve_general", "choose_auto_cores",
      "make_mask", "decode_model_id", "beta_names_by_model", "full_parameter_order",
      "uses_gamma_penalty_model", "gamma_grid_for_model",
      "pack_spatial_theta", "unpack_spatial_u", "unpack_spatial_theta", "spatial_param_names", "build_multistart_theta",
      "logdet_univ", "logdet_kron2", "standardize_W_rows", "safe_numeric_matrix_from_csv",
      "read_inputs", "read_gamma_values", "build_period_data",
      "fit_univariate_gns", "build_common_inputs", "profile_beta_sigma", "fit_mstr_model",
      "sigma_se_table", "which_min_finite", "theta_hessian_profile", "calc_deff", "attach_profile_deff", "assemble_parameter_table", "selected_fit_metrics",
      "fit_over_gamma_grid", "create_output_workbook", "write_replace_sheet",
      "fill_workbook", "run_one_model", "worker_fun",
      "GAMMA_VALUES", "GAMMA_CSV_FILE", "GAMMA_SOURCE", "GAMMA_PENALTY_MODEL_IDS", "W_ROW_STANDARDIZE", "W_ROW_TOL",
      "common_inputs", "out_dir"
    ),
    envir = environment()
  )
  results <- parallel::parLapply(cl, model_ids, worker_fun)
  names(results) <- model_ids
  failures <- vapply(results, function(x) !isTRUE(x$ok), logical(1))
  if (any(failures) && isTRUE(retry_failed_sequentially)) {
    failed_ids <- names(results)[failures]
    message("Retrying failed models sequentially: ", paste(failed_ids, collapse = ", "))
    retry_res <- lapply(failed_ids, worker_fun)
    names(retry_res) <- failed_ids
    for (id in failed_ids) if (isTRUE(retry_res[[id]]$ok)) results[[id]] <- retry_res[[id]]
    failures <- vapply(results, function(x) !isTRUE(x$ok), logical(1))
  }
  if (any(failures)) warning("Some models failed: ", paste(names(results)[failures], collapse = ", "))
  invisible(results)
}

# ------------------------------------------------------------
# Main execution block
# ------------------------------------------------------------
# This block validates input files, initializes shared objects, runs all model
# specifications, writes gamma-path and selected-model outputs, and then creates
# the final summary table comparing the pAIC and pBIC selections.

if (!file.exists(DATA_FILE)) {
  alt_data <- file.path(DATA_DIR, sprintf("simulated_data_%s_n%d_T%d.csv", DATA_ID, N_VALUE, T_VALUE))
  if (file.exists(alt_data)) DATA_FILE <- alt_data
}
if (!file.exists(WEIGHTS_FILE)) {
  alt_w <- file.path(DATA_DIR, sprintf("spatial_weights_n%d.csv", N_VALUE))
  if (file.exists(alt_w)) WEIGHTS_FILE <- alt_w
}
if (!file.exists(DATA_FILE)) stop("Data file not found: ", normalizePath(DATA_FILE, mustWork = FALSE))
if (!file.exists(WEIGHTS_FILE)) stop("Weights file not found: ", normalizePath(WEIGHTS_FILE, mustWork = FALSE))
OUT_DIR <- ensure_writable_dir(OUT_DIR, label = "OUT_DIR")

message("========================================")
message("mstr_mgns_full11_commented_complete.r")
message("N_VALUE      : ", N_VALUE)
message("DATA_FILE    : ", normalizePath(DATA_FILE, mustWork = TRUE))
message("WEIGHTS_FILE : ", normalizePath(WEIGHTS_FILE, mustWork = TRUE))
message("TEMPLATE_FILE: ", if (file.exists(TEMPLATE_FILE)) normalizePath(TEMPLATE_FILE, mustWork = TRUE) else "not used / not found")
message("GAMMA_CSV    : ", if (file.exists(GAMMA_CSV_FILE)) normalizePath(GAMMA_CSV_FILE, mustWork = TRUE) else "not used / not found")
message("OUT_DIR      : ", normalizePath(OUT_DIR, mustWork = FALSE))
message("TARGET_TIME  : ", TARGET_TIME)
message("DEFF_MODE    : ", DEFF_EVALUATION_MODE)
message("WRITE_XLSX   : ", WRITE_XLSX, "; USE_TEMPLATE_FOR_XLSX: ", USE_TEMPLATE_FOR_XLSX)
message("========================================")

common_inputs <- build_common_inputs(
  data_file = DATA_FILE,
  weights_file = WEIGHTS_FILE,
  template_file = TEMPLATE_FILE,
  target_time = TARGET_TIME
)

message("Gamma grid source: ", common_inputs$gamma_source,
        "; points: ", length(common_inputs$gamma_values),
        "; range: ", min(common_inputs$gamma_values), " - ", max(common_inputs$gamma_values))
message("W row-standardization: ", common_inputs$W_stats$row_standardized,
        "; max row-sum deviation before = ", common_inputs$W_stats$max_abs_deviation_before,
        "; after = ", common_inputs$W_stats$max_abs_deviation_after)
safe_write_csv(common_inputs$W_stats,
          file.path(OUT_DIR, "weight_matrix_diagnostics.csv"), row.names = FALSE)
safe_write_csv(data.frame(gamma = common_inputs$gamma_values),
          file.path(OUT_DIR, "gamma_grid_used.csv"), row.names = FALSE)

results_all <- run_all_models(
  common_inputs = common_inputs,
  model_ids = MODEL_IDS,
  out_dir = OUT_DIR,
  parallel = USE_PARALLEL,
  n_cores = N_CORES,
  retry_failed_sequentially = RETRY_FAILED_SEQUENTIALLY
)

#' Extract a selected scalar value from one model result for the final summary.
#'
#' The helper is used to build run_summary_n<N>.csv for both pAIC and pBIC
#' selections without duplicating indexing logic.
summary_selected_value <- function(res, crit = c("pAIC", "pBIC"), field) {
  crit <- match.arg(crit)
  if (!isTRUE(res$ok) || is.null(res$fit) || is.null(res$fit$gamma_tbl)) return(NA_real_)
  tbl <- res$fit$gamma_tbl
  idx <- if (crit == "pAIC") which_min_finite(tbl$pAIC) else which_min_finite(tbl$pBIC)
  if (!(field %in% names(tbl))) return(NA_real_)
  as.numeric(tbl[[field]][idx])
}

run_summary <- data.frame(
  model_id = names(results_all),
  ok = vapply(results_all, function(x) isTRUE(x$ok), logical(1)),
  error = vapply(results_all, function(x) if (is.null(x$error)) "" else as.character(x$error), character(1)),
  pAIC = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "pAIC"),
  pAIC_gamma = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "gamma"),
  pAIC_deff = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "deff"),
  pAIC_logLik = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "logLik"),
  pAIC_free_spatial_parameter_count = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "free_spatial_param_count"),
  pAIC_beta_sigma_count = vapply(results_all, summary_selected_value, numeric(1), crit = "pAIC", field = "beta_sigma_count"),
  pBIC = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "pBIC"),
  pBIC_gamma = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "gamma"),
  pBIC_deff = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "deff"),
  pBIC_logLik = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "logLik"),
  pBIC_free_spatial_parameter_count = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "free_spatial_param_count"),
  pBIC_beta_sigma_count = vapply(results_all, summary_selected_value, numeric(1), crit = "pBIC", field = "beta_sigma_count"),
  stringsAsFactors = FALSE
)

summary_file <- file.path(OUT_DIR, sprintf("run_summary_n%d.csv", N_VALUE))
safe_write_csv(run_summary, summary_file, row.names = FALSE)
message("Saved summary: ", normalizePath(summary_file, mustWork = FALSE))
print(run_summary, row.names = FALSE)

# Workbook output is postponed until the CSV summary has already been saved and printed.
write_model_outputs(results_all, common_inputs, out_dir = OUT_DIR)

if ("1111" %in% names(results_all) && isTRUE(results_all[["1111"]]$ok)) {
  build_true_evaluation_1111(results_all[["1111"]]$fit, out_dir = OUT_DIR)
} else {
  warning("Model 1111 did not finish successfully; true-value evaluation was not created.")
}

# ------------------------------------------------------------
# Save whole-program elapsed time
# ------------------------------------------------------------
# This block is intentionally placed at the very end of the script, after all
# estimation and result-writing routines.  It measures the total elapsed time
# from PROGRAM_START_TIME to successful completion of the program and saves the
# result as a CSV file in OUT_DIR.
PROGRAM_END_TIME <- Sys.time()
PROGRAM_END_ELAPSED <- unname(proc.time()[["elapsed"]])
PROGRAM_ELAPSED_SECONDS <- PROGRAM_END_ELAPSED - PROGRAM_START_ELAPSED

elapsed_time_report <- data.frame(
  script_name = "implement-B-with-elapsed-time.r",
  n_value = N_VALUE,
  target_time = TARGET_TIME,
  start_time = format(PROGRAM_START_TIME, "%Y-%m-%d %H:%M:%S %Z"),
  end_time = format(PROGRAM_END_TIME, "%Y-%m-%d %H:%M:%S %Z"),
  elapsed_seconds = round(PROGRAM_ELAPSED_SECONDS, 3),
  elapsed_minutes = round(PROGRAM_ELAPSED_SECONDS / 60, 3),
  elapsed_hours = round(PROGRAM_ELAPSED_SECONDS / 3600, 3),
  completion_status = "completed",
  stringsAsFactors = FALSE
)

elapsed_time_file <- file.path(
  OUT_DIR,
  sprintf("elapsed_time_n%d.csv", N_VALUE)
)
safe_write_csv(elapsed_time_report, elapsed_time_file, row.names = FALSE)

message("========================================")
message("Program completed successfully.")
message(sprintf("Elapsed time: %.3f seconds (%.3f minutes / %.3f hours)",
                PROGRAM_ELAPSED_SECONDS,
                PROGRAM_ELAPSED_SECONDS / 60,
                PROGRAM_ELAPSED_SECONDS / 3600))
message("Saved elapsed-time report: ",
        normalizePath(elapsed_time_file, mustWork = FALSE))
message("========================================")

