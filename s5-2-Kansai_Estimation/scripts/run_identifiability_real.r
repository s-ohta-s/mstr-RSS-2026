# run_identifiability_real.r
#
# Numerical check of the two sufficient identifiability conditions
# (Propositions 1 and 2, §3.2 of the paper) on the Kansai real data (n = 198, T = 2).
#
# Real data have no true parameters, so Proposition 1 is evaluated at the
# FITTED mean vector mu_hat instead of mu0:
#
#     mu_hat = (I - R_hat x W)^{-1} X beta_hat      (MSAR / MGNS)
#     mu_hat = X beta_hat                           (MSEM / VARX / OLS, R = 0)
#
# No model is re-estimated. beta_hat and R_hat are read back from
#   output/comparison_table_pAIC.csv  and  output/comparison_table_pBIC.csv,
# and the reconstructed mu_hat is cross-checked against the y*_pred columns of
#   output/fitted_residuals_<criterion>/fitted_residuals_<model>.csv
# which hold the same quantity. The script aborts if they disagree, since that
# would mean the coefficient ordering was misread.
#
# Proposition 2 depends on W only and is therefore model-independent.
#
# Usage (run via source() from the RStudio console, or with Rscript):
#   source("scripts/run_identifiability_real.r")
#
# Output:
#   output/identifiability_diag_real.csv

cat("\n")
cat("Identifiability diagnostics on the real data (n = 198)\n")

# Phase 0: configuration

# Resolve paths relative to the script (independent of the working directory)
.get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
SCRIPT_DIR <- .get_script_dir()
ROOT <- dirname(SCRIPT_DIR)

DATA_FILE   <- file.path(ROOT, "data", "real", "transformed_data198_for_R.csv")
WEIGHT_FILE <- file.path(ROOT, "data", "real", "W_row_standardized_198_for_R.csv")
RESULT_DIR  <- file.path(ROOT, "output")
OUT_FILE    <- file.path(RESULT_DIR, "identifiability_diag_real.csv")

# Estimation settings -- must match run_real_data.r
Y_VARS <- c("y1", "y2")
X_VARS <- list(
  y1 = c("x_common1", "x_common2", "x_specific1_1"),
  y2 = c("x_common1", "x_common2", "x_specific2_1")
)
TIME_VAR   <- "time"
TIME_POINT <- 2
REGION_VAR <- "region"

# Tolerance for the mu_hat cross-check against the exported fitted values
MU_CHECK_TOL <- 1e-8

stopifnot(file.exists(DATA_FILE), file.exists(WEIGHT_FILE), dir.exists(RESULT_DIR))

if (!requireNamespace("spdep", quietly = TRUE)) {
  stop("the package spdep is required: install.packages('spdep')")
}

# Source only what is needed (no estimation code is loaded)
for (f in c("multivar_common.r", "spatial_utilities.r", "identifiability_diag.r")) {
  source(file.path(ROOT, "R", f))
}

# Parameter names, in the same order as the columns of the design matrix X
# (build_design_matrix_extended: y1 block, y2 block, then the y_lag block with
#  A in row-major order)
BETA_BASE_NAMES <- c(
  "beta_intercept_y1", "beta_x_common1_y1", "beta_x_common2_y1", "beta_x_specific1_1_y1",
  "beta_intercept_y2", "beta_x_common1_y2", "beta_x_common2_y2", "beta_x_specific2_1_y2"
)
A_NAMES <- c("A[1,1]", "A[1,2]", "A[2,1]", "A[2,2]")
R_NAMES <- c("R[1,1]", "R[1,2]", "R[2,1]", "R[2,2]")

# Phase 1: data preparation (X and W, one variant per time-lag setting)

cat("\n### Phase 1: data preparation ###\n")

.data_cache <- list()
get_data_list <- function(include_time_lag) {
  key <- as.character(include_time_lag)
  if (is.null(.data_cache[[key]])) {
    .data_cache[[key]] <<- prepare_data_extended(
      data_file = DATA_FILE, weight_file = WEIGHT_FILE,
      y_vars = Y_VARS, x_vars = X_VARS,
      time_var = TIME_VAR, time_point = TIME_POINT, region_var = REGION_VAR,
      include_intercept = TRUE, include_time_lag = include_time_lag,
      verbose = FALSE
    )
    dl <- .data_cache[[key]]
    cat(sprintf("  include_time_lag = %-5s : X is %d x %d, n = %d, k = %d\n",
                include_time_lag, nrow(dl$X), ncol(dl$X), dl$n, dl$k))
  }
  .data_cache[[key]]
}

dl0 <- get_data_list(TRUE)
N <- dl0$n
K <- dl0$k
W <- dl0$W

# Phase 2: Proposition 2 (depends on W only)

cat("\n### Phase 2: Proposition 2 ###\n")

p2 <- prop2_gram_diag(W)
cat(sprintf("  min eigenvalue of the normalized Gram matrix : %.6g\n", p2$min_eig))
cat(sprintf("  corresponding min singular value            : %.6g\n", p2$min_sv))
cat(sprintf("  W is %s\n", ifelse(p2$W_symmetric, "SYMMETRIC (!)", "nonsymmetric")))
if (p2$min_eig <= 0) {
  cat("  ! WARNING: Proposition 2 is NOT satisfied\n")
}

# Phase 3: Proposition 1 at mu_hat, for every model and both criteria

cat("\n### Phase 3: Proposition 1 at mu_hat ###\n")

# Pull one parameter block out of a comparison table as a numeric vector
get_block <- function(tab, mid, pnames) {
  col <- paste0("estimate.", mid)
  if (!col %in% colnames(tab)) return(rep(NA_real_, length(pnames)))
  v <- suppressWarnings(as.numeric(tab[[col]]))
  v[match(pnames, tab$parameter)]
}

rows <- list()
mu_check_failed <- character(0)

for (crit in c("pAIC", "pBIC")) {

  tab_file <- file.path(RESULT_DIR, sprintf("comparison_table_%s.csv", crit))
  if (!file.exists(tab_file)) {
    cat(sprintf("\n  %s: %s not found -- skipped\n", crit, basename(tab_file)))
    next
  }
  tab <- read.csv(tab_file, stringsAsFactors = FALSE, check.names = FALSE)

  model_ids <- sub("^estimate\\.", "",
                   grep("^estimate\\.", colnames(tab), value = TRUE))

  cat(sprintf("\n  --- %s (%d models) ---\n", crit, length(model_ids)))

  for (mid in model_ids) {

    beta_base <- get_block(tab, mid, BETA_BASE_NAMES)
    if (any(is.na(beta_base))) {
      cat(sprintf("    %-5s: regression coefficients incomplete -- skipped\n", mid))
      next
    }

    # A present <=> the model was fitted with include_time_lag = TRUE.
    # For the diagonal models the off-diagonal entries of A are structurally
    # absent and are treated as zero.
    A_vals <- get_block(tab, mid, A_NAMES)
    has_lag <- any(!is.na(A_vals))
    beta_hat <- if (has_lag) c(beta_base, ifelse(is.na(A_vals), 0, A_vals)) else beta_base

    # R absent (MSEM / VARX / OLS) -> zero matrix; off-diagonals of a diagonal
    # R are likewise treated as zero
    R_vals <- get_block(tab, mid, R_NAMES)
    R_hat <- if (all(is.na(R_vals))) matrix(0, K, K) else
      matrix(ifelse(is.na(R_vals), 0, R_vals), K, K, byrow = TRUE)

    dl <- get_data_list(has_lag)
    if (length(beta_hat) != ncol(dl$X)) {
      cat(sprintf("    %-5s: beta has %d entries but X has %d columns -- skipped\n",
                  mid, length(beta_hat), ncol(dl$X)))
      next
    }

    mu_hat <- compute_mu(R_hat, beta_hat, dl$X, dl$W, K, N)

    # Cross-check against the exported fitted values (same quantity)
    fr_file <- file.path(RESULT_DIR, paste0("fitted_residuals_", crit),
                         sprintf("fitted_residuals_%s.csv", mid))
    mu_diff <- NA_real_
    if (file.exists(fr_file)) {
      fr <- read.csv(fr_file, stringsAsFactors = FALSE)
      fr <- fr[order(fr[[REGION_VAR]]), ]
      pred <- c(fr$y1_pred, fr$y2_pred)
      if (length(pred) == length(mu_hat)) {
        mu_diff <- max(abs(mu_hat - pred))
        if (!is.na(mu_diff) && mu_diff > MU_CHECK_TOL) {
          mu_check_failed <- c(mu_check_failed, sprintf("%s/%s (%.3g)", crit, mid, mu_diff))
        }
      }
    }

    p1 <- prop1_rank_margin_from_mu(mu_hat, dl$X, dl$W, K, N)

    rows[[length(rows) + 1]] <- data.frame(
      scope = "model", criterion = crit, model = mid,
      include_time_lag = has_lag,
      n = N, k = K, p = ncol(dl$X), required_rank = p1$required,
      rank = p1$rank, rank_ok = p1$rank_ok,
      min_sv = p1$min_sv, max_sv = p1$max_sv, cond = p1$cond,
      mu_check_max_abs_diff = mu_diff,
      prop2_min_eig = NA_real_, prop2_min_sv = NA_real_, W_symmetric = NA,
      stringsAsFactors = FALSE
    )

    cat(sprintf("    %-5s: rank %2d/%2d %-3s  min_sv = %.6g  cond = %.4g  (mu check %.2g)\n",
                mid, p1$rank, p1$required, ifelse(p1$rank_ok, "OK", "NG"),
                p1$min_sv, p1$cond, mu_diff))
  }
}

# Phase 4: assemble and export

cat("\n### Phase 4: export ###\n")

if (length(mu_check_failed) > 0) {
  stop(sprintf(paste0("mu_hat does not match the exported fitted values for: %s\n",
                      "  The coefficient ordering read from the comparison table is wrong; ",
                      "no output was written."),
               paste(mu_check_failed, collapse = ", ")))
}

w_row <- data.frame(
  scope = "W", criterion = NA_character_, model = NA_character_,
  include_time_lag = NA, n = N, k = K, p = NA_integer_, required_rank = NA_integer_,
  rank = NA_integer_, rank_ok = NA, min_sv = NA_real_, max_sv = NA_real_, cond = NA_real_,
  mu_check_max_abs_diff = NA_real_,
  prop2_min_eig = p2$min_eig, prop2_min_sv = p2$min_sv, W_symmetric = p2$W_symmetric,
  stringsAsFactors = FALSE
)

out <- rbind(w_row, do.call(rbind, rows))
write.csv(out, OUT_FILE, row.names = FALSE)
cat(sprintf("  written: %s (%d rows)\n", basename(OUT_FILE), nrow(out)))

# Phase 5: summary for the manuscript

cat("\n### Summary ###\n")

mrows <- out[out$scope == "model", ]
cat(sprintf("  Proposition 1: full column rank in %d of %d model/criterion combinations\n",
            sum(mrows$rank_ok), nrow(mrows)))
cat(sprintf("  Proposition 1: smallest singular value over all combinations = %.6g (%s, %s)\n",
            min(mrows$min_sv),
            mrows$model[which.min(mrows$min_sv)],
            mrows$criterion[which.min(mrows$min_sv)]))
cat(sprintf("  Proposition 1: largest singular value over all combinations  = %.6g\n",
            max(mrows$min_sv)))
cat(sprintf("  Proposition 2: min eigenvalue = %.6g  (min singular value = %.6g)\n",
            p2$min_eig, p2$min_sv))
cat(sprintf("  mu_hat cross-check: max abs difference = %.3g (tolerance %.0e)\n",
            max(mrows$mu_check_max_abs_diff, na.rm = TRUE), MU_CHECK_TOL))

cat("\nDone\n")
