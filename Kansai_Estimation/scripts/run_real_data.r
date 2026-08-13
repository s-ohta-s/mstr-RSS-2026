# run_real_data.r
#
# Real-data estimation runner script (Kansai-region panel data, n=198 municipalities × T=2 time points)
#
# Estimate 11 models on the Kansai-region panel data and
# output γ search, significance tests, comparison tables, and fitted values/residuals.
#
# Data (data/real/):
#   transformed_data198_for_R.csv    — transformed panel data
#   W_row_standardized_198_for_R.csv — row-standardized spatial weight matrix
#   region_id_mapping198.csv         — municipality-ID mapping table (for reference)
#
# Model structure:
#   y1 equation: std taxable income       = f(ln aging rate, ln population density, ln establishments per capita)
#   y2 equation: std population change ratio = f(ln aging rate, ln population density, ln medical facilities per capita)
#
# Variable mapping:
#   y1             = std taxable income        (standardized ln(taxable income))
#   y2             = std population change ratio (standardized ln(1 + population change rate/100))
#   x_common1      = ln aging rate
#   x_common2      = ln population density
#   x_specific1_1  = ln establishments per capita
#   x_specific2_1  = ln medical facilities per capita
#   time=1 -> year 2016, time=2 -> year 2021
#   region=1..198 -> see region_id_mapping198.csv
#
# Usage (run via source() from the RStudio console):
#   source("scripts/run_real_data.r")
#
#   - Paths are resolved relative to the script, so the working directory does not matter.

cat("\n")
cat("Real-data estimation: n = 198 municipalities × T = 2 time points\n")
total_start_time <- Sys.time()

# Phase 0: configuration

cat("\n### Phase 0: configuration ###\n")

# Output verbosity (TRUE shows per-function detailed logs and timing; default is quiet)
VERBOSE <- FALSE

# Resolve paths relative to the script (independent of the working directory; do not use setwd())
# This script lives under scripts/, so its parent is taken as the project root.
.get_script_dir <- function() {
  # When run with Rscript: get the script path from the --file= argument
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  # When run with source(): walk up the call stack to find the ofile set by source()
  # (in RStudio etc., sys.frames()[[1]] is not necessarily the source frame, so search all)
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) {
      return(dirname(normalizePath(of)))
    }
  }
  getwd()
}
SCRIPT_DIR <- .get_script_dir()
ROOT <- dirname(SCRIPT_DIR)   # project root (parent of scripts/)

N <- 198
DATA_FILE   <- file.path(ROOT, "data", "real", "transformed_data198_for_R.csv")
WEIGHT_FILE <- file.path(ROOT, "data", "real", "W_row_standardized_198_for_R.csv")
OUTPUT_DIR  <- file.path(ROOT, "output")

# γ-search grid (two-stage search: coarse -> fine)
# Stage 1: log10 step 0.5, upper bound 10^6 (to detect boundary solutions)
# Stage 2: refine the neighborhood of the optimal γ* in log10 steps of 0.05
GAMMAS_COARSE <- c(0, 10^seq(-2, 6, by = 0.5))

Y_VARS <- c("y1", "y2")
X_VARS <- list(
  y1 = c("x_common1", "x_common2", "x_specific1_1"),
  y2 = c("x_common1", "x_common2", "x_specific2_1")
)
TIME_VAR    <- "time"
TIME_POINT  <- 2        # estimate at year 2021 (time=2), use year 2016 (time=1) as the time lag
REGION_VAR  <- "region"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("  Loading source files...\n")
# Packages imported directly (Matrix is a transitive dependency of spdep/spatialreg, so not listed)
required_packages <- c("spdep", "spatialreg", "numDeriv", "systemfit")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(sprintf("the package %s is required: install.packages('%s')", pkg, pkg))
  }
}

# Load the function library (in dependency order)
core_files <- c(
  "multivar_common.r", "msar_likelihood.r", "spatial_utilities.r",
  "initial_estimation_extended.r", "msem_mgns_likelihood.r", "model_methods.r",
  "build_result_object.r", "penalized_spatial.r", "fit_indsar.r", "fit_indsem.r",
  "fit_indgns.r", "fit_varx.r", "fit_indreg.r", "full_hessian_inference.r",
  "inference_and_export.r", "experiment_output_functions.r"
)
for (f in core_files) source(file.path(ROOT, "R", f))

stopifnot(file.exists(DATA_FILE))
stopifnot(file.exists(WEIGHT_FILE))

cat(sprintf("  Data: %s\n", DATA_FILE))
cat(sprintf("  Weight matrix: %s\n", WEIGHT_FILE))
cat(sprintf("  Output: %s/\n", OUTPUT_DIR))
cat(sprintf("  γ grid (coarse): %d points (%.4g ~ %.4g)\n",
            length(GAMMAS_COARSE), min(GAMMAS_COARSE), max(GAMMAS_COARSE)))

# Phase 1: γ search (6 full models)

cat("\n")
cat("### Phase 1: γ search ###\n")

phase1_start <- Sys.time()

# γ-search configuration: model_id -> (fit_func, include_time_lag)
gamma_search_configs <- list(
  "1111" = list(func = fit_mgns_penalized, time_lag = TRUE,  label = "MGNS(1111)"),
  "0111" = list(func = fit_msem_penalized,  time_lag = TRUE,  label = "MSEM(0111)"),
  "1011" = list(func = fit_msar_penalized,  time_lag = TRUE,  label = "MSAR(1011)"),
  "1101" = list(func = fit_mgns_penalized, time_lag = FALSE, label = "MGNS(1101)"),
  "0101" = list(func = fit_msem_penalized,  time_lag = FALSE, label = "MSEM(0101)"),
  "1001" = list(func = fit_msar_penalized,  time_lag = FALSE, label = "MSAR(1001)")
)

gamma_searches <- list()

for (mid in names(gamma_search_configs)) {
  cfg <- gamma_search_configs[[mid]]

  cat(sprintf("\n--- γ search for %s ---\n", cfg$label))

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
    cat(sprintf("  x Error: %s\n", e$message))
    NULL
  })

  if (!is.null(gs)) {
    gamma_searches[[mid]] <- gs

    # Output 1: γ-search CSV
    export_gamma_search_csv(
      gs, mid,
      output_file = file.path(OUTPUT_DIR, sprintf("gamma_search_%s.csv", mid)),
      verbose = VERBOSE
    )
  }
}

# Output 1 summary: optimal γ across models
gamma_summary <- build_gamma_optimal_summary(
  gamma_searches,
  output_file = file.path(OUTPUT_DIR, "gamma_optimal_summary.csv"),
  verbose = VERBOSE
)

phase1_time <- difftime(Sys.time(), phase1_start, units = "mins")
if (VERBOSE) cat(sprintf("\nPhase 1 complete: %.1f min\n", as.numeric(phase1_time)))

# Phase 2: estimate all 11 models

cat("\n")
cat("### Phase 2: estimate all 11 models ###\n")

phase2_start <- Sys.time()

# 2a: full models -- take the optimal-γ result from the γ-search results

results_pAIC <- list()   # pAIC criterion
results_pBIC <- list()   # pBIC criterion
gamma_info_pAIC <- list() # γ info (pAIC criterion)
gamma_info_pBIC <- list() # γ info (pBIC criterion)

for (mid in names(gamma_searches)) {
  gs <- gamma_searches[[mid]]
  s <- gs$summary

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

# 2b: diagonal models
cat("\n--- Estimating the diagonal models ---\n")

common_diag_args <- list(
  data_file = DATA_FILE, weight_file = WEIGHT_FILE,
  y_vars = Y_VARS, x_vars = X_VARS,
  time_var = TIME_VAR, time_point = TIME_POINT,
  region_var = REGION_VAR, include_time_lag = TRUE,
  verbose = VERBOSE
)

diag_models <- list(
  "d0dd" = list(func = fit_indsar,  label = "IndSAR"),
  "0ddd" = list(func = fit_indsem,  label = "IndSEM"),
  "dddd" = list(func = fit_indgns, label = "IndGNS")
)

for (mid in names(diag_models)) {
  dm <- diag_models[[mid]]
  cat(sprintf("\n--- %s (%s) ---\n", mid, dm$label))

  result <- tryCatch({
    do.call(dm$func, common_diag_args)
  }, error = function(e) {
    cat(sprintf("  x Error: %s\n", e$message))
    NULL
  })

  if (!is.null(result)) {
    results_pAIC[[mid]] <- result
    results_pBIC[[mid]] <- result
    gamma_info_pAIC[[mid]] <- list(gamma = 0, pAIC = result$fit$AIC,
                                    pBIC = result$fit$BIC, d_eff = result$fit$num_params)
    gamma_info_pBIC[[mid]] <- gamma_info_pAIC[[mid]]
    cat(sprintf("  %s complete (AIC=%.2f)\n", mid, result$fit$AIC))
  }
}

# 2c: VARX
cat("\n--- VARX (0011) ---\n")
result_varx <- tryCatch({
  fit_varx(
    data_file = DATA_FILE,
    y_vars = Y_VARS, x_vars = X_VARS,
    time_var = TIME_VAR, time_point = TIME_POINT,
    region_var = REGION_VAR, verbose = VERBOSE
  )
}, error = function(e) {
  cat(sprintf("  x Error: %s\n", e$message))
  NULL
})

if (!is.null(result_varx)) {
  results_pAIC[["0011"]] <- result_varx
  results_pBIC[["0011"]] <- result_varx
  gamma_info_pAIC[["0011"]] <- list(gamma = 0, pAIC = result_varx$fit$AIC,
                                     pBIC = result_varx$fit$BIC, d_eff = result_varx$fit$num_params)
  gamma_info_pBIC[["0011"]] <- gamma_info_pAIC[["0011"]]
  cat(sprintf("  0011 complete (AIC=%.2f)\n", result_varx$fit$AIC))
}

# 2d: OLS diagonal
cat("\n--- OLS diagonal (000d) ---\n")
result_ols <- tryCatch({
  fit_indreg(
    data_file = DATA_FILE,
    y_vars = Y_VARS, x_vars = X_VARS,
    time_var = TIME_VAR, time_point = TIME_POINT,
    region_var = REGION_VAR, verbose = VERBOSE
  )
}, error = function(e) {
  cat(sprintf("  x Error: %s\n", e$message))
  NULL
})

if (!is.null(result_ols)) {
  results_pAIC[["000d"]] <- result_ols
  results_pBIC[["000d"]] <- result_ols
  gamma_info_pAIC[["000d"]] <- list(gamma = 0, pAIC = result_ols$fit$AIC,
                                     pBIC = result_ols$fit$BIC, d_eff = result_ols$fit$num_params)
  gamma_info_pBIC[["000d"]] <- gamma_info_pAIC[["000d"]]
  cat(sprintf("  000d complete (AIC=%.2f)\n", result_ols$fit$AIC))
}

phase2_time <- difftime(Sys.time(), phase2_start, units = "mins")
if (VERBOSE) cat(sprintf("\nPhase 2 complete: %.1f min\n", as.numeric(phase2_time)))
cat(sprintf("  Number of estimated models: pAIC=%d, pBIC=%d\n",
            length(results_pAIC), length(results_pBIC)))

# Phase 3: significance tests

cat("\n")
cat("### Phase 3: significance tests ###\n")

phase3_start <- Sys.time()

full_model_ids <- c("1111", "0111", "1011", "1101", "0101", "1001")

# pAIC-criterion version
cat("\n--- Significance tests (pAIC-criterion version) ---\n")
for (mid in full_model_ids) {
  if (is.null(results_pAIC[[mid]])) next
  g <- gamma_info_pAIC[[mid]]$gamma
  cat(sprintf("\n  %s (γ=%.4g): ", mid, g))

  # Ψ-based SE + profile Hessian SE
  results_pAIC[[mid]] <- tryCatch({
    add_inference(results_pAIC[[mid]], compute_spatial_se = TRUE, gamma = g, verbose = VERBOSE)
  }, error = function(e) {
    cat(sprintf("add_inference error: %s ", e$message))
    results_pAIC[[mid]]
  })

  # Full Hessian sandwich SE
  results_pAIC[[mid]] <- tryCatch({
    add_full_inference(results_pAIC[[mid]], gamma = g, verbose = VERBOSE)
  }, error = function(e) {
    cat(sprintf("add_full_inference error: %s ", e$message))
    results_pAIC[[mid]]
  })

  cat("done\n")
}

# pBIC-criterion version
cat("\n--- Significance tests (pBIC-criterion version) ---\n")
for (mid in full_model_ids) {
  if (is.null(results_pBIC[[mid]])) next
  g <- gamma_info_pBIC[[mid]]$gamma

  # Skip if γ equals the pAIC one (same object)
  if (!is.null(gamma_info_pAIC[[mid]]) &&
      g == gamma_info_pAIC[[mid]]$gamma) {
    results_pBIC[[mid]] <- results_pAIC[[mid]]
    cat(sprintf("  %s: γ*_pBIC == γ*_pAIC -> reuse the pAIC version\n", mid))
    next
  }

  cat(sprintf("\n  %s (γ=%.4g): ", mid, g))

  results_pBIC[[mid]] <- tryCatch({
    add_inference(results_pBIC[[mid]], compute_spatial_se = TRUE, gamma = g, verbose = VERBOSE)
  }, error = function(e) {
    cat(sprintf("add_inference error: %s ", e$message))
    results_pBIC[[mid]]
  })

  results_pBIC[[mid]] <- tryCatch({
    add_full_inference(results_pBIC[[mid]], gamma = g, verbose = VERBOSE)
  }, error = function(e) {
    cat(sprintf("add_full_inference error: %s ", e$message))
    results_pBIC[[mid]]
  })

  cat("done\n")
}

# add_inference for the diagonal models and VARX
cat("\n--- Significance tests for non-full models ---\n")
non_full_ids <- c("d0dd", "0ddd", "dddd", "0011", "000d")
for (mid in non_full_ids) {
  if (is.null(results_pAIC[[mid]])) next
  cat(sprintf("  %s: ", mid))
  results_pAIC[[mid]] <- tryCatch({
    add_inference(results_pAIC[[mid]], compute_spatial_se = TRUE, gamma = 0, verbose = VERBOSE)
  }, error = function(e) {
    cat(sprintf("Error: %s ", e$message))
    results_pAIC[[mid]]
  })
  results_pBIC[[mid]] <- results_pAIC[[mid]]
  cat("done\n")
}

phase3_time <- difftime(Sys.time(), phase3_start, units = "mins")
if (VERBOSE) cat(sprintf("\nPhase 3 complete: %.1f min\n", as.numeric(phase3_time)))

# Phase 4: generate output tables

cat("\n")
cat("### Phase 4: generate output tables ###\n")

phase4_start <- Sys.time()

cat("\n--- Output 2: parameter comparison table (pAIC criterion) ---\n")
build_comparison_table(
  results_list = results_pAIC,
  gamma_info = gamma_info_pAIC,
  output_file = file.path(OUTPUT_DIR, "comparison_table_pAIC.csv"),
  verbose = VERBOSE
)

cat("\n--- Output 3: parameter comparison table (pBIC criterion) ---\n")
build_comparison_table(
  results_list = results_pBIC,
  gamma_info = gamma_info_pBIC,
  output_file = file.path(OUTPUT_DIR, "comparison_table_pBIC.csv"),
  verbose = VERBOSE
)

cat("\n--- Output 4: model-selection summary table ---\n")
build_model_selection_table(
  results_pAIC = results_pAIC,
  results_pBIC = results_pBIC,
  gamma_info_pAIC = gamma_info_pAIC,
  gamma_info_pBIC = gamma_info_pBIC,
  output_file = file.path(OUTPUT_DIR, "model_selection_summary.csv"),
  verbose = VERBOSE
)

# Since this is real data, the true parameters are unknown; no bias evaluation is performed.
cat("\n--- Output 5: comparison table against true values -> skipped (no true values for real data) ---\n")

cat("\n--- Output 6: Ψ vs Hessian SE comparison table (pAIC criterion) ---\n")
full_results_pAIC <- results_pAIC[intersect(full_model_ids, names(results_pAIC))]
build_se_comparison_table(
  results_list = full_results_pAIC,
  gamma_info = gamma_info_pAIC,
  criterion = "pAIC",
  output_file = file.path(OUTPUT_DIR, "se_comparison_full_models_pAIC.csv"),
  verbose = VERBOSE
)

cat("\n--- Output 6: Ψ vs Hessian SE comparison table (pBIC criterion) ---\n")
full_results_pBIC <- results_pBIC[intersect(full_model_ids, names(results_pBIC))]
build_se_comparison_table(
  results_list = full_results_pBIC,
  gamma_info = gamma_info_pBIC,
  criterion = "pBIC",
  output_file = file.path(OUTPUT_DIR, "se_comparison_full_models_pBIC.csv"),
  verbose = VERBOSE
)

# pred = trend prediction ŷ (d=0): MSAR/MGNS ŷ=(I−R̂⊗W)⁻¹Xβ̂, MSEM/VARX/OLS ŷ=Xβ̂
# resid = innovation residual ε̂, resid_reduced = obs − pred
cat("\n--- Output 7: fitted values / residuals (pAIC criterion, all models) ---\n")
export_fitted_residuals_all_models(
  results_list = results_pAIC,
  output_dir = file.path(OUTPUT_DIR, "fitted_residuals_pAIC"),
  verbose = VERBOSE
)

cat("\n--- Output 7: fitted values / residuals (pBIC criterion, all models) ---\n")
export_fitted_residuals_all_models(
  results_list = results_pBIC,
  output_dir = file.path(OUTPUT_DIR, "fitted_residuals_pBIC"),
  verbose = VERBOSE
)

phase4_time <- difftime(Sys.time(), phase4_start, units = "mins")
if (VERBOSE) cat(sprintf("\nPhase 4 complete: %.1f min\n", as.numeric(phase4_time)))

# Phase 5: output check and completion report

cat("\n")
cat("### Phase 5: output check and completion report ###\n")

cat("\nList of output files:\n")
files <- list.files(OUTPUT_DIR, full.names = TRUE)
for (f in files) {
  finfo <- file.info(f)
  cat(sprintf("  %s (%.1f KB)\n", basename(f), finfo$size / 1024))
}

total_time <- difftime(Sys.time(), total_start_time, units = "mins")
if (VERBOSE) cat(sprintf("\nExecution-time summary:\n"))
if (VERBOSE) cat(sprintf("  Phase 1 (γ search):    %.1f min\n", as.numeric(phase1_time)))
if (VERBOSE) cat(sprintf("  Phase 2 (estimation):  %.1f min\n", as.numeric(phase2_time)))
if (VERBOSE) cat(sprintf("  Phase 3 (significance): %.1f min\n", as.numeric(phase3_time)))
if (VERBOSE) cat(sprintf("  Phase 4 (tables):      %.1f min\n", as.numeric(phase4_time)))
if (VERBOSE) cat(sprintf("  Total:                 %.1f min\n", as.numeric(total_time)))

cat("\n")
cat(sprintf("Real-data estimation n=%d complete\n", N))

# Variable-name reference notes (consult when interpreting results)
#
# Parameter names in the CSV output -> actual variable meanings:
#
#   beta_intercept_y1     -> intercept of std taxable income
#                           (original scale: β × σ₁ = β × 0.1242 gives the change in ln(taxable income))
#   beta_intercept_y2     -> intercept of std population change ratio
#                           (original scale: β × σ₂' = β × 0.01081 gives the change in ln(1+x/100))
#   beta_x_common1_y1     -> std taxable income <- ln aging rate
#   beta_x_common1_y2     -> std population change ratio <- ln aging rate
#   beta_x_common2_y1     -> std taxable income <- ln population density
#   beta_x_common2_y2     -> std population change ratio <- ln population density
#   beta_x_specific1_1_y1 -> std taxable income <- ln establishments per capita
#   beta_x_specific2_1_y2 -> std population change ratio <- ln medical facilities per capita
#   R[i,j]                -> spatial lag matrix (spatial spillover to y)
#   Lambda[i,j]           -> spatial error matrix (spatial correlation of the errors)
#   A[i,j]                -> time-lag matrix (AR(1) coefficients)
#   Sigma[i,j]            -> SUR error covariance matrix
#
# ---- Standardization parameters ----
#
# y1: μ₁  = 7.999338875981922,  σ₁  = 0.124198460868784
# y2: μ₂' = -0.00974936746851443, σ₂' = 0.01081453957123827
#
# ---- Recovery to the original scale ----
#
# Recovering y1:
#   ln(taxable income) = std_y1 × 0.1242 + 7.9993
#   taxable income     = exp(std_y1 × 0.1242 + 7.9993)
#
# Recovering y2:
#   ln(1 + population change rate/100) = std_y2 × 0.01081 + (-0.00975)
#   population change rate (%)         = 100 × (exp(std_y2 × 0.01081 - 0.00975) − 1)
#
# ---- Recovering β coefficients to the original scale ----
#
# y1 equation:
#   Δ(ln taxable income)   = β₁ × σ₁  = β₁ × 0.1242
#
# y2 equation:
#   Δ(continuous growth rate) = β₂ × σ₂' = β₂ × 0.01081
#   Δ(population change rate %) ≈ 100 × β₂ × σ₂' = β₂ × 1.0815
#   (approximation valid for small changes, since ln(1+x/100) ≈ x/100)
#
# ---- Recovering per-capita variables ----
#
#   establishment density (per area) = exp(ln establishments per capita + ln population density)
#   medical facility density (per area) = exp(ln medical facilities per capita + ln population density) − 0.001
