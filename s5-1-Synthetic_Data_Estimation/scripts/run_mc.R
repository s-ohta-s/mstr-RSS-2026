# run_mc.R — Monte Carlo driver
#
# Runs experiments A (estimator validity: true model only) and B (model
# selection: all 11 candidates) over cells (dgp × paramset × n) × replications,
# with replications parallelized via future.apply (multisession).
#
# Usage:
#   Rscript scripts/run_mc.R                        # quick wiring check
#   Rscript scripts/run_mc.R --reps=300 --full      # production run
#
# Options (all optional):
#   --reps=N        replications per cell        (default 5)
#   --full          all 3 dgp × 4 paramsets × 3 n (default: subset)
#   --dgps=...      comma list, e.g. --dgps=1111,0111
#   --paramsets=... e.g. --paramsets=weak,realistic
#   --ns=...        e.g. --ns=100,400
#   --workers=N     parallel workers               (default 10)
#   --noB / --noA   skip experiment B / A
#   --serial        run sequentially (debugging)
#   --outdir=PATH   output directory               (default output/sim)
#   --nocheckpoint  disable per-replication checkpointing
#
# Checkpointing (ON by default): every completed replication is saved to
# <outdir>/checkpoints/ as an .rds file. Re-running the SAME command after an
# interruption skips completed replications and resumes. Progress at any time:
#   (Get-ChildItem <outdir>\checkpoints).Count   # completed / total tasks
# After a fully successful run the checkpoint directory is deleted.
#
# Outputs (written under --outdir):
#   results_A.csv            tidy per-parameter rows
#   results_B.csv            one selection row per rep
#   identifiability_diag.csv prop2 per W + prop1 aggregates per cell
#   mc_errors.csv            failed replications (cell, rep, message)
#   mc_run_info.txt          configuration + timing record

## ---------------------------------------------------------------- options
.args <- commandArgs(trailingOnly = TRUE)
optval <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), .args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^--", name, "="), "", hit[1])
}
optflag <- function(name) any(.args == paste0("--", name))

N_REPS     <- as.integer(optval("reps", "5"))
N_WORKERS  <- as.integer(optval("workers", "10"))
RUN_A      <- !optflag("noA")
RUN_B      <- !optflag("noB")
SERIAL     <- optflag("serial")
CHECKPOINT <- !optflag("nocheckpoint")
BASE_SEED  <- 20260711L

DGPS      <- strsplit(optval("dgps",      if (optflag("full")) "1111,0111,1011" else "1111,0111"), ",")[[1]]
PARAMSETS <- strsplit(optval("paramsets", if (optflag("full")) "weak,mid,strong,realistic" else "weak,realistic"), ",")[[1]]
NS        <- as.integer(strsplit(optval("ns", if (optflag("full")) "100,400,900" else "100"), ",")[[1]])

## ---------------------------------------------------------------- paths
.get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  getwd()
}
ROOT <- dirname(.get_script_dir())
OUTDIR <- optval("outdir", file.path(ROOT, "output", "sim"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MC_SOURCE_FILES <- c(
  "multivar_common.r", "msar_likelihood.r", "spatial_utilities.r",
  "initial_estimation_extended.r", "msem_mgns_likelihood.r", "model_methods.r",
  "build_result_object.r", "penalized_spatial.r", "fit_indsar.r", "fit_indsem.r",
  "fit_indgns.r", "fit_varx.r", "fit_indreg.r", "full_hessian_inference.r",
  "inference_and_export.r", "experiment_output_functions.r",
  "data_prep_fast.r", "stationarity_diag.r", "identifiability_diag.r",
  "sim_parameters.r", "generate_dataset.r", "mc_experiment.r"
)
source_all <- function(root) {
  for (f in MC_SOURCE_FILES) source(file.path(root, "R", f))
}
suppressMessages(source_all(ROOT))

## ---------------------------------------------------------------- seeds
## strides must exceed max rep count (300), so
## seed = BASE + 1e7*dgp_idx + 1e6*paramset_idx + 1e5*n_idx + rep_id
DGP_IDX      <- c("1111" = 1L, "0111" = 2L, "1011" = 3L)
PARAMSET_IDX <- c(weak = 1L, mid = 2L, strong = 3L, realistic = 4L)
N_IDX        <- c("100" = 1L, "400" = 2L, "900" = 3L)
mc_seed <- function(dgp, ps, n, rep) {
  BASE_SEED + 1e7L * DGP_IDX[[dgp]] + 1e6L * PARAMSET_IDX[[ps]] +
    1e5L * N_IDX[[as.character(n)]] + rep
}

## ---------------------------------------------------------------- plan
cells <- expand.grid(dgp = DGPS, paramset = PARAMSETS, n = NS,
                     stringsAsFactors = FALSE)
cat(sprintf("MC driver: %d cells x %d reps | A=%s B=%s | workers=%d%s\n",
            nrow(cells), N_REPS, RUN_A, RUN_B, N_WORKERS,
            ifelse(SERIAL, " (serial)", "")))
cat(sprintf("Cells: dgp={%s} paramset={%s} n={%s}\n",
            paste(DGPS, collapse = ","), paste(PARAMSETS, collapse = ","),
            paste(NS, collapse = ",")))
t_start <- Sys.time()

## static quantities + prop2 diagnostics per n (computed once)
statics <- list()
prop2_rows <- list()
for (n in NS) {
  W <- get_or_create_W(n, dir = file.path(ROOT, "data", "simulated"))
  statics[[as.character(n)]] <- prepare_static(W)
  p2 <- prop2_gram_diag(W)
  prop2_rows[[length(prop2_rows) + 1]] <- data.frame(
    scope = "W", n = n, dgp = NA, paramset = NA,
    prop1_pass_rate = NA, prop1_min_sv_min = NA, prop1_min_sv_median = NA,
    prop2_min_eig = p2$min_eig, W_symmetric = p2$W_symmetric,
    stringsAsFactors = FALSE
  )
  cat(sprintf("  n=%d: prop2 min_eig = %.6f (W %s)\n", n, p2$min_eig,
              ifelse(p2$W_symmetric, "symmetric (!)", "nonsymmetric")))
}

## save true-parameter tables per (dgp, paramset) cell
for (ci in seq_len(nrow(cells))) {
  dgp <- cells$dgp[ci]; ps <- cells$paramset[ci]
  f <- file.path(OUTDIR, sprintf("true_parameters_%s_%s.csv", dgp, ps))
  if (!file.exists(f)) {
    write.csv(params_to_true_df(get_dgp_params(ps, dgp)), f, row.names = FALSE)
  }
}

## ---------------------------------------------------------------- run
## one task per (cell, rep); workers source the R files once and cache statics
tasks <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) {
  data.frame(ci = ci, rep = seq_len(N_REPS))
}))

CKPT_DIR <- file.path(OUTDIR, "checkpoints")
if (CHECKPOINT) {
  dir.create(CKPT_DIR, showWarnings = FALSE, recursive = TRUE)
  n_done <- length(list.files(CKPT_DIR, pattern = "\\.rds$"))
  cat(sprintf("checkpointing ON: %s (%d/%d tasks already done)\n",
              CKPT_DIR, n_done, nrow(tasks)))
}

run_task <- function(ti) {
  ci <- tasks$ci[ti]; rep <- tasks$rep[ti]
  dgp <- cells$dgp[ci]; ps <- cells$paramset[ci]; n <- cells$n[ci]

  ## resume: return the stored result if this replication already completed
  ckpt <- if (CHECKPOINT) {
    file.path(CKPT_DIR, sprintf("ckpt_%s_%s_n%d_rep%03d.rds", dgp, ps, n, rep))
  } else NULL
  if (!is.null(ckpt) && file.exists(ckpt)) {
    prev <- tryCatch(readRDS(ckpt), error = function(e) NULL)
    if (!is.null(prev)) return(prev)   # corrupt file -> recompute below
  }

  if (!exists(".mc_worker_ready", envir = globalenv())) {
    suppressMessages(source_all(ROOT))
    assign(".mc_worker_ready", TRUE, envir = globalenv())
    assign(".mc_worker_statics", list(), envir = globalenv())
    assign(".mc_worker_solvers", list(), envir = globalenv())
  }
  st <- get(".mc_worker_statics", envir = globalenv())
  key <- as.character(n)
  if (is.null(st[[key]])) {
    W <- get_or_create_W(n, dir = file.path(ROOT, "data", "simulated"))
    st[[key]] <- prepare_static(W)
    assign(".mc_worker_statics", st, envir = globalenv())
  }
  static <- st[[key]]

  ## per-cell solver cache (spatial systems fixed within a cell)
  sv <- get(".mc_worker_solvers", envir = globalenv())
  skey <- paste(dgp, ps, n, sep = "|")
  if (is.null(sv[[skey]])) {
    sv[[skey]] <- make_spatial_solvers(get_dgp_params(ps, dgp), static)
    assign(".mc_worker_solvers", sv, envir = globalenv())
  }

  res <- mc_run_rep(dgp, ps, as.character(n), static, rep,
                    seed = mc_seed(dgp, ps, n, rep),
                    run_A = RUN_A, run_B = RUN_B, solvers = sv[[skey]])
  res$cell <- list(dgp = dgp, paramset = ps, n = n, rep = rep)

  ## checkpoint: atomic write (tmp + rename) so an interrupted save never
  ## leaves a truncated .rds that would be mistaken for a completed task
  if (!is.null(ckpt) && is.null(res$error)) {
    tmp <- paste0(ckpt, ".tmp")
    ok <- tryCatch({ saveRDS(res, tmp); file.rename(tmp, ckpt) },
                   error = function(e) FALSE)
    if (!isTRUE(ok)) unlink(tmp)
  }
  res
}

if (SERIAL) {
  all_res <- lapply(seq_len(nrow(tasks)), run_task)
} else {
  suppressMessages({
    library(future)
    library(future.apply)
  })
  plan(multisession, workers = N_WORKERS)
  on.exit(plan(sequential), add = TRUE)
  all_res <- future_lapply(seq_len(nrow(tasks)), run_task,
                           future.seed = TRUE,
                           future.globals = c("tasks", "cells", "ROOT",
                                              "MC_SOURCE_FILES", "source_all",
                                              "RUN_A", "RUN_B", "BASE_SEED",
                                              "DGP_IDX", "PARAMSET_IDX", "N_IDX",
                                              "mc_seed", "CHECKPOINT", "CKPT_DIR"),
                           future.packages = c("spdep", "spatialreg",
                                               "numDeriv", "Matrix"))
}

## ---------------------------------------------------------------- collect
A_list <- list(); B_list <- list(); err_list <- list(); time_list <- list()
for (r in all_res) {
  if (!is.null(r$error)) {
    err_list[[length(err_list) + 1]] <- data.frame(
      dgp = r$cell$dgp, paramset = r$cell$paramset, n = r$cell$n,
      rep = r$cell$rep, error = r$error, stringsAsFactors = FALSE)
  }
  if (!is.null(r$A)) A_list[[length(A_list) + 1]] <- r$A
  if (!is.null(r$B)) B_list[[length(B_list) + 1]] <- r$B
  if (!is.null(r$timing)) time_list[[length(time_list) + 1]] <- r$timing
}

if (length(time_list) > 0) {
  timings <- do.call(rbind, time_list)
  write.csv(timings, file.path(OUTDIR, "mc_timings.csv"), row.names = FALSE)
  tm <- aggregate(cbind(sec_gen, sec_A, sec_B) ~ dgp + paramset + n,
                  data = timings, FUN = function(x) round(mean(x, na.rm = TRUE), 1),
                  na.action = na.pass)
  cat("mean per-rep seconds (gen / A / B):\n")
  print(tm, row.names = FALSE)
}

if (length(A_list) > 0) {
  results_A <- do.call(rbind, A_list)
  write.csv(results_A, file.path(OUTDIR, "results_A.csv"), row.names = FALSE)
  cat(sprintf("results_A.csv: %d rows (%d reps)\n",
              nrow(results_A), length(A_list)))

  ## prop1 aggregates per cell -> identifiability_diag.csv
  key_cols <- c("dgp", "paramset", "n")
  reps1 <- results_A[!duplicated(results_A[, c(key_cols, "rep")]), ]
  agg <- aggregate(cbind(prop1_min_sv) ~ dgp + paramset + n, data = reps1,
                   FUN = function(x) c(min = min(x), median = median(x)))
  pass <- aggregate(prop1_rank_ok ~ dgp + paramset + n, data = reps1, FUN = mean)
  for (i in seq_len(nrow(agg))) {
    prop2_rows[[length(prop2_rows) + 1]] <- data.frame(
      scope = "cell", n = agg$n[i], dgp = agg$dgp[i], paramset = agg$paramset[i],
      prop1_pass_rate = pass$prop1_rank_ok[i],
      prop1_min_sv_min = agg$prop1_min_sv[i, "min"],
      prop1_min_sv_median = agg$prop1_min_sv[i, "median"],
      prop2_min_eig = NA, W_symmetric = NA, stringsAsFactors = FALSE)
  }
}

if (length(B_list) > 0) {
  results_B <- do.call(rbind, B_list)
  write.csv(results_B, file.path(OUTDIR, "results_B.csv"), row.names = FALSE)
  cat(sprintf("results_B.csv: %d rows\n", nrow(results_B)))
}

write.csv(do.call(rbind, prop2_rows),
          file.path(OUTDIR, "identifiability_diag.csv"), row.names = FALSE)

if (length(err_list) > 0) {
  errs <- do.call(rbind, err_list)
  write.csv(errs, file.path(OUTDIR, "mc_errors.csv"), row.names = FALSE)
  cat(sprintf("WARNING: %d replication(s) failed -> mc_errors.csv\n", nrow(errs)))
  if (CHECKPOINT) cat("checkpoints kept for resume:", CKPT_DIR, "\n")
} else {
  cat("no replication errors\n")
  if (CHECKPOINT && dir.exists(CKPT_DIR)) {
    unlink(CKPT_DIR, recursive = TRUE)   # outputs written; checkpoints redundant
    cat("checkpoints removed (run fully successful)\n")
  }
}

elapsed <- difftime(Sys.time(), t_start, units = "mins")
info <- c(
  sprintf("timestamp: %s", format(Sys.time())),
  sprintf("reps=%d workers=%d serial=%s A=%s B=%s", N_REPS, N_WORKERS, SERIAL, RUN_A, RUN_B),
  sprintf("dgps=%s paramsets=%s ns=%s", paste(DGPS, collapse = ","),
          paste(PARAMSETS, collapse = ","), paste(NS, collapse = ",")),
  sprintf("base_seed=%d", BASE_SEED),
  sprintf("elapsed_min=%.2f", as.numeric(elapsed)),
  sprintf("failed=%d/%d", length(err_list), nrow(tasks))
)
writeLines(info, file.path(OUTDIR, "mc_run_info.txt"))
cat(sprintf("done: %.1f min\n", as.numeric(elapsed)))
