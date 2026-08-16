# mc_experiment.r
# Monte Carlo experiment core
# Both experiments share the same generated data within a replication.

## Standard variable configuration (fixed across the whole study)
MC_Y_VARS <- c("y1", "y2")
MC_X_VARS <- list(
  y1 = c("x_common1", "x_common2", "x_specific1_1"),
  y2 = c("x_common1", "x_common2", "x_specific2_1")
)

## Full (penalized) candidate models: fit function, likelihood type, time lag
MC_FULL_MODELS <- list(
  "1111" = list(fit = "fit_mgns_penalized", type = "MGNS", time_lag = TRUE),
  "0111" = list(fit = "fit_msem_penalized", type = "MSEM", time_lag = TRUE),
  "1011" = list(fit = "fit_msar_penalized", type = "MSAR", time_lag = TRUE),
  "1101" = list(fit = "fit_mgns_penalized", type = "MGNS", time_lag = FALSE),
  "0101" = list(fit = "fit_msem_penalized", type = "MSEM", time_lag = FALSE),
  "1001" = list(fit = "fit_msar_penalized", type = "MSAR", time_lag = FALSE)
)

## Simple candidate models (γ-free): fit function + whether W is needed
MC_SIMPLE_MODELS <- list(
  "0011" = list(fit = "fit_varx",   needs_W = FALSE),
  "000d" = list(fit = "fit_indreg", needs_W = FALSE),
  "d0dd" = list(fit = "fit_indsar", needs_W = TRUE),
  "0ddd" = list(fit = "fit_indsem", needs_W = TRUE),
  "dddd" = list(fit = "fit_indgns", needs_W = TRUE)
)

## Default single-stage γ grid for the MC
MC_GAMMAS <- c(0, 10^seq(-2, 3, by = 0.5))

#' Quiet single-stage γ search with GIC at every grid point
#'
#' pAIC/pBIC selection requires df_eff(γ) — i.e. the profile Hessian — at each
#' candidate γ, exactly as in compare_gamma(), but without console output and
#' returning structured results. The Hessian at γ = 0 doubles as the local-
#' identification diagnostic (min eigenvalue).
#'
#' @param model_id one of names(MC_FULL_MODELS)
#' @param data_list prebuilt data_list (attach_data), lag config must match
#' @param gammas γ grid
#' @return list(summary, best_pAIC, best_pBIC, hess_min_eig_gamma0) or NULL if
#'         every fit failed. best_* = list(gamma, res, pAIC, pBIC, df_eff)
gamma_search_mc <- function(model_id, data_list, gammas = MC_GAMMAS) {

  cfg <- MC_FULL_MODELS[[model_id]]
  fit_fun <- get(cfg$fit)
  k <- data_list$k

  k1 <- if (cfg$type == "MGNS") 2 * k^2 else k^2

  # The S1 initial estimation (spatialreg per-variable fits) does not depend on
  # γ — at n = 900 it is ~86% of a single fit's cost — so run it ONCE and pass
  # the same inits to every γ point (numerically identical to the internal S1).
  init_args <- tryCatch({
    if (cfg$type == "MSAR") {
      ii <- initial_estimation_sly_extended(data_list, verbose = 0)
      list(R_init = ii$R_init, Sigma_init = ii$Sigma_init)
    } else if (cfg$type == "MSEM") {
      ii <- initial_estimation_sem_extended(data_list, verbose = 0)
      list(T_init = ii$T_init, Sigma_init = ii$Sigma_init)
    } else {
      ii <- initial_estimation_sdem_extended(data_list, verbose = 0)
      list(R_init = ii$R_init, T_init = ii$T_init, Sigma_init = ii$Sigma_init)
    }
  }, error = function(e) list())   # empty -> each fit falls back to internal S1

  rows <- list()
  results <- list()
  hess0 <- NA_real_

  for (g in gammas) {
    res <- tryCatch(
      do.call(fit_fun, c(
        list(data_file = NULL, weight_file = NULL,
             y_vars = MC_Y_VARS, x_vars = MC_X_VARS,
             gamma = g, verbose = FALSE, data_list = data_list),
        init_args)),
      error = function(e) NULL)
    if (is.null(res)) next

    k2 <- res$fit$num_params - k1

    H <- tryCatch({
      if (cfg$type == "MSAR") {
        compute_profile_hessian_msar(
          R_hat = res$coefficients$R, beta_hat = res$coefficients$beta,
          Sigma_hat = res$coefficients$Sigma, data_list = data_list)
      } else if (cfg$type == "MSEM") {
        compute_profile_hessian_msem(
          Lambda_hat = res$coefficients$Lambda, beta_hat = res$coefficients$beta,
          Sigma_hat = res$coefficients$Sigma, data_list = data_list)
      } else {
        compute_profile_hessian_mgns(
          R_hat = res$coefficients$R, Lambda_hat = res$coefficients$Lambda,
          beta_hat = res$coefficients$beta, Sigma_hat = res$coefficients$Sigma,
          data_list = data_list)
      }
    }, error = function(e) NULL)

    if (!is.null(H)) {
      gic <- tryCatch(
        compute_gic(H = H, gamma = g, loglik = res$fit$loglik,
                    k2 = k2, n_obs = res$fit$num_obs),
        error = function(e) NULL)
    } else gic <- NULL

    if (g == 0 && !is.null(H)) {
      hess0 <- min(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values)
    }

    df_eff <- if (!is.null(gic)) gic$df_eff else res$fit$num_params
    paic   <- if (!is.null(gic)) gic$GIC_AIC else res$fit$AIC
    pbic   <- if (!is.null(gic)) gic$GIC_BIC else res$fit$BIC

    results[[as.character(g)]] <- res
    rows[[length(rows) + 1]] <- data.frame(
      gamma = g, loglik = res$fit$loglik, AIC = res$fit$AIC, BIC = res$fit$BIC,
      df_eff = df_eff, GIC_AIC = paic, GIC_BIC = pbic,
      converged = isTRUE(res$convergence$converged),
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) return(NULL)
  s <- do.call(rbind, rows)

  pick <- function(crit) {
    idx <- which.min(s[[crit]])
    list(gamma = s$gamma[idx],
         res = results[[as.character(s$gamma[idx])]],
         pAIC = s$GIC_AIC[idx], pBIC = s$GIC_BIC[idx],
         df_eff = s$df_eff[idx])
  }

  list(summary = s,
       best_pAIC = pick("GIC_AIC"),
       best_pBIC = pick("GIC_BIC"),
       hess_min_eig_gamma0 = hess0)
}

#' Experiment A: fit the true model, return tidy per-parameter rows
#'
#' @param dgp_id "1111" | "0111" | "1011" (also the fitted model)
#' @param params get_dgp_params() output (true values)
#' @param data_list attach_data() output (include_time_lag = TRUE)
#' @param gammas γ grid
#' @param check_temporal FALSE for conditional (realistic) cells
#' @param se_full also run add_full_inference (full-Hessian sandwich SE);
#'        expensive — OFF by default, covered95_hess is then NA
#' @return data.frame, one row per parameter (24 rows for K = 2)
mc_experiment_A <- function(dgp_id, params, data_list, gammas = MC_GAMMAS,
                            check_temporal = TRUE, se_full = FALSE,
                            gs = NULL) {

  # gs: precomputed gamma_search_mc() output (shared with experiment B when
  # both run on the same replication — avoids repeating the true model's search)
  if (is.null(gs)) gs <- gamma_search_mc(dgp_id, data_list, gammas)
  if (is.null(gs)) stop(sprintf("experiment A: all fits failed for %s", dgp_id))

  best <- gs$best_pAIC
  res <- best$res

  # SEs: β/α via Ψ (paper Eq. (18)), spatial via the
  # penalized numerical Hessian — both provided by add_inference(gamma = γ*)
  res <- tryCatch(
    add_inference(res, compute_spatial_se = TRUE, gamma = best$gamma, verbose = FALSE),
    error = function(e) res)
  if (se_full) {
    res <- tryCatch(add_full_inference(res, gamma = best$gamma, verbose = FALSE),
                    error = function(e) res)
  }

  tab <- extract_params_uniform(res, dgp_id, gamma = best$gamma,
                                pAIC = best$pAIC, pBIC = best$pBIC,
                                d_eff = best$df_eff)

  true_df <- params_to_true_df(params)
  out <- merge(true_df, tab[, c("parameter", "estimate", "se_psi", "se_hessian")],
               by = "parameter", sort = FALSE)
  # merge() may reorder; restore the canonical parameter order
  out <- out[match(intersect(true_df$parameter, out$parameter), out$parameter), ]
  # drop parameters the fitted model does not estimate (e.g. R rows under an
  # MSEM fit, Λ rows under an MSAR fit — extract emits them as NA)
  out <- out[!is.na(out$estimate), ]

  cover <- function(true, est, se) {
    ifelse(is.na(se) | se <= 0, NA_integer_,
           as.integer(abs(est - true) <= qnorm(0.975) * se))
  }
  out$covered95_psi  <- cover(out$true, out$estimate, out$se_psi)
  out$covered95_hess <- cover(out$true, out$estimate, out$se_hessian)
  # main criterion: the se_psi column holds Ψ-SE for β/α rows and
  # the penalized-profile-Hessian SE for R/Λ rows -> use it for all rows
  out$covered95_main <- out$covered95_psi

  # Replication-level diagnostics (repeated on every row of this rep)
  vd <- check_estimate_validity(res, data_list, check_temporal = check_temporal)
  beta_true <- params_to_beta_vector(params, include_time_lag = TRUE)
  p1 <- tryCatch(
    prop1_rank_margin(params$R, beta_true, data_list$X,
                      if (!is.null(data_list$W_sp)) data_list$W_sp else data_list$W,
                      k = data_list$k, n = data_list$n),
    error = function(e) list(rank_ok = NA, min_sv = NA_real_, cond = NA_real_))

  out$gamma      <- best$gamma
  out$gamma_pBIC <- gs$best_pBIC$gamma

  # pBIC 選択 γ でのフィットから推定値のみ抽出して併記する
  # (SE/カバレッジは pAIC 側と同様に add_inference が必要なため出力しない)
  tab_bic <- extract_params_uniform(gs$best_pBIC$res, dgp_id,
                                    gamma = gs$best_pBIC$gamma,
                                    pAIC = gs$best_pBIC$pAIC,
                                    pBIC = gs$best_pBIC$pBIC,
                                    d_eff = gs$best_pBIC$df_eff)
  out$estimate_pBIC <- tab_bic$estimate[match(out$parameter, tab_bic$parameter)]
  out$df_eff     <- best$df_eff
  out$loglik     <- res$fit$loglik
  out$converged  <- isTRUE(res$convergence$converged)
  out$stationarity_ok <- vd$stationarity_ok
  out$rho_dyn    <- vd$rho_dyn
  out$min_det_R      <- vd$min_det_R
  out$min_det_Lambda <- vd$min_det_Lambda
  out$prop1_rank_ok  <- p1$rank_ok
  out$prop1_min_sv   <- p1$min_sv
  out$hess_min_eig_gamma0 <- gs$hess_min_eig_gamma0

  rownames(out) <- NULL
  out
}

#' Experiment B: fit all 11 candidate models, return one selection row
#'
#' @param dgp_id true DGP ("1111" | "0111" | "1011")
#' @param data_list_lag attach_data() with include_time_lag = TRUE
#' @param data_list_nolag attach_data() with include_time_lag = FALSE
#' @param data_long long-format data.frame (rbind(df_lag, df_current)) for the
#'        simple models, W dense weight matrix
#' @param W dense spatial weight matrix
#' @param gammas γ grid for the penalized models
#' @return one-row data.frame (selection + per-model pAIC/pBIC values)
mc_experiment_B <- function(dgp_id, data_list_lag, data_list_nolag,
                            data_long, W, gammas = MC_GAMMAS,
                            keep_search_for = NULL) {

  model_ids <- c(names(MC_FULL_MODELS), names(MC_SIMPLE_MODELS))
  paic <- setNames(rep(NA_real_, length(model_ids)), model_ids)
  pbic <- paic
  kept_search <- NULL

  for (mid in names(MC_FULL_MODELS)) {
    dl <- if (MC_FULL_MODELS[[mid]]$time_lag) data_list_lag else data_list_nolag
    gs <- tryCatch(gamma_search_mc(mid, dl, gammas), error = function(e) NULL)
    if (!is.null(gs)) {
      paic[mid] <- gs$best_pAIC$pAIC
      pbic[mid] <- gs$best_pBIC$pBIC
      if (identical(mid, keep_search_for)) kept_search <- gs
    }
  }

  for (mid in names(MC_SIMPLE_MODELS)) {
    cfg <- MC_SIMPLE_MODELS[[mid]]
    fit_fun <- get(cfg$fit)
    args <- list(data_file = NULL, y_vars = MC_Y_VARS, x_vars = MC_X_VARS,
                 time_var = "time", time_point = 2, region_var = "region",
                 verbose = FALSE, data = data_long)
    if (cfg$needs_W) {
      args$weight_file <- NULL
      args$W_matrix <- W
      args$include_time_lag <- TRUE
    }
    res <- tryCatch(do.call(fit_fun, args), error = function(e) NULL)
    if (!is.null(res)) {
      paic[mid] <- res$fit$AIC   # γ = 0, df_eff = num_params -> pAIC = AIC
      pbic[mid] <- res$fit$BIC
    }
  }

  sel_paic <- if (all(is.na(paic))) NA_character_ else names(which.min(paic))
  sel_pbic <- if (all(is.na(pbic))) NA_character_ else names(which.min(pbic))

  row <- data.frame(
    selected_pAIC = sel_paic,
    selected_pBIC = sel_pbic,
    true_model = dgp_id,
    hit_pAIC = as.integer(identical(sel_paic, dgp_id)),
    hit_pBIC = as.integer(identical(sel_pbic, dgp_id)),
    n_models_ok = sum(!is.na(paic)),
    stringsAsFactors = FALSE
  )
  for (mid in model_ids) {
    row[[paste0("pAIC_", mid)]] <- paic[mid]
    row[[paste0("pBIC_", mid)]] <- pbic[mid]
  }
  list(row = row, kept_search = kept_search)
}

#' One Monte Carlo replication: generate data, run experiment A and/or B
#'
#' @param dgp_id,paramset_id,n_label cell identifiers
#' @param static prepare_static() output for this n
#' @param rep_id replication number (1-based)
#' @param seed RNG seed for the generator
#' @param run_A,run_B which experiments to run
#' @param gammas γ grid
#' @param solvers optional cached make_spatial_solvers() output
#' @return list(A = data.frame|NULL, B = data.frame|NULL, error = NULL|string)
mc_run_rep <- function(dgp_id, paramset_id, n_label, static, rep_id, seed,
                       run_A = TRUE, run_B = TRUE, gammas = MC_GAMMAS,
                       solvers = NULL) {

  out <- list(A = NULL, B = NULL, error = NULL, timing = NULL)

  tryCatch({
    params <- get_dgp_params(paramset_id, dgp_id)
    conditional <- identical(params$generation, "conditional")

    t_gen <- system.time(
      gen <- generate_cell_dataset(params, static, seed, solvers = solvers)
    )["elapsed"]

    dl_lag <- attach_data(gen$df_current, gen$df_lag, static,
                          y_vars = MC_Y_VARS, x_vars = MC_X_VARS,
                          include_time_lag = TRUE, time_point = 2, time_lag = 1)

    id_cols <- data.frame(dgp = dgp_id, paramset = paramset_id, n = static$n,
                          rep = rep_id, seed = seed, stringsAsFactors = FALSE)

    t_A <- NA_real_; t_B <- NA_real_
    gs_true <- NULL

    # B first when both run: its search of the true model is reused by A
    if (run_B) {
      t_B <- system.time({
        dl_nolag <- attach_data(gen$df_current, NULL, static,
                                y_vars = MC_Y_VARS, x_vars = MC_X_VARS,
                                include_time_lag = FALSE, time_point = 2)
        data_long <- rbind(gen$df_lag, gen$df_current)
        b <- mc_experiment_B(dgp_id, dl_lag, dl_nolag, data_long, static$W,
                             gammas = gammas,
                             keep_search_for = if (run_A) dgp_id else NULL)
      })["elapsed"]
      out$B <- cbind(id_cols, b$row)
      gs_true <- b$kept_search
    }

    if (run_A) {
      t_A <- system.time({
        a <- mc_experiment_A(dgp_id, params, dl_lag, gammas = gammas,
                             check_temporal = !conditional, gs = gs_true)
      })["elapsed"]
      out$A <- cbind(id_cols[rep(1, nrow(a)), ], a)
      rownames(out$A) <- NULL
    }

    out$timing <- cbind(id_cols, data.frame(
      sec_gen = as.numeric(t_gen), sec_A = as.numeric(t_A),
      sec_B = as.numeric(t_B), stringsAsFactors = FALSE))
  }, error = function(e) {
    out$error <<- conditionMessage(e)
  })

  out
}
