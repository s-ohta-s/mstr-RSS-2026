# data_prep_fast.r
# Two-stage data preparation for the Monte Carlo driver

#' Static spatial quantities for a given W (compute once per n)
#'
#' @param W n×n row-standardized spatial weight matrix (dense), OR
#'          a path to a weight-matrix CSV (same format as prepare_data_extended)
#' @param verbose verbose output
#' @return list(W, W_sp, eigen_W, W_listw, n)
prepare_static <- function(W, verbose = FALSE) {

  if (is.character(W)) {
    W_raw <- read.csv(W, header = TRUE, stringsAsFactors = FALSE)
    first_col <- W_raw[, 1]
    is_id_column <- !is.numeric(first_col) || all(first_col == 1:nrow(W_raw))
    W <- if (is_id_column) as.matrix(W_raw[, -1]) else as.matrix(W_raw)
  }
  dimnames(W) <- NULL
  n <- nrow(W)

  row_sums <- rowSums(W)
  if (!all(abs(row_sums - 1) < 1e-6)) {
    warning("The spatial weight matrix may not be row-standardized")
  }

  if (verbose) cat(sprintf("prepare_static: n = %d, computing eigen(W) / listw ...\n", n))

  list(
    W       = W,
    W_sp    = if (requireNamespace("Matrix", quietly = TRUE)) Matrix::Matrix(W, sparse = TRUE) else NULL,
    eigen_W = compute_eigen_W(W),
    W_listw = spdep::mat2listw(W, style = "W"),
    n       = n
  )
}

#' Assemble a data_list from in-memory data frames (per replication)
#'
#' Produces the same structure as prepare_data_extended(), reusing the static
#' quantities from prepare_static(). df_current / df_lag must each contain one
#' row per region, with the region variable identifying the ordering.
#'
#' @param df_current data.frame at the estimation time point (n rows)
#' @param df_lag data.frame at the lag time point (n rows); NULL if
#'        include_time_lag = FALSE
#' @param static return value of prepare_static()
#' @param y_vars vector of dependent variable names (e.g., c("y1","y2"))
#' @param x_vars named list of regressor names per response
#' @param region_var region id variable name
#' @param include_intercept include an intercept column per equation
#' @param include_time_lag include lagged responses as regressors
#' @param time_point,time_lag labels stored in data_info (documentation only)
#' @return data_list compatible with all fit_* functions
attach_data <- function(
  df_current,
  df_lag,
  static,
  y_vars,
  x_vars,
  region_var = "region",
  include_intercept = TRUE,
  include_time_lag = TRUE,
  time_point = NULL,
  time_lag = NULL
) {

  k <- length(y_vars)
  n <- static$n

  regions <- sort(unique(df_current[[region_var]]))
  if (length(regions) != n || nrow(df_current) != n) {
    stop(sprintf("attach_data: df_current must have one row per region (%d rows, n=%d)",
                 nrow(df_current), n))
  }
  df_current <- df_current[order(df_current[[region_var]]), ]

  y_lag <- NULL
  if (include_time_lag) {
    if (is.null(df_lag)) stop("attach_data: df_lag is required when include_time_lag = TRUE")
    if (nrow(df_lag) != n) stop("attach_data: df_lag must have one row per region")
    df_lag <- df_lag[order(df_lag[[region_var]]), ]
    y_lag <- matrix(NA_real_, nrow = n, ncol = k)
    for (i in 1:k) y_lag[, i] <- df_lag[[y_vars[i]]]
  }

  # Stacked response vector (kn×1)
  y <- numeric(k * n)
  for (i in 1:k) y[((i - 1) * n + 1):(i * n)] <- df_current[[y_vars[i]]]

  # Design matrix (identical construction to prepare_data_extended)
  X <- build_design_matrix_extended(
    data = df_current,
    y_vars = y_vars,
    x_vars = x_vars,
    y_lag = y_lag,
    include_intercept = include_intercept,
    include_time_lag = include_time_lag,
    n = n,
    k = k,
    verbose = FALSE
  )

  p0 <- if (include_time_lag) ncol(X) - k^2 else ncol(X)

  list(
    y = y,
    X = X,
    W = static$W,
    W_sp = static$W_sp,
    W_listw = static$W_listw,
    y_lag = y_lag,
    eigen_W = static$eigen_W,
    n = n,
    k = k,
    regions = regions,
    p0 = p0,
    data_info = list(
      y_vars = y_vars,
      x_vars = x_vars,
      time_point = time_point,
      time_lag = time_lag,
      region_var = region_var,
      time_var = "time",
      include_intercept = include_intercept,
      include_time_lag = include_time_lag
    )
  )
}
