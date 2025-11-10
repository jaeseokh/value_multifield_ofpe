# ============================================================
# functions_for_ensembles.R — ensemble training utilities
#   - Group split: leave-one-location-out by farm_field
#   - Monotone XGB (+1 for n_rate), smoother RF
#   - Curve grid at representative medians (vary n only)
#   - Roughness metric (sum of squared 2nd differences)
#   - Optional SCAM baseline for shape-constrained reference
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(stringr)
  library(ranger)
  library(xgboost)
  suppressWarnings(suppressMessages(requireNamespace("scam", quietly = TRUE)))
  library(Metrics)
})

`%||%` <- function(a, b) if (is.null(a)) b else a
is_num  <- function(x) is.numeric(x) || is.integer(x)

# ----------------------- IDs & Splits -----------------------

parse_farm_field <- function(ffy_id, sep = "_") {
  parts <- str_split(ffy_id, fixed(sep), simplify = TRUE)
  if (ncol(parts) < 3) stop("ffy_id must look like 'farm_field_year' e.g. '8_1_2023'.")
  paste(parts[,1], parts[,2], sep = "_")
}

add_farm_field <- function(df, ffy_col = "ffy_id", out_col = "farm_field", sep = "_") {
  if (!ffy_col %in% names(df)) stop("Missing ffy_id column: ", ffy_col)
  mutate(df, "{out_col}" := parse_farm_field(.data[[ffy_col]], sep = sep))
}

make_lolo_index <- function(df, location_col = "farm_field", ffy_col = "ffy_id") {
  tests <- df %>% distinct(!!sym(location_col), !!sym(ffy_col)) %>% arrange(!!sym(location_col), !!sym(ffy_col))
  list(locations = unique(tests[[location_col]]), tests = tests)
}

# ------------------ Representative row & N grid --------------

representative_row <- function(df, x_vars, exclude = NULL) {
  keep <- setdiff(x_vars, exclude %||% character(0))
  vals <- lapply(keep, function(v) {
    x <- df[[v]]
    if (is.factor(x)) factor(names(sort(table(x), decreasing = TRUE))[1], levels = levels(x))
    else if (is.character(x)) names(sort(table(x), decreasing = TRUE))[1]
    else stats::median(x, na.rm = TRUE)
  })
  names(vals) <- keep
  as_tibble(vals)
}

build_n_seq <- function(df, n_col, n_points = 100) {
  r <- range(df[[n_col]], na.rm = TRUE)
  if (!all(is.finite(r))) stop("Non-finite n_col range.")
  if (r[1] == r[2]) r <- c(r[1] - 1e-6, r[2] + 1e-6)
  seq(r[1], r[2], length.out = n_points)
}

# ------------------------- Roughness -------------------------

roughness_ssd <- function(x, y) {
  ord <- order(x)
  y2  <- diff(y[ord], differences = 2)
  sum(y2^2, na.rm = TRUE)
}

# ---------------------- Model matrices ----------------------

.mm_build <- function(df, x_vars) model.matrix(~ . - 1, data = df[, x_vars, drop = FALSE])
.mm_align <- function(df, x_vars, ref_cols) {
  mm <- .mm_build(df, x_vars)
  out <- matrix(0, nrow = nrow(mm), ncol = length(ref_cols))
  colnames(out) <- ref_cols
  common <- intersect(colnames(mm), ref_cols)
  out[, common] <- mm[, common, drop = FALSE]
  out
}

# ---------------------- Fit base learners -------------------

fit_rf <- function(train_df, y_var, x_vars, rf = list()) {
  fml <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))
  ranger::ranger(
    formula         = fml,
    data            = train_df[, c(y_var, x_vars), drop = FALSE],
    num.trees       = rf$num.trees %||% 800,
    mtry            = rf$mtry %||% max(1, floor(length(x_vars)/4)),
    min.node.size   = rf$min.node.size %||% 40,
    max.depth       = rf$max.depth %||% 12,
    sample.fraction = rf$sample.fraction %||% 0.7,
    seed            = rf$seed %||% 123
  )
}

fit_xgb <- function(train_df, y_var, x_vars, n_col, xgb = list()) {
  mm  <- .mm_build(train_df, x_vars)
  dtr <- xgboost::xgb.DMatrix(data = mm, label = train_df[[y_var]])

  mono_vec <- if ((xgb$monotone %||% TRUE)) {
    as.numeric(colnames(mm) == n_col)
  } else rep(0, ncol(mm))

  params <- list(
    objective = "reg:squarederror",
    eta = xgb$eta %||% 0.03,
    max_depth = xgb$max_depth %||% 3,
    min_child_weight = xgb$min_child_weight %||% 20,
    gamma = xgb$gamma %||% 2,
    subsample = xgb$subsample %||% 0.7,
    colsample_bytree = xgb$colsample_bytree %||% 0.7,
    lambda = xgb$lambda %||% 2.0,
    alpha  = xgb$alpha %||% 0.0,
    monotone_constraints = paste0("(", paste(mono_vec, collapse = ","), ")")
  )

  bst <- xgboost::xgb.train(
    params = params, data = dtr,
    nrounds = xgb$nrounds %||% 1200,
    watchlist = list(train = dtr), verbose = 0
  )
  list(fit = bst, refcols = colnames(mm))
}

fit_bundle <- function(train_df, y_var, x_vars, n_col,
                       use_rf = TRUE, use_xgb = TRUE,
                       rf = list(), xgb = list()) {

  td <- train_df %>%
    mutate(across(all_of(x_vars), ~ if (is.character(.)) factor(.) else .)) %>%
    tidyr::drop_na(all_of(c(y_var, x_vars)))

  fits <- list(); meta <- list()

  if (isTRUE(use_rf))  fits$RF  <- fit_rf(td, y_var, x_vars, rf)
  if (isTRUE(use_xgb)) {
    xgb_fit <- fit_xgb(td, y_var, x_vars, n_col, xgb)
    fits$XGB <- xgb_fit$fit
    meta$XGB_refcols <- xgb_fit$refcols
  }

  list(fits = fits, meta = meta, x_vars = x_vars, y_var = y_var, n_col = n_col, train_df = td)
}

# -------------------------- Predict -------------------------

predict_models <- function(bundle, newdata) {
  xvars <- bundle$x_vars

  present <- intersect(xvars, names(newdata))
  missing <- setdiff(xvars, names(newdata))
  if (length(missing)) warning("newdata is missing predictors: ", paste(missing, collapse = ", "),
                               ". Proceeding with present predictors only.")
  if (!length(present)) stop("No predictor columns available in newdata.")

  nd <- newdata %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(present),
                                ~ if (is.character(.)) factor(.) else .))

  out <- list()

  if (!is.null(bundle$fits$RF)) {
    out$RF <- predict(bundle$fits$RF, nd[, present, drop = FALSE])$predictions
  }

  if (!is.null(bundle$fits$XGB)) {
    mm <- .mm_align(nd[, present, drop = FALSE], present, bundle$meta$XGB_refcols)
    out$XGB <- predict(bundle$fits$XGB, xgboost::xgb.DMatrix(mm))
  }

  tibble::as_tibble(out) %>%
    dplyr::mutate(ENSEMBLE = rowMeans(dplyr::across(everything()), na.rm = TRUE))
}

predict_curve_models <- function(bundle, x_vars, n_col, n_seq) {
  rep_row <- representative_row(bundle$train_df, x_vars, exclude = n_col)
  grid <- tibble(!!n_col := n_seq) %>%
    tidyr::crossing(rep_row) %>%
    dplyr::relocate(dplyr::all_of(n_col))
  preds <- predict_models(bundle, grid)
  dplyr::bind_cols(grid, preds)
}

# ------------------------- SCAM (optional) ------------------

build_scam_formula <- function(y_var, x_vars, n_col, k_all = 6) {
  n_term <- sprintf("s(%s, k=%d, bs='mpi')", n_col, k_all) # monotone ↑ in N
  rhs <- c(n_term, setdiff(x_vars, n_col))
  as.formula(paste(y_var, "~", paste(rhs, collapse = " + ")))
}

fit_scam <- function(df, y_var, x_vars, n_col, k_all = 6) {
  if (!requireNamespace("scam", quietly = TRUE)) stop("Package 'scam' required for SCAM fit.")
  dd <- df %>% mutate(across(all_of(x_vars), ~ if (is.character(.)) factor(.) else .)) %>%
    drop_na(all_of(c(y_var, x_vars)))
  scam::scam(build_scam_formula(y_var, x_vars, n_col, k_all), data = dd)
}

predict_curve_scam <- function(scam_fit, ref_df, x_vars, n_col, n_seq) {
  rep_row <- representative_row(ref_df, x_vars, exclude = n_col)
  grid <- tibble(!!n_col := n_seq) %>% tidyr::crossing(rep_row) %>% dplyr::relocate(dplyr::all_of(n_col))
  yh <- as.numeric(predict(scam_fit, newdata = grid, type = "response"))
  dplyr::bind_cols(grid, tibble(SCAM = yh))
}

# ---------------------- LOLO evaluation ---------------------

eval_holdout <- function(y_true, preds_df) {
  mods <- colnames(preds_df)
  dplyr::bind_rows(lapply(mods, function(m) {
    tibble(model = m,
           RMSE = Metrics::rmse(y_true, preds_df[[m]]),
           MAE  = Metrics::mae(y_true,  preds_df[[m]]))
  }))
}

lolo_run_once <- function(df, test_ffy_id, y_var, x_vars, n_col,
                          use_rf, use_xgb, rf, xgb,
                          n_points = 100, use_scam = TRUE, scam_k = 6) {

  test_loc <- df %>% dplyr::filter(ffy_id == test_ffy_id) %>% dplyr::distinct(farm_field) %>% dplyr::pull()
  if (length(test_loc) != 1) stop("Unique farm_field not found for test ffy_id: ", test_ffy_id)

  train_df <- df %>% dplyr::filter(farm_field != test_loc)
  test_df  <- df %>% dplyr::filter(ffy_id == test_ffy_id)

  bundle <- fit_bundle(train_df, y_var, x_vars, n_col, use_rf, use_xgb, rf, xgb)

  # Holdout predictions on real subplot rows
  ph <- predict_models(bundle, test_df)

  model_cols_metrics <- intersect(c("RF","XGB","ENSEMBLE"), names(ph))
  metrics_models <- eval_holdout(test_df[[y_var]], ph %>% dplyr::select(dplyr::all_of(model_cols_metrics)))

  # Curve & roughness (median covariates)
  n_seq   <- build_n_seq(train_df, n_col, n_points)
  curves  <- predict_curve_models(bundle, x_vars, n_col, n_seq)

  model_cols_curves <- intersect(c("RF","XGB","ENSEMBLE"), names(curves))
  curve_long <- curves %>%
    tidyr::pivot_longer(cols = dplyr::all_of(model_cols_curves),
                        names_to = "model", values_to = "yhat")

  rough_tbl <- curve_long %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(roughness = roughness_ssd(.data[[n_col]], yhat), .groups = "drop")

  # Optional SCAM baseline
  scam_bits <- NULL
  if (isTRUE(use_scam)) {
    sf <- fit_scam(test_df, y_var, x_vars, n_col, k_all = scam_k)
    scam_hold <- tibble(SCAM = as.numeric(predict(sf, newdata = test_df, type = "response")))
    met_scam  <- eval_holdout(test_df[[y_var]], scam_hold) %>% dplyr::mutate(model = "SCAM")
    curves_scam <- predict_curve_scam(sf, ref_df = train_df, x_vars, n_col, n_seq)
    rough_scam  <- tibble(model = "SCAM",
                          roughness = roughness_ssd(curves_scam[[n_col]], curves_scam$SCAM))
    scam_bits <- list(sf = sf, curves = curves_scam, metrics = met_scam, rough = rough_scam)
  }

  list(
    bundle = bundle,
    test_ffy_id = test_ffy_id,
    test_loc = test_loc,
    holdout_metrics = metrics_models,
    curve_table = curves,
    curve_roughness = rough_tbl,
    scam = scam_bits
  )
}
