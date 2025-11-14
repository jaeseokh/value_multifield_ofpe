# ============================================================
# functions_for_trees.R — tree-model utilities for LOLO
#   - Group split: leave-one-location-out by farm_field
#   - Monotone XGB (+1 for n_rate), RF baseline
#   - Curve grid at TEST medians (vary N only)
#   - Roughness metric (sum of squared 2nd differences)
#   - Optional mgcv GAM baseline (fit on TEST field only)
#   - EONR via grid profit search across scenarios
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(tibble); library(stringr)
  library(ranger)
  library(xgboost)
  library(mgcv)     # <— use mgcv instead of scam
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
  td <- train_df[, c(y_var, x_vars), drop = FALSE]
  td <- td %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(x_vars), ~ if (is.character(.)) factor(.) else .)) %>%
    tidyr::drop_na(dplyr::all_of(c(y_var, x_vars)))

  n <- nrow(td)
  if (n < 10) {
    warning(sprintf("RF skipped: too few training rows after drop_na (%d).", n))
    return(NULL)
  }

  # Safe defaults for tiny n
  mtry_default  <- max(1, floor(length(x_vars) / 4))
  min_node_safe <- min(rf$min.node.size %||% 40, max(2, floor(n / 5)))
  samp_frac     <- max(rf$sample.fraction %||% 0.7, min(1, 1 / n + 1e-8))  # ensure >= 1 obs
  replace_arg   <- rf$replace %||% TRUE

  fml <- as.formula(paste(y_var, "~", paste(x_vars, collapse = " + ")))

  out <- tryCatch(
    ranger::ranger(
      formula         = fml,
      data            = td,
      num.trees       = rf$num.trees %||% 800,
      mtry            = rf$mtry %||% mtry_default,
      min.node.size   = min_node_safe,
      max.depth       = rf$max.depth %||% 12,
      sample.fraction = samp_frac,
      replace         = replace_arg,
      importance      = rf$importance %||% "impurity",
      seed            = rf$seed %||% 123
    ),
    error = function(e) {
      warning(sprintf("RF retry with sample.fraction=1.0 & min.node.size=2 due to: %s", e$message))
      tryCatch(
        ranger::ranger(
          formula         = fml,
          data            = td,
          num.trees       = rf$num.trees %||% 800,
          mtry            = rf$mtry %||% mtry_default,
          min.node.size   = 2,
          max.depth       = rf$max.depth %||% 12,
          sample.fraction = 1.0,
          replace         = TRUE,
          importance      = rf$importance %||% "impurity",
          seed            = rf$seed %||% 123
        ),
        error = function(e2) {
          warning(sprintf("RF failed after retry: %s. Returning NULL.", e2$message))
          NULL
        }
      )
    }
  )

  out
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

  tibble::as_tibble(out)   # RF/XGB columns only
}

# NOTE: medians from a provided reference df (e.g., TEST field)
predict_curve_models <- function(bundle, x_vars, n_col, n_seq, ref_df) {
  rep_row <- representative_row(ref_df, x_vars, exclude = n_col)
  grid <- tibble(!!n_col := n_seq) %>%
    tidyr::crossing(rep_row) %>%
    dplyr::relocate(dplyr::all_of(n_col))
  preds <- predict_models(bundle, grid)
  dplyr::bind_cols(grid, preds)
}

# ------------------------- mgcv GAM baseline ----------------

build_mgcv_formula <- function(y_var, x_vars, n_col, k_all = 6, smooth_others = FALSE) {
  # Smooth N with cubic regression spline; others linear by default for speed
  n_term <- sprintf("s(%s, k=%d, bs='cr')", n_col, k_all)
  if (isTRUE(smooth_others)) {
    other <- setdiff(x_vars, n_col)
    other_terms <- if (length(other)) sprintf("s(%s, k=4, bs='cr')", other) else character(0)
    rhs <- c(n_term, other_terms)
  } else {
    rhs <- c(n_term, setdiff(x_vars, n_col))
  }
  as.formula(paste(y_var, "~", paste(rhs, collapse = " + ")))
}

fit_mgcv <- function(df, y_var, x_vars, n_col, k_all = 6, smooth_others = FALSE) {
  dd <- df %>%
    mutate(across(all_of(x_vars), ~ if (is.character(.)) factor(.) else .)) %>%
    tidyr::drop_na(all_of(c(y_var, x_vars)))
  fml <- build_mgcv_formula(y_var, x_vars, n_col, k_all = k_all, smooth_others = smooth_others)
  mgcv::gam(fml, data = dd, method = "REML", select = TRUE)
}

predict_curve_mgcv <- function(mgcv_fit, ref_df, x_vars, n_col, n_seq) {
  rep_row <- representative_row(ref_df, x_vars, exclude = n_col)
  grid <- tibble(!!n_col := n_seq) %>%
    tidyr::crossing(rep_row) %>%
    dplyr::relocate(dplyr::all_of(n_col))
  yh <- as.numeric(predict(mgcv_fit, newdata = grid, type = "response"))
  dplyr::bind_cols(grid, tibble(GAM = yh))
}

# ---- Holdout metrics helper (RMSE/MAE) ----
eval_holdout <- function(y_true, preds_df) {
  # guardrails
  if (is.null(preds_df) || !nrow(preds_df)) {
    return(tibble::tibble(model = character(), RMSE = numeric(), MAE = numeric()))
  }
  if (!requireNamespace("Metrics", quietly = TRUE)) {
    stop("Package 'Metrics' is required for RMSE/MAE. install.packages('Metrics')")
  }
  mods <- names(preds_df)
  purrr::map_dfr(mods, function(m) {
    p <- preds_df[[m]]
    ok <- is.finite(y_true) & is.finite(p)
    if (!any(ok)) {
      tibble::tibble(model = m, RMSE = NA_real_, MAE = NA_real_)
    } else {
      tibble::tibble(
        model = m,
        RMSE  = Metrics::rmse(y_true[ok], p[ok]),
        MAE   = Metrics::mae (y_true[ok], p[ok])
      )
    }
  })
}

# ---------------------- LOLO main ---------------------------

lolo_run_once <- function(df, test_ffy_id, y_var, x_vars, n_col,
                          use_rf, use_xgb, rf, xgb,
                          n_points = 100, use_gam = TRUE, gam_k = 6, gam_smooth_others = FALSE) {

  # Identify held-out location from the requested field-year
  test_loc <- df %>% dplyr::filter(ffy_id == test_ffy_id) %>% dplyr::distinct(farm_field) %>% dplyr::pull()
  if (length(test_loc) != 1) stop("Unique farm_field not found for test ffy_id: ", test_ffy_id)

  # Split
  train_df <- df %>% dplyr::filter(farm_field != test_loc)
  test_df  <- df %>% dplyr::filter(ffy_id == test_ffy_id)

  # Fit bundle on TRAIN locations only
  bundle <- fit_bundle(train_df, y_var, x_vars, n_col, use_rf, use_xgb, rf, xgb)

  # Holdout predictions on real subplot rows
  ph <- predict_models(bundle, test_df)
  model_cols_metrics <- intersect(c("RF","XGB"), names(ph))
  metrics_models <- eval_holdout(test_df[[y_var]], ph %>% dplyr::select(dplyr::all_of(model_cols_metrics)))

  # N grid from TEST, covariates fixed at TEST medians
  n_seq  <- build_n_seq(test_df, n_col, n_points)
  curves <- predict_curve_models(bundle, x_vars, n_col, n_seq, ref_df = test_df)

  # Roughness per model (RF/XGB curves)
  model_cols_curves <- intersect(c("RF","XGB"), names(curves))
  curve_long <- curves %>%
    tidyr::pivot_longer(cols = dplyr::all_of(model_cols_curves),
                        names_to = "model", values_to = "yhat")

  rough_tbl <- curve_long %>%
    dplyr::group_by(model) %>%
    dplyr::summarise(roughness = roughness_ssd(.data[[n_col]], yhat), .groups = "drop")

  # Optional mgcv baseline (fit on TEST DF only to reflect that field's shape)
  gam_bits <- NULL
  if (isTRUE(use_gam)) {
    gam_fit   <- fit_mgcv(test_df, y_var, x_vars, n_col, k_all = gam_k, smooth_others = gam_smooth_others)
    gam_hold  <- tibble::tibble(GAM = as.numeric(predict(gam_fit, newdata = test_df, type = "response")))
    met_gam   <- eval_holdout(test_df[[y_var]], gam_hold) %>% dplyr::mutate(model = "GAM")
    curves_gam <- predict_curve_mgcv(gam_fit, ref_df = test_df, x_vars, n_col, n_seq)
    rough_gam  <- tibble::tibble(model = "GAM",
                                 roughness = roughness_ssd(curves_gam[[n_col]], curves_gam$GAM))
    gam_bits <- list(gam = gam_fit, curves = curves_gam, metrics = met_gam, rough = rough_gam)
  }

  list(
    bundle = bundle,
    test_ffy_id = test_ffy_id,
    test_loc = test_loc,
    holdout_metrics = metrics_models,
    curve_table = curves,          # RF/XGB
    curve_roughness = rough_tbl,   # RF/XGB
    gam = gam_bits                 # GAM curves/metrics in $gam
  )
}
