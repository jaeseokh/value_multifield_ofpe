# ============================================================
# functions_for_bayesian.R — curve-first Bayesian attribution
#   - Extract curve features per field-year (plateau, steepness, curvature, EONR)
#   - Build "clue" predictors from timing-weather & a few soil/topo interactions
#   - Fit sparse (horseshoe) Bayesian regressions of features on clues (brms)
#   - Produce N- vs Water-side driver scores with uncertainty
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr)
  library(stringr); library(tibble)
  library(ggplot2)
  library(brms)     # Bayesian regression
  library(matrixStats)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ------------------------------------------------------------
# 1) Curve features (computed from an N–yield grid per ffy_id)
#    Expected input: a data.frame with columns: ffy_id, n_rate, and
#    one or more model columns (e.g., GAM, RF, XGB).
#    We pick ONE column as the "observed curve" for feature extraction.
# ------------------------------------------------------------

finite_diff <- function(x, y) {
  # first derivative on equally/unequally spaced x via central differences
  n <- length(x); dy <- rep(NA_real_, n)
  ord <- order(x); x <- x[ord]; y <- y[ord]
  if (n >= 3) {
    for (i in 2:(n-1)) {
      dy[i] <- (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    }
  }
  dy[ord] <- dy
  dy
}

second_diff <- function(x, y) {
  # second derivative via discrete second difference on ordered x
  ord <- order(x); x <- x[ord]; y <- y[ord]
  n <- length(x); d2 <- rep(NA_real_, n)
  if (n >= 3) {
    for (i in 2:(n-1)) {
      d2[i] <- ( (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]) ) /
               (0.5*((x[i+1]-x[i])+(x[i]-x[i-1])))
    }
  }
  d2[ord] <- d2
  d2
}

roughness_ssd <- function(x, y) {
  ord <- order(x); y2 <- diff(y[ord], differences = 2)
  sum(y2^2, na.rm = TRUE)
}

feature_from_curve <- function(df_curve, n_col = "n_rate", y_col = "GAM") {
  # df_curve: rows for a single ffy_id
  stopifnot(all(c(n_col, y_col) %in% names(df_curve)))
  x <- df_curve[[n_col]]
  y <- df_curve[[y_col]]
  # Plateau (max over grid)
  plateau <- max(y, na.rm = TRUE)
  # Steepness at low-N (25th percentile of N)
  n_q25 <- as.numeric(quantile(x, 0.25, na.rm = TRUE))
  d1 <- finite_diff(x, y)
  steep_low <- median(d1[which.min(abs(x - n_q25)) + (-1:1)], na.rm = TRUE)
  if (!is.finite(steep_low)) steep_low <- median(d1, na.rm = TRUE)
  # Curvature index (larger = more wiggly; invert to "concavity")
  curv <- roughness_ssd(x, y)   # we will use log(1+curv)
  # ΔY over +50 kg N from 25th to 75th percentile anchors
  n_q75 <- as.numeric(quantile(x, 0.75, na.rm = TRUE))
  y_q25 <- approx(x, y, xout = n_q25)$y
  y_q75 <- approx(x, y, xout = n_q75)$y
  deltaY_50 <- (y_q75 - y_q25) / max(1, (n_q75 - n_q25)) * 50
  tibble(
    plateau = plateau,
    steep_low = steep_low,
    curv_log = log1p(curv),
    dY50 = as.numeric(deltaY_50)
  )
}

# Optional: rudimentary SEs via small moving-window smoothing (stability proxy)
feature_se_proxy <- function(df_curve, n_col = "n_rate", y_col = "GAM") {
  # compute features on slightly smoothed variants and take sd
  x <- df_curve[[n_col]]
  y <- df_curve[[y_col]]
  if (length(unique(x)) < 6) {
    return(tibble(plateau_se = NA_real_, steep_low_se = NA_real_,
                  curv_log_se = NA_real_, dY50_se = NA_real_))
  }
  variants <- list(
    y, stats::filter(y, rep(1/3,3), sides = 2),
    stats::filter(y, rep(1/5,5), sides = 2)
  )
  feats <- purrr::map_dfr(variants, ~ feature_from_curve(
    tibble(!!n_col := x, !!y_col := as.numeric(.x)), n_col, y_col))
  sds <- sapply(feats, sd, na.rm = TRUE)
  tibble(
    plateau_se  = sds[["plateau"]]  %||% NA_real_,
    steep_low_se= sds[["steep_low"]]%||% NA_real_,
    curv_log_se = sds[["curv_log"]] %||% NA_real_,
    dY50_se     = sds[["dY50"]]     %||% NA_real_
  )
}

# Wrapper to compute features (+ SE proxy) per ffy_id
build_feature_table <- function(curves_csv_dir,
                                prefer = c("GAM","RF","XGB"),
                                n_col = "n_rate") {
  files <- list.files(curves_csv_dir, pattern = "n_curves\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No n_curves.csv files found under: ", curves_csv_dir)
  out <- purrr::map_dfr(files, function(fp) {
    df <- readr::read_csv(fp, show_col_types = FALSE)
    ffy <- basename(dirname(fp))
    # pick first available model column in priority order
    y_col <- (prefer[prefer %in% names(df)])[1]
    if (is.na(y_col)) return(NULL)
    feat  <- feature_from_curve(df, n_col = n_col, y_col = y_col)
    seprx <- feature_se_proxy  (df, n_col = n_col, y_col = y_col)
    bind_cols(tibble(ffy_id = ffy, model_used = y_col), feat, seprx)
  })
  out
}

# ------------------------------------------------------------
# 2) Build "clue" predictors from stacked analysis table
#    Input: subplot-level table; we aggregate to field-year means
# ------------------------------------------------------------

default_clue_spec <- function() {
  list(
    N_side = c(
      "heavy_rain_days_postN_d15",
      "precip_postN_d15",
      "precip_N_to_yield"
    ),
    W_side = c(
      "gdd_N_to_yield",
      "edd_postN_d15",
      "dry_days_N_to_yield",
      "max_dry_spell_N_to_yield"
    ),
    soil_topo = c("sand","water_storage","slope","TPI")
  )
}

# interactions we pre-approve to avoid a zoo
default_interactions <- function() {
  list(
    c("sand","heavy_rain_days_postN_d15"),   # N loss on sand + heavy rain
    c("sand","precip_postN_d15"),
    c("water_storage","precip_postN_d15"),   # buffering
    c("slope","heavy_rain_days_postN_d15")   # runoff proxy
  )
}

zscore <- function(x) as.numeric(scale(x))

build_clue_table <- function(stacked_df, clue_spec = default_clue_spec(),
                             interactions = default_interactions()) {
  # aggregate to field-year
  need <- unique(c(unlist(clue_spec), unlist(interactions)))
  stopifnot(all(c("ffy_id", need) %in% names(stacked_df)))
  agg <- stacked_df %>%
    group_by(ffy_id) %>%
    summarise(across(all_of(need), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  # z-score
  agg <- agg %>% mutate(across(all_of(need), zscore))

  # build interactions
  for (pair in interactions) {
    nm <- paste0(pair[1], "_x_", pair[2])
    if (all(pair %in% names(agg))) {
      agg[[nm]] <- agg[[pair[1]]] * agg[[pair[2]]]
    }
  }

  # label groups for partitioned effects
  N_side <- intersect(clue_spec$N_side, names(agg))
  W_side <- intersect(clue_spec$W_side, names(agg))
  soil_topo <- intersect(clue_spec$soil_topo, names(agg))
  ints <- setdiff(names(agg), c("ffy_id", N_side, W_side, soil_topo))

  list(
    table = agg,
    groups = list(N_side = N_side, W_side = W_side, soil_topo = soil_topo, interactions = ints)
  )
}

# ------------------------------------------------------------
# 3) Fit sparse Bayesian models (one per feature) with brms
#    We use horseshoe prior; allow sign-anchored priors on a few coefficients
# ------------------------------------------------------------

make_brms_data <- function(features_tbl, clues_list) {
  df <- clues_list$table %>% left_join(features_tbl, by = "ffy_id")
  df
}

# Prior helpers
# --- version-safe horseshoe prior for brms ---
# --- version-safe horseshoe prior for brms (no do.call in stored call) ---
hs_prior <- function(tau0 = 0.2, par_ratio = 0.2, use_slab = TRUE) {
  # which arguments does this brms build support?
  hs_formals <- names(formals(brms::horseshoe))
  # build a call like horseshoe(par_ratio=..., scale_global=..., slab_df=..., slab_scale=...)
  call_elems <- list(as.name("horseshoe"))
  if ("par_ratio"    %in% hs_formals) call_elems$par_ratio    <- par_ratio
  if ("scale_global" %in% hs_formals) call_elems$scale_global <- tau0
  if (use_slab && "slab_df"    %in% hs_formals) call_elems$slab_df    <- 4
  if (use_slab && "slab_scale" %in% hs_formals) call_elems$slab_scale <- 1
  hs_call <- as.call(call_elems)
  # evaluate the call inside brms namespace so the object’s recorded call is horseshoe(...)
  hs_obj <- eval(hs_call, envir = asNamespace("brms"))
  brms::prior(hs_obj, class = "b")
}

# (unchanged) weakly-informative extras
default_extras <- function() {
  c(
    brms::prior(normal(0, 5), class = "Intercept"),
    brms::prior(exponential(1), class = "sigma")
  )
}

# (unchanged) model wrapper
fit_feature_model <- function(df, y_name, x_names,
                              family = gaussian(),
                              core_priors = hs_prior(tau0 = 0.2),
                              extra_priors = NULL,
                              iter = 4000, chains = 4, seed = 123,
                              adapt_delta = 0.9, max_treedepth = 12) {
  fml <- stats::as.formula(paste(y_name, "~ 1 +", paste(x_names, collapse = " + ")))
  pri <- c(core_priors, default_extras())
  if (!is.null(extra_priors) && length(extra_priors)) pri <- c(pri, extra_priors)
  brms::brm(
    formula = fml,
    data = df,
    family = family,
    prior = pri,
    iter = iter, chains = chains, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 0
  )
}


# weakly-informative extras
default_extras <- function() {
  c(
    brms::prior(normal(0, 5), class = "Intercept"),
    brms::prior(exponential(1), class = "sigma")
  )
}


sign_anchor_priors <- function() {
  # Example: sand × heavy rain postN -> negative effect on steepness
  list(
    steep_low = c(prior(normal(-0.2, 0.2), class = "b",
                        coef = "sand_x_heavy_rain_days_postN_d15")),
    curv_log  = c(),
    plateau   = c(),
    dY50      = c()
  )
}

fit_feature_model <- function(df, y_name, x_names,
                              family = gaussian(),
                              core_priors = hs_prior(tau0 = 0.2),
                              extra_priors = NULL,
                              iter = 4000, chains = 4, seed = 123,
                              adapt_delta = 0.9, max_treedepth = 12) {
  fml <- stats::as.formula(paste(y_name, "~ 1 +", paste(x_names, collapse = " + ")))
  pri <- c(core_priors, default_extras())
  if (!is.null(extra_priors) && length(extra_priors)) pri <- c(pri, extra_priors)
  brms::brm(
    formula = fml,
    data = df,
    family = family,
    prior = pri,
    iter = iter, chains = chains, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = 0
  )
}


# ------------------------------------------------------------
# 4) Driver scores (N vs Water) per field-year with uncertainty
#    Given posterior draws of betas and the clue matrix Z
# ------------------------------------------------------------

posterior_driver_scores <- function(fit, df, group_lists, y_name) {
  # Extract draws for coefficients
  dr <- as_draws_df(fit)
  # Map column names b_var to var
  bn <- grep("^b_", names(dr), value = TRUE)
  coef_names <- sub("^b_", "", bn)
  B <- as.matrix(dr[, bn, drop = FALSE])
  # Predictor design
  X <- as.matrix(df[, coef_names, drop = FALSE])
  # Contributions (draws x obs): B %*% t(X) => transpose result
  contrib <- X %*% t(B)  # obs x draws
  # Split per group
  out <- list()
  for (g in c("N_side","W_side")) {
    vars_g <- intersect(group_lists[[g]], coef_names)
    if (!length(vars_g)) next
    idx <- match(vars_g, coef_names)
    Xg <- as.matrix(df[, vars_g, drop = FALSE])
    Bg <- t(B[, idx, drop = FALSE])  # vars x draws
    Cg <- Xg %*% Bg                   # obs x draws
    out[[g]] <- Cg
  }
  # Summaries
  summ <- tibble(ffy_id = df$ffy_id,
                 feature = y_name,
                 S_N_mean = rowMeans(out$N_side, na.rm = TRUE) %||% NA_real_,
                 S_W_mean = rowMeans(out$W_side, na.rm = TRUE) %||% NA_real_,
                 S_N_lo = rowQuantiles(out$N_side, probs = 0.1, na.rm = TRUE) %||% NA_real_,
                 S_N_hi = rowQuantiles(out$N_side, probs = 0.9, na.rm = TRUE) %||% NA_real_,
                 S_W_lo = rowQuantiles(out$W_side, probs = 0.1, na.rm = TRUE) %||% NA_real_,
                 S_W_hi = rowQuantiles(out$W_side, probs = 0.9, na.rm = TRUE) %||% NA_real_)
  list(scores = summ, raw = out)
}

# Convenience plot
plot_driver_quadrant <- function(scores_tbl, title = "Driver scores") {
  ggplot(scores_tbl, aes(x = S_N_mean, y = S_W_mean)) +
    geom_hline(yintercept = 0, linetype = 3) +
    geom_vline(xintercept = 0, linetype = 3) +
    geom_point(alpha = 0.7) +
    geom_errorbar(aes(ymin = S_W_lo, ymax = S_W_hi), width = 0) +
    geom_errorbarh(aes(xmin = S_N_lo, xmax = S_N_hi), height = 0) +
    labs(x = "N-side contribution (posterior mean, 80% PI)",
         y = "Water-side contribution (posterior mean, 80% PI)",
         title = title) +
    theme_minimal(base_size = 12)
}
