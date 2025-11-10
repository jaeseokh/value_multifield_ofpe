# ============================================================
# modeling_bayes.R -- Bayesian hierarchical model + updates
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(readr); library(ggplot2)
  library(brms); library(posterior)
})

# ---- base formula (edit as you like) ----
bayes_formula_default <- function() {
  # varying intercepts & N_rate slopes by field
  bf(yield ~ 1 + N_rate + gdd_early + rain_early + (1 + N_rate | field_id))
}

bayes_priors_default <- function() {
  c(
    prior(normal(0, 10), class = "Intercept"),
    prior(normal(0, 5),  class = "b"),
    prior(exponential(1), class = "sd"),
    prior(exponential(1), class = "sigma")
  )
}

# ---- fit model ----
fit_bayes_hier <- function(df, formula = bayes_formula_default(),
                           priors  = bayes_priors_default(),
                           family  = gaussian(),
                           iter = 4000, adapt_delta = 0.95, seed = 123) {
  set.seed(seed)
  brm(
    formula = formula,
    data    = df,
    family  = family,
    prior   = priors,
    cores   = max(1, parallel::detectCores() - 1),
    iter    = iter,
    control = list(adapt_delta = adapt_delta),
    silent = 2, refresh = 0
  )
}

# ---- summaries & export ----
extract_posterior_summary <- function(fit, path_csv = NULL, pattern = c("^b_","^sd_")) {
  summ <- posterior_summary(fit, pars = pattern) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("parameter")
  if (!is.null(path_csv)) readr::write_csv(summ, path_csv)
  summ
}

# ---- counterfactual curves over N_rate ----
posterior_curves_over_N <- function(fit, df, n_grid = 50) {
  rng  <- range(df$N_rate, na.rm = TRUE)
  grid <- tibble::tibble(N_rate = seq(rng[1], rng[2], length.out = n_grid))
  # pick typical values for weather covariates
  grid$gdd_early  <- median(df$gdd_early,  na.rm = TRUE)
  grid$rain_early <- median(df$rain_early, na.rm = TRUE)
  # predict at population-level (marginal over random effects)
  ep <- fitted(fit, newdata = grid, re_formula = NA, summary = TRUE)
  cbind(grid, as.data.frame(ep))
}

plot_bayes_curves <- function(fit, df, save_path = NULL) {
  ce <- posterior_curves_over_N(fit, df)
  p <- ggplot(ce, aes(N_rate, Estimate)) +
    geom_line(linewidth = 0.9) +
    geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), alpha = 0.2) +
    labs(y = "Predicted yield", title = "Posterior mean yield–N curve") +
    theme_minimal()
  if (!is.null(save_path)) ggsave(save_path, p, width = 6, height = 4)
  invisible(p)
}

# ---- sequential updating (simple prior-from-posterior) ----
build_informative_priors <- function(fit_prior, shrink = 1.5) {
  draws <- as_draws_df(fit_prior)
  pri <- list()
  # Intercept
  if ("b_Intercept" %in% names(draws)) {
    pri <- c(pri, prior(normal(mean(draws$b_Intercept), sd(draws$b_Intercept)/shrink), class = "Intercept"))
  }
  # Fixed effects b_*
  b_cols <- grep("^b_", names(draws), value = TRUE)
  for (bc in setdiff(b_cols, "b_Intercept")) {
    nm <- sub("^b_", "", bc)
    pri <- c(pri, prior(normal(mean(draws[[bc]]), sd(draws[[bc]])/shrink), class = "b", coef = nm))
  }
  # Random effects scales, residual sigma
  pri <- c(pri, prior(exponential(1), class = "sd"), prior(exponential(1), class = "sigma"))
  pri
}

update_bayes_with_prior <- function(fit_prior, new_data,
                                    formula = bayes_formula_default(),
                                    family = gaussian(),
                                    iter = 4000, adapt_delta = 0.98, seed = 123) {
  set.seed(seed)
  pri <- build_informative_priors(fit_prior)
  brm(
    formula = formula,
    data    = new_data,
    family  = family,
    prior   = pri,
    cores   = max(1, parallel::detectCores() - 1),
    iter    = iter,
    control = list(adapt_delta = adapt_delta),
    silent = 2, refresh = 0
  )
}

# ---- disk IO helpers for stage paths ----
save_bayes_artifacts <- function(fit, out_tbl_dir, out_fig_dir, df_for_plot) {
  dir.create(out_tbl_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_fig_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(fit, file.path(out_tbl_dir, "bayes_fit.rds"))
  extract_posterior_summary(fit, file.path(out_tbl_dir, "bayes_coef_summary.csv"))
  plot_bayes_curves(fit, df_for_plot, file.path(out_fig_dir, "bayes_curves.png"))
}
