# Code/src/functions_apsim_compare.R
# ----------------------------------------------------------
# Helpers for APSIM–OFPE comparison (prior vs OFPE likelihood)
#   - mukey intersection for OFPE fields
#   - Region-specific N support (North / Central / South)
#   - APSIM 30-year prior curves on region-consistent N range
#   - OFPE empirical yield–N curves (per ffy_id, no clustering)
#   - N-level alignment and APSIM–OFPE discrepancy
# ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(data.table)
  library(sf)
  library(ggplot2)
  library(ggnewscale)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ----------------------------------------------------------
# 1) mukey intersection for a single OFPE field-year
#    (requires functions_for_processing.R for read_boundary()
#     and get_ssurgo_for_boundary())
# ----------------------------------------------------------

compute_mukey_for_merged <- function(
    ffy_id,
    base_dir_raw    = here::here("Data","Raw"),
    base_dir_merged = here::here("Data","Processed","Analysis_ready"),
    top_depth = 0,
    bottom_depth = 150
) {

  message("→ mukey for ", ffy_id)

  out_path <- file.path(base_dir_merged, paste0(ffy_id, "_merged_mukey.rds"))

  if (file.exists(out_path)) {
    message("✓ Already exists → skipping: ", out_path)
    merged_mukey_sf <- readRDS(out_path)
    return(list(
      ffy_id          = ffy_id,
      merged_mukey_sf = merged_mukey_sf
    ))
  }

  merged_path <- file.path(base_dir_merged, paste0(ffy_id, "_merged_data.rds"))
  if (!file.exists(merged_path)) {
    warning("Merged file not found for ", ffy_id)
    return(NULL)
  }

  merged_sf <- readRDS(merged_path)

  if (!inherits(merged_sf, "sf")) stop("Merged object is not sf.")
  if (!"obs_id" %in% names(merged_sf)) stop("merged_data has no obs_id.")

  bdry <- read_boundary(ffy_id, base_dir = base_dir_raw)

  ssurgo_sf <- tryCatch(
    get_ssurgo_for_boundary(
      boundary_sf  = bdry,
      vars         = c("sandtotal_r","silttotal_r","claytotal_r","awc_r"),
      top_depth    = top_depth,
      bottom_depth = bottom_depth
    ),
    error = function(e) {
      warning("SSURGO query failed for ", ffy_id, ": ", e$message)
      NULL
    }
  )

  if (is.null(ssurgo_sf) || nrow(ssurgo_sf) == 0) {
    warning("No SSURGO polygons for ", ffy_id, " → saving no-mukey version.")

    merged_mukey_sf <- merged_sf %>%
      mutate(
        mukey_dom  = NA_character_,
        musym_dom  = NA_character_,
        muname_dom = NA_character_
      )

    saveRDS(merged_mukey_sf, out_path)
    return(list(
      ffy_id          = ffy_id,
      merged_mukey_sf = merged_mukey_sf
    ))
  }

  merged_sf2 <- st_transform(merged_sf, st_crs(ssurgo_sf))

  inter <- suppressWarnings(
    st_intersection(
      dplyr::select(merged_sf2, obs_id),
      ssurgo_sf
    )
  )

  if (nrow(inter) == 0) {
    warning("Empty intersection for ", ffy_id, " → saving NAs for mukey.")

    merged_mukey_sf <- merged_sf %>%
      mutate(
        mukey_dom  = NA_character_,
        musym_dom  = NA_character_,
        muname_dom = NA_character_
      )

    saveRDS(merged_mukey_sf, out_path)
    return(list(
      ffy_id          = ffy_id,
      merged_mukey_sf = merged_mukey_sf
    ))
  }

  inter$area <- as.numeric(st_area(inter))

  inter_dt <- inter |>
    st_drop_geometry() |>
    data.table::as.data.table()

  inter_dt[, area_pct := area / sum(area, na.rm = TRUE), by = obs_id]

  obs_dom <- inter_dt[order(obs_id, -area_pct)][, .SD[1], by = obs_id]

  obs_mukey_dt <- tibble::tibble(
    obs_id    = obs_dom$obs_id,
    mukey_dom = as.character(obs_dom$mukey),
    musym_dom = obs_dom$musym,
    muname_dom = obs_dom$muname
  )

  merged_mukey_sf <- merged_sf %>%
    dplyr::left_join(obs_mukey_dt, by = "obs_id")

  saveRDS(merged_mukey_sf, out_path)

  return(list(
    ffy_id          = ffy_id,
    merged_mukey_sf = merged_mukey_sf
  ))
}

# ----------------------------------------------------------
# 2) Region-specific N support for APSIM priors (corn–soy baseline)
#    These correspond to the MRTN-extended bands we discussed:
#      North:   130–200
#      Central: 140–210
#      South:   150–220
# ----------------------------------------------------------

region_n_support_tbl <- function() {
  tibble::tribble(
    ~region,    ~N_min, ~N_max,
    "North",     130,    200,
    "Central",   140,    210,
    "South",     150,    220
  )
}

# ----------------------------------------------------------
# 3) APSIM 30-year prior curves on region-consistent N range
#    Input:
#      - yield_curve_field_dt: APSIM long table (N_fert, Y_corn, region, id_10, ...)
#      - ofpe_il_cells: mapping ffy_id -> region, id_10 (only IL fields)
#    Output:
#      - apsim_prior_curves_il: APSIM prior curves on truncated N range
#        columns: region, id_10, N_lbs_ac, n_years, Y_apsim_mean, Y_apsim_var, Y_apsim_sd
# ----------------------------------------------------------

build_apsim_prior_curves <- function(yield_curve_field_dt,
                                     ofpe_il_cells,
                                     support_tbl = region_n_support_tbl()) {

  # APSIM full curves (30-year) in lbs/ac and bu/ac
  apsim_cell_curves <- yield_curve_field_dt %>%
    mutate(
      N_lbs_ac = N_fert * 0.8921,         # kg/ha -> lb/ac
      Y_bu_ac  = Y_corn / 62.77           # kg/ha -> bu/ac
    ) %>%
    group_by(region, id_10, N_lbs_ac) %>%
    summarise(
      n_years      = n(),
      Y_apsim_mean = mean(Y_bu_ac, na.rm = TRUE),
      Y_apsim_var  = var(Y_bu_ac, na.rm = TRUE),
      .groups      = "drop"
    ) %>%
    mutate(
      Y_apsim_sd = sqrt(pmax(Y_apsim_var, 0))
    )

  # Only APSIM cells that intersect OFPE IL fields
  apsim_cell_curves_il <- apsim_cell_curves %>%
    inner_join(ofpe_il_cells %>% distinct(region, id_10),
               by = c("region","id_10"))

  # Apply region-specific N truncation
  apsim_prior_curves_il <- apsim_cell_curves_il %>%
    left_join(support_tbl, by = "region") %>%
    filter(
      !is.na(N_min),
      N_lbs_ac >= N_min,
      N_lbs_ac <= N_max
    ) %>%
    arrange(region, id_10, N_lbs_ac)

  apsim_prior_curves_il
}

# ----------------------------------------------------------
# 4) OFPE empirical yield–N curves per field-year (no clustering)
#    Reads *_merged_data.rds and averages yield by exact N rate.
# ----------------------------------------------------------

build_ofpe_curve_from_merged <- function(ffy_id,
                                         analysis_dir,
                                         n_col = "n_rate",
                                         y_col = "yield") {

  merged_path <- file.path(analysis_dir, paste0(ffy_id, "_merged_data.rds"))
  if (!file.exists(merged_path)) {
    warning("No merged_data for ", ffy_id)
    return(NULL)
  }

  merged_sf <- readRDS(merged_path)

  if (!inherits(merged_sf, "sf")) {
    stop("merged_data for ", ffy_id, " is not an sf object.")
  }

  df <- sf::st_drop_geometry(merged_sf)

  if (!n_col %in% names(df)) {
    stop("Column ", n_col, " not found in merged data for ", ffy_id)
  }
  if (!y_col %in% names(df)) {
    stop("Column ", y_col, " not found in merged data for ", ffy_id)
  }

  df %>%
    dplyr::filter(
      is.finite(.data[[n_col]]),
      is.finite(.data[[y_col]])
    ) %>%
    dplyr::group_by(!!rlang::sym(n_col)) %>%
    dplyr::summarise(
      Y_ofpe = mean(.data[[y_col]], na.rm = TRUE),
      n_obs  = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::rename(N_fert = !!rlang::sym(n_col)) %>%
    dplyr::mutate(
      ffy_id = ffy_id,
      .before = 1
    ) %>%
    dplyr::arrange(N_fert)
}



# ----------------------------------------------------------
# 5) Align OFPE N levels to APSIM prior N levels
#    (nearest neighbor, with tolerance in lb/ac)
# ----------------------------------------------------------

align_n_levels <- function(ofpe_curve, apsim_curve, tol = 5) {
  ofpe_curve <- ofpe_curve %>% arrange(N_fert)
  apsim_curve <- apsim_curve %>% arrange(N_lbs_ac)

  match_idx <- sapply(ofpe_curve$N_fert, function(n) {
    diffs <- abs(apsim_curve$N_lbs_ac - n)
    if (all(is.na(diffs))) return(NA_integer_)
    j <- which.min(diffs)
    if (diffs[j] > tol) return(NA_integer_)
    j
  })

  apsim_matched <- apsim_curve[match_idx, , drop = FALSE]

  ofpe_curve %>%
    mutate(
      N_apsim      = apsim_matched$N_lbs_ac,
      Y_apsim_mean = apsim_matched$Y_apsim_mean
    )
}

# ----------------------------------------------------------
# 6) Build APSIM–OFPE discrepancy for a single field-year
#    Uses:
#      - apsim_prior_curves: output of build_apsim_prior_curves()
#      - ofpe_il_cells: mapping ffy_id -> region, id_10
# ----------------------------------------------------------

build_discrepancy_for_field_prior <- function(ffy_id,
                                              apsim_prior_curves,
                                              ofpe_il_cells,
                                              analysis_dir,
                                              n_col = "n_rate",
                                              y_col = "yield_bu_ac",
                                              tol = 5) {

  # 1. OFPE curve for this ffy_id
  ofpe_curve <- build_ofpe_curve_from_merged(
    ffy_id       = ffy_id,
    analysis_dir = analysis_dir,
    n_col        = n_col,
    y_col        = y_col
  )

  if (is.null(ofpe_curve) || !nrow(ofpe_curve)) {
    warning("No OFPE curve for ", ffy_id)
    return(NULL)
  }

  # 2. Find its APSIM cell (region + id_10)
  loc_row <- ofpe_il_cells %>%
    filter(ffy_id == !!ffy_id) %>%
    slice(1)

  if (!nrow(loc_row)) {
    warning("No APSIM cell mapping for ", ffy_id)
    return(NULL)
  }

  apsim_curve <- apsim_prior_curves %>%
    filter(region == loc_row$region,
           id_10  == loc_row$id_10)

  if (!nrow(apsim_curve)) {
    warning("No APSIM prior curve rows for ", ffy_id)
    return(NULL)
  }

  # 3. Align N levels and compute discrepancies
  joined <- align_n_levels(
    ofpe_curve,
    apsim_curve,
    tol = tol
  )

  joined %>%
    mutate(
      D_Y   = Y_ofpe - Y_apsim_mean,
      D_N   = N_apsim - N_fert,
      ffy_id = ffy_id,
      region = loc_row$region,
      id_10  = loc_row$id_10,
      .before = 1
    )
}

# ----------------------------------------------------------
# 7) Plot helpers: APSIM prior vs OFPE OFPE curves, by year, with county yield line
#    - combined_panel_df must have:
#         ffy_id, year, farm_field_label, region,
#         N, Y, type ("APSIM"/"OFPE"), Y_sd
#    - region_cols: named vector, e.g. c(North="forestgreen", Central="steelblue", South="purple")
# ----------------------------------------------------------

make_year_plot <- function(combined_panel_df,
                           year_sel,
                           region_cols,
                           z_band = 1.28,
                           county_lines_df = NULL) {

  df_year <- combined_panel_df %>%
    dplyr::filter(year == year_sel)

  # county yield (one value per panel)
  county_lines <- NULL
  if (!is.null(county_lines_df)) {
    county_lines <- county_lines_df %>%
      dplyr::filter(year == year_sel) %>%
      dplyr::distinct(farm_field_label, Y_county_curr)
  }

  ggplot() +
    # APSIM band
    geom_ribbon(
      data = df_year %>% dplyr::filter(type == "APSIM"),
      aes(
        x    = N,
        ymin = Y - z_band * Y_sd,
        ymax = Y + z_band * Y_sd
      ),
      fill  = "#ffd6d6",
      alpha = 0.5
    ) +
    # APSIM mean (red)
    geom_line(
      data = df_year %>% dplyr::filter(type == "APSIM"),
      aes(x = N, y = Y),
      colour   = "red",
      linewidth = 0.7
    ) +
    # county yield (horizontal dashed grey line, per panel)
    {
      if (!is.null(county_lines)) {
        geom_hline(
          data        = county_lines,
          aes(yintercept = Y_county_curr),
          colour      = "grey30",
          linetype    = "dashed",
          linewidth   = 0.5,
          inherit.aes = FALSE
        )
      }
    } +
    # OFPE points (by region)
    geom_point(
      data = df_year %>% dplyr::filter(type == "OFPE"),
      aes(x = N, y = Y, colour = region),
      alpha = 0.35,
      size  = 0.6
    ) +
    scale_colour_manual(
      values = region_cols,
      name   = "Region"
    ) +
    ggnewscale::new_scale_colour() +
    # OFPE smooth (black, economics-style concave GAM)
    geom_smooth(
      data = df_year %>% dplyr::filter(type == "OFPE"),
      aes(x = N, y = Y, group = farm_field_label, colour = "OFPE smooth"),
      method = "gam",
      formula = y ~ s(x, k = 4),
      se      = FALSE,
      linewidth = 0.7
    ) +
    scale_colour_manual(
      name   = "Curve",
      values = c("OFPE smooth" = "black")
    ) +
    facet_wrap(~ farm_field_label, scales = "free_y") +
    labs(
      x = "N rate (lb/ac)",
      y = "Yield (bu/ac)",
      title = paste0(
        "APSIM prior (red, band = mean ± ", z_band, " sd) vs OFPE (points, black smooth), year ",
        year_sel
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "bottom",
      strip.text      = element_text(size = 8)
    )
}


# 8 align_field

align_field <- function(ffy_id_sel, ofpe_df, apsim_df, max_tol_N = 30) {
  df_o <- ofpe_df %>% filter(ffy_id == ffy_id_sel)
  df_a <- apsim_df %>% filter(ffy_id == ffy_id_sel)

  if (nrow(df_o) == 0 || nrow(df_a) == 0) {
    return(tibble())
  }

  purrr::map_dfr(seq_len(nrow(df_o)), function(i) {
    N_obs <- df_o$N[i]
    Y_obs <- df_o$Y[i]

    diffs <- abs(df_a$N - N_obs)
    j     <- which.min(diffs)

    if (length(j) == 0 || diffs[j] > max_tol_N) {
      return(tibble())
    }

    tibble(
      ffy_id  = ffy_id_sel,
      farm    = df_o$farm[i],
      field   = df_o$field[i],
      year    = df_o$year[i],
      region  = df_o$region[i],
      id_10   = df_o$id_10[i],
      N_ofpe  = N_obs,
      Y_ofpe  = Y_obs,
      N_apsim = df_a$N[j],
      Y_apsim = df_a$Y[j]
    )
  })
}



# 9. helper: compute raw + central moments for one field-year
compute_ofpe_moments_field <- function(ffy_id,
                                       analysis_dir,
                                       n_col = "n_rate",
                                       y_col = "yield") {

  merged_path <- file.path(analysis_dir, paste0(ffy_id, "_merged_mukey.rds"))
  if (!file.exists(merged_path)) {
    warning("No merged_data for ", ffy_id)
    return(NULL)
  }

  merged_sf <- readRDS(merged_path)

  if (!inherits(merged_sf, "sf")) {
    stop("merged_data for ", ffy_id, " is not an sf object.")
  }

  df <- sf::st_drop_geometry(merged_sf)

  if (!n_col %in% names(df)) {
    stop("Column ", n_col, " not found in merged data for ", ffy_id)
  }
  if (!y_col %in% names(df)) {
    stop("Column ", y_col, " not found in merged data for ", ffy_id)
  }

  df %>%
    dplyr::filter(
      is.finite(.data[[n_col]]),
      is.finite(.data[[y_col]])
    ) %>%
    dplyr::group_by(N_fert = .data[[n_col]]) %>%
    dplyr::summarise(
      n_obs = dplyr::n(),
      m1    = mean(.data[[y_col]], na.rm = TRUE),
      m2    = mean((.data[[y_col]])^2, na.rm = TRUE),
      m3    = mean((.data[[y_col]])^3, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      mu1   = m1,
      mu2   = pmax(m2 - m1^2, 0),                            # variance
      mu3   = m3 - 3 * m2 * m1 + 2 * m1^3,                  # 3rd central
      sd    = sqrt(mu2),
      skew  = dplyr::if_else(mu2 > 0, mu3 / (mu2^(3/2)), NA_real_),
      ffy_id = ffy_id,
      .before = 1
    ) %>%
    dplyr::arrange(N_fert)
}
