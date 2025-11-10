# ============================================================
# functions_for_processing.R — Modular data-processing utilities for OFPE
# Goal: produce an analysis-ready table per field-year that contains
#        experimental vars (yield, n_rate, s_rate, IDs),
#        SSURGO soils (area-weighted), Topography/DEM (weighted),
#        timing-specific weather features (Daymet), and geometry (sf, EPSG:4326).
# ============================================================

suppressPackageStartupMessages({
  # core
  library(dplyr); library(tidyr); library(purrr); library(data.table); library(stringr)
  library(readr); library(lubridate); library(here)
  # spatial
  library(sf); library(terra); library(raster)
  # soils
  library(soilDB)
  # elevation
  library(elevatr)
  # weather
  library(daymetr)
})

# ----------------------------
# General utilities
# ----------------------------
log_message <- function(...) cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", ..., "\n")

st_set_4326 <- function(x) {
  if (is.na(sf::st_crs(x))) {
    x <- sf::st_set_crs(x, 4326)
    message("⚠️ set missing CRS to EPSG:4326")
  }
  x
}

ensure_obs_id <- function(exp_sf) {
  if (!"obs_id" %in% names(exp_sf)) exp_sf$obs_id <- seq_len(nrow(exp_sf))
  exp_sf
}

safe_union_cols <- function(dfs) {
  all_cols <- Reduce(union, lapply(dfs, names))
  lapply(dfs, function(d) {
    miss <- setdiff(all_cols, names(d))
    if (length(miss)) for (m in miss) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  })
}

save_processed <- function(object, fname, out_dir = here("Data","Processed","Analysis_ready")) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file.path(out_dir, fname))
  invisible(TRUE)
}

read_processed <- function(fname, in_dir = here("Data","Processed","Analysis_ready")) {
  readRDS(file.path(in_dir, fname))
}

# ----------------------------
# I/O helpers (experimental + boundary)
# ----------------------------
read_trial_lists <- function(base_dir = here("Data","Raw")) {
  data_ids <- list.files(file.path(base_dir, "exp_tb_data")) |>
    stringr::str_subset("_tb\\.rds$") |>
    stringr::str_remove("_tb\\.rds$")
  bdry_ids <- list.files(file.path(base_dir, "exp_bdry_data")) |>
    stringr::str_subset("_bdry\\.rds$") |>
    stringr::str_remove("_bdry\\.rds$")
  list(data_ids = data_ids, bdry_ids = bdry_ids)
}

read_exp_sf <- function(ffy_id, base_dir = here("Data","Raw")) {
  tb_path <- file.path(base_dir, "exp_tb_data",  paste0(ffy_id, "_tb.rds"))
  stopifnot(file.exists(tb_path))
  exp_tb <- readRDS(tb_path)

  # Handle geometry column robustly
  if ("geom" %in% names(exp_tb)) {
    g <- exp_tb$geom; exp_tb$geom <- NULL
    exp_sf <- sf::st_sf(exp_tb, geometry = g)
  } else if ("geometry" %in% names(exp_tb) && inherits(exp_tb$geometry, "sfc")) {
    g <- exp_tb$geometry; exp_tb$geometry <- NULL
    exp_sf <- sf::st_sf(exp_tb, geometry = g)
  } else {
    stop("Experimental table lacks an sf geometry column ('geom' or 'geometry').")
  }

  exp_sf <- st_set_4326(exp_sf)
  ensure_obs_id(exp_sf)
}

read_boundary <- function(ffy_id, base_dir = here("Data","Raw")) {
  bd_path <- file.path(base_dir, "exp_bdry_data", paste0(ffy_id, "_bdry.rds"))
  stopifnot(file.exists(bd_path))
  bd <- readRDS(bd_path)
  st_set_4326(bd)
}

# ----------------------------
# SOIL PROCESS
# ----------------------------
# Fetch SSURGO mapunit polygons within boundary and attach selected properties
get_ssurgo_for_boundary <- function(boundary_sf,
                                    vars = c("sandtotal_r","silttotal_r","claytotal_r","awc_r"),
                                    top_depth = 0, bottom_depth = 150) {
  boundary_sp <- as(boundary_sf, "Spatial")
  res <- SDA_spatialQuery(boundary_sp, what = "mupolygon", db = "SSURGO", geomIntersection = TRUE)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  ssurgo_geom <- sf::st_as_sf(res) |> st_set_4326()

  mukeydata <- tryCatch({
    get_SDA_property(property = vars, method = "Weighted Average",
                     mukeys = ssurgo_geom$mukey,
                     top_depth = top_depth, bottom_depth = bottom_depth)
  }, error = function(e) NULL)

  if (is.null(mukeydata) || nrow(mukeydata) == 0) return(NULL)

  ssurgo_geom |>
    dplyr::left_join(mukeydata, by = "mukey") |>
    dplyr::select(mukey, musym, muname,
                  sand = sandtotal_r, silt = silttotal_r, clay = claytotal_r,
                  water_storage = awc_r)
}

# Area-weight SSURGO properties for each experimental polygon
area_weight_soils_to_exp <- function(exp_sf, ssurgo_sf) {
  n <- nrow(exp_sf)
  if (is.null(ssurgo_sf) || nrow(ssurgo_sf) == 0) {
    return(data.frame(
      clay = NA_real_, sand = NA_real_, silt = NA_real_, water_storage = NA_real_
    )[rep(1, n), , drop = FALSE])
  }
  exp_sf <- st_transform(exp_sf, st_crs(ssurgo_sf))
  inter <- suppressWarnings(sf::st_intersection(dplyr::select(exp_sf, obs_id), ssurgo_sf))
  if (nrow(inter) == 0) {
    return(data.frame(
      clay = NA_real_, sand = NA_real_, silt = NA_real_, water_storage = NA_real_
    )[rep(1, n), , drop = FALSE])
  }
  inter$area <- as.numeric(sf::st_area(inter))
  dt <- inter |>
    sf::st_drop_geometry() |>
    data.table::as.data.table()
  if (!"obs_id" %in% names(dt)) dt[, obs_id := .I]
  dt[, area_pct := area / sum(area), by = obs_id]
  agg <- dt[, .(
    clay = weighted.mean(clay, w = area_pct, na.rm = TRUE),
    sand = weighted.mean(sand, w = area_pct, na.rm = TRUE),
    silt = weighted.mean(silt, w = area_pct, na.rm = TRUE),
    water_storage = weighted.mean(water_storage, w = area_pct, na.rm = TRUE)
  ), by = obs_id][order(obs_id)]
  idx <- match(exp_sf$obs_id, agg$obs_id)
  out <- as.data.frame(agg[idx, .(clay, sand, silt, water_storage)])
  if (nrow(out) != n) {
    out_full <- data.frame(clay = NA_real_, sand = NA_real_, silt = NA_real_, water_storage = NA_real_)[rep(1, n), , drop = FALSE]
    if (!is.null(idx)) out_full[!is.na(idx), ] <- out[!is.na(idx), ]
    out <- out_full
  }
  out
}

# ----------------------------
# TOPOGRAPHY PROCESS
# ----------------------------
# Build elev/slope/aspect/TPI SpatRaster for boundary bbox
get_topography_stack <- function(boundary_sf, z = 14) {
  bb_sf <- boundary_sf |> st_transform(4326) |> st_bbox() |> st_as_sfc() |> st_as_sf()
  r_rast <- elevatr::get_elev_raster(bb_sf, clip = "locations", z = z)
  names(r_rast) <- "elev"
  r <- terra::rast(r_rast)
  slope  <- terra::terrain(r$elev, v = "slope",  unit = "degrees")
  aspect <- terra::terrain(r$elev, v = "aspect", unit = "degrees")
  tpi    <- terra::terrain(r$elev, v = "TPI")
  topo <- c(r$elev, slope, aspect, tpi)
  names(topo) <- c("elev","slope","aspect","TPI")
  topo
}

# Weighted extract topo to experimental polygons
weighted_extract_topo <- function(exp_sf, topo_r) {
  # project to raster CRS for extraction
  v <- terra::vect(sf::st_transform(exp_sf, terra::crs(topo_r)))
  vals <- terra::extract(topo_r, v, exact = TRUE, weights = TRUE)
  dt <- data.table::as.data.table(vals)
  keep <- intersect(c("ID","elev","slope","aspect","TPI","weight"), names(dt))
  dt <- dt[, ..keep]
  if (!nrow(dt)) {
    return(data.frame(elev = NA_real_, slope = NA_real_, aspect = NA_real_, TPI = NA_real_)[rep(1, nrow(exp_sf)), ])
  }
  out <- dt[, .(
    elev   = weighted.mean(elev,   w = weight, na.rm = TRUE),
    slope  = weighted.mean(slope,  w = weight, na.rm = TRUE),
    aspect = weighted.mean(aspect, w = weight, na.rm = TRUE),
    TPI    = weighted.mean(TPI,    w = weight, na.rm = TRUE)
  ), by = ID][order(ID)]
  out$ID <- NULL
  as.data.frame(out)
}

# ----------------------------
# WEATHER PROCESS (Daymet + timing features)
# ----------------------------
# Degree-day helpers
gdd_day <- function(tmin, tmax, base = 10, cap = 30) pmax(pmin((tmin + tmax)/2, cap) - base, 0)
edd_day <- function(tmax, thresh = 30) pmax(tmax - thresh, 0)

# Download Daymet for a point and year-range
download_daymet_point <- function(lat, lon, start_year, end_year = start_year) {
  out <- tryCatch({
    daymetr::download_daymet(lat = lat, lon = lon, start = start_year, end = end_year, internal = TRUE)$data
  }, error = function(e) NULL)
  if (is.null(out)) stop("Daymet download failed for lat=", lat, ", lon=", lon)
  data.table::as.data.table(out)
}

compute_postN_features <- function(daily_dt, n_time, y_time) {
  nm <- c(
    "precip_N_to_yield","gdd_N_to_yield","edd_N_to_yield",
    "dry_days_N_to_yield","max_dry_spell_N_to_yield","heavy_rain_days_N_to_yield",
    paste0("precip_postN_d", c(7,15,30)),
    paste0("heavy_rain_days_postN_d", c(7,15,30)),
    paste0("gdd_postN_d", c(7,15,30)),
    paste0("edd_postN_d", c(7,15,30))
  )
  if (is.na(n_time) || is.na(y_time)) return(as.list(setNames(rep(NA_real_, length(nm)), nm)))

  dt <- daily_dt[, .(yday, prcp = `prcp..mm.day.`, tmax = `tmax..deg.c.`, tmin = `tmin..deg.c.`)]
  dt[, gdd := gdd_day(tmin, tmax)]
  dt[, edd := edd_day(tmax)]
  nd <- yday(n_time); yd <- yday(y_time)
  seg <- dt[yday > nd & yday <= yd]
  dry_run <- rle(as.integer(seg$prcp < 1e-6))
  roll_sum <- function(var, k) sum(dt[yday > nd & yday <= (nd + k)][[var]], na.rm = TRUE)

  list(
    precip_N_to_yield        = sum(seg$prcp, na.rm = TRUE),
    gdd_N_to_yield           = sum(seg$gdd,  na.rm = TRUE),
    edd_N_to_yield           = sum(seg$edd,  na.rm = TRUE),
    dry_days_N_to_yield      = sum(seg$prcp < 1e-6, na.rm = TRUE),
    max_dry_spell_N_to_yield = ifelse(any(dry_run$values == 1), max(dry_run$lengths[dry_run$values == 1]), 0),
    heavy_rain_days_N_to_yield = sum(seg$prcp >= 10, na.rm = TRUE),

    precip_postN_d7   = roll_sum("prcp", 7),
    precip_postN_d15  = roll_sum("prcp", 15),
    precip_postN_d30  = roll_sum("prcp", 30),
    heavy_rain_days_postN_d7  = sum(dt[yday > nd & yday <= (nd + 7)]$prcp  >= 10, na.rm = TRUE),
    heavy_rain_days_postN_d15 = sum(dt[yday > nd & yday <= (nd + 15)]$prcp >= 10, na.rm = TRUE),
    heavy_rain_days_postN_d30 = sum(dt[yday > nd & yday <= (nd + 30)]$prcp >= 10, na.rm = TRUE),
    gdd_postN_d7  = sum(dt[yday > nd & yday <= (nd + 7)]$gdd, na.rm = TRUE),
    gdd_postN_d15 = sum(dt[yday > nd & yday <= (nd + 15)]$gdd, na.rm = TRUE),
    gdd_postN_d30 = sum(dt[yday > nd & yday <= (nd + 30)]$gdd, na.rm = TRUE),
    edd_postN_d7  = sum(dt[yday > nd & yday <= (nd + 7)]$edd, na.rm = TRUE),
    edd_postN_d15 = sum(dt[yday > nd & yday <= (nd + 15)]$edd, na.rm = TRUE),
    edd_postN_d30 = sum(dt[yday > nd & yday <= (nd + 30)]$edd, na.rm = TRUE)
  )
}

compute_stage_features <- function(daily_dt, s_time, n_time, y_time) {
  nm <- as.vector(outer(c("precip","dry_days","max_dry_spell","heavy_rain_days","gdd","edd"), paste0("_S",1:4), paste0))
  if (any(is.na(c(s_time, n_time, y_time)))) {
    out <- as.list(setNames(rep(NA_real_, length(nm)), nm))
    out$n_app_stage <- NA_character_
    out$precip_15_post_N <- NA_real_
    out$days_to_N_app <- NA_real_
    return(out)
  }
  dt <- daily_dt[, .(yday, prcp = `prcp..mm.day.`, tmax = `tmax..deg.c.`, tmin = `tmin..deg.c.`)]
  dt[, gdd := gdd_day(tmin, tmax)]
  dt[, edd := edd_day(tmax)]
  ds <- yday(s_time); dn <- yday(n_time); dy <- yday(y_time)
  seg <- list(
    S1 = dt[yday >= ds & yday <  (ds + 15)],
    S2 = dt[yday >= (ds + 15) & yday < dn],
    S3 = dt[yday >= dn & yday < (dn + 15)],
    S4 = dt[yday >= (dn + 15) & yday <= dy]
  )
  agg <- function(x) {
    dr <- rle(as.integer(x$prcp < 1e-6))
    list(
      precip = sum(x$prcp, na.rm = TRUE),
      dry_days = sum(x$prcp < 1e-6, na.rm = TRUE),
      max_dry_spell = ifelse(any(dr$values == 1), max(dr$lengths[dr$values == 1]), 0),
      heavy_rain_days = sum(x$prcp >= 10, na.rm = TRUE),
      gdd = sum(x$gdd, na.rm = TRUE),
      edd = sum(x$edd, na.rm = TRUE)
    )
  }
  A <- lapply(seg, agg)
  out <- list(
    precip_S1 = A$S1$precip, dry_days_S1 = A$S1$dry_days, max_dry_spell_S1 = A$S1$max_dry_spell, heavy_rain_days_S1 = A$S1$heavy_rain_days, gdd_S1 = A$S1$gdd, edd_S1 = A$S1$edd,
    precip_S2 = A$S2$precip, dry_days_S2 = A$S2$dry_days, max_dry_spell_S2 = A$S2$max_dry_spell, heavy_rain_days_S2 = A$S2$heavy_rain_days, gdd_S2 = A$S2$gdd, edd_S2 = A$S2$edd,
    precip_S3 = A$S3$precip, dry_days_S3 = A$S3$dry_days, max_dry_spell_S3 = A$S3$max_dry_spell, heavy_rain_days_S3 = A$S3$heavy_rain_days, gdd_S3 = A$S3$gdd, edd_S3 = A$S3$edd,
    precip_S4 = A$S4$precip, dry_days_S4 = A$S4$dry_days, max_dry_spell_S4 = A$S4$max_dry_spell, heavy_rain_days_S4 = A$S4$heavy_rain_days, gdd_S4 = A$S4$gdd, edd_S4 = A$S4$edd
  )
  out$n_app_stage      <- dplyr::case_when(dn < (ds + 15) ~ "S1", TRUE ~ "S3_or_S4")
  out$precip_15_post_N <- sum(dt[yday > dn & yday <= (dn + 15)]$prcp, na.rm = TRUE)
  out$days_to_N_app    <- dn - ds
  out
}

build_weather_features <- function(ffy_id, boundary_sf, s_time, n_time, y_time,
                                   start_pad_years = 30, apr_to_sep_only = TRUE) {
  # centroid point
  cent <- boundary_sf |> st_centroid() |> st_coordinates()
  lat <- cent[1, "Y"]; lon <- cent[1, "X"]
  ffy_year <- as.numeric(sub(".*_(\\d{4})$", "\\1", ffy_id))

  dm_all <- download_daymet_point(lat, lon, ffy_year - start_pad_years, ffy_year)
  dm_year <- dm_all |> dplyr::filter(year == ffy_year) |>
    dplyr::rename(prcp = `prcp..mm.day.`, tmax = `tmax..deg.c.`, tmin = `tmin..deg.c.`) |>
    dplyr::mutate(gdd = gdd_day(tmin,tmax), edd = edd_day(tmax)) |>
    dplyr::select(prcp, tmax, tmin, yday, gdd, edd)

  if (apr_to_sep_only) {
    day_to_month <- function(yday) {
      bounds <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366)
      labs   <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
      labs[findInterval(yday, bounds)]
    }
    dm_year <- dm_year |> mutate(month = day_to_month(yday)) |> filter(month %in% c("Apr","May","Jun","Jul","Aug","Sep"))
  }

  # seasonal totals and climatology means
  year_sum <- dm_year |> summarise(
    prcp_t = round(sum(prcp, na.rm=TRUE),1),
    gdd_t  = round(sum(gdd,  na.rm=TRUE),1),
    edd_t  = round(sum(edd,  na.rm=TRUE),1)
  )

  avgN <- function(n) dm_all |>
    dplyr::filter(year >= (ffy_year - (n-1)) & year <= ffy_year) |>
    dplyr::rename(prcp = `prcp..mm.day.`, tmax = `tmax..deg.c.`, tmin = `tmin..deg.c.`) |>
    dplyr::mutate(gdd = gdd_day(tmin,tmax), edd = edd_day(tmax)) |>
    summarise(prcp = round(mean(prcp, na.rm=TRUE),1), gdd = round(mean(gdd, na.rm=TRUE),1), edd = round(mean(edd, na.rm=TRUE),1))

  mean5  <- avgN(5);  names(mean5)  <- paste0(c("prcp","gdd","edd"),"_5")
  mean30 <- avgN(30); names(mean30) <- paste0(c("prcp","gdd","edd"),"_30")

  postN  <- compute_postN_features(as.data.table(dm_all[year == ffy_year]), n_time, y_time)
  stages <- compute_stage_features(as.data.table(dm_all[year == ffy_year]), s_time, n_time, y_time)

  c(
    list(ffy_id = ffy_id, s_time = s_time, n_time = n_time, yield_time = y_time),
    as.list(year_sum), as.list(mean5), as.list(mean30), postN, stages
  ) |> as.data.frame()
}

# ----------------------------
# ASSEMBLY per field-year
# ----------------------------
# date_manifest must contain ffy_id, s_time, n_time, yield_time
build_field_dataset <- function(ffy_id,
                                base_dir = here("Data","Raw"),
                                date_manifest = NULL,
                                cache_dir = here("Data","Processed","Analysis_ready"),
                                elev_z = 14,
                                crs_out = 4326) {
  log_message("• Building field:", ffy_id)

  exp_sf <- read_exp_sf(ffy_id, base_dir)
  boundary_sf <- read_boundary(ffy_id, base_dir)

  # --- SOILS ---
  ssurgo_sf <- tryCatch(get_ssurgo_for_boundary(boundary_sf), error = function(e) NULL)
  soils_df  <- tryCatch(area_weight_soils_to_exp(exp_sf, ssurgo_sf), error = function(e) NULL)
  if (is.null(soils_df)) soils_df <- data.frame(clay=NA_real_, sand=NA_real_, silt=NA_real_, water_storage=NA_real_)[rep(1, nrow(exp_sf)), ]

  # --- TOPO ---
  topo_r <- tryCatch(get_topography_stack(boundary_sf, z = elev_z), error = function(e) NULL)
  topo_df <- if (!is.null(topo_r)) tryCatch(weighted_extract_topo(exp_sf, topo_r), error = function(e) NULL) else NULL
  if (is.null(topo_df) || nrow(topo_df) == 0) topo_df <- data.frame(elev=NA_real_, slope=NA_real_, aspect=NA_real_, TPI=NA_real_)[rep(1, nrow(exp_sf)), ]

  # --- WEATHER ---
  if (is.null(date_manifest)) stop("date_manifest (s_time, n_time, yield_time) is required for weather features.")
  di <- date_manifest |> dplyr::filter(ffy_id == {{ffy_id}})
  if (!nrow(di)) stop("No date manifest row for ", ffy_id)
  s_time <- di$s_time[[1]]; n_time <- di$n_time[[1]]; y_time <- di$yield_time[[1]]
  weather_df <- tryCatch(build_weather_features(ffy_id, boundary_sf, s_time, n_time, y_time), error = function(e) NULL)
  if (is.null(weather_df)) weather_df <- data.frame(ffy_id = ffy_id)

  # --- COMBINE ---
  exp_tbl <- sf::st_drop_geometry(exp_sf)
  combined <- cbind(exp_tbl, topo_df, soils_df)
  out_sf <- sf::st_sf(combined, geometry = sf::st_geometry(exp_sf)) |> sf::st_transform(crs_out)

  # attach field-level weather by ffy_id (broadcast to rows)
  out_sf$ffy_id <- ffy_id
  weather_df$ffy_id <- ffy_id
  out_sf <- out_sf |> dplyr::left_join(weather_df, by = "ffy_id")

  # cache per-field
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(out_sf, file.path(cache_dir, paste0(ffy_id, "_merged_data.rds")))

  out_sf
}

# ----------------------------
# STACK across fields
# ----------------------------
stack_fields <- function(sf_list) {
  drop_geo <- lapply(sf_list, sf::st_drop_geometry)
  aligned  <- safe_union_cols(drop_geo)
  stacked  <- do.call(rbind, aligned)
  # geometries are not stacked here; return stacked data.frame
  stacked
}
