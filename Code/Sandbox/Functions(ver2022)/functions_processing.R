# ---
# SCRIPT: functions_for_process.R
#
# PURPOSE: 
#   This script contains a collection of helper functions used during the 
#   data processing stage (executed in '1.Data_Process.Rmd'). These functions 
#   are designed to handle specific tasks such as fetching, cleaning, and
#   integrating various data sources. The primary goal is to streamline the 
#   preparation of a clean, analysis-ready dataset.
#
# KEY FUNCTIONS:

#   - get_non_exp_data: Retrieves and processes external geospatial data, 
#     including topographic variables from a Digital Elevation Model (DEM) and 
#     soil properties from the SSURGO database.

#   - calc_stage_specific_gdd_edd: Downloads and calculates weather metrics 
#     like Growing Degree Days (GDD) and Extreme Degree Days (EDD) from 
#     Daymet data for specific in-season periods.

#   - convert_N_unit: Standardizes different nitrogen fertilizer application
#     rates and units into a common nitrogen equivalent ('N_equiv').

#   - Geospatial helpers (stars_to_stack, st_set_4326): Utility functions
#     for working with spatial data objects.
#
# AUTHOR: 
#   Jaeseok Hwang
#
# ---

# Revised function to calculate stage-specific GDD and EDD with proper NA handling

compute_stage_features <- function(daily_weather, s_time, n_time, y_time,
                                   heavy_mm = 25.4) {
  # daily_weather columns expected from daymet: yday, prcp..mm.day., tmax..deg.c., tmin..deg.c.
  stopifnot(!is.na(s_time), !is.na(y_time))
  s_doy <- yday(s_time); y_doy <- yday(y_time)
  if (y_doy <= s_doy) {
    # Handle year wrap if needed (rare for corn, but guard anyway)
    # Here we just bail out to NAs
    out <- as.list(rep(NA_real_, 6*4 + 3))
    names(out) <- c(
      as.vector(outer(c("precip","dry_days","max_dry_spell","heavy_rain_days","gdd","edd"),
                      paste0("_S",1:4), paste0)),
      c("n_app_stage","precip_15_post_N","days_to_N_app")
    )
    return(out)
  }

  dw <- as.data.table(daily_weather)[yday >= s_doy & yday <= y_doy]
  if (nrow(dw) == 0L) {
    out <- as.list(rep(NA_real_, 6*4 + 3))
    names(out) <- c(
      as.vector(outer(c("precip","dry_days","max_dry_spell","heavy_rain_days","gdd","edd"),
                      paste0("_S",1:4), paste0)),
      c("n_app_stage","precip_15_post_N","days_to_N_app")
    )
    return(out)
  }

  # Rename and compute daily GDD/EDD
  setnames(dw,
           c("prcp..mm.day.","tmax..deg.c.","tmin..deg.c."),
           c("prcp","tmax","tmin"), skip_absent = TRUE)
  dw[, gdd := pmax(0, ((tmax + tmin)/2) - 10)]
  dw[, edd := pmax(0, tmax - 30)]

  # Split into 4 equal stages by DOY
  N <- y_doy - s_doy + 1
  cuts <- floor(seq(s_doy, y_doy, length.out = 5))
  dw[, stage := fifelse(yday <= cuts[2], 1L,
                 fifelse(yday <= cuts[3], 2L,
                 fifelse(yday <= cuts[4], 3L, 4L)))]


  # thresholds to tune
dry_thresh  <- 0.5   # mm/day to count as "dry"
heavy_mm    <- 25.4  # mm/day to count as "heavy rain"



  # Per-stage aggregations
  stage_agg <- dw[, {
  # logical vector for "dry" within this stage/group
  flags <- (prcp <= dry_thresh) | is.na(prcp)

  # longest consecutive run of TRUEs
  r <- rle(flags)
  max_dry <- if (!length(r$lengths) || !any(r$values)) 0L else max(r$lengths[r$values])

  .(
    precip          = sum(prcp, na.rm = TRUE),
    dry_days        = sum(flags, na.rm = TRUE),
    max_dry_spell   = as.integer(max_dry),
    heavy_rain_days = sum(prcp >= heavy_mm, na.rm = TRUE),
    gdd             = sum(gdd, na.rm = TRUE),
    edd             = sum(edd, na.rm = TRUE)
  )
}, by = stage][order(stage)]

  # Ensure all stages exist
  full <- data.table(stage=1:4)[stage_agg, on="stage"]
  for (col in c("precip","dry_days","max_dry_spell","heavy_rain_days","gdd","edd")) {
    set(full, i=which(is.na(full[[col]])), j=col, value=0)
  }

  # Flatten names: precip_S1, ..., edd_S4
  out <- as.list(unlist(lapply(1:4, function(s) {
    c(precip = full[stage==s, precip],
      dry_days = full[stage==s, dry_days],
      max_dry_spell = full[stage==s, max_dry_spell],
      heavy_rain_days = full[stage==s, heavy_rain_days],
      gdd = full[stage==s, gdd],
      edd = full[stage==s, edd])
  })))
  names(out) <- as.vector(outer(c("precip","dry_days","max_dry_spell","heavy_rain_days","gdd","edd"),
                                paste0("_S",1:4), paste0))

  # Nitrogen post-application metrics (if n_time present)
  n_app_stage <- NA_character_
  precip_15_post_N <- NA_real_
  days_to_N_app <- NA_real_

  if (!is.na(n_time)) {
    n_doy <- yday(n_time)
    if (n_doy >= s_doy && n_doy <= y_doy) {
      # what stage was N applied?
      n_app_stage_int <- if (n_doy <= cuts[2]) 1L else if (n_doy <= cuts[3]) 2L else if (n_doy <= cuts[4]) 3L else 4L
      n_app_stage <- paste0("S", n_app_stage_int)
      # precip 15 days after N
      precip_15_post_N <- dw[yday > n_doy & yday <= (n_doy + 15), sum(prcp, na.rm=TRUE)]
      # days from sowing to N
      days_to_N_app <- n_doy - s_doy
    }
  }

  c(out, list(
    n_app_stage = n_app_stage,
    precip_15_post_N = precip_15_post_N,
    days_to_N_app = days_to_N_app
  ))
}



# ---- Post-N-only weather features (no s_time required) ----
# daily_weather must contain at least: yday, prcp..mm.day., tmax..deg.c., tmin..deg.c.
compute_postN_features <- function(daily_weather, n_time, y_time,
                                   dry_thresh = 0.5,   # mm to call a "dry" day
                                   heavy_mm   = 25.4,  # mm for "heavy rain"
                                   gdd_base   = 10,    # °C
                                   edd_thresh = 30,    # °C
                                   post_windows = c(7, 15, 30),
                                   pre_windows  = c(7, 15, 30)) {
  out_names <- c(
    # whole-window N -> harvest
    "precip_N_to_yield","gdd_N_to_yield","edd_N_to_yield",
    "dry_days_N_to_yield","max_dry_spell_N_to_yield","heavy_rain_days_N_to_yield",
    # rolling post-N windows
    paste0("precip_postN_d", post_windows),
    paste0("heavy_rain_days_postN_d", post_windows),
    paste0("gdd_postN_d", post_windows),
    paste0("edd_postN_d", post_windows),
    # antecedent windows before N
    paste0("precip_preN_d", pre_windows),
    paste0("dry_days_preN_d", pre_windows),
    paste0("heavy_rain_days_preN_d", pre_windows)
  )
  out <- as.list(rep(NA_real_, length(out_names)))
  names(out) <- out_names

  if (is.na(n_time) || is.na(y_time)) return(out)

  n_doy <- yday(n_time); y_doy <- yday(y_time)
  if (y_doy <= n_doy) return(out)  # guard against wrap; clip if you prefer

  dw <- as.data.table(daily_weather)
  setnames(dw,
           c("prcp..mm.day.","tmax..deg.c.","tmin..deg.c."),
           c("prcp","tmax","tmin"), skip_absent = TRUE)
  if (!all(c("yday","prcp","tmax","tmin") %in% names(dw))) return(out)

  dw[, gdd := pmax(0, ((tmax + tmin)/2) - gdd_base)]
  dw[, edd := pmax(0, tmax - edd_thresh)]

  # ---- N -> harvest window aggregates ----
  seg <- dw[yday >= n_doy & yday <= y_doy]
  if (nrow(seg)) {
    flags <- (seg$prcp <= dry_thresh) | is.na(seg$prcp)
    r <- rle(flags)
    max_dry <- if (!length(r$lengths) || !any(r$values)) 0L else max(r$lengths[r$values])

    out$precip_N_to_yield           <- sum(seg$prcp, na.rm = TRUE)
    out$gdd_N_to_yield              <- sum(seg$gdd,  na.rm = TRUE)
    out$edd_N_to_yield              <- sum(seg$edd,  na.rm = TRUE)
    out$dry_days_N_to_yield         <- sum(flags,    na.rm = TRUE)
    out$max_dry_spell_N_to_yield    <- as.numeric(max_dry)
    out$heavy_rain_days_N_to_yield  <- sum(seg$prcp >= heavy_mm, na.rm = TRUE)
  }

  # ---- Fixed post-N windows (d=7/15/30) ----
  for (w in post_windows) {
    win <- dw[yday > n_doy & yday <= (n_doy + w)]
    if (nrow(win)) {
      out[[paste0("precip_postN_d", w)]]            <- sum(win$prcp, na.rm = TRUE)
      out[[paste0("heavy_rain_days_postN_d", w)]]   <- sum(win$prcp >= heavy_mm, na.rm = TRUE)
      out[[paste0("gdd_postN_d", w)]]               <- sum(win$gdd,  na.rm = TRUE)
      out[[paste0("edd_postN_d", w)]]               <- sum(win$edd,  na.rm = TRUE)
    }
  }

  # ---- Antecedent pre-N windows (d=7/15/30) ----
  for (w in pre_windows) {
    win <- dw[yday >= max(1, n_doy - w) & yday < n_doy]
    if (nrow(win)) {
      flags <- (win$prcp <= dry_thresh) | is.na(win$prcp)
      out[[paste0("precip_preN_d", w)]]             <- sum(win$prcp, na.rm = TRUE)
      out[[paste0("dry_days_preN_d", w)]]           <- sum(flags,     na.rm = TRUE)
      out[[paste0("heavy_rain_days_preN_d", w)]]    <- sum(win$prcp >= heavy_mm, na.rm = TRUE)
    }
  }

  out
}




day_to_month <- function(yday) {
  # Define boundaries for each month
  month_boundaries <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 366)
  months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  # Find the month index
  month_index <- findInterval(yday, month_boundaries)
  
  # Return the corresponding month
  return(months[month_index])
}

st_set_4326 <- function(data_sf) {
  if (is.na(st_crs(data_sf))) {
    data_sf <- st_set_crs(data_sf, 4326)
    cat("Warning: valid crs was not set for this data. Check carefully if this has caused any problems below.")
  }

  return(data_sf)
}

convert_N_unit <- function(
  form,
  unit,
  rate,
  reporting_unit,
  conversion_type = "to_n_equiv"
) {
  
  conv_table <- 
  fromJSON(
    here("Data", "Raw", "nitrogen_conversion.json"), 
    flatten = TRUE
  ) %>%
  data.table() %>%
  .[, conv_factor := as.numeric(conv_factor)] %>%
  .[, form_unit := paste(type, unit, sep = "_")] %>%
  as.data.frame()

  if (form == "N_equiv") {
    conv_factor_n <- 1
  } else {
    conv_factor_n <- which(conv_table[, "form_unit"] %in% paste(form, unit, sep = "_")) %>%
      conv_table[., "conv_factor"]
  }

  if (reporting_unit == "metric") {
    conv_factor_n <- conv_factor_n * conv_unit(1, "lbs", "kg") * conv_unit(1, "hectare", "acre")
  }

  if (conversion_type == "to_n_equiv") {
    converted_rate <- (conv_factor_n)*rate
  } else {
    converted_rate <- (1/conv_factor_n)*rate
  }

  return(as.numeric(converted_rate))
}



get_base_rate <- function(input_data, input_type){
  if(input_type %in% c("NH3", "urea", "uan32", "uan28", "1_2_1(36)", "LAN(26)", "MAP", "1_0_0", "1_0_1", "2_3_2(22)",
                       "15_10_6", "3_0_1", "2_3_4(32)", "4_3_4(33)", "5_1_5", "Sp", "N_equiv", "24-0-0-3 UAN","chicken_manure")){
    is_base <- "base" %in% input_data[, strategy]
    
    if (is_base) {
      base_rate <- input_data[strategy == "base", ] %>% 
        rowwise() %>% 
        mutate(
          n_equiv_rate = convert_N_unit(
            form = form, 
            unit = unit, 
            rate = rate, 
            reporting_unit = w_field_data$reporting_unit
          ) 
        ) %>% 
        data.table() %>% 
        .[, sum(n_equiv_rate)]
    } else {
      base_rate <- 0
    }
  }else{
    base_rate = 0
  }
  
  return(base_rate)
}


get_gc_rate <- function(gc_rate, input_type, form, unit, convert, base_rate){
  if((input_type %in% c("NH3", "urea", "uan32", "uan28", "1_2_1(36)", "LAN(26)", "MAP", "1_0_0", "1_0_1", "2_3_2(22)",
                       "15_10_6", "3_0_1", "2_3_4(32)", "4_3_4(33)", "5_1_5", "Sp", "N_equiv", "24-0-0-3 UAN","chicken_manure")) 
     & (convert == TRUE)){
    if (!is.numeric(gc_rate)) {
      Rx_file <- file.path(
        here("Data/Growers", ffy, "Raw"), 
        paste0(gc_rate_n, ".shp")
      )
      
      if (file.exists(Rx_file)){
        #--- if the Rx file exists ---#
        gc_type <- "Rx"
        gc_rate <- Rx_file_n
        
      }
    } else {
      gc_rate <- data.table(gc_rate = gc_rate, 
                            form = form,
                            unit = unit) %>%
        rowwise() %>%
        mutate(gc_rate :=  convert_N_unit(form, unit, gc_rate, reporting_unit) + base_rate) %>%
        ungroup() %>%
        as.data.frame() %>%
        pull("gc_rate")
      
      gc_type <- "uniform"
    }
  }else{
    gc_rate = gc_rate
  }
  
  return(gc_rate)
}



stars_to_stack <- function(stars) {
  stack <- lapply(1:length(stars), function(x) as(stars[x], "Raster")) %>%
    do.call(raster::stack, .)
  return(stack)
}



get_ssurgo_props <- function(field_bound, vars, summarize = FALSE) {

  # Get SSURGO mukeys for polygon intersection
  ssurgo_geom <-
    SDA_spatialQuery(
      field_bound,
      what = 'geom',
      db = 'SSURGO',
      geomIntersection = T
    ) %>%
    st_as_sf() %>%
    mutate(
      area = as.numeric(st_area(.)),
      area_weight = area / sum(area)
    )
  # Get soil properties for each mukey
  mukeydata <-
    get_SDA_property(
      property = vars,
      method = 'Weighted Average',
      mukeys = ssurgo_geom$mukey,
      top_depth = 0,
      bottom_depth = 150
    )
  ssurgo_data <- left_join(ssurgo_geom, mukeydata, by = 'mukey')
  if (summarize == TRUE) {
    ssurgo_data_sum <-
      ssurgo_data %>%
      data.table() %>%
      .[,
        lapply(.SD, weighted.mean, w = area_weight),
        .SDcols = vars
      ]
    return(ssurgo_data_sum)
  } else {
    return(ssurgo_data)
  }
}


get_non_exp_data <- function(ffy_id) {
  # Initialize boundary_file and td_file
  boundary_sf <- readRDS(here("Data","Raw","exp_bdry_data",paste0(ffy_id,"_bdry.rds")))%>%
      st_set_4326() %>%
      st_transform(4326)%>%
      st_bbox() %>%
      st_as_sfc() %>%
      st_as_sf()
  # Topographic Data Retrieval
  # Get elevation raster
  dem <- get_elev_raster(
    boundary_sf,
    clip = "locations",
    z = 14
  )
  names(dem) <- 'elev'
  
  # Calculate terrain attributes
  dem_slope <- terrain(dem$elev, "slope", unit = "degrees")
  dem_aspect <- terrain(dem$elev, "aspect", unit = "degrees")
  dem_rast <- rast(dem)
  dem_curv <- spatialEco::curvature(dem_rast, type = "mcnab")
  names(dem_curv) <- "curv"
  dem_tpi <- terrain(dem, "tpi", unit = "m")
  
  # Stack all the rasters together
  topo_dem <- stack(dem, dem_slope, dem_aspect, dem_curv, dem_tpi) %>%
    st_as_stars() %>%
    split(3)

  # SSURGO Data Retrieval
  # Convert boundary_sf to Spatial
  boundary_sp <- as(boundary_sf, "Spatial")
  
  # Define variables to retrieve from SSURGO
  vars <- c("sandtotal_r", "silttotal_r", "claytotal_r", "awc_r", "om_r", "dbovendry_r")
  
  # Download SSURGO properties
  ssurgo <- get_ssurgo_props(boundary_sp, vars = vars)
  
  # Select and rename columns
  ssurgo <- ssurgo %>%
    dplyr::select(c("mukey", "musym", "muname", "sandtotal_r", "silttotal_r", "claytotal_r", "awc_r")) %>%
    dplyr::rename(
      clay = "claytotal_r",
      sand = "sandtotal_r",
      silt = "silttotal_r",
      water_storage = "awc_r"
    )

  return(list(topo_dem = topo_dem, ssurgo = ssurgo))  # Return both datasets as a list
}
