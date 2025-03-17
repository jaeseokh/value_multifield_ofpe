

# Revised function to calculate stage-specific GDD and EDD with proper NA handling

calc_stage_specific_gdd_edd <- function(ffy_id) {
  # Extract centroid of the field boundary
  ffy_year <- temp_daymet <- daymet_t <- daymet_30 <- boundary_sf <- centroid <- NULL
  
  boundary_sf <- readRDS(here("Data", "Raw", "exp_bdry_data", paste0(ffy_id, "_bdry.rds")))
  
  centroid <- boundary_sf %>% 
    st_centroid() %>% 
    st_coordinates()

  # Extract trial year from the ffy string
  ffy_year <- as.numeric(sub(".*_(\\d{4})$", "\\1", ffy_id))

  # Download 30 years of Daymet data for the field centroid
  temp_daymet <- download_daymet(
    lat = centroid[1, "Y"],
    lon = centroid[1, "X"],
    start = ffy_year - 30,
    end = ffy_year
  ) %>% 
    .$data %>% 
    data.table()

  # Extract trial year data
  daymet_t <- temp_daymet %>%
    filter(year == ffy_year) %>%
    rename(prcp = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c.) %>%
    mutate(
      gdd = ifelse(tmax > 10, (tmax + tmin) * 0.5 - 10, 0),
      gdd = pmax(gdd, 0),
      edd = ifelse(tmax > 30, tmax - 30, 0)
    ) %>%
    dplyr::select(prcp, tmax, tmin, yday, gdd, edd) %>%
    mutate(month = day_to_month(yday)) %>%
    filter(month %in% c('Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep')) %>%
    summarize(
      prcp_t = round(sum(prcp, na.rm = TRUE), 1),
      gdd_t = round(sum(gdd, na.rm = TRUE), 1),
      edd_t = round(sum(edd, na.rm = TRUE), 1)
    )

  # Calculate 5-year averages (including ffy_year)
  daymet_5 <- temp_daymet %>%
    filter(year >= (ffy_year - 4) & year <= ffy_year) %>%
    rename(prcp = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c.) %>%
    mutate(
      gdd = ifelse(tmax > 10, (tmax + tmin) * 0.5 - 10, 0),
      gdd = pmax(gdd, 0),
      edd = ifelse(tmax > 30, tmax - 30, 0)
    ) %>%
    summarize(
      prcp_5 = round(mean(prcp, na.rm = TRUE), 1),
      gdd_5 = round(mean(gdd, na.rm = TRUE), 1),
      edd_5 = round(mean(edd, na.rm = TRUE), 1)
    )

  # Calculate 30-year averages (including ffy_year)
  daymet_30 <- temp_daymet %>%
    filter(year >= (ffy_year - 30) & year <= ffy_year) %>%
    rename(prcp = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c.) %>%
    mutate(
      gdd = ifelse(tmax > 10, (tmax + tmin) * 0.5 - 10, 0),
      gdd = pmax(gdd, 0),
      edd = ifelse(tmax > 30, tmax - 30, 0)
    ) %>%
    summarize(
      prcp_30 = round(mean(prcp, na.rm = TRUE), 1),
      gdd_30 = round(mean(gdd, na.rm = TRUE), 1),
      edd_30 = round(mean(edd, na.rm = TRUE), 1)
    )

  # Get the planting (s_time), nitrogen application (n_time), and harvest time (yield_time) from date_info
  s_time <- date_info %>% filter(ffy_id_in == ffy_id) %>% pull(s_time) %>% first()
  n_time <- date_info %>% filter(ffy_id_in == ffy_id) %>% pull(n_time) %>% first()
  yield_time <- date_info %>% filter(ffy_id_in == ffy_id) %>% pull(yield_time) %>% first()

  # Calculate GDD_t and EDD_t for the period between planting and harvest (or nitrogen to harvest)
  # In-season GDD and EDD from planting (s_time) to harvest (yield_time) or nitrogen (n_time)
  stage_gdd_edd <- temp_daymet %>%
    filter(year == ffy_year) %>%
    rename(prcp = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c.) %>%
    mutate(
      # Calculate GDD and EDD only between planting (s_time) and harvest (yield_time)
      in_season = (yday >= as.numeric(format(s_time, "%j"))) & (yday <= as.numeric(format(yield_time, "%j"))),
      gdd_stage = ifelse(tmax > 10, (tmax + tmin) * 0.5 - 10, 0),
      gdd_stage = pmax(gdd_stage, 0),
      edd_stage = ifelse(tmax > 30, tmax - 30, 0)
    ) %>%
    filter(in_season) %>%
    summarize(
      stage_gdd_t = sum(gdd_stage, na.rm = TRUE),
      stage_edd_t = sum(edd_stage, na.rm = TRUE)
    )

  # Calculate GDD and EDD for nitrogen to harvest (n_time to yield_time)
  stage_gdd_edd_nitrogen <- temp_daymet %>%
    filter(year == ffy_year) %>%
    rename(prcp = prcp..mm.day., tmax = tmax..deg.c., tmin = tmin..deg.c.) %>%
    mutate(
      # Calculate GDD and EDD for the period between nitrogen application and harvest
      in_season_nitrogen = (yday >= as.numeric(format(n_time, "%j"))) & (yday <= as.numeric(format(yield_time, "%j"))),
      gdd_stage_nitrogen = ifelse(tmax > 10, (tmax + tmin) * 0.5 - 10, 0),
      gdd_stage_nitrogen = pmax(gdd_stage_nitrogen, 0),
      edd_stage_nitrogen = ifelse(tmax > 30, tmax - 30, 0)
    ) %>%
    filter(in_season_nitrogen) %>%
    summarize(
      stage_gdd_nitrogen = sum(gdd_stage_nitrogen, na.rm = TRUE),
      stage_edd_nitrogen = sum(edd_stage_nitrogen, na.rm = TRUE)
    )

  # Return the results as a list including the new stage-specific GDD and EDD
  return(list(
    prcp_t = daymet_t$prcp_t, 
    gdd_t = daymet_t$gdd_t, 
    edd_t = daymet_t$edd_t,
    prcp_5 = daymet_5$prcp_5, 
    gdd_5 = daymet_5$gdd_5, 
    edd_5 = daymet_5$edd_5,
    prcp_30 = daymet_30$prcp_30, 
    gdd_30 = daymet_30$gdd_30, 
    edd_30 = daymet_30$edd_30,
    stage_gdd_t = stage_gdd_edd$stage_gdd_t,
    stage_edd_t = stage_gdd_edd$stage_edd_t,
    stage_gdd_nitrogen = stage_gdd_edd_nitrogen$stage_gdd_nitrogen,
    stage_edd_nitrogen = stage_gdd_edd_nitrogen$stage_edd_nitrogen,
    s_time = s_time,
    n_time = n_time,
    yield_time = yield_time
  ))
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



