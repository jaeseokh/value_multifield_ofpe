# Code/src/get_il_yield_nass.R
# -------------------------------------------------------------------
# Build county-level yield benchmarks for IL corn from NASS (2016–2023),
# attach them to OFPE field-years, and save:
#   - Data/External/il_county_yield.rds       (table)
#   - Data/External/il_county_yield_sf.rds    (sf polygons)
#   - Data/Processed/APSIP_compare_ready/field_benchmarks.rds
# -------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(sf)
  library(tigris)
  library(rnassqs)
  library(here)
  library(purrr)
})

options(tigris_use_cache = TRUE)

# -------------------------------------------------------------------
# 0. Parameters and folders
# -------------------------------------------------------------------

years_target <- 2010:2024
state_alpha  <- "IL"

nass_dir <- here("Data","Processed","NASS")
apsim_ready_dir <- here("Data", "Processed", "APSIP_compare_ready")
analysis_dir <- here("Data", "Processed", "Analysis_ready")

# -------------------------------------------------------------------
# 1. NASS Quick Stats: county-level corn yields for IL
# -------------------------------------------------------------------
#  NASS API key as an env var, e.g.
 Sys.setenv(NASS_KEY = "882975A4-231A-36E1-99AA-9B3BA4355E33")

if (is.na(Sys.getenv("NASS_KEY", unset = NA))) {
  stop("Please set Sys.setenv(NASS_KEY = '...') before running this script.")
}


# Check API key status
nassqs_auth(Sys.getenv("NASS_KEY"))
options(nassqs.auth = list(api_key = "882975A4-231A-36E1-99AA-9B3BA4355E33"))

getOption("nassqs.auth")

# helper: query one year at a time with the *known working* combo
get_year_county <- function(yy) {
  message("→ Downloading NASS county yields for ", yy, " ...")
  nassqs(list(
    commodity_desc    = "CORN",
    year              = yy,
    state_alpha       = state_alpha,
    statisticcat_desc = "YIELD",
    agg_level_desc    = "COUNTY"
  ))
}

raw_list <- purrr::map(years_target, get_year_county)
raw_county <- bind_rows(raw_list)


il_county_yield <- raw_county %>%
  dplyr::mutate(
    year        = as.integer(year),
    fips        = paste0("17", county_code),
    yield_bu_ac = as.numeric(gsub(",", "", Value))
  ) %>%
  dplyr::filter(
    util_practice_desc == "GRAIN",
    unit_desc          == "BU / ACRE",
    year %in% years_target
  ) %>%
  dplyr::select(year, fips, county_name, county_code, yield_bu_ac)


# wide table: one row per county, columns Y_county_2010 ... Y_county_2024
il_county_yield_wide <- il_county_yield %>%
  tidyr::pivot_wider(
    id_cols     = c(fips, county_name, county_code),
    names_from  = year,
    values_from = yield_bu_ac,
    names_prefix = "Y_county_"
  )


# contemporaneous (same year) as a separate helper
county_current <- il_county_yield %>%
  dplyr::select(
    year,
    fips,
    Y_county_curr = yield_bu_ac
  )

# -------------------------------------------------------------------
# 2. Make an sf object for mapping / spatial joins (optional but useful)
# -------------------------------------------------------------------

message("→ Building sf polygons for IL counties ...")

il_counties_sf <- tigris::counties(state = state_alpha, cb = TRUE, year = 2020) %>%
  st_transform(4326) %>%
  mutate(fips = GEOID)

il_county_yield_sf <- il_counties_sf %>%
  left_join(il_county_yield, by = "fips")

saveRDS(il_county_yield_sf,
        file = file.path(nass_dir, "il_county_yield_sf.rds"))

# -------------------------------------------------------------------
# 3. Attach county yield benchmarks to OFPE field-years
# -------------------------------------------------------------------

message("→ Attaching county yields to field_loc2 ...")

field_loc2 <- readRDS(file.path(analysis_dir, "field_loc2.rds"))

il_ffy_ids <- field_loc2 %>%
  dplyr::filter(stusps == "IL") %>%
  dplyr::pull(ffy_id) %>%
  unique()

# helper: for a single ffy_id, find its county via centroid and counties_sf
compute_field_county <- function(ffy_id,
                                 counties_sf,
                                 base_dir_raw = here::here("Data", "Raw")) {

  bdry <- read_boundary(ffy_id, base_dir = base_dir_raw)

  if (is.null(bdry) || nrow(bdry) == 0) {
    warning("No boundary for ", ffy_id)
    return(NULL)
  }

  # match CRS
  bdry <- sf::st_transform(bdry, sf::st_crs(counties_sf))

  # use centroid/point-on-surface to avoid multi-polygon headaches
  cen <- sf::st_point_on_surface(bdry)

  hit <- sf::st_join(cen, counties_sf["fips"], left = TRUE)

  if (nrow(hit) == 0 || is.na(hit$fips[1])) {
    warning("No county match for ", ffy_id)
    return(NULL)
  }

  tibble::tibble(
    ffy_id = ffy_id,
    fips   = hit$fips[1]
  )
}

field_counties_il <- purrr::map_dfr(
  il_ffy_ids,
  compute_field_county,
  counties_sf = il_counties_sf
)

saveRDS(field_counties_il,
        file.path(apsim_ready_dir, "field_counties_il.rds"))


field_benchmarks <- field_loc2 %>%
  dplyr::filter(stusps == "IL") %>%
  dplyr::select(-stusps, -fips) %>%   # drop any old fips / state code
  # attach county fips for each field
  dplyr::left_join(field_counties_il, by = "ffy_id") %>%
  # attach full county yield history (Y_county_2010...Y_county_2024)
  dplyr::left_join(il_county_yield_wide, by = "fips") %>%
  # attach contemporaneous county yield for that year
  dplyr::left_join(county_current, by = c("year", "fips"))


saveRDS(
  field_benchmarks,
  file = file.path(nass_dir, "field_benchmarks.rds")
)
