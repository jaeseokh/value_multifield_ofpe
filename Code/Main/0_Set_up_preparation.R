library(here)
library(data.table)
library(magrittr)
library(dplyr)
library(sf)

### Variable name uniformity  ##


# Function to safely convert to POSIXct using mdy_hms if necessary
convert_to_datetime <- function(time_value) {
  # Check if time_value is already in Date or POSIXct format
  if (inherits(time_value, "Date") | inherits(time_value, "POSIXct")) {
    return(time_value)  # Return as is if already in the correct format
  }
  
  # If time_value is logical or not in the correct format, use mdy_hms to convert
  if (is.logical(time_value) | is.character(time_value)) {
    tryCatch({
      return(mdy_hms(time_value))  # Attempt to convert using mdy_hms
    }, error = function(e) {
      return(NA)  # If conversion fails, return NA
    })
  }
  
  return(NA)  # Return NA if the type is unrecognized
}



# 1. choose all type of seed and nitrogen rate variable names
vars_to_select <- c("obs_id","yield", "seed_rate", "s_rate", "n_rate", "uan32_rate", "nh3_rate", "uan28_rate", "urea_rate")

# 2. sort out only the experimental  input variables 
vars_to_check <-  c( "seed_rate", "s_rate", "n_rate", "uan32_rate", "nh3_rate", "uan28_rate", "urea_rate")

# 3. Field specific variables to be included in the Regression ##
field_reg_vars <- c('elev','slope','aspect','tpi', 'clay', 'sand', 
                         'silt', 'water_storage')
 

### Price information (extract and calculate crop and input prices information from the raw data)

# 1. Extract the corn price-received by year
 corn_price_raw <- fread(here("Data","Raw","corn_price_received_by_year.csv"))
 
   corn_price_tab <- corn_price_raw[Period == "MARKETING YEAR",.(Year, Value)]
     
     setnames(corn_price_tab, "Value", "corn")
# 2. Extract the nitrogen price (N equivalent) by year

  n_price_raw <- fread(here("Data","Raw","nitrogen_fertilizer_prices.csv"))
 
  n_price_tab <- n_price_raw[!(Year %in% c("2014", "2015")), 
                      .(Year = as.numeric(Year), 
                        nitrogen = rowMeans(.SD, na.rm = TRUE)), 
                      .SDcols = 2:4]

                    

# 3. Merge the corn price and seed cost information by year
price_tab <- merge(corn_price_tab, n_price_tab, , by = "Year")

price_tab <- price_tab %>% data.table() %>%
  setnames(tolower(names(price_tab)))


  
### Extract USDA reported annual average corn seed rates by state and year

# 1. read the US county map shape file
  us_map_sf <- st_read(here("Data","Raw","us_state_2023.shp"))  %>%
         st_transform(4326)
 
  # 2. sort out the state fips code for the states where experimental(ofpe) field are located
   ofpe_fips <- c(17, 18, 19, 20, 21, 26, 27, 29, 31, 38, 39, 46, 55)

    ofpe_sf <- us_map_sf %>%
               setnames(names(.), tolower(names(.))) %>%
               filter(statefp %in% ofpe_fips) %>% 
               mutate(fips = as.numeric(statefp)) %>%
               dplyr::select(fips,stusps)


# 3. add the USDA reported seed rate information by state and year

  seed_usda_raw <- fread(here("Data","Raw","corn_annual_plant_pop.csv"))

   seed_usda <-seed_usda_raw %>%
    setnames(names(.), tolower(names(.))) %>%        
     setnames("state ansi", "fips") %>%
     .[period=='YEAR',.(year, fips, s_rate =as.numeric(gsub("[^0-9.]", "", value))/1000)] 



