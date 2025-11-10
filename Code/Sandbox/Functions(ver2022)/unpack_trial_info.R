

###############################################################
### Read the information of field-specific parameters ####
## (To make experimental data have their field specific information)

#--- get field parameters for the field-year ---#

w_field_data <- field_data[field_year == ffy_id, ]

# Get the crop specific Information (unit and machine width)
ffy_id <- w_field_data[, field_year]
crop <- tolower(w_field_data[, crop])
crop_unit <- w_field_data[, crop_unit] 
land_unit <- w_field_data[, land_unit] 
reporting_unit <- w_field_data[, reporting_unit] 
harvester_width <- w_field_data[, h_width][[1]]
crop_unit <- w_field_data[, crop_unit] 


input_data <- dplyr::select(w_field_data, starts_with("input.")) %>% 
    lapply(., function(x) x[[1]]) %>% 
    rbindlist(fill = TRUE)

## trial_info table 
# (to be used, as it is merged with processed exp data)

trial_type <- input_data %>%
  filter(strategy == "trial") %>%
  pull(form) %>%
  paste0(., collapse = "_")

trial_info <- 
  tibble(
    ffy_id = ffy_id,
    crop = crop, 
    input_type = input_data$form,
    form = input_data$form,
    unit = input_data$unit,
    process = ifelse(input_data$strategy == "trial", TRUE, FALSE),
    use_td =  input_data$use_target_rate_instead,
    gc_rate = input_data$sq_rate,
    gc_type = ifelse(is.numeric(input_data$sq_rate), "uniform", "Rx")
  ) %>% 
  filter(process == TRUE) %>%
  rowwise() %>%
  mutate(base_rate = get_base_rate(input_data, input_type)) %>%
  mutate(gc_rate_conv = get_gc_rate(gc_rate, input_type, form, unit, convert="T", base_rate)) 

