

process_info_table <- function(input_tb_binded) {
  # Step 2: Calculate the quadratic fit and yield maximizing seeding rate for each field
  for(i in 1:length(input_tb_binded$ffy_id)){
    
    eval_s <- readRDS(here("Data", "processed", "Analysis_results", paste0(input_tb_binded$ffy_id[i], "_eval_tb.rds")))
    
    # Apply quadratic fit to each element in eval_s
    quad_fit <- lm(yield_hat ~ s_rate + I(s_rate^2), data = eval_s)
    
    # Check the coefficients of the quadratic fit
    coef_2nd <- quad_fit$coef[3]
    
    # Update input_tb_binded with the quadratic fit coefficient
    input_tb_binded[i, `:=`(quad_fit = coef_2nd)]
  }
  
  # Step 3: Generate response type variable based on the conditions
  input_tb_binded <- input_tb_binded[, resp_type := fcase(
    dif_s < 0 & quad_fit < 0 & ymsr == eosr, "A1",  # under_seed & corner sol
    dif_s < 0 & quad_fit < 0 & ymsr != eosr, "A2",
    dif_s > 0 & quad_fit < 0 & ymsr == eosr, "B1",  # over_seed & corner sol
    dif_s > 0 & quad_fit < 0 & ymsr != eosr, "B2",  # over_seed & corner sol
    quad_fit > 0, "C",
    default = NA_character_  # Default value if no conditions match
  )]

  # Filter out rows with NA response type
  input_tb_resp <- input_tb_binded[!is.na(resp_type)]
  
  # Convert resp_type to a factor with specified levels
  input_tb_resp <- input_tb_resp[, resp_type := factor(resp_type, levels = c("A1", "A2", "B1", "B2", "C"))]

  return(input_tb_resp)
}
