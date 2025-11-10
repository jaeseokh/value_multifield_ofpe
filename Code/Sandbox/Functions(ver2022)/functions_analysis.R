# ---
# SCRIPT: functions_for_analysis.R


# PURPOSE: 
#   This script provides a suite of functions dedicated to the data analysis 
#   phase of the pipeline (executed in '2.Data_Analysis.Rmd'). These functions 
#   support the core modeling tasks, including model specification, evaluation,
#   and the final economic analysis.
#

# KEY FUNCTIONS:
#   - generate_gam_formulas: Systematically creates a comprehensive set of 
#     Generalized Additive Model (GAM) formulas to test different functional 
#     forms (linear vs. splines) for predictor variables.
#   - run_and_evaluate_gams: Iterates through the generated formulas, fits a 
#     GAM for each, and identifies the best-performing model based on the 
#     Generalized Cross-Validation (GCV) score.
#   - predict_yield_range: Uses a fitted model to predict yield outcomes across
#     a range of input rates.
#   - estimate_profit_fcn: Calculates the estimated profit based on predicted
#     yields and input/crop price data, which is essential for determining the
#     Economically Optimal Nitrogen Rate (EONR).
#

# AUTHOR: 
#   Jaeseok Hwang
#

# ---

# Function to generate all possible GAM formulas
#  for a given set of input variables 
#  and field regression variables
#  each input and set of field variables can be either linear or spline
#  for the spline function, input variables can choose k=3 or k=4
#  in case of field variables, all of them will be either linear or spline with k=3
#  the function returns a list of formulas


generate_gam_formulas <- function(input_vars, field_reg_vars) {
  # Define possible spline options for input_vars
  input_var_options <- list(
    "linear" = input_vars,
    "spline_k3" = paste0("s(", input_vars, ",k=3)"),
    "spline_k4" = paste0("s(", input_vars, ",k=4)")
  )
  
  # Possible choices for field_reg_vars (either all linear or all splines with k=3)
  field_reg_linear <- paste(field_reg_vars, collapse = "+")
  field_reg_spline <- paste0("s(", field_reg_vars, ",k=3)", collapse = "+")
  
  field_var_options <- list("linear" = field_reg_linear, "spline_k3" = field_reg_spline)
  
  # Function to generate formula for given choices of input and field variables
  create_formula <- function(input_var_comb, field_var_choice) {
    input_part <- paste(input_var_comb, collapse = "+")
    formula <- paste("yield ~", input_part, "+", field_var_choice)
    return(formula)
  }
  
  # Generate all combinations of input_vars and field_reg_vars
  formulas <- list()
  
  # If input_vars is of length 1
  if (length(input_vars) == 1) {
    input_combinations <- expand.grid(
      input_1 = names(input_var_options),
      field = names(field_var_options),
      stringsAsFactors = FALSE
    )
    
    # Generate formula for each combination
    for (i in seq_len(nrow(input_combinations))) {
      input_1_choice <- input_var_options[[input_combinations$input_1[i]]]
      field_choice <- field_var_options[[input_combinations$field[i]]]
      formulas[[i]] <- create_formula(input_1_choice, field_choice)
    }
  }
  
  # If input_vars is of length 2
  if (length(input_vars) == 2) {
    input_combinations <- expand.grid(
      input_1 = names(input_var_options),
      input_2 = names(input_var_options),
      field = names(field_var_options),
      stringsAsFactors = FALSE
    )
    
    # Generate formula for each combination
    for (i in seq_len(nrow(input_combinations))) {
      input_1_choice <- input_var_options[[input_combinations$input_1[i]]][1]
      input_2_choice <- input_var_options[[input_combinations$input_2[i]]][2]
      field_choice <- field_var_options[[input_combinations$field[i]]]
      formulas[[i]] <- create_formula(c(input_1_choice, input_2_choice), field_choice)
    }
  }
  
  return(formulas)
}


# Function to run and evaluate GAM models for a given set of formulas and a data set.
# The function fits GAM models using REML and GCV.Cp methods and by each formula,
# evaluates them based on REML and GCV scores.
# The function returns the two best models based on REML and GCV scores,respectively.


run_and_evaluate_gams <- function(gam_formulas, data) {
  results <- data.frame()
  
  # Loop through each formula
  for (i in seq_along(gam_formulas)) {
    formula <- as.formula(gam_formulas[[i]])
    
    # Fit the GAM model
    gam_model_gcv <- gam(formula, data = data, method = "GCV.Cp")
    
    # Extract evaluation metrics
    gcv_score <- gam_model_gcv$gcv.ubre  # GCV score
    
    # Store the results in a data frame
    results <- rbind(results, data.frame(
      formula = paste(gam_formulas[[i]]), 
      gcv = gcv_score
    ))

    rownames(results) <- NULL
  }

  best_gcv <- as.formula(results[which.min(results$gcv),"formula"])
  
    gam_gcv <- gam(best_gcv, data = data, method = "GCV.Cp")
    

 return(list(gam_best = gam_gcv))
}


     predict_yield_range <- function( formula,data_for_evaluation) {

  #--- predict yield ---#
  yield_prediction <- predict(formula, newdata =  data_for_evaluation, se = TRUE)

 evaluation_data <- data_for_evaluation[, .(s_rate)][, `:=`(
    yield_hat = yield_prediction$fit,
    yield_hat_se = yield_prediction$se.fit
     )]

  return(evaluation_data)

}

  estimate_profit_fcn <- function(info_table, evaluation_data, price_table) {
  # Get prices for the given year from price_tab
  corn_price_year <- price_table[Year == info_table$year, corn_price]
  seed_price_year <- price_table[Year == info_table$year, seed_price]

  # Calculate profit estimates in eval_data
  evaluation_table <- evaluation_data %>%
    .[, profit_hat_year := corn_price_year * yield_hat - seed_price_year * s_rate] %>%
    .[, profit_hat_low := corn_price_low * yield_hat - seed_price_low * s_rate] %>%
    .[, profit_hat_high := corn_price_high * yield_hat - seed_price_high * s_rate] %>%
    .[, profit_hat_se := corn_price_year * yield_hat_se]

  return(evaluation_table)
}
