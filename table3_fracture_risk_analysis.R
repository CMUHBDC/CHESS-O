#===============================================================================
# CHESS-O: Fracture Risk Analysis (Cause-Specific Cox Proportional Hazards Models)
# Table 3: Risk of incident nontraumatic major osteoporotic fracture
#          associated with CHESS-O-estimated osteoporosis
#===============================================================================
#
# Description:
#   This script reproduces the cause-specific Cox proportional hazards analysis
#   of incident nontraumatic major osteoporotic fracture risk associated with
#   CHESS-O-estimated osteoporosis status (categorical) and BMD score
#   (continuous), across four follow-up windows (1-year, 2-year, 5-year, and
#   overall follow-up).
#
#   Cause-specific Cox models were used, with death treated as a competing
#   event: participants who died during follow-up without experiencing a
#   fracture were censored at the time of death.
#
#   For each follow-up window, four nested Cox models were fitted:
#     - Crude:   exposure only
#     - Model 1: + Age, Sex
#     - Model 2: + Age, Sex, Cancer, RA, COPD
#     - Model 3: + Age, Sex, Cancer, RA, COPD, Steroid use
#
#   Fracture incidence rates (per 1,000 person-years) were also calculated
#   for each follow-up window.
#
# Data availability:
#   Individual-level data are not included because of institutional
#   data-governance and participant-privacy restrictions.
#
#   The object `analytic_data` is a placeholder for the analytic dataset
#   (study population for clinical outcome validation, n = 76,051). 
#
# Required variables:
#   - Osteoporosis:      Binary/categorical exposure estimated from the
#                        CHESS-O model's predicted classification 
#                        (0 = non-osteoporosis, 1 = osteoporosis)
#   - bmd_score:         Continuous BMD-based risk score estimated by the
#                        CHESS-O model
#   - bmd_score_neg_01:  -bmd_score / 0.1 (This rescaling expresses the hazard
#                         ratio per 0.1-unit decrease in bmd_score.)
#   - Age, Sex:          Demographic covariates
#   - Cancer, RA, COPD:  Comorbidity covariates
#   - Steroid:           Medication covariate
#   - FT_Fracture_1yr / Fracture_1yr: 1-year follow-up time / event indicator
#   - FT_Fracture_2yr / Fracture_2yr: 2-year follow-up time / event indicator
#   - FT_Fracture_5yr / Fracture_5yr: 5-year follow-up time / event indicator
#   - FT_Fracture / Fracture:         Overall follow-up time / event indicator
#
# Note: format_values() rounds and formats numeric values for table display
#       with precision that scales with magnitude.
#
#===============================================================================


#-------------------------------------------------------------------------------
# 1. Environment and packages
#-------------------------------------------------------------------------------
library(survival)
library(dplyr)
library(tidyr)


#-------------------------------------------------------------------------------
# 2. Function: fit nested Cox proportional hazards models
#-------------------------------------------------------------------------------
create_cox_models_fracture <- function(data, time_var, event_var, exposure_var = "Osteoporosis") {
  
  if(!is.factor(data[[exposure_var]]) && !is.numeric(data[[exposure_var]])) {
    data[[exposure_var]] <- factor(data[[exposure_var]])
  }
  
  models <- list()
  
  tryCatch({
    models$crude <- coxph(
      as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", exposure_var)),
      data = data
    )
  }, error = function(e) cat("Crude model error:", e$message, "\n"))
  
  if(all(c("Age", "Sex") %in% names(data))) {
    tryCatch({
      models$model1 <- coxph(
        as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                          exposure_var, " + Age + Sex")),
        data = data
      )
    }, error = function(e) cat("Model 1 error:", e$message, "\n"))
  }
  
  if(all(c("Age", "Sex", "Cancer", "RA", "COPD") %in% names(data))) {
    tryCatch({
      models$model2 <- coxph(
        as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                          exposure_var, " + Age + Sex + Cancer + RA + COPD")),
        data = data
      )
    }, error = function(e) cat("Model 2 error:", e$message, "\n"))
  }
  
  if(all(c("Age", "Sex", "Cancer", "RA", "COPD", "Steroid") %in% names(data))) {
    tryCatch({
      models$model3 <- coxph(
        as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ",
                          exposure_var, " + Age + Sex + Cancer + RA + COPD + Steroid")),
        data = data
      )
    }, error = function(e) cat("Model 3 error:", e$message, "\n"))
  }
  
  attr(models, "exposure_var") <- exposure_var
  return(models)
}


#-------------------------------------------------------------------------------
# 3. Function: extract hazard ratios and 95% CIs from fitted models
#-------------------------------------------------------------------------------
extract_cox_results_fracture <- function(models) {
  exposure_var <- attr(models, "exposure_var")
  results_list <- list()
  
  for(model_name in names(models)) {
    model <- models[[model_name]]
    if(is.null(model)) next
    
    summary_model <- summary(model)
    coef_names    <- rownames(summary_model$coefficients)
    target_coef   <- grep(exposure_var, coef_names, value = TRUE)[1]
    
    if(length(target_coef) > 0 && !is.na(target_coef) && target_coef %in% coef_names) {
      hr       <- summary_model$coefficients[target_coef, "exp(coef)"]
      ci_lower <- summary_model$conf.int[target_coef, "lower .95"]
      ci_upper <- summary_model$conf.int[target_coef, "upper .95"]
      p_value  <- summary_model$coefficients[target_coef, "Pr(>|z|)"]
      
      results_list[[model_name]] <- data.frame(
        model    = model_name,
        exposure = exposure_var,
        hr       = hr,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        p_value  = p_value,
        n        = model$n,
        n_events = model$nevent,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if(length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    combined_results$hr_ci <- paste0(
      format_values(combined_results$hr), " (",
      format_values(combined_results$ci_lower), ", ",
      format_values(combined_results$ci_upper), ")"
    )
    return(combined_results)
  }
  return(data.frame())
}


#-------------------------------------------------------------------------------
# 4. Functions: fracture incidence rate (per 1,000 person-years)
#-------------------------------------------------------------------------------
#
# calculate_fracture_incidence(): stratified by exposure group (categorical
#   exposure). The first group (by row order) is flagged as the reference
#   category for hazard ratio reporting.
# calculate_fracture_incidence_continuous(): overall cohort incidence
#   (continuous exposure; no stratification).
#-------------------------------------------------------------------------------
calculate_fracture_incidence <- function(data, time_var, event_var, group_var) {
  data %>%
    group_by(.data[[group_var]]) %>%
    summarise(
      Events              = sum(.data[[event_var]] == 1, na.rm = TRUE),
      Non_events          = sum(.data[[event_var]] == 0, na.rm = TRUE),
      person_years        = sum(.data[[time_var]], na.rm = TRUE),
      incidence_per1000PY = Events / person_years * 1000,
      .groups = "drop"
    ) %>%
    mutate(is_reference = row_number() == 1)
}

calculate_fracture_incidence_continuous <- function(data, time_var, event_var) {
  data %>%
    summarise(
      Events              = sum(.data[[event_var]] == 1, na.rm = TRUE),
      Non_events          = sum(.data[[event_var]] == 0, na.rm = TRUE),
      person_years        = sum(.data[[time_var]], na.rm = TRUE),
      incidence_per1000PY = Events / person_years * 1000
    )
}


#-------------------------------------------------------------------------------
# 5. Main analysis wrapper: loop over follow-up windows
#-------------------------------------------------------------------------------
run_fracture_analysis_comprehensive <- function(data, exposure_var = "Osteoporosis", continuous = FALSE) {
  
  fracture_outcomes <- list(
    "1-year"  = list(time = "FT_Fracture_1yr", event = "Fracture_1yr"),
    "2-year"  = list(time = "FT_Fracture_2yr", event = "Fracture_2yr"),
    "5-year"  = list(time = "FT_Fracture_5yr", event = "Fracture_5yr"),
    "Overall" = list(time = "FT_Fracture",    event = "Fracture")
  )
  
  all_results   <- list()
  all_incidence <- list()
  
  for(outcome_name in names(fracture_outcomes)) {
    time_var  <- fracture_outcomes[[outcome_name]]$time
    event_var <- fracture_outcomes[[outcome_name]]$event
    
    if (!all(c(time_var, event_var) %in% names(data))) {
      next
    }
    
    tryCatch({
      models  <- create_cox_models_fracture(data, time_var, event_var, exposure_var)
      results <- extract_cox_results_fracture(models)
      if(nrow(results) > 0) {
        results$outcome <- outcome_name
        all_results[[outcome_name]] <- results
      }
    }, error = function(e) cat("Cox model error:", e$message, "\n"))
    
    tryCatch({
      incidence <- if(continuous) {
        calculate_fracture_incidence_continuous(data, time_var, event_var)
      } else {
        calculate_fracture_incidence(data, time_var, event_var, exposure_var)
      }
      incidence$outcome <- outcome_name
      all_incidence[[outcome_name]] <- incidence
    }, error = function(e) cat("Incidence calculation error:", e$message, "\n"))
  }
  
  cox_results       <- if(length(all_results)   > 0) do.call(rbind, all_results)   else data.frame()
  incidence_results <- if(length(all_incidence) > 0) do.call(rbind, all_incidence) else data.frame()
  
  return(list(cox_results = cox_results, incidence = incidence_results))
}


#-------------------------------------------------------------------------------
# 6. Function: assemble summary table (incidence + hazard ratios by model)
#-------------------------------------------------------------------------------
create_fracture_table <- function(cox_results, incidence_results, exposure_label, continuous = FALSE) {
  
  model_n <- cox_results %>%
    group_by(model) %>% slice(1) %>% ungroup() %>%
    select(model, n) %>%
    mutate(n_label = paste0("N = ", n))
  
  if(continuous) {
    table_data <- incidence_results %>%
      left_join(
        cox_results %>%
          select(outcome, model, hr_ci) %>%
          pivot_wider(names_from = model, values_from = hr_ci),
        by = "outcome"
      ) %>%
      mutate(
        Events_NonEvents    = paste0(Events, " / ", Non_events),
        person_years        = format_values(person_years),
        incidence_per1000PY = format_values(incidence_per1000PY)
      ) %>%
      select(outcome, Events_NonEvents, person_years, incidence_per1000PY,
             crude, model1, model2, model3)
  } else {
    table_data <- incidence_results %>%
      left_join(
        cox_results %>%
          select(outcome, model, hr_ci) %>%
          pivot_wider(names_from = model, values_from = hr_ci),
        by = "outcome"
      ) %>%
      mutate(
        Events_NonEvents    = paste0(Events, " / ", Non_events),
        person_years        = format_values(person_years),
        incidence_per1000PY = format_values(incidence_per1000PY),
        crude  = if_else(is_reference, "Ref.", crude),
        model1 = if_else(is_reference, "Ref.", model1),
        model2 = if_else(is_reference, "Ref.", model2),
        model3 = if_else(is_reference, "Ref.", model3)
      ) %>%
      select(outcome, Events_NonEvents, person_years, incidence_per1000PY,
             crude, model1, model2, model3)
  }
  
  n_row <- tibble(
    outcome             = "N",
    Events_NonEvents    = "",
    person_years        = "",
    incidence_per1000PY = "",
    crude  = model_n %>% filter(model == "crude")  %>% pull(n_label) %>% first(),
    model1 = model_n %>% filter(model == "model1") %>% pull(n_label) %>% first(),
    model2 = model_n %>% filter(model == "model2") %>% pull(n_label) %>% first(),
    model3 = model_n %>% filter(model == "model3") %>% pull(n_label) %>% first()
  )
  
  bind_rows(n_row, table_data)
}


#-------------------------------------------------------------------------------
# 7. Run analysis
#-------------------------------------------------------------------------------
#
# Exposure 1: Osteoporosis (categorical)
# Exposure 2: bmd_score_neg_01 (continuous; hazard ratio per 0.1-unit
#             decrease in bmd_score)
#-------------------------------------------------------------------------------
results_osteoporosis <- run_fracture_analysis_comprehensive(
  data         = analytic_data,
  exposure_var = "Osteoporosis",
  continuous   = FALSE
)

results_bmd <- run_fracture_analysis_comprehensive(
  data         = analytic_data,
  exposure_var = "bmd_score_neg_01",
  continuous   = TRUE
)

table_osteoporosis <- create_fracture_table(
  results_osteoporosis$cox_results,
  results_osteoporosis$incidence,
  exposure_label = "Osteoporosis",
  continuous     = FALSE
)

table_bmd <- create_fracture_table(
  results_bmd$cox_results,
  results_bmd$incidence,
  exposure_label = "bmd_score_neg_01",
  continuous     = TRUE
)
