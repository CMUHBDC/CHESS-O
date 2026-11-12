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
#   Absolute risk increases (ARIs, per 1,000 persons; last column of Table 3)
#   were estimated for the 1-, 2-, and 5-year windows by standardization
#   with death as a competing risk:
#     1. Two cause-specific Cox models with Model 3 covariates were fitted
#        within each window: fracture (death censored; identical to Model 3
#        above) and death without prior fracture (fracture censored).
#     2. The two models were combined to predict each participant's
#        cumulative incidence of fracture (CIF) at t years.
#     3. Categorical exposure: CIFs were averaged over the whole cohort with
#        everyone set to osteoporosis (AR_exposed) and to non-osteoporosis
#        (AR_reference); ARI = AR_exposed - AR_reference.
#        Continuous exposure: average CIF difference if every participant's
#        bmd_score were 0.1 unit lower than observed.
#     4. 95% CIs: nonparametric percentile bootstrap.
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
#   - Death_competing:   1 = died during overall follow-up without a prior
#                        fracture (time of death = FT_Fracture), 0 = otherwise
#                        (used only for the absolute risk analysis)
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
library(purrr)
library(prodlim)         # Hist() for competing-risk outcomes
library(riskRegression)  # CSC(), predictRisk()

# Absolute risk analysis settings
COVARS_MODEL3 <- c("Age", "Sex", "Cancer", "RA", "COPD", "Steroid")
ARI_PER       <- 1000   # absolute risks expressed per 1,000 persons
B_BOOT        <- 1000   # number of bootstrap replicates
SEED          <- 123    # replicate b uses seed SEED + b (reproducible in parallel)
N_CORES       <- 1      # > 1 runs the bootstrap in parallel


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
# 6. Functions: absolute risk increase (ARI) 
#    (cause-specific Cox models, death as a competing risk, Model 3 covariates)
#-------------------------------------------------------------------------------
ari_windows <- tibble(
  outcome   = c("1-year", "2-year", "5-year"),
  t         = c(1, 2, 5),
  time_var  = c("FT_Fracture_1yr", "FT_Fracture_2yr", "FT_Fracture_5yr"),
  event_var = c("Fracture_1yr", "Fracture_2yr", "Fracture_5yr")
)

# Competing-risk outcome for a t-year window:
#   status_t = 1: fracture within t years (window event indicator)
#   status_t = 2: death without prior fracture within t years
#   status_t = 0: censored
make_ari_window <- function(data, w, exposure_var, covars = COVARS_MODEL3) {
  data %>%
    mutate(
      time_t   = pmin(.data[[w$time_var]], w$t),
      status_t = case_when(
        .data[[w$event_var]] == 1 & .data[[w$time_var]] <= w$t ~ 1L,
        Death_competing == 1 & FT_Fracture <= w$t             ~ 2L,
        TRUE                                                   ~ 0L
      )
    ) %>%
    select(time_t, status_t, all_of(c(exposure_var, covars)),
           all_of(c(w$time_var, w$event_var))) %>%
    drop_na()
}

csc_formula <- function(exposure_var, covars = COVARS_MODEL3) {
  as.formula(paste("Hist(time_t, status_t) ~",
                   paste(c(exposure_var, covars), collapse = " + ")))
}

# Counterfactual exposure assignment
#   categorical: everyone set to the reference (level 1 / 0) or exposed (level 2 / 1)
#   continuous:  bmd_score_neg_01 + 1 = bmd_score 0.1 unit lower than observed
set_exposure <- function(d, exposure_var, exposed, continuous) {
  x <- d[[exposure_var]]
  if(continuous) {
    if(exposed) d[[exposure_var]] <- x + 1
  } else if(is.factor(x)) {
    lv <- levels(x)
    d[[exposure_var]] <- factor(lv[if(exposed) 2 else 1], levels = lv)
  } else {
    d[[exposure_var]] <- if(exposed) 1 else 0
  }
  d
}

predict_cif <- function(fit, newdata, t, vars) {
  key <- do.call(paste, c(lapply(newdata[vars], as.character), sep = "\r"))
  u   <- !duplicated(key)
  r   <- as.numeric(predictRisk(fit, newdata = newdata[u, , drop = FALSE],
                                times = t, cause = 1))
  r[match(key, key[u])]
}

# Standardized absolute risks (probability scale)
g_comp_ari <- function(d, exposure_var, t, continuous, covars = COVARS_MODEL3) {
  fit  <- CSC(csc_formula(exposure_var, covars), data = d, cause = 1)
  vars <- c(exposure_var, covars)
  r0   <- predict_cif(fit, set_exposure(d, exposure_var, FALSE, continuous), t, vars)
  r1   <- predict_cif(fit, set_exposure(d, exposure_var, TRUE,  continuous), t, vars)
  c(AR_reference = mean(r0), AR_exposed = mean(r1), ARI = mean(r1 - r0))
}

boot_g_comp_ari <- function(d, exposure_var, t, continuous, B = B_BOOT) {
  one_rep <- function(b) {
    set.seed(SEED + b)
    db <- d[sample.int(nrow(d), replace = TRUE), , drop = FALSE]
    tryCatch(g_comp_ari(db, exposure_var, t, continuous),
             error = function(e) c(AR_reference = NA, AR_exposed = NA, ARI = NA))
  }
  if(N_CORES > 1) {
    cl <- parallel::makeCluster(N_CORES)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, {
      library(dplyr); library(survival); library(prodlim); library(riskRegression)
    })
    parallel::clusterExport(cl, c("g_comp_ari", "csc_formula", "set_exposure",
                                  "predict_cif", "COVARS_MODEL3", "SEED"),
                            envir = globalenv())
    res <- parallel::parLapply(cl, seq_len(B), one_rep)
  } else {
    res <- lapply(seq_len(B), one_rep)
  }
  do.call(rbind, res)
}

# ARI (95% CI) per 1,000 persons for each fixed-horizon window
run_ari_analysis <- function(data, exposure_var = "Osteoporosis", continuous = FALSE) {

  if(!continuous && !is.factor(data[[exposure_var]]) && !is.numeric(data[[exposure_var]])) {
    data[[exposure_var]] <- factor(data[[exposure_var]])
  }

  map_dfr(seq_len(nrow(ari_windows)), function(i) {
    w   <- ari_windows[i, ]
    d   <- make_ari_window(data, w, exposure_var)
    est <- g_comp_ari(d, exposure_var, w$t, continuous)
    bt  <- boot_g_comp_ari(d, exposure_var, w$t, continuous)
    ci  <- apply(bt, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)

    tibble(
      outcome   = w$outcome,
      exposure  = exposure_var,
      t         = w$t,
      n         = nrow(d),
      fractures = sum(d$status_t == 1),
      deaths    = sum(d$status_t == 2),
      AR_reference       = est["AR_reference"] * ARI_PER,
      AR_reference_lower = ci[1, "AR_reference"] * ARI_PER,
      AR_reference_upper = ci[2, "AR_reference"] * ARI_PER,
      AR_exposed         = est["AR_exposed"] * ARI_PER,
      AR_exposed_lower   = ci[1, "AR_exposed"] * ARI_PER,
      AR_exposed_upper   = ci[2, "AR_exposed"] * ARI_PER,
      ARI       = est["ARI"] * ARI_PER,
      ARI_lower = ci[1, "ARI"] * ARI_PER,
      ARI_upper = ci[2, "ARI"] * ARI_PER,
      boot_failed = sum(is.na(bt[, "ARI"]))
    ) %>%
      mutate(ari_ci = paste0(format_values(ARI), " (",
                             format_values(ARI_lower), ", ",
                             format_values(ARI_upper), ")"))
  })
}


#-------------------------------------------------------------------------------
# 7. Function: assemble summary table (incidence + hazard ratios by model
#    + ARI per 1,000 persons, Model 3)
#-------------------------------------------------------------------------------
create_fracture_table <- function(cox_results, incidence_results, exposure_label,
                                  continuous = FALSE, ari_results = NULL) {

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
      )
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
      )
  }

  table_cols <- c("outcome", "Events_NonEvents", "person_years", "incidence_per1000PY",
                  "crude", "model1", "model2", "model3")

  if(!is.null(ari_results)) {
    table_data <- table_data %>%
      left_join(ari_results %>% select(outcome, ari_ci), by = "outcome") %>%
      mutate(
        ARI_model3 = coalesce(ari_ci, "—"),
        ARI_model3 = if(continuous) ARI_model3 else if_else(is_reference, "Ref.", ARI_model3)
      )
    table_cols <- c(table_cols, "ARI_model3")
  }

  table_data <- table_data %>% select(all_of(table_cols))

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
  if(!is.null(ari_results)) n_row$ARI_model3 <- ""

  bind_rows(n_row, table_data)
}


#-------------------------------------------------------------------------------
# 8. Run analysis
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


ari_osteoporosis <- run_ari_analysis(
  data         = analytic_data,
  exposure_var = "Osteoporosis",
  continuous   = FALSE
)

ari_bmd <- run_ari_analysis(
  data         = analytic_data,
  exposure_var = "bmd_score_neg_01",
  continuous   = TRUE
)

table_osteoporosis <- create_fracture_table(
  results_osteoporosis$cox_results,
  results_osteoporosis$incidence,
  exposure_label = "Osteoporosis",
  continuous     = FALSE,
  ari_results    = ari_osteoporosis
)

table_bmd <- create_fracture_table(
  results_bmd$cox_results,
  results_bmd$incidence,
  exposure_label = "bmd_score_neg_01",
  continuous     = TRUE,
  ari_results    = ari_bmd
)
