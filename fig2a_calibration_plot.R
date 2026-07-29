#===============================================================================
# CHESS-O: Decile-Based Calibration Plot and Brier Scores
# Figure 2A — Calibration Analysis Across Validation Cohorts
#===============================================================================
#
# Description:
#   Within each validation cohort, participants were divided into 10 groups 
#   using cohort-specific deciles of CHESS-O-predicted osteoporosis probability.
#
#   Cohort-specific Brier scores were also calculated. The Brier score is the
#   mean squared difference between the predicted probability and the binary
#   observed outcome. Lower values indicate better overall probabilistic
#   accuracy. 
#
# Data availability:
#   Individual-level validation data are not included because of institutional
#   data-governance and participant-privacy restrictions.
#
#   The objects `cmuh_file`, `auh_file`, and `segmed_file` are placeholders
#   for the corresponding analytic datasets.
#
# Input datasets:
#   - CMUH validation cohort (internal held-out):   n = 3,794
#   - AUH validation cohort (external):             n = 6,427
#   - Segmed validation cohort (external):          n = 400
#
# Required variables after harmonization:
#   - prob_col:
#       CHESS-O-predicted probability of osteoporosis, ranging from 0 to 1
#
#   - event_col:
#       Binary reference-standard osteoporosis classification
#       1 = osteoporosis
#       0 = non-osteoporosis
#
# Calibration grouping:
#   Risk groups are defined separately within each validation cohort rather
#   than using pooled cut points across cohorts.
#
#===============================================================================


#-------------------------------------------------------------------------------
# 1. Environment and packages
#-------------------------------------------------------------------------------
rm(list = ls())
library("tidyverse")
library("data.table")
library("ggtext")
options(scipen = 999)


#-------------------------------------------------------------------------------
# 2. Construct cohort-specific quantile-based calibration groups
#-------------------------------------------------------------------------------
get_calib = function(df, prob_col, event_col, model_label){
  
  df = df %>% mutate(pred = .[[prob_col]], event = .[[event_col]]) 
  
  q = quantile(df$pred, probs = seq(0, 1, 0.1), na.rm = TRUE, type = 8)
  
  df = df %>% 
    mutate(decile = case_when(pred <= q[2] ~ 1, 
                              pred > q[2] & pred <= q[3] ~ 2,
                              pred > q[3] & pred <= q[4] ~ 3,
                              pred > q[4] & pred <= q[5] ~ 4,
                              pred > q[5] & pred <= q[6] ~ 5,
                              pred > q[6] & pred <= q[7] ~ 6,
                              pred > q[7] & pred <= q[8] ~ 7,
                              pred > q[8] & pred <= q[9] ~ 8,
                              pred > q[9] & pred <= q[10] ~ 9,
                              pred > q[10] ~ 10))
  
  setDT(df)
  
  tabN = df[, .N, by = decile][order(decile)]
  tabOutcome = df[event == 1, .N, by = decile][order(decile)]
  tabOutcome = merge(tabN, tabOutcome, by = "decile", all.x = TRUE)
  setnames(tabOutcome, c("decile", "N", "Outcome"))
  tabOutcome[is.na(Outcome), Outcome := 0]
  
  tabPred = df[, .(Pred = mean(pred, na.rm = TRUE)), by = decile][order(decile)]
  
  tab = merge(tabOutcome, tabPred, by = "decile")
  tab[, Obs := Outcome / N]
  tab[, Model := model_label]
  
  tab = tab %>% rename(Group = decile) %>% select(Model, Group, N, Outcome, Pred, Obs)
  
  return(tab)
}

calib_table = map_dfr(names(datasets), function(dataset_name) {
  get_calib(datasets[[dataset_name]], prob_col, event_col, dataset_name)
})

color_values = c("indianred1", "royalblue", "gray55")
shape_values = c(15, 16, 17)

calib_data <- calib_table %>%
  select(Model, Group, Pred, Obs) %>%
  mutate(Model = factor(Model, levels = c("CMUH", "AUH", "Segmed")))

#-------------------------------------------------------------------------------
# 3. Brier score
#-------------------------------------------------------------------------------
brier_score <- function(df, prob_col, event_col) {
  pred  <- df[[prob_col]]
  event <- df[[event_col]]
  mean((pred - event)^2, na.rm = TRUE)
}

brier_results <- map_dfr(names(datasets), function(dataset_name) {
  tibble(
    Model      = dataset_name,
    BrierScore = round(brier_score(datasets[[dataset_name]], prob_col, event_col), 4)
  )
})

group_label_bs <- brier_results %>%
  mutate(label = paste0("**", Model, "**", " (BS = ", BrierScore, ")")) %>%
  pull(label)


#-------------------------------------------------------------------------------
# 4. Generate calibration plot
#-------------------------------------------------------------------------------
ggplot(data = calib_data, 
       aes(x = Pred, y = Obs, group = Model, colour = Model, shape = Model)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2, color = "gray55") +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 3.5) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  xlab("Predicted Probability") + ylab("Observed Probability") +
  theme_bw() +
  theme(
    text = element_text(family = "Calibri"),
    legend.position       = c(0.78, 0.10),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.title          = element_blank(),
    legend.text       = element_markdown(size = 14),
    legend.key.height = unit(1.4, "lines"),  
    legend.spacing.y  = unit(0.5, "lines"),        
    legend.key.size       = unit(0.9, "lines"),
    legend.key.width      = unit(2.0, "lines"),
    panel.grid.major      = element_blank(),
    panel.grid.minor      = element_blank(),
    axis.title.x          = element_text(size = 14, face = "bold",          
                                         margin = margin(t = 5, unit = "mm")),
    axis.title.y          = element_text(size = 14, face = "bold",          
                                         margin = margin(r = 5, unit = "mm")),
    axis.text             = element_text(size = 13, colour = "black")) +
  scale_colour_manual(values = color_values, labels = group_label_bs) + 
  scale_shape_manual(values  = shape_values, labels = group_label_bs) +
  guides(colour = guide_legend(ncol = 1),
         shape  = guide_legend(ncol = 1))


