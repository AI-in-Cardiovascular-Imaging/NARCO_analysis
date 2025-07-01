library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)
library(readxl)
library(pROC)
library(yardstick)
library(writexl)
library(caret)
library(cvms)
library(psych)
library(irr)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# invasive <- read_excel("C:/WorkingData/Documents/3_Research/IVUS_data.xlsx")

# Load the data
# baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))
baseline <- readRDS("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/baseline.rds")

baseline <- baseline %>% 
  filter(record_id != 109, patient_id != "NARCO_250") %>%
  filter(!is.na(funct_spect_date) | !is.na(funct_pet_date)) %>%
  filter(!is.na(inv_ffrdobu)) %>%
  filter(caa_course___0 == "yes")

baseline[11, ]$funct_spect_date <- NA

# remove duplicates of patient_id
baseline <- baseline %>% distinct(patient_id, .keep_all = TRUE)

baseline_copy <- baseline

intra <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/patient_ids.xlsx")
intra <- intra %>% mutate(
  elliptic_ratio = round(major_axis / minor_axis, 1),
  MLN = round(100 - (MLA / distal_area) * 100, 0),
  MLN = ifelse(MLN < 0, 0, MLN),
  # MLA = round(MLA, 0),
  # OLA = round(OLA, 0),
  elliptic_ratio_ostial = round(ostial_major_axis / ostial_minor_axis, 1),
  MLN_ostial = round(100 - (OLA / distal_area) * 100, 0),
  MLN_ostial = ifelse(MLN_ostial < 0, 0, MLN_ostial),
)
intra <- intra %>% 
  rename(
    "patient_id" = `Study ID`,
    ccta_mla_a = MLA,
    ccta_mla_elliptic = elliptic_ratio,
    ccta_mla_ln = MLN,
    ccta_mla_h = major_axis,
    ccta_mla_w = minor_axis,
    ccta_ostial_a = OLA,
    ccta_ostial_elliptic = elliptic_ratio_ostial,
    ccta_ostial_pn = MLN_ostial,
    ccta_ostial_w = ostial_minor_axis
  )

intra_copy <- intra

intra <- intra %>%
  rename(
    ccta_ostial_w_intra = ccta_ostial_w,
    ccta_ostial_a_intra = ccta_ostial_a,
    ccta_ostial_elliptic_intra = ccta_ostial_elliptic
  ) %>%
  select(patient_id, ccta_ostial_w_intra, ccta_ostial_a_intra, ccta_ostial_elliptic_intra)

baseline <- baseline %>%
  left_join(intra, by = "patient_id") %>%
  mutate(
    ccta_ostial_w = ccta_ostial_w_intra,
    ccta_ostial_a = ccta_ostial_a_intra,
    ccta_ostial_elliptic = ccta_ostial_elliptic_intra
  )

baseline <- baseline %>% select(-ccta_ostial_w_intra, -ccta_ostial_a_intra, -ccta_ostial_elliptic_intra)

test <- baseline %>%
  rename(
    ccta_ostial_w_baseline = ccta_ostial_w,
    ccta_ostial_a_baseline = ccta_ostial_a,
    ccta_ostial_elliptic_baseline = ccta_ostial_elliptic
  ) %>%
  select(patient_id, ccta_ostial_w_baseline, ccta_ostial_a_baseline, ccta_ostial_elliptic_baseline)

intra <- intra_copy %>%
  left_join(test, by = "patient_id") %>%
  mutate(
    ccta_ostial_w = ccta_ostial_w_baseline,
    ccta_ostial_a = ccta_ostial_a_baseline,
    ccta_ostial_elliptic = ccta_ostial_elliptic_baseline
  )

intra <- intra %>% select(-ccta_ostial_w_baseline, -ccta_ostial_a_baseline, -ccta_ostial_elliptic_baseline)

baseline <- baseline %>% mutate(
  funct_pos = case_when(
    funct_pet_mismatchcaa == "yes" | funct_spect_mismatchcaa == "yes" ~ 1,
    funct_pet_mismatchcaa == "no" | funct_spect_mismatchcaa == "no" ~ 0,
    TRUE ~ 0
  ),
)

baseline <- baseline %>%
  mutate(percent_stenosis = cut(ccta_mla_ln, 
                        breaks = c(-Inf, 50, 70, 90, Inf), 
                        labels = c("<50", "50-70", "70-90", ">90")))
# Map the bins to specific sizes
size_values <- c("<50" = 2, "50-70" = 5, "70-90" = 10, ">90" = 20)


baseline <- baseline %>% mutate(
  ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0),
  ffr_0.75 = ifelse(inv_ffrdobu < 0.75, 1, 0),
)

# for all funct_pet_mis_aha___0 to funct_pet_mis_aha___16 unfactor
for (i in 0:16) {
  baseline <- baseline %>% mutate(
    !!paste0("funct_pet_mis_aha___", i) := ifelse(!!sym(paste0("funct_pet_mis_aha___", i)) == "yes", 1, 0)
  )
}

baseline <- baseline %>% mutate(ergo_pos = case_when(
  funct_ergo_findings == "clinical positive, electrical negative" | 
  funct_ergo_findings == "clinical negative, electrical positive" |
  funct_ergo_findings == "clinical and electrical positive" ~ 1,
  funct_ergo_findings == "clinical and electrical negative" ~ 0,
  TRUE ~ NA_real_
),
funct_flowres_diff = (funct_pet_cx_flowres + funct_pet_lad_flowres) / 2 - funct_pet_rca_flowres,
funct_pet_score = rowSums(select(baseline, starts_with("funct_pet_mis_aha___"))),
funct_pet_mismatchcaa = case_when(
  funct_pet_mismatchcaa == "yes" ~ 1,
  funct_pet_mismatchcaa == "no" ~ 0,
  TRUE ~ 0
),
funct_spect_mismatchcaa = case_when(
  funct_spect_mismatchcaa == "yes" ~ 1,
  funct_spect_mismatchcaa == "no" ~ 0,
  TRUE ~ 0
),
funct_spect_hr_diff = funct_spect_maxhr - inv_dobu_hr,
funct_spect_aosys_diff = funct_spect_bp_sys_stress - inv_dobu_aosys,
funct_pet_hr_diff = funct_pet_maxhr - inv_dobu_hr,
funct_pet_aosys_diff = funct_pet_bp_sys_stress - inv_dobu_aosys,
funct_rest_hr = case_when(
  !is.na(funct_spect_date) ~ funct_spect_minhr,
  !is.na(funct_pet_date) ~ funct_pet_minhr,
  TRUE ~ NA_real_
),
funct_rest_aosys = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_sys_rest,
  !is.na(funct_pet_date) ~ funct_pet_bp_sys_rest,
  TRUE ~ NA_real_
),
funct_rest_aodia = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_dia_rest,
  !is.na(funct_pet_date) ~ funct_pet_bp_dia_rest,
  TRUE ~ NA_real_
),
funct_rest_aomean = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_mean_rest,
  !is.na(funct_pet_date) ~ funct_pet_bp_mean_rest,
  TRUE ~ NA_real_
),
funct_stress_hr = case_when(
  !is.na(funct_spect_date) ~ funct_spect_maxhr,
  !is.na(funct_pet_date) ~ funct_pet_maxhr,
  TRUE ~ NA_real_
),
funct_stress_aosys = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_sys_stress,
  !is.na(funct_pet_date) ~ funct_pet_bp_sys_stress,
  TRUE ~ NA_real_
),
funct_stress_aodia = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_dia_stress,
  !is.na(funct_pet_date) ~ funct_pet_bp_dia_stress,
  TRUE ~ NA_real_
),
funct_stress_aomean = case_when(
  !is.na(funct_spect_date) ~ funct_spect_bp_mean_stress,
  !is.na(funct_pet_date) ~ funct_pet_bp_mean_stress,
  TRUE ~ NA_real_
),
funct_stress_hrper = case_when(
  !is.na(funct_spect_date) ~ funct_spect_hrper,
  !is.na(funct_pet_date) ~ funct_pet_hrper,
  TRUE ~ NA_real_
)
)

baseline <- baseline %>% mutate(
  inv_hr_change = inv_dobu_hr - inv_rest_hr,
  inv_aosys_change = inv_dobu_aosys - inv_rest_aosys,
  inv_aodia_change = inv_dobu_aodia - inv_rest_aodia,
  inv_aomean_change = inv_dobu_aomean - inv_rest_aomean,
  funct_hr_change = funct_stress_hr - funct_rest_hr,
  funct_aosys_change = funct_stress_aosys - funct_rest_aosys,
  funct_aodia_change = funct_stress_aodia - funct_rest_aodia,
  funct_aomean_change = funct_stress_aomean - funct_rest_aomean,
  pet_hr_change = funct_pet_maxhr - funct_pet_minhr,
  pet_aosys_change = funct_pet_bp_sys_stress - funct_pet_bp_sys_rest,
  pet_aodia_change = funct_pet_bp_dia_stress - funct_pet_bp_dia_rest,
  pet_aomean_change = funct_pet_bp_mean_stress - funct_pet_bp_mean_rest,
  spect_hr_change = funct_spect_maxhr - funct_spect_minhr, 
  spect_aosys_change = funct_spect_bp_sys_stress - funct_spect_bp_sys_rest,
  spect_aodia_change = funct_spect_bp_dia_stress - funct_spect_bp_dia_rest,
  spect_aomean_change = funct_spect_bp_mean_stress - funct_spect_bp_mean_rest
)

baseline <- baseline %>% mutate(
  ccta_pos = ifelse(ccta_mla_a <= 5.59, 1, 0)
)

mean(baseline$inv_dobu_hr, na.rm = T)
mean(baseline$funct_spect_maxhr, na.rm = T)
mean(baseline$funct_pet_maxhr, na.rm = T)
mean(baseline$inv_dobu_aosys, na.rm = T)
mean(baseline$funct_spect_bp_sys_stress, na.rm = T)
mean(baseline$funct_pet_bp_sys_stress, na.rm = T)

pet <- baseline %>% filter(!is.na(funct_pet_date))
spect <- baseline %>% filter(!is.na(funct_spect_date))

baseline %>% select(patient_id, inv_ffrdobu, ccta_ostial_a, funct_pos, funct_spect_mismatchcaa, funct_pet_mismatchcaa) %>% View()

baseline %>% 
  filter(funct_pos == 0, inv_ffrdobu > 0.8) %>% 
  select(patient_id, inv_ffrdobu, ccta_ostial_w) %>% 
  arrange(desc(ccta_ostial_w))

#########################################################################################################################################
# Functions
plot_ffr_data <- function(baseline, var1, var2, method1_name = "Method 1", method2_name = "Method 2", scale_breaks = seq(0.4, 1, 0.1)) {
  n <- nrow(baseline)
  inv_var1 <- pull(baseline, {{var1}})
  inv_var2 <- pull(baseline, {{var2}})
  method_1 <- rep(method1_name, n)
  method_2 <- rep(method2_name, n)
  ffr_df <- data.frame(method = c(method_1, method_2), 
                       value = c(inv_var1, inv_var2))
  colnames(ffr_df) <- c("method", "value")
  
  plot_ffr_change <- ggpaired(ffr_df, x = "method", y = "value", color = "method", 
                              line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b")) + 
    stat_compare_means(paired = TRUE) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    scale_y_continuous(breaks = scale_breaks)
  
  scatter_ffr <- ggplot(ffr_df, aes(x = method, y = value)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    geom_boxplot(color = c("#193bac", "#ebc90b")) +
    geom_jitter(width = 0.05, color = "#bebebec7") +
    theme_classic() +
    scale_y_continuous(breaks = scale_breaks) + 
    ylab("Value")
  
  return(list(plot_ffr_change, scatter_ffr))
}

plot_with_patho_ffr <- function(baseline, var1, var2, method1_name = "Method 1", method2_name = "Method 2", scale_breaks = seq(0.4, 1, 0.1)) {
  n <- nrow(baseline)
  inv_var1 <- pull(baseline, {{var1}})
  inv_var2 <- pull(baseline, {{var2}})
  inv_ffrado <- pull(baseline, inv_ffrado)
  inv_ffrado <- ifelse(inv_ffrado <= 0.8, 1, 0)
  inv_ffrdobu <- pull(baseline, inv_ffrdobu)
  inv_ffrdobu <- ifelse(inv_ffrdobu <= 0.8, 1, 0)
  method_1 <- rep(method1_name, n)
  method_2 <- rep(method2_name, n)
  ivus_df <- data.frame(method = c(method_1, method_2), 
                       value = c(inv_var1, inv_var2),
                       ffr_value = c(inv_ffrado, inv_ffrdobu))
  ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
  ivus_df$method <- factor(ivus_df$method, levels = c(method1_name, method2_name), ordered = TRUE)
  colnames(ivus_df) <- c("method", "value", "ffr_value")
  
  plot_ffr_change <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#ebc90b", "#193bac", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = scale_breaks)

  ivus_df$ivus_method <- factor(ivus_df$method, levels = c(method1_name, method2_name), ordered = TRUE)
  scatter_imla <- ggplot(ivus_df, aes(x = method, y = value)) +
    geom_boxplot(aes(color = method)) +
    geom_jitter(width = 0.05, color = "#bebebec7") +
    theme_classic() +
    scale_y_continuous(breaks = scale_breaks) + 
    scale_color_manual(values = c(var1 = "#193bac", var2 = "#ebc90b"))
  
  return(list(plot_ffr_change, scatter_imla))
}

create_boxplot <- function(data, cols, y_breaks) {
  # Reshape the data to long format
  df_long <- pivot_longer(data, cols = all_of(cols), 
                          names_to = "variable", values_to = "value")
  
  # Specify the order of variable levels
  df_long$variable <- factor(df_long$variable, levels = cols)
  
  # Plotting
  plot <- ggplot(df_long, aes(x = record_id, y = value, fill = variable)) +
    geom_boxplot() +
    labs(title = paste("Change from", cols[1], "to", cols[2], "to", cols[3]),
         x = "Record ID",
         y = "Value",
         fill = "Variable") +
    scale_fill_manual(values = c("#193bac", "#ebc90b", "red")) +
    scale_y_continuous(breaks = y_breaks) +
    theme_classic()
  
  return(plot)
}

# perform_logistic_regression <- function(response_var, data, funct_vars) {
#   logistic_results <- list()
  
#   for (i in 1:length(funct_vars)) {
#     funct_var <- funct_vars[i]
    
#     p_value <- NA
#     odds_ratio <- NA
#     ci_lower <- NA
#     ci_upper <- NA
    
#     tryCatch({
#       ivus_glm <- glm(as.formula(paste(response_var, "~", funct_var)), family = "binomial", data = data)
      
#       p_value <- round(summary(ivus_glm)$coefficients[2, 4], 3)
#       odds_ratio <- round(exp(coef(ivus_glm)[2]), 2)
#       ci_lower <- round(exp(confint(ivus_glm)[2, 1]), 2)
#       ci_upper <- round(exp(confint(ivus_glm)[2, 2]), 2)
#     }, error = function(e) {
#       # If error occurs, leave p_value and odds_ratio as NA
#     })
    
#     logistic_results[[i]] <- data.frame(
#       Variable = funct_var,
#       P_Value = p_value,
#       Odds_Ratio = odds_ratio,
#       CI_Lower = ci_lower,
#       CI_Upper = ci_upper
#     )
#   }
  
#   logistic_results_df <- do.call(rbind, logistic_results)
  
#   return(logistic_results_df)
# }

perform_logistic_regression <- function(response_var, data, funct_vars) {
  logistic_results <- list()
  
  for (i in 1:length(funct_vars)) {
    funct_var <- funct_vars[i]
    
    p_value <- NA
    odds_ratio <- NA
    ci_lower <- NA
    ci_upper <- NA
    
    tryCatch({
      ivus_glm <- glm(as.formula(paste(response_var, "~", funct_var)), family = "binomial", data = data)
      
      p_value <- summary(ivus_glm)$coefficients[2, 4]
      odds_ratio <- exp(coef(ivus_glm)[2])
      ci_lower <- exp(confint(ivus_glm)[2, 1])
      ci_upper <- exp(confint(ivus_glm)[2, 2])
    }, error = function(e) {
      # If error occurs, leave p_value and odds_ratio as NA
    })
    
    logistic_results[[i]] <- data.frame(
      Variable = funct_var,
      P_Value = p_value,
      Odds_Ratio = odds_ratio,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    )
  }
  
  logistic_results_df <- do.call(rbind, logistic_results)
  
  return(logistic_results_df)
}

perform_roc_analysis <- function(response_var, data, funct_vars) {
  roc_results <- list()
  
  for (i in 1:length(funct_vars)) {
    funct_var <- funct_vars[i]
    
    auc <- NA
    sensitivity <- NA
    specificity <- NA
    threshold <- NA
    
    tryCatch({
      roc_model <- glm(as.formula(paste(response_var, "~", funct_var)), data = data, family = binomial)
      
      roc_pred <- predict(roc_model, data, type = "response")
      
      roc_curve <- roc(data[[response_var]], roc_pred)
      
      auc <- round(auc(roc_curve), 2)
      
      roc_threshold <- roc(data[[response_var]], data[[funct_var]])
      opt_threshold <- coords(roc_threshold, "best", ret = "threshold")
      threshold <- round(opt_threshold, 2)
      
      sens_spec <- coords(roc_curve, "best", ret = c("sensitivity", "specificity"))
      sensitivity <- round(sens_spec["sensitivity"], 2)
      specificity <- round(sens_spec["specificity"], 2)
    }, error = function(e) {
      # If error occurs, leave variables as NA
    })    

    roc_results[[i]] <- data.frame(
      Variable = funct_var,
      AUC = auc,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Threshold = threshold
    )
  }

  roc_results_df <- do.call(rbind, roc_results)
  
  return(roc_results_df)
}

simple_logistic_regression <- function(baseline_data, response, explanatory, family = "binomial") {
    data <- baseline_data %>% 
        select(!!sym(response), !!sym(explanatory)) %>% 
        drop_na()
    formula <- as.formula(paste(response, "~", explanatory))
    mdl <- glm(formula, data = data, family = family)

    max <- max(data[[explanatory]], na.rm = TRUE)
    min <- min(data[[explanatory]], na.rm = TRUE)

    if (min<0) {
        min <- 0
    }

    explanatory_data <- tibble(
        !!sym(explanatory) := seq(min, max, length.out = 1000)
    )
    
    prediction_data <- explanatory_data %>% 
        mutate(
            response = predict(mdl, explanatory_data, type = "response"),
            most_likely_outcome = round(response),
            odds_ratio = response / (1 - response),
            log_odds_ratio = log(odds_ratio),
            log_odds_ratio2 = predict(mdl, explanatory_data)
        )

    # because of bug rename response here
    names(prediction_data)[2] <- response
    
    plot_log <- ggplot(data, aes(x = !!sym(explanatory), y = !!sym(response))) +
        geom_point() + 
        geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
        scale_x_continuous(breaks = seq(round(min, 0), round(max), round(max / 10, 1)))
    
    plot_odds <- ggplot(prediction_data, aes(x = !!sym(explanatory), y = odds_ratio)) + 
        geom_point() + 
        geom_line() +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_x_continuous(breaks = seq(round(min, 0), round(max), round(max / 10, 1))) +
        scale_y_continuous(breaks = seq(0, round(max(prediction_data$odds_ratio)), 1))

    confusion <- NULL
    plot_acc <- NULL
    tryCatch(
        {
            actual_response <- data[[response]]
            predicted_response <- round(fitted(mdl))
            outcome <- table(predicted_response, actual_response)
            confusion <- conf_mat(outcome)
            plot_acc <- autoplot(confusion)
        },
        error = function(e) {
            print("Check model! No binary response")
        }
    )
    return(list(mdl, prediction_data, plot_log, plot_odds, confusion, plot_acc))
}

create_bland_altman_plot <- function(df, col_A, col_B, ylim_lower= -10, ylim_upper=10) {
  
  # Create new column for average measurement
  df$avg <- rowMeans(df[, c(col_A, col_B)])
  
  # Create new column for difference in measurements
  df$diff <- df[[col_A]] - df[[col_B]]
  
  # Calculate mean difference and standard deviation
  mean_diff <- mean(df$diff)
  lower <- mean_diff - 1.96 * sd(df$diff)
  upper <- mean_diff + 1.96 * sd(df$diff)
  
  # Create the Bland-Altman plot
  plot <- ggplot(df, aes(x = avg, y = diff)) +
    geom_point(size = 2) +
    geom_hline(yintercept = mean_diff, color = "blue", linetype = "solid") +
    geom_hline(yintercept = lower, color = "red", linetype = "dashed") +
    geom_hline(yintercept = upper, color = "red", linetype = "dashed") +
    ggtitle("Bland-Altman Plot") +
    ylab("Difference Between Measurements") +
    xlab("Average Measurement") +
    ylim(c(ylim_lower, ylim_upper)) +
    theme_minimal()
  
  # Return the plot
  return(plot)
}

###### Functional Analysis ######################################################################################################
table <- baseline %>%
  select(funct_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_nuclear <- conf_mat(table, truth = ffr_0.8, estimate = funct_pos)
result_nuclear <- summary(confusion_nuclear, event_level = "second")
# confusion matrix color
data <- data.frame(
  "target" = baseline$ffr_0.8,
  "prediction" = baseline$funct_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

# conf_mat_all <- plot_confusion_matrix(eval, add_sums = TRUE, rm_zero_percentages = FALSE, rm_zero_text = FALSE, add_zero_shading = FALSE)
conf_mat_all <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

table <- pet %>%
  select(funct_pet_mismatchcaa, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_pet <- conf_mat(table, truth = ffr_0.8, estimate = funct_pet_mismatchcaa)
result_pet <- summary(confusion_pet, event_level = "second")

data <- data.frame(
  "target" = pet$ffr_0.8,
  "prediction" = pet$funct_pet_mismatchcaa,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_pet <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

table <- spect %>%
  select(funct_spect_mismatchcaa, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_spect <- conf_mat(table, truth = ffr_0.8, estimate = funct_spect_mismatchcaa)
result_spect <- summary(confusion_spect, event_level = "second")

data <- data.frame(
  "target" = spect$ffr_0.8,
  "prediction" = spect$funct_spect_mismatchcaa,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_spect <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

baseline <- baseline %>% mutate(
  ccta_pos = ifelse(ccta_mla_a < 3.91, 1, 0), # 3.91 derived from the ROC curve
  combined_pos = funct_pos + ccta_pos,
  combined_pos = case_when(
    combined_pos == 2 ~ 1,
    combined_pos == 1 ~ 1,
    TRUE ~ 0
  )
)

table <- baseline %>%
  select(combined_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_combined_pos <- conf_mat(table)
result_combined_pos <- summary(confusion_combined_pos, event_level = "second")

data <- data.frame(
  "target" = baseline$ffr_0.8,
  "prediction" = baseline$combined_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_combined_pos <- plot_confusion_matrix(eval, add_sums = TRUE, rm_zero_percentages = FALSE, rm_zero_text = FALSE, add_zero_shading = FALSE)

# Subanalysis 0.75
table <- pet %>%
  select(funct_pet_mismatchcaa, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_pet <- conf_mat(table, truth = ffr_0.8, estimate = funct_pet_mismatchcaa)
result_pet <- summary(confusion_pet, event_level = "second")

data <- data.frame(
  "target" = pet$ffr_0.8,
  "prediction" = pet$funct_pet_mismatchcaa,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_pet <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_matrix_nuclear.png", conf_mat_all, width = 3, height = 3)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_matrix_pet.png", conf_mat_pet, width = 3, height = 3)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_matrix_spect.png", conf_mat_spect, width = 3, height = 3)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_matrix_combined_pos.png", conf_mat_combined_pos, width = 6, height = 5)

# save all results in a excel sheet
results_list <- list(
  "Nuclear" = result_nuclear,
  "PET" = result_pet,
  "SPECT" = result_spect,
  "Combined" = result_combined_pos
)

write_xlsx(results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/functional_analysis_results.xlsx")

# baseline IVUS variables to predict FFR change
# get all variables that start with funct_
funct_vars <- colnames(baseline)[grep("^funct_", colnames(baseline))]

logistic_results_nuclear <- perform_logistic_regression("ffr_0.8", baseline, "funct_pos")
logistic_results_pet <- perform_logistic_regression("ffr_0.8", pet, "funct_pet_mismatchcaa")
logistic_results_spect <- perform_logistic_regression("ffr_0.8", spect, "funct_spect_mismatchcaa")
logistic_results_all <- perform_logistic_regression("ffr_0.8", baseline, funct_vars)

######### ROC analysis for IVUS variables ###############################################################
# Perform ROC analysis for all the specified response variables
mdl_funct_pos <- simple_logistic_regression(baseline, "ffr_0.8", "funct_pos")[[1]]
prediction_data_funct_pos <- simple_logistic_regression(baseline, "ffr_0.8", "funct_pos")[[2]]
data_funct <- baseline %>% select(ffr_0.8, funct_pos) %>% drop_na()
roc_funct <- roc(data_funct$ffr_0.8, mdl_funct_pos$fitted.values)
mdl_pet_pos <- simple_logistic_regression(pet, "ffr_0.8", "funct_pet_mismatchcaa")[[1]]
prediction_data_pet_pos <- simple_logistic_regression(pet, "ffr_0.8", "funct_pet_mismatchcaa")[[2]]
data_pet <- pet %>% select(ffr_0.8, funct_pet_mismatchcaa) %>% drop_na()
roc_pet <- roc(data_pet$ffr_0.8, mdl_pet_pos$fitted.values)
mdl_spect_pos <- simple_logistic_regression(spect, "ffr_0.8", "funct_spect_mismatchcaa")[[1]]
prediction_data_spect_pos <- simple_logistic_regression(spect, "ffr_0.8", "funct_spect_mismatchcaa")[[2]]
data_spect <- spect %>% select(ffr_0.8, funct_spect_mismatchcaa) %>% drop_na()
roc_spect <- roc(data_spect$ffr_0.8, mdl_spect_pos$fitted.values)
ggroc(roc_funct)
ggroc(roc_pet)
ggroc(roc_spect)


###### CCTA Analysis ######################################################################################################
# linear relationships
mla_ffr <- ggplot(baseline, aes(x = ccta_mla_a, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = percent_stenosis, alpha = 0.7, color = ccta_mla_elliptic)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 10, 0.5)) +
  scale_color_gradient(low = "#070775", high = "#ff9100") +  # Gradient color scale
  scale_size_manual(values = size_values) +
  labs(size = "Percent stenosis [%]") +
  xlab("Minimal lumen area [mm²]") +
  ylab("FFR")

baseline <- baseline %>% mutate(
  area_categories = ifelse(ccta_mla_a < 2, "<2", 
                           ifelse(ccta_mla_a < 4, "2-4", 
                                  ifelse(ccta_mla_a < 6, "4-6", 
                                         ifelse(ccta_mla_a < 8, "6-8", ">8"))))
  )

area_sizes <- c("<2" = 0.5, "2-4" = 2, "4-6" = 4, "6-8" = 7, ">8" = 10)

mln_ffr <- ggplot(baseline, aes(x = ccta_mla_ln, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = area_categories, alpha = 0.7, color = ccta_mla_elliptic)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_gradient(low = "#070775", high = "#ff9100") +
  scale_size_manual(values = area_sizes) +
  labs(size = "MLA categories [mm²]") +
  xlab("Maximal lumen narrowing [%]") +
  ylab("FFR")

mla_ellip_ffr <- ggplot(baseline, aes(x = ccta_mla_elliptic, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = area_categories, alpha = 0.7, color = ccta_mla_ln)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 7, 0.5)) +
  scale_color_gradient(low = "#070775", high = "#ff9100") +
  scale_size_manual(values = area_sizes) +
  labs(size = "MLA categories [mm²]") +
  xlab("Minimal lumen elliptic ratio") +
  ylab("FFR")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ccta_mla_ffr.png", mla_ffr, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ccta_mln_ffr.png", mln_ffr, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ccta_mla_ellip_ffr.png", mla_ellip_ffr, width = 6, height = 5)

# Define the variables for comparison
rest_vars <- c("inv_ffrado", 
                "inv_rest_hr", 
                "inv_dobu_hrpred",
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean",
               "funct_rest_hr",
               "funct_stress_hrper",
               "funct_rest_aosys",
               "funct_rest_aodia",
               "funct_rest_aomean")

dobu_vars <- c("inv_ffrdobu", 
                "inv_dobu_hr",
                "inv_dobu_hrpred",
               "inv_dobu_aosys", 
               "inv_dobu_aodia", 
               "inv_dobu_aomean",
               "funct_stress_hr",
               "funct_stress_hrper",
               "funct_stress_aosys",
               "funct_stress_aodia",
               "funct_stress_aomean")

# Initialize a list to store results
results <- list()

# Loop over variables to perform the tests and compute statistics
for (i in 1:length(rest_vars)) {
  rest_var <- rest_vars[i]
  dobu_var <- dobu_vars[i]
  
  # Initialize variables to store test results
  rest_summary <- ""
  dobu_summary <- ""
  test_method <- ""
  p_value <- NA
  
  # Try performing Shapiro-Wilk test
  tryCatch({
    rest_shapiro <- shapiro.test(as.numeric(baseline[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(baseline[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(baseline[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(baseline[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(baseline[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(baseline[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(baseline[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(baseline[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(baseline[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(baseline[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(baseline[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(baseline[[dobu_var]], 0.75, na.rm = TRUE)
      rest_summary <- sprintf("%.2f (%.2f-%.2f)", rest_median, rest_q1, rest_q3)
      dobu_summary <- sprintf("%.2f (%.2f-%.2f)", dobu_median, dobu_q1, dobu_q3)
      test_method <- "Wilcoxon"
    }
    p_value <- test_result$p.value
  }, error = function(e) {
    # Fill in blank line if error occurs
    rest_var <- ""
    dobu_var <- ""
  })
  
  # Store results
  results[[i]] <- data.frame(
    Variable = rest_var,
    Rest_Summary = rest_summary,
    Dobu_Summary = dobu_summary,
    Test = test_method,
    P_Value = round(p_value,3)
  )
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Print the results
print(results_df)

inv_changes <- c("inv_hr_change",
"inv_aosys_change",
"inv_aodia_change",
"inv_aomean_change")

funct_changes <- c("funct_hr_change",
"funct_aosys_change",
"funct_aodia_change",
"funct_aomean_change")

p_values <- data.frame(
  Variable = c("Heart rate", "Aortic systolic pressure", "Aortic diastolic pressure", "Aortic mean pressure"),
  Test = NA,
  P_Value = NA
)

for (i in 1:length(inv_changes)) {
  if (shapiro.test(baseline[[inv_changes[i]]])$p.value > 0.05 & shapiro.test(baseline[[funct_changes[i]]])$p.value > 0.05) {
    p_value <- t.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "t-test"
    p_values$P_Value[i] <- p_value
  } else {
    p_value <- wilcox.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "Wilcoxon"
    p_values$P_Value[i] <- p_value
  }
}

results_list <- list(
  "Functional vs. Invasive Changes" = results_df,
  "P-Values" = p_values
)
# Optionally, save the results to a CSV file
write_xlsx(results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/functional_vs_invasive_changes.xlsx")

# baseline IVUS variables to predict FFR change
funct_vars <- c("ccta_mla_a", 
                "ccta_mla_w",
               "ccta_mla_elliptic",
               "ccta_mla_ln",
               "ccta_ostial_a",
                "ccta_ostial_w",
                "ccta_ostial_elliptic",
                "ccta_ostial_pn", 
               "funct_pos")

logistic_results_ffr_0.8 <- perform_logistic_regression("ffr_0.8", baseline, funct_vars)

results_list <- list(
  "ffr_0.8" = logistic_results_ffr_0.8
)

# Write the list of data frames to an Excel file
write_xlsx(results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/ccta_logistic_regression_results.xlsx")

######### ROC analysis for IVUS variables ###############################################################
# Perform ROC analysis for all the specified response variables
roc_results_ffr_0.8 <- perform_roc_analysis("ffr_0.8", baseline, funct_vars)

# Create a list of data frames to be written to the Excel file
roc_results_list <- list(
  "ffr_0.8" = roc_results_ffr_0.8
)

# Write the list of data frames to an Excel file
write_xlsx(roc_results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/funct_roc_analysis_results.xlsx")

############################################################################################################
# best cutoffs
# Fit models and get prediction data
mdl_mla_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")[[1]]
mdl_mla_w_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_w")[[1]]
mdl_mla_ellip_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_ln")[[1]]
mdl_ostial_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a")[[1]]
mdl_ostial_w_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_w")[[1]]
mdl_ostial_ellip_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")[[1]]
mdl_ostial_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")[[1]]

prediction_data_mla <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")[[2]]
prediction_data_mla_w <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_w")[[2]]
prediction_data_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")[[2]]
prediction_data_mla_ln <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_ln")[[2]]
prediction_data_ostial <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a")[[2]]
prediction_data_ostial_w <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_w")[[2]]
prediction_data_ostial_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")[[2]]
prediction_data_ostial_ln <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")[[2]]

data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_a) %>% drop_na()
roc_mla <- roc(data_ivus$ffr_0.8, mdl_mla_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_w) %>% drop_na()
roc_mla_w <- roc(data_ivus$ffr_0.8, mdl_mla_w_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_elliptic) %>% drop_na()
roc_mla_ellip <- roc(data_ivus$ffr_0.8, mdl_mla_ellip_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_ostial_a) %>% drop_na()
roc_ostial <- roc(data_ivus$ffr_0.8, mdl_ostial_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_ostial_w) %>% drop_na()
roc_ostial_w <- roc(data_ivus$ffr_0.8, mdl_ostial_w_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_ostial_elliptic) %>% drop_na()
roc_ostial_ellip <- roc(data_ivus$ffr_0.8, mdl_ostial_ellip_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_ostial_pn) %>% drop_na()
roc_ostial_ln <- roc(data_ivus$ffr_0.8, mdl_ostial_ln_ffr$fitted.values)

# Define test sets
test_set <- list(
  coords(roc_mla), coords(roc_mla_w), coords(roc_mla_ellip),
  coords(roc_mla_ln),
  coords(roc_ostial), coords(roc_ostial_w), coords(roc_ostial_ellip),
  coords(roc_ostial_ln)
)

prediction_data <- list(
  prediction_data_mla, prediction_data_mla_w, prediction_data_ellip,
  prediction_data_mla_ln,
  prediction_data_ostial, prediction_data_ostial_w, prediction_data_ostial_ellip,
  prediction_data_ostial_ln
)

variables <- c(
  "CCTA MLA", "CCTA MLA minor axis", "CCTA elliptic ratio",
  "CCTA MLA luminal narrowing",
  "CCTA Ostial", "CCTA Ostial minor axis", "CCTA ostial elliptic ratio",
  "CCTA Ostial luminal narrowing"
)

# Extract thresholds
thresholds <- list()

for (i in 1:length(test_set)) {
  coordinates <- test_set[[i]]
  coordinates <- coordinates[-c(1, nrow(coordinates)), ]
  coordinates <- as_tibble(coordinates)
  
  result <- coordinates %>%
    group_by(sensitivity) %>%
    filter(specificity == max(specificity)) %>%
    select(threshold) %>%
    distinct()
  
  threshold_values <- result$threshold
  thresholds[[variables[i]]] <- threshold_values
}

# Generate confusion matrices and calculate metrics
cutoff_tables <- list()

thresholds_testing <- thresholds[[1]]
prediction <- prediction_data[[1]]
name <- names(prediction)[1]
threshold <- thresholds_testing[1]
for (i in 1:length(thresholds_testing)) {
  threshold <- thresholds_testing[i]
  closest_index <- which.min(abs(prediction$ffr_0.8 - threshold))
  value <- prediction[[name]][closest_index]
}

for (i in 1:length(thresholds)) {
  thresholds_testing <- thresholds[[i]]
  prediction <- prediction_data[[i]]
  name <- names(prediction)[1]
  
  table_thresholds <- tibble()
  for (j in 1:length(thresholds_testing)) {
    threshold <- thresholds_testing[j]
    closest_index <- which.min(abs(prediction$ffr_0.8 - threshold))
    value <- prediction[[name]][closest_index]
    
    if (name %in% c("ccta_mla_a", "ccta_mla_w", "ccta_ostial_a", "ccta_ostial_w")) {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) > value, 0, 1)) %>%
        select(pred_label, ffr_0.8) %>%
        drop_na() %>%
        table()
    } else {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) < value, 0, 1)) %>%
        select(pred_label, ffr_0.8) %>%
        drop_na() %>%
        table()
    }
    
    confusion <- conf_mat(table)
    result <- summary(confusion, event_level = "second")
    
    sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
    spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
    ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
    npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
    acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()
    
    df <- tibble(
      value = as.numeric(round(value, 2)),
      sensitivity = as.numeric(round(sens * 100, 0)),
      specificity = as.numeric(round(spec * 100, 0)),
      ppv = as.numeric(round(ppv * 100, 0)),
      npv = as.numeric(round(npv * 100, 0)),
      accuracy = as.numeric(round(acc * 100, 0)),
      number_of_true_negatives = as.numeric(table[1, 1]),
      number_of_false_negatives = as.numeric(table[1, 2]),
      number_of_false_positives = as.numeric(table[2, 1]),
      number_of_true_positives = as.numeric(table[2, 2])
    )
    names(df)[1] <- name
    table_thresholds <- bind_rows(table_thresholds, df)
  }
  cutoff_tables[[variables[i]]] <- table_thresholds
}

# Save all as CSV
path <- "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/"
for (i in 1:length(cutoff_tables)) {
  name <- variables[i]
  table <- cutoff_tables[[i]]
  write_csv(table, paste0(path, name, ".csv"))
}

roc_list_mla <- list(roc_mla, roc_mla_ln, roc_mla_ellip)
ggroc_mla <- ggroc(roc_list_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue", "darkgreen")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC MLN =", round(roc_mla_ln$auc, 2)), color = "darkblue") +
    annotate("text", x = 0.7, y = 0.1, label = paste("AUC elliptic ratio =", round(roc_mla_ellip$auc, 2)), color = "darkgreen")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/roc_ccta.png", ggroc_mla, width = 6, height = 5)

############################################################################################################
# ROC Cross-Validation for FFR 0.8
funct_vars <- c("ccta_mla_a", 
                "ccta_mla_w",
               "ccta_mla_elliptic",
               "ccta_mla_ln",
               "ccta_ostial_a",
                "ccta_ostial_w",
                "ccta_ostial_elliptic",
                "ccta_ostial_pn"
                )

cross <- baseline %>% drop_na(ffr_0.8, all_of(funct_vars))
cross$ffr_0.8 <- factor(cross$ffr_0.8, levels = c(0, 1), labels = c("Class0", "Class1"))

# Create a list of seeds
set.seed(69)
seeds <- sample(1:1000, 100)

results_list <- list()

# Function to train and evaluate models
train_and_evaluate <- function(seed) {
  set.seed(seed)
  
  # Split the data into training (78%) and testing (22%) sets --> 44/12 patients
  train_indices <- createDataPartition(cross$ffr_0.8, p = 0.78, list = FALSE) # automatically stratified
  train_data <- cross[train_indices, ]
  test_data <- cross[-train_indices, ]
  
  # Check if both classes are present in training and test sets
  if (length(unique(train_data$ffr_0.8)) < 2 || length(unique(test_data$ffr_0.8)) < 2) {
    message(sprintf("Seed %d: Training or Test set does not contain both classes.", seed))
    return(NULL)
  }

  # Initialize a list to store results for this seed
  seed_results <- list()
  
  # Train logistic regression models and evaluate ROC for each IVUS variable
  for (funct_var in funct_vars) {
    formula <- as.formula(paste("ffr_0.8 ~", funct_var))
    
    # Model fitting with error handling
    model <- tryCatch({
      glm(formula, data = train_data, family = "binomial")
    }, error = function(e) {
      message(sprintf("Error fitting model for variable %s: %s", funct_var, e))
      return(NULL)
    })

    if (is.null(model)) next

    # Predictions with error handling
    predictions <- tryCatch({
      predict(model, test_data, type = "response")
    }, error = function(e) {
      message(sprintf("Error predicting for variable %s: %s", funct_var, e))
      return(NULL)
    })

    if (is.null(predictions)) next

    # ROC curve calculation with error handling
    roc_curve <- tryCatch({
      roc(test_data$ffr_0.8, predictions)
    }, error = function(e) {
      message(sprintf("Error calculating ROC for variable %s: %s", funct_var, e))
      return(NULL)
    })

    if (is.null(roc_curve)) next

    best_coords <- tryCatch({
      coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
    }, error = function(e) {
      message(sprintf("Error in coords for variable %s: %s", funct_var, e))
      return(data.frame(threshold = NA, sensitivity = NA, specificity = NA))
    })
    
    if (nrow(best_coords) > 1) {
      best_coords <- best_coords[1, ]
    }

    # Calculate the optimal thresholds for sensitivity and specificity
    coords_list <- coords(roc_curve, "all", ret = c("threshold", "sensitivity", "specificity"))
    coords_list <- coords_list[-c(1, nrow(coords_list)), ]
    best_coords_sens <- coords_list[order(coords_list$sensitivity, decreasing = TRUE),][1, ]
    best_coords_spec <- coords_list[order(coords_list$specificity, decreasing = TRUE),][1, ]

    thresholds <- c(as.numeric(best_coords$threshold), as.numeric(best_coords_sens$threshold), as.numeric(best_coords_spec$threshold))
    sensitivities <- c(as.numeric(best_coords$sensitivity), as.numeric(best_coords_sens$sensitivity), as.numeric(best_coords_spec$specificity))
    specificities <- c(as.numeric(best_coords$specificity), as.numeric(best_coords_sens$specificity), as.numeric(best_coords_spec$sensitivity))

    for (i in 1:length(thresholds)) {
      threshold <- thresholds[i]
      sens_other <- sensitivities[i]
      spec_other <- specificities[i]
      
      # Find the threshold closest to the predicted values
      closest_index <- which.min(abs(predictions - threshold))
      value <- test_data[[funct_var]][closest_index]

      table <- if (funct_var %in% c("ccta_mla_a", "ccta_ostial_a", "ccta_mla_w", "ccta_ostial_w")) {
        test_data %>%
          mutate(pred_label = ifelse(!!sym(funct_var) > value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      } else {
        test_data %>%
          mutate(pred_label = ifelse(!!sym(funct_var) < value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      }

      # Confusion matrix with error handling
      if (nrow(table) != 2 || ncol(table) != 2) {
        message(sprintf("Confusion matrix not square for variable %s at seed %d", funct_var, seed))
        next
      }

      confusion <- conf_mat(table)
      result <- summary(confusion, event_level = "second")

      sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
      spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
      ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
      npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
      acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()

      # Naming convention funct_var_1, funct_var_2, funct_var_3
      seed_results[[paste(funct_var, i, sep = "_")]] <- list(
        AUC = as.numeric(auc(roc_curve)),
        Threshold = as.numeric(value),
        Sensitivity = as.numeric(sens),
        Specificity = as.numeric(spec),
        PPV = as.numeric(ppv),
        NPV = as.numeric(npv),
        Accuracy = as.numeric(acc),
        TrueNegatives = as.numeric(table[1, 1]),
        FalseNegatives = as.numeric(table[1, 2]),
        FalsePositives = as.numeric(table[2, 1]),
        TruePositives = as.numeric(table[2, 2]),
        Sensitivity_other = as.numeric(sens_other),
        Specificity_other = as.numeric(spec_other)
      )
    }
  }
  
  return(seed_results)
}

# Run the training and evaluation for each seed and store results
for (seed in seeds) {
  seed_results <- tryCatch({
    train_and_evaluate(seed)
  }, error = function(e) {
    message(sprintf("Error in seed %d: %s", seed, e))
    return(NULL)
  })

  if (!is.null(seed_results)) {
    results_list[[as.character(seed)]] <- seed_results
  }
}

# Aggregate the results across different seeds
aggregate_results <- function(results_list) {
  aggregate_list <- list()
  
  for (funct_var in funct_vars) {
    for (i in 1:3) {
      var_name <- paste(funct_var, i, sep = "_")
      aucs <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$AUC else NA)
      thresholds <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Threshold else NA)
      sensitivities <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Sensitivity else NA)
      specificities <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Specificity else NA)
      ppvs <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$PPV else NA)
      npvs <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$NPV else NA)
      accuracies <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Accuracy else NA)
      true_negatives <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$TrueNegatives else NA)
      false_negatives <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$FalseNegatives else NA)
      false_positives <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$FalsePositives else NA)
      true_positives <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$TruePositives else NA)
      sensitivity_others <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Sensitivity_other else NA)
      specificity_others <- sapply(results_list, function(x) if (!is.null(x[[var_name]])) x[[var_name]]$Specificity_other else NA)
      
      aggregate_list[[var_name]] <- data.frame(
        Seed = seeds,
        AUC = aucs,
        Threshold = thresholds,
        Sensitivity = sensitivities,
        Specificity = specificities,
        PPV = ppvs,
        NPV = npvs,
        Accuracy = accuracies,
        TrueNegatives = true_negatives,
        FalseNegatives = false_negatives,
        FalsePositives = false_positives,
        TruePositives = true_positives,
        Sensitivity_other = sensitivity_others,
        Specificity_other = specificity_others
      )
    }
  }
  
  return(aggregate_list)
}

final_results <- aggregate_results(results_list)

# Convert list of lists to a data frame
roc_summary_df_ffr0.8 <- do.call(rbind, lapply(names(final_results), function(var_name) {
  res <- final_results[[var_name]]
  data.frame(
    Variable = var_name,
    Mean_AUC = mean(res$AUC, na.rm = TRUE),
    SD_AUC = sd(res$AUC, na.rm = TRUE),
    Mean_Threshold = mean(res$Threshold, na.rm = TRUE),
    SD_Threshold = sd(res$Threshold, na.rm = TRUE),
    Mean_Sensitivity = mean(res$Sensitivity, na.rm = TRUE),
    SD_Sensitivity = sd(res$Sensitivity, na.rm = TRUE),
    Mean_Specificity = mean(res$Specificity, na.rm = TRUE),
    SD_Specificity = sd(res$Specificity, na.rm = TRUE),
    Mean_PPV = mean(res$PPV, na.rm = TRUE),
    SD_PPV = sd(res$PPV, na.rm = TRUE),
    Mean_NPV = mean(res$NPV, na.rm = TRUE),
    SD_NPV = sd(res$NPV, na.rm = TRUE),
    Mean_Accuracy = mean(res$Accuracy, na.rm = TRUE),
    SD_Accuracy = sd(res$Accuracy, na.rm = TRUE),
    Mean_TrueNegatives = mean(res$TrueNegatives, na.rm = TRUE),
    SD_TrueNegatives = sd(res$TrueNegatives, na.rm = TRUE),
    Mean_FalseNegatives = mean(res$FalseNegatives, na.rm = TRUE),
    SD_FalseNegatives = sd(res$FalseNegatives, na.rm = TRUE),
    Mean_FalsePositives = mean(res$FalsePositives, na.rm = TRUE),
    SD_FalsePositives = sd(res$FalsePositives, na.rm = TRUE),
    Mean_TruePositives = mean(res$TruePositives, na.rm = TRUE),
    SD_TruePositives = sd(res$TruePositives, na.rm = TRUE),
    Mean_Sensitivity_other = mean(res$Sensitivity_other, na.rm = TRUE),
    SD_Sensitivity_other = sd(res$Sensitivity_other, na.rm = TRUE),
    Mean_Specificity_other = mean(res$Specificity_other, na.rm = TRUE),
    SD_Specificity_other = sd(res$Specificity_other, na.rm = TRUE)
  )
}))

# Create a new summary with rounded values for better readability
roc_summary_df_ffr0.8_rounded <- roc_summary_df_ffr0.8 %>%
  mutate(
    AUC = paste(round(Mean_AUC, 2), "±", round(SD_AUC, 2)),
    Threshold = paste(round(Mean_Threshold, 2), "±", round(SD_Threshold, 2)),
    Sensitivity = paste(round(Mean_Sensitivity * 100, 0), "±", round(SD_Sensitivity * 100, 0)),
    Specificity = paste(round(Mean_Specificity * 100, 0), "±", round(SD_Specificity * 100, 0)),
    PPV = paste(round(Mean_PPV * 100, 0), "±", round(SD_PPV * 100, 0)),
    NPV = paste(round(Mean_NPV * 100, 0), "±", round(SD_NPV * 100, 0)),
    Accuracy = paste(round(Mean_Accuracy * 100, 0), "±", round(SD_Accuracy * 100, 0)),
    TrueNegatives = paste(round(Mean_TrueNegatives, 0), "±", round(SD_TrueNegatives, 0)),
    FalseNegatives = paste(round(Mean_FalseNegatives, 0), "±", round(SD_FalseNegatives, 0)),
    FalsePositives = paste(round(Mean_FalsePositives, 0), "±", round(SD_FalsePositives, 0)),
    TruePositives = paste(round(Mean_TruePositives, 0), "±", round(SD_TruePositives, 0)),
    Sensitivity_other = paste(round(Mean_Sensitivity_other * 100, 0), "±", round(SD_Sensitivity_other * 100, 0)),
    Specificity_other = paste(round(Mean_Specificity_other * 100, 0), "±", round(SD_Specificity_other * 100, 0))
  ) %>%
  select(Variable, AUC, Threshold, Sensitivity, Specificity, PPV, NPV, Accuracy, TrueNegatives, FalseNegatives, FalsePositives, TruePositives, Sensitivity_other, Specificity_other)

# Print the rounded summary dataframe
print(roc_summary_df_ffr0.8_rounded)

# Save the rounded summary dataframe to an Excel file
write_xlsx(roc_summary_df_ffr0.8_rounded, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/roc_summary_ffr0.8_rounded_summary.xlsx")

############################################################################################################
# forest plot
# for hemodynamic and anatomic relevance as predictor
forest_df_bsa <- data.frame(
    var = c(
            "Minimal lumen area",
            "Minimal lumen elliptic ratio",
            "Maximal lumen narrowing [%]"
            ),
    odds_ratio = c(
            exp(mdl_mla_ffr$coefficients[2]),
            exp(mdl_mla_ellip_ffr$coefficients[2]),
            exp(mdl_mla_ln_ffr$coefficients[2])),
    lower_ci = c(
            exp(confint(mdl_mla_ffr))[2],
            exp(confint(mdl_mla_ellip_ffr))[2],
            exp(confint(mdl_mla_ln_ffr))[2]), 
    upper_ci = c(
            exp(confint(mdl_mla_ffr))[4],
            exp(confint(mdl_mla_ellip_ffr))[4],
            exp(confint(mdl_mla_ln_ffr))[4])
)

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_ffr_0.8 <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 6, 0.5), labels = parse(text = seq(0, 6, 0.5))) +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/figures/forest_ffr_0.8.png", forest_ffr_0.8, width = 5, height = 2)

graphical_mla <- ggroc(roc_mla, legacy.axes = TRUE, color = "#23AD13") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen area [mm²]") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC =", round(roc_mla$auc, 2)), color = "#23AD13")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/graphical_mla.png", graphical_mla, width = 2.5, height = 2.5)

# get all columns starting with caa
caa_vars <- names(baseline)[grep("^caa", names(baseline))]

treatment <- baseline %>% filter(synp_surgery_yn == "yes")
non_treatment <- baseline %>% anti_join(treatment, by = "record_id") %>% filter(inv_ffrdobu <= 0.8)

shapiro.test(treatment$inv_ffrdobu)
median(treatment$inv_ffrdobu)
min(treatment$inv_ffrdobu)
max(treatment$inv_ffrdobu)
shapiro.test(non_treatment$inv_ffrdobu)
median(non_treatment$inv_ffrdobu)
min(non_treatment$inv_ffrdobu)
max(non_treatment$inv_ffrdobu)

wilcox.test(treatment$inv_ffrdobu, non_treatment$inv_ffrdobu)

rest_vars <- c("inv_rest_hr", 
                "inv_dobu_hrpred",
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean",
               "funct_pet_minhr",
               "funct_pet_hrper",
               "funct_pet_bp_sys_rest",
               "funct_pet_bp_dia_rest",
               "funct_pet_bp_mean_rest")

dobu_vars <- c("inv_dobu_hr",
                "inv_dobu_hrpred",
               "inv_dobu_aosys", 
               "inv_dobu_aodia", 
               "inv_dobu_aomean",
               "funct_pet_maxhr",
               "funct_pet_hrper",
               "funct_pet_bp_sys_stress",
               "funct_pet_bp_dia_stress",
               "funct_pet_bp_mean_stress")

# Initialize a list to store results
results <- list()

# Loop over variables to perform the tests and compute statistics
for (i in 1:length(rest_vars)) {
  rest_var <- rest_vars[i]
  dobu_var <- dobu_vars[i]
  
  # Initialize variables to store test results
  rest_summary <- ""
  dobu_summary <- ""
  test_method <- ""
  p_value <- NA
  
  # Try performing Shapiro-Wilk test
  tryCatch({
    rest_shapiro <- shapiro.test(as.numeric(baseline[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(baseline[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(baseline[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(baseline[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(baseline[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(baseline[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(baseline[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(baseline[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(baseline[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(baseline[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(baseline[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(baseline[[dobu_var]], 0.75, na.rm = TRUE)
      rest_summary <- sprintf("%.2f (%.2f-%.2f)", rest_median, rest_q1, rest_q3)
      dobu_summary <- sprintf("%.2f (%.2f-%.2f)", dobu_median, dobu_q1, dobu_q3)
      test_method <- "Wilcoxon"
    }
    p_value <- test_result$p.value
  }, error = function(e) {
    # Fill in blank line if error occurs
    rest_var <- ""
    dobu_var <- ""
  })
  
  # Store results
  results[[i]] <- data.frame(
    Variable = rest_var,
    Rest_Summary = rest_summary,
    Dobu_Summary = dobu_summary,
    Test = test_method,
    P_Value = round(p_value,3)
  )
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Print the results
print(results_df)

inv_changes <- c("inv_hr_change",
"inv_aosys_change",
"inv_aodia_change",
"inv_aomean_change")

funct_changes <- c("pet_hr_change",
"pet_aosys_change",
"pet_aodia_change",
"pet_aomean_change")

p_values <- data.frame(
  Variable = c("Heart rate", "Aortic systolic pressure", "Aortic diastolic pressure", "Aortic mean pressure"),
  Test = NA,
  P_Value = NA
)

for (i in 1:length(inv_changes)) {
  if (shapiro.test(baseline[[inv_changes[i]]])$p.value > 0.05 & shapiro.test(baseline[[funct_changes[i]]])$p.value > 0.05) {
    p_value <- t.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "t-test"
    p_values$P_Value[i] <- p_value
  } else {
    p_value <- wilcox.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "Wilcoxon"
    p_values$P_Value[i] <- p_value
  }
}

results_list <- list(
  "Functional vs. Invasive Changes" = results_df,
  "P-Values" = p_values
)
# Optionally, save the results to a CSV file
write_xlsx(results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/pet_vs_invasive_changes.xlsx")

test <- baseline %>% filter(!is.na(funct_pet_hrper))
chisq.test(test$inv_dobu_hrpred, test$funct_pet_hrper)

rest_vars <- c("inv_ffrado", 
                "inv_rest_hr", 
                "inv_dobu_hrpred",
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean",
               "funct_spect_minhr",
               "funct_spect_hrper",
               "funct_spect_bp_sys_rest",
               "funct_spect_bp_dia_rest",
               "funct_spect_bp_mean_rest")

dobu_vars <- c("inv_ffrdobu", 
                "inv_dobu_hr",
                "inv_dobu_hrpred",
               "inv_dobu_aosys", 
               "inv_dobu_aodia", 
               "inv_dobu_aomean",
               "funct_spect_maxhr",
               "funct_spect_hrper",
               "funct_spect_bp_sys_stress",
               "funct_spect_bp_dia_stress",
               "funct_spect_bp_mean_stress")

# Initialize a list to store results
results <- list()

# Loop over variables to perform the tests and compute statistics
for (i in 1:length(rest_vars)) {
  rest_var <- rest_vars[i]
  dobu_var <- dobu_vars[i]
  
  # Initialize variables to store test results
  rest_summary <- ""
  dobu_summary <- ""
  test_method <- ""
  p_value <- NA
  
  # Try performing Shapiro-Wilk test
  tryCatch({
    rest_shapiro <- shapiro.test(as.numeric(baseline[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(baseline[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(baseline[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(baseline[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(baseline[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(baseline[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(baseline[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(baseline[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(baseline[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(baseline[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(baseline[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(baseline[[dobu_var]], 0.75, na.rm = TRUE)
      rest_summary <- sprintf("%.2f (%.2f-%.2f)", rest_median, rest_q1, rest_q3)
      dobu_summary <- sprintf("%.2f (%.2f-%.2f)", dobu_median, dobu_q1, dobu_q3)
      test_method <- "Wilcoxon"
    }
    p_value <- test_result$p.value
  }, error = function(e) {
    # Fill in blank line if error occurs
    rest_var <- ""
    dobu_var <- ""
  })
  
  # Store results
  results[[i]] <- data.frame(
    Variable = rest_var,
    Rest_Summary = rest_summary,
    Dobu_Summary = dobu_summary,
    Test = test_method,
    P_Value = round(p_value,3)
  )
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Print the results
print(results_df)

inv_changes <- c("inv_hr_change",
"inv_aosys_change",
"inv_aodia_change",
"inv_aomean_change")

funct_changes <- c("spect_hr_change",
"spect_aosys_change",
"spect_aodia_change",
"spect_aomean_change")

p_values <- data.frame(
  Variable = c("Heart rate", "Aortic systolic pressure", "Aortic diastolic pressure", "Aortic mean pressure"),
  Test = NA,
  P_Value = NA
)

for (i in 1:length(inv_changes)) {
  if (shapiro.test(baseline[[inv_changes[i]]])$p.value > 0.05 & shapiro.test(baseline[[funct_changes[i]]])$p.value > 0.05) {
    p_value <- t.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "t-test"
    p_values$P_Value[i] <- p_value
  } else {
    p_value <- wilcox.test(baseline[[inv_changes[i]]], baseline[[funct_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "Wilcoxon"
    p_values$P_Value[i] <- p_value
  }
}

results_list <- list(
  "Functional vs. Invasive Changes" = results_df,
  "P-Values" = p_values
)
# Optionally, save the results to a CSV file
write_xlsx(results_list, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Marius_data/spect_vs_invasive_changes.xlsx")

test <- baseline %>% filter(!is.na(funct_spect_hrper))
chisq.test(test$inv_dobu_hrpred, test$funct_spect_hrper)

roc_list_mla <- list(roc_mla, roc_funct)
ggroc_mla_funct <- ggroc(roc_list_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC nuclear imaging =", round(roc_funct$auc, 2)), color = "darkblue")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/roc_mla_funct.png", ggroc_mla_funct, width = 6, height = 5)

ggroc_mla <- ggroc(roc_mla, legacy.axes = TRUE, color = "darkred") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkred")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/roc_mla.png", ggroc_mla, width = 6, height = 5)

############################################################################################################
table <- baseline %>%
  select(ccta_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_ccta <- conf_mat(table)
result_ccta <- summary(confusion_ccta, event_level = "second")

data <- data.frame(
  "target" = baseline$ffr_0.8,
  "prediction" = baseline$ccta_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_ccta_all <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

table <- pet %>%
  select(ccta_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_ccta <- conf_mat(table)
result_ccta <- summary(confusion_ccta, event_level = "second")

data <- data.frame(
  "target" = pet$ffr_0.8,
  "prediction" = pet$ccta_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_ccta_pet <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

table <- spect %>%
  select(ccta_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_ccta <- conf_mat(table)
result_ccta <- summary(confusion_ccta, event_level = "second")

data <- data.frame(
  "target" = spect$ffr_0.8,
  "prediction" = spect$ccta_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_ccta_spect <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_ccta_all.png", conf_mat_ccta_all, width = 3, height = 3)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_ccta_pet.png", conf_mat_ccta_pet, width = 3, height = 3)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_ccta_spect.png", conf_mat_ccta_spect, width = 3, height = 3)


############################################################################################################
plot_ffr_funct <- ggplot(baseline, aes(x = inv_ffrdobu, y = ccta_mla_a)) +
  geom_smooth(method = "gam", formula = y ~ poly(x, 3), se = FALSE, linetype = "dashed", aes(color = "#d8d6d6"), alpha=0.5) +
  geom_point(aes(color = as.factor(funct_pos))) +
  scale_color_manual(values = c("0" = "darkblue", "1" = "orange")) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  annotate("text", x = 0.6, y = 5.5, label = "MLA >5.6 mm²", color = "darkgreen", vjust = -1) +
  geom_hline(yintercept = 5.6, linetype = "dashed", color = "darkgreen") +
  labs(x = "Fractional Flow Reserve", y = "Minimal lumen area [mm²]") +
  theme_classic()

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/plot_ffr_funct.png", plot_ffr_funct, width = 5, height = 4)

# load patient list
patient_ids <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/Aktuelle_Patientenliste_NARCO_Retro_und_Prosp.xlsx")
patient_ids <- patient_ids %>% select(`Study ID`, Anrede, Name, Vorname, Birthdate, PID)
patient_ids$`Study ID` <- paste0("NARCO_", patient_ids$`Study ID`)

# keep only the columns in patient_ids where Study ID matches the cells in patient_id in baseline
patient_ids <- patient_ids %>% filter(`Study ID` %in% baseline$patient_id)
patient_ids <- patient_ids %>%
  mutate(Birthdate = as.Date(as.numeric(Birthdate), origin = "1899-12-30"))
# save as csv
# write_csv(patient_ids, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/patient_ids.csv")

############################################################################################################
###################################### ICC analysis ####################################################### 
############################################################################################################
baseline <- baseline_copy %>% 
  mutate(
        # ccta_mla_a = round(ccta_mla_a, 0),
        # ccta_ostial_a = round(ccta_ostial_a, 0), 
        ccta_mla_elliptic = round(ccta_mla_elliptic, 1),
        ccta_ostial_elliptic = round(ccta_ostial_elliptic, 1),
        ccta_mla_ln = round(ccta_mla_ln, 0),
        ccta_ostial_pn = round(ccta_ostial_pn, 0))

intra <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/patient_ids.xlsx")
intra <- intra %>% mutate(
  elliptic_ratio = round(major_axis / minor_axis, 1),
  MLN = round(100 - (MLA / distal_area) * 100, 0),
  MLN = ifelse(MLN < 0, 0, MLN),
  # MLA = round(MLA, 0),
  # OLA = round(OLA, 0),
  elliptic_ratio_ostial = round(ostial_major_axis / ostial_minor_axis, 1),
  MLN_ostial = round(100 - (OLA / distal_area) * 100, 0),
  MLN_ostial = ifelse(MLN_ostial < 0, 0, MLN_ostial),
)
intra <- intra %>% 
  rename(
    "patient_id" = `Study ID`,
    ccta_mla_a = MLA,
    ccta_mla_elliptic = elliptic_ratio,
    ccta_mla_ln = MLN,
    ccta_mla_h = major_axis,
    ccta_mla_w = minor_axis,
    ccta_ostial_a = OLA,
    ccta_ostial_elliptic = elliptic_ratio_ostial,
    ccta_ostial_pn = MLN_ostial,
    ccta_ostial_w = ostial_minor_axis
  )

inter <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/SPECTPETpaperCCTAanalysisMRB.xlsx")
inter <- inter %>% rename(
  "patient_id" = `Study ID`,
  ccta_ostial_a = `MLA ostium mm2`,
  ccta_ostial_h = `Major ostium cm`,
  ccta_ostial_w = `Minor ostium cm`,
  ccta_mla_a = `MLA mm2`,
  ccta_mla_h = `Major MLA cm`,
  ccta_mla_w = `Minor MLA cm`,
  ccta_distal_a = `Distal reference mm2`,
)

inter_og <- inter %>% mutate(
  ccta_mla_elliptic = round(as.numeric(ccta_mla_h) / as.numeric(ccta_mla_w), 1),
  ccta_mla_ln = round(100 - (as.numeric(ccta_mla_a) / as.numeric(ccta_distal_a)) * 100, 0),
  ccta_mla_ln = ifelse(ccta_mla_ln < 0, 0, ccta_mla_ln),
  ccta_mla_h = as.numeric(ccta_mla_h) * 10,
  ccta_mla_w = as.numeric(ccta_mla_w) * 10,
  ccta_ostial_elliptic = round(as.numeric(ccta_ostial_h) / as.numeric(ccta_ostial_w), 1),
  ccta_ostial_pn = round(100 - (as.numeric(ccta_ostial_a) / as.numeric(ccta_distal_a)) * 100, 0),
  ccta_ostial_pn = ifelse(ccta_ostial_pn < 0, 0, ccta_ostial_pn),
  ccta_ostial_h = as.numeric(ccta_ostial_h) * 10,
  ccta_ostial_w = as.numeric(ccta_ostial_w) * 10,
  # ccta_mla_a = round(ccta_mla_a, 0),
  # ccta_ostial_a = round(ccta_ostial_a, 0),
)

intra <- intra %>% select(patient_id, ccta_mla_a, ccta_mla_elliptic, ccta_mla_ln, ccta_mla_h, ccta_mla_w, ccta_ostial_a, ccta_ostial_elliptic, ccta_ostial_pn, ccta_ostial_w)
inter <- inter_og %>% select(patient_id, ccta_mla_a, ccta_mla_elliptic, ccta_mla_ln, ccta_mla_h, ccta_mla_w, ccta_ostial_a, ccta_ostial_elliptic, ccta_ostial_pn, ccta_ostial_w)

baseline_icc <- baseline %>% 
  select(patient_id, ccta_mla_a, ccta_mla_elliptic, ccta_mla_ln, ccta_mla_h, ccta_mla_w, ccta_ostial_a, ccta_ostial_elliptic, ccta_ostial_pn, ccta_ostial_w, ccta_quali, ccta_photon) %>%
  left_join(intra, by = "patient_id")

baseline_icc <- baseline_icc %>%
  left_join(inter, by = "patient_id")

baseline_icc <- baseline_icc %>% drop_na()

baseline_icc$ccta_quali <- as.numeric(as.character(baseline_icc$ccta_quali))
baseline_icc <- baseline_icc %>% mutate(ccta_quali = ifelse(ccta_photon == "yes", ccta_quali + 1, ccta_quali))

# copy every row n-times with n being ccta_quali
baseline_icc_no_weight <- baseline_icc
baseline_icc <- baseline_icc %>% 
  slice(rep(row_number(), baseline_icc$ccta_quali * 2)) %>%
  select(-ccta_quali, -ccta_photon)

area_intra <- ICC(cbind(baseline_icc$ccta_mla_a.x, baseline_icc$ccta_mla_a.y))$results
area_inter <- ICC(cbind(baseline_icc$ccta_mla_a.x, baseline_icc$ccta_mla_a))$results
icc(cbind(baseline_icc$ccta_mla_a.x, baseline_icc$ccta_mla_a.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc$ccta_mla_a.x, baseline_icc$ccta_mla_a.y, baseline_icc_no_weight$ccta_mla_a), model = "twoway", type = "agreement", unit = "single")

area <- baseline_icc %>% select(patient_id, ccta_mla_a.x, ccta_mla_a.y, ccta_mla_a) %>% drop_na()
ba_area_intra <- create_bland_altman_plot(area, "ccta_mla_a.x", "ccta_mla_a.y", ylim_lower = -6, ylim_upper = 6)
ba_area_inter <- create_bland_altman_plot(area, "ccta_mla_a.x", "ccta_mla_a", ylim_lower = -6, ylim_upper = 6)

elliptic_intra <- ICC(cbind(baseline_icc$ccta_mla_elliptic.x, baseline_icc$ccta_mla_elliptic.y))$results
elliptic_inter <- ICC(cbind(baseline_icc$ccta_mla_elliptic.x, baseline_icc$ccta_mla_elliptic))$results
icc(cbind(baseline_icc$ccta_mla_elliptic.x, baseline_icc$ccta_mla_elliptic.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_no_weight$ccta_mla_elliptic.x, baseline_icc_no_weight$ccta_mla_elliptic.y), model = "twoway", type = "agreement", unit = "single")

elliptic <- baseline_icc %>% select(patient_id, ccta_mla_elliptic.x, ccta_mla_elliptic.y, ccta_mla_elliptic) %>% drop_na()
ba_elliptic_intra <- create_bland_altman_plot(elliptic, "ccta_mla_elliptic.x", "ccta_mla_elliptic.y", ylim_lower = -2.5, ylim_upper = 2.5)
ba_elliptic_inter <- create_bland_altman_plot(elliptic, "ccta_mla_elliptic.x", "ccta_mla_elliptic", ylim_lower = -2.5, ylim_upper = 2.5)

mln_intra <- ICC(cbind(baseline_icc$ccta_mla_ln.x, baseline_icc$ccta_mla_ln.y))$results
mln_inter <- ICC(cbind(baseline_icc$ccta_mla_ln.x, baseline_icc$ccta_mla_ln))$results
icc(cbind(baseline_icc$ccta_mla_ln.x, baseline_icc$ccta_mla_ln.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_no_weight$ccta_mla_ln.x, baseline_icc_no_weight$ccta_mla_ln.y), model = "twoway", type = "agreement", unit = "single")

ln <- baseline_icc %>% select(patient_id, ccta_mla_ln.x, ccta_mla_ln.y, ccta_mla_ln) %>% drop_na()
ba_ln_intra <- create_bland_altman_plot(ln, "ccta_mla_ln.x", "ccta_mla_ln.y", ylim_lower = -40, ylim_upper = 40)
ba_ln_inter <- create_bland_altman_plot(ln, "ccta_mla_ln.x", "ccta_mla_ln", ylim_lower = -40, ylim_upper = 40)

minor_intra <- ICC(cbind(baseline_icc$ccta_mla_w.x, baseline_icc$ccta_mla_w.y))$results
minor_inter <- ICC(cbind(baseline_icc$ccta_mla_w.x, baseline_icc$ccta_mla_w))$results
icc(cbind(baseline_icc$ccta_mla_w.x, baseline_icc$ccta_mla_w.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_no_weight$ccta_mla_w.x, baseline_icc_no_weight$ccta_mla_w.y), model = "twoway", type = "agreement", unit = "single")

w <- baseline_icc %>% select(patient_id, ccta_mla_w.x, ccta_mla_w.y, ccta_mla_w) %>% drop_na()
ba_minor_intra <- create_bland_altman_plot(w, "ccta_mla_w.x", "ccta_mla_w.y", -1.5, 1.5)
ba_minor_inter <- create_bland_altman_plot(w, "ccta_mla_w.x", "ccta_mla_w", -1.5, 1.5)

area_ostial_intra <- ICC(cbind(baseline_icc$ccta_ostial_a.x, baseline_icc$ccta_ostial_a.y))$results
ellitpic_ostial_intra <- ICC(cbind(baseline_icc$ccta_ostial_elliptic.x, baseline_icc$ccta_ostial_elliptic.y))$results
mln_ostial_intra <- ICC(cbind(baseline_icc$ccta_ostial_pn.x, baseline_icc$ccta_ostial_pn.y))$results
w_ostial_intra <- ICC(cbind(baseline_icc$ccta_ostial_w.x, baseline_icc$ccta_ostial_w.y))$results

area_ostial <- baseline_icc %>% select(patient_id, ccta_ostial_a.x, ccta_ostial_a.y, ccta_ostial_a) %>% drop_na()
ba_area_ostial_intra <- create_bland_altman_plot(area_ostial, "ccta_ostial_a.x", "ccta_ostial_a.y", ylim_lower = -6, ylim_upper = 6)
ellitpic_ostial <- baseline_icc %>% select(patient_id, ccta_ostial_elliptic.x, ccta_ostial_elliptic.y, ccta_ostial_elliptic) %>% drop_na()
mln_ostial <- baseline_icc %>% select(patient_id, ccta_ostial_pn.x, ccta_ostial_pn.y, ccta_ostial_pn) %>% drop_na()
ba_elliptic_ostial_intra <- create_bland_altman_plot(ellitpic_ostial, "ccta_ostial_elliptic.x", "ccta_ostial_elliptic.y", ylim_lower = -2.5, ylim_upper = 2.5)
ba_mln_ostial_intra <- create_bland_altman_plot(mln_ostial, "ccta_ostial_pn.x", "ccta_ostial_pn.y", ylim_lower = -40, ylim_upper = 40)
w_ostial <- baseline_icc %>% select(patient_id, ccta_ostial_w.x, ccta_ostial_w.y, ccta_ostial_w) %>% drop_na()
ba_minor_ostial_intra <- create_bland_altman_plot(w_ostial, "ccta_ostial_w.x", "ccta_ostial_w.y", -1.5, 1.5)


ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_area_intra.png", ba_area_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_elliptic_intra.png", ba_elliptic_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ln_intra.png", ba_ln_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_minor_intra.png", ba_minor_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_area_inter.png", ba_area_inter, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_elliptic_inter.png", ba_elliptic_inter, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ln_inter.png", ba_ln_inter, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_minor_inter.png", ba_minor_inter, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_area_ostial_intra.png", ba_area_ostial_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_elliptic_ostial_intra.png", ba_elliptic_ostial_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_mln_ostial_intra.png", ba_mln_ostial_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_minor_ostial_intra.png", ba_minor_ostial_intra, width = 6, height = 5)

intra_data_frames <- list(
  "Area" = area_intra,
  "Elliptic" = elliptic_intra,
  "MLN" = mln_intra,
  "Minor" = minor_intra,
  "Area_Ostial" = area_ostial_intra,
  "Elliptic_Ostial" = ellitpic_ostial_intra, 
  "MLN_Ostial" = mln_ostial_intra,
  "Minor_Ostial" = w_ostial_intra
)

inter_data_frames <- list(
  "Area" = area_inter,
  "Elliptic" = elliptic_inter,
  "MLN" = mln_inter,
  "Minor" = minor_inter,
)
# get index as column
intra_data_frames <- lapply(intra_data_frames, function(x) {
  x$Variable <- rownames(x)
  rownames(x) <- NULL
  x
})

inter_data_frames <- lapply(inter_data_frames, function(x) {
  x$Variable <- rownames(x)
  rownames(x) <- NULL
  x
})

# round all values to 3 decimals if numeric
intra_data_frames <- lapply(intra_data_frames, function(x) {
  x[] <- lapply(x, function(y) if(is.numeric(y)) round(y, 3) else y)
  x
})

inter_data_frames <- lapply(inter_data_frames, function(x) {
  x[] <- lapply(x, function(y) if(is.numeric(y)) round(y, 3) else y)
  x
})

# save all in one excel file
write_xlsx(intra_data_frames, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/intra_icc.xlsx")
write_xlsx(inter_data_frames, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/inter_icc.xlsx")

############################################################################################################
# test for ostium
ostium <- inter_og %>% select(patient_id, ccta_ostial_a, ccta_ostial_elliptic, ccta_ostial_pn, ccta_ostial_h, ccta_ostial_w)

baseline_icc_ostium <- baseline %>% 
  select(patient_id, ccta_ostial_a, ccta_ostial_elliptic, ccta_ostial_pn, ccta_ostial_h, ccta_ostial_w, ccta_quali, ccta_photon) %>%
  left_join(ostium, by = "patient_id")

baseline_icc_ostium <- baseline_icc_ostium %>% drop_na()

baseline_icc_ostium$ccta_quali <- as.numeric(as.character(baseline_icc_ostium$ccta_quali))
baseline_icc_ostium <- baseline_icc_ostium %>% mutate(ccta_quali = ifelse(ccta_photon == "yes", ccta_quali + 1, ccta_quali))

baseline_icc_ostium_no_weight <- baseline_icc_ostium
# copy every row n-times with n being ccta_quali
baseline_icc_ostium <- baseline_icc_ostium %>% 
  slice(rep(row_number(), baseline_icc_ostium$ccta_quali * 2)) %>%
  select(-ccta_quali, -ccta_photon)

ostial_area_intra <- ICC(cbind(baseline_icc_ostium$ccta_ostial_a.x, baseline_icc_ostium$ccta_ostial_a.y))$results
icc(cbind(baseline_icc_ostium$ccta_ostial_a.x, baseline_icc_ostium$ccta_ostial_a.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_ostium_no_weight$ccta_ostial_a.x, baseline_icc_ostium_no_weight$ccta_ostial_a.y), model = "twoway", type = "agreement", unit = "single")

area <- baseline_icc_ostium %>% select(patient_id, ccta_ostial_a.x, ccta_ostial_a.y)
ba_ostial_area_intra <- create_bland_altman_plot(area, "ccta_ostial_a.x", "ccta_ostial_a.y")

ostial_elliptic_intra <- ICC(cbind(baseline_icc_ostium$ccta_ostial_elliptic.x, baseline_icc_ostium$ccta_ostial_elliptic.y))$results
icc(cbind(baseline_icc_ostium$ccta_ostial_elliptic.x, baseline_icc_ostium$ccta_ostial_elliptic.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_ostium_no_weight$ccta_ostial_elliptic.x, baseline_icc_ostium_no_weight$ccta_ostial_elliptic.y), model = "twoway", type = "agreement", unit = "single")

elliptic <- baseline_icc_ostium %>% select(patient_id, ccta_ostial_elliptic.x, ccta_ostial_elliptic.y)
ba_ostial_elliptic_intra <- create_bland_altman_plot(elliptic, "ccta_ostial_elliptic.x", "ccta_ostial_elliptic.y")

ostial_ln_intra <- ICC(cbind(baseline_icc_ostium$ccta_ostial_pn.x, baseline_icc_ostium$ccta_ostial_pn.y))$results
icc(cbind(baseline_icc_ostium$ccta_ostial_pn.x, baseline_icc_ostium$ccta_ostial_pn.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_ostium_no_weight$ccta_ostial_pn.x, baseline_icc_ostium_no_weight$ccta_ostial_pn.y), model = "twoway", type = "agreement", unit = "single")

ln <- baseline_icc_ostium %>% select(patient_id, ccta_ostial_pn.x, ccta_ostial_pn.y)
ba_ostial_ln_intra <- create_bland_altman_plot(ln, "ccta_ostial_pn.x", "ccta_ostial_pn.y")

ostial_minor_intra <- ICC(cbind(baseline_icc_ostium$ccta_ostial_w.x, baseline_icc_ostium$ccta_ostial_w.y))$results
icc(cbind(baseline_icc_ostium$ccta_ostial_w.x, baseline_icc_ostium$ccta_ostial_w.y), model = "twoway", type = "agreement", unit = "single")
icc(cbind(baseline_icc_ostium_no_weight$ccta_ostial_w.x, baseline_icc_ostium_no_weight$ccta_ostial_w.y), model = "twoway", type = "agreement", unit = "single")

w <- baseline_icc_ostium %>% select(patient_id, ccta_ostial_w.x, ccta_ostial_w.y)
ba_ostial_minor_inter <- create_bland_altman_plot(w, "ccta_ostial_w.x", "ccta_ostial_w.y", -1.5, 1.5)

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ostial_area_intra.png", ba_ostial_area_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ostial_elliptic_intra.png", ba_ostial_elliptic_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ostial_ln_intra.png", ba_ostial_ln_intra, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ba_ostial_minor_inter.png", ba_ostial_minor_inter, width = 6, height = 5)

ostium_data_frames <- list(
  "Area" = ostial_area_intra,
  "Elliptic" = ostial_elliptic_intra,
  "MLN" = ostial_ln_intra,
  "Minor" = ostial_minor_intra
)

# get index as column
ostium_data_frames <- lapply(ostium_data_frames, function(x) {
  x$Variable <- rownames(x)
  rownames(x) <- NULL
  x
})

# round all values to 3 decimals if numeric
ostium_data_frames <- lapply(ostium_data_frames, function(x) {
  x[] <- lapply(x, function(y) if(is.numeric(y)) round(y, 3) else y)
  x
})

# save all in one excel file
write_xlsx(ostium_data_frames, "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/ostium_icc.xlsx")

ostial_w_plot <- ggroc(roc_ostial_w, legacy.axes = TRUE, color = "darkred") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC ostial =", round(roc_ostial_w$auc, 2)), color = "darkred")

intra <- intra %>% select(patient_id, ccta_ostial_w)
test <- baseline %>% select(patient_id, ffr_0.8) %>%
  left_join(intra, by = "patient_id")

roc_ostial_w_2 <- roc(test$ffr_0.8, test$ccta_ostial_w)
ostial_w_2_plot <- ggroc(roc_ostial_w_2, legacy.axes = TRUE, color = "darkred") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC ostial =", round(roc_ostial_w_2$auc, 2)), color = "darkred")

test_2 <- baseline %>% select(patient_id, ffr_0.8, ccta_ostial_w) %>%
  left_join(intra, by = "patient_id")
test_2 <- test_2 %>% mutate(
  ccta_ostial_w = (ccta_ostial_w.x + ccta_ostial_w.y) / 2)

roc_ostial_w_3 <- roc(test_2$ffr_0.8, test_2$ccta_ostial_w)
ostial_w_3_plot <- ggroc(roc_ostial_w_3, legacy.axes = TRUE, color = "darkred") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC ostial =", round(roc_ostial_w_3$auc, 2)), color = "darkred")

outliers <- ggplot(data = test_2, aes(x = ccta_ostial_w.x, y = ccta_ostial_w.y)) +
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_abline(intercept = 0.37, slope = 1, linetype = "dotted", color = "red") +
  geom_abline(intercept = -0.37, slope = 1, linetype = "dotted", color = "red") +
  theme_classic()

mdl <- lm(ccta_ostial_w.y ~ ccta_ostial_w.x, data = test_2)
summary(mdl)

ostial <- ostium %>% select(patient_id, ccta_ostial_w)
test_4 <- baseline %>% select(patient_id, ffr_0.8) %>%
  left_join(ostial, by = "patient_id")

roc_ostial_w_4 <- roc(test_4$ffr_0.8, test_4$ccta_ostial_w)

ostial_w_4 <- ggroc(roc_ostial_w_4, legacy.axes = TRUE, color = "darkred") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC ostial =", round(roc_ostial_w_4$auc, 2)), color = "darkred")

# for test_5 randomly take test_2 and randomly pick either value from ccta_ostial_w.x or ccta_ostial_w.y to create ccta_ostial_w
# Function to generate ROC curve and extract AUC
generate_roc <- function(data) {
  data <- data %>% mutate(
    ccta_ostial_w = ifelse(runif(nrow(data)) > 0.5,
                           data[["ccta_ostial_w.x"]],
                           data[["ccta_ostial_w.y"]])
  )
  roc_curve <- roc(data$ffr_0.8, data$ccta_ostial_w)
  return(roc_curve)
}

# Run the procedure 50 times
roc_curves <- replicate(50, generate_roc(test_2), simplify = FALSE)

# Extract AUC values
auc_values <- sapply(roc_curves, function(x) x$auc)
mean_roc <- roc(roc_curves[[1]]$response,
                rowMeans(sapply(roc_curves, function(x) x$predictor)))
# Calculate mean and CI for ROC curve
ci_se_roc <- ci.se(mean_roc)
ci_sp_roc <- ci.sp(mean_roc)

# Create a data frame for sensitivity CI
ci_se_df <- data.frame(
  x = x_ci,
  ymin = ymin_ci_se,
  ymax = ymax_ci_se
)

# Create a data frame for specificity CI
ci_sp_df <- data.frame(
  x = x_ci,
  ymin = ymin_ci_sp,
  ymax = ymax_ci_sp
)

# Add identifier columns to the data frames
ci_se_df$type <- "Sensitivity"
ci_sp_df$type <- "Specificity"

# Combine the two data frames
ci_combined_df <- rbind(ci_se_df, ci_sp_df)

# # Plot with combined data
# ostial_w_averaged <- ggroc(mean_roc, legacy.axes = FALSE, color = "darkred") +
#   geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
#   ggtitle("Mean ROC Curve with 95% CI") +
#   annotate("text", x = 0.3, y = 0.3,
#            label = paste("Mean AUC =", round(mean(auc_values), 2)),
#            color = "darkred") +
#   theme_classic() +
#   geom_ribbon(data = ci_combined_df, aes(x = x, ymin = ymin, ymax = ymax, fill = type),
#               alpha = 0.2) +
#   scale_fill_manual(values = c("Sensitivity" = "blue", "Specificity" = "green"))


ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_plot.png", ostial_w_plot, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_2_plot.png", ostial_w_2_plot, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_3_plot.png", ostial_w_3_plot, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/outliers.png", outliers, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_4.png", ostial_w_4, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_averaged.png", ostial_w_averaged, width = 6, height = 5)

roc_highest <- data.frame(
  specificity = roc_ostial_w$specificities,
  sensitivity = roc_ostial_w$sensitivities
)

# order specificity ascending and sensitivity descending tidyverse
roc_highest <- roc_highest %>% arrange(specificity, desc(sensitivity))

roc_lowest <- data.frame(
  specificity = roc_ostial_w_2$specificities,
  sensitivity = roc_ostial_w_2$sensitivities
)

# order specificity ascending and sensitivity descending tidyverse
roc_lowest <- roc_lowest %>% arrange(specificity, desc(sensitivity))

# Plot with combined data
ostial_w_averaged_combined <- ggroc(mean_roc, legacy.axes = FALSE, color = "darkred") +
  geom_abline(intercept = 1, slope = 1, linetype = "dashed") +
  ggtitle("Mean ROC Curve with 95% CI") +
  annotate("text", x = 0.3, y = 0.3,
           label = paste("Mean AUC =", round(mean(auc_values), 2)),
           color = "darkred") +
  theme_classic() +         
  geom_ribbon(data = ci_combined_df, aes(x = x, ymin = ymin, ymax = ymax, fill = type), 
              alpha = 0.2) +
  scale_fill_manual(values = c("Sensitivity" = "blue", "Specificity" = "green")) +
  geom_line(data = roc_highest, aes(x = specificity, y = sensitivity), color = "blue", linetype = "dotted") +
  geom_line(data = roc_lowest, aes(x = specificity, y = sensitivity), color = "blue", linetype = "dotted")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ostial_w_averaged_combined.png", ostial_w_averaged_combined, width = 6, height = 5)

mla_vs_minor <- ggplot(baseline, aes(ccta_mla_a, ccta_mla_w)) +
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) +
  ggtitle("MLA versus MLA minor axis") +
  xlab("MLA [mm²]") +
  ylab("MLA minor axis [mm]") +
  theme_classic()

ola_vs_minor <- ggplot(baseline, aes(ccta_ostial_a, ccta_ostial_w)) +
  geom_point() +
  geom_smooth(method ="lm", se = FALSE) +
  ggtitle("Ostial area versus ostial minor axis") +
  xlab("Ostial area [mm²]") +
  ylab("OLA minor axis [mm]") +
  theme_classic()

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/mla_vs_minor.png", mla_vs_minor, width = 6, height = 5)
ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/ola_vs_minor.png", ola_vs_minor, width = 6, height = 5)

# confusion matrix ostial area
# cut-off ostial area: 7.66 mm²
baseline$ccta_ostial_a_pos <- ifelse(baseline$ccta_ostial_a > 7.66, 0, 1)

table <- baseline %>%
  select(ccta_ostial_a_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_ccta <- conf_mat(table)
result_ccta <- summary(confusion_ccta, event_level = "second")

data <- data.frame(
  "target" = baseline$ffr_0.8,
  "prediction" = baseline$ccta_ostial_a_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_ccta_ostial_a <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_ccta_ostial_a.png", conf_mat_ccta_ostial_a, width = 3, height = 3)

# confusion matrix ostial area
# cut-off ostial area: 1.75mm
baseline$ccta_ostial_w_pos <- ifelse(baseline$ccta_ostial_w > 1.75, 0, 1)

table <- baseline %>%
  select(ccta_ostial_w_pos, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion_ccta <- conf_mat(table)
result_ccta <- summary(confusion_ccta, event_level = "second")

data <- data.frame(
  "target" = baseline$ffr_0.8,
  "prediction" = baseline$ccta_ostial_w_pos,
  stringsAsFactors = FALSE
)

# Evaluate predictions and create confusion matrix
eval <- evaluate(
  data = data,
  target_col = "target",
  prediction_cols = "prediction",
  type = "binomial"
)

conf_mat_ccta_ostial_w <- plot_confusion_matrix(
  eval,
  add_sums = FALSE,
  add_row_percentages = FALSE,
  add_col_percentages = FALSE,
  rm_zero_percentages = FALSE,
  rm_zero_text = FALSE,
  counts_on_top = TRUE,
  palette = "Blues",
  intensity_by = "counts",
  digits = 1,
  darkness = 0.8
)

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/confusion_ccta_ostial_w.png", conf_mat_ccta_ostial_w, width = 3, height = 3)

ccta_variables <- c("ccta_mla_a", "ccta_mla_elliptic", "ccta_mla_ln", "ccta_mla_w", "ccta_ostial_a", "ccta_ostial_elliptic", "ccta_ostial_pn", "ccta_ostial_w")

# for all create a df with mean, sd, median and iqr per variable (tidyverse)
summary_ccta <- baseline %>%
  select(all_of(ccta_variables)) %>%
  gather(key = "variable", value = "value") %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    iqr = IQR(value, na.rm = TRUE),
    p_value = shapiro.test(value)$p.value
  )

# roc curve for roc_ostial
plot_roc_ostial <- ggroc(roc_ostial, legacy.axes = TRUE, color = "darkred") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("Minimal lumen") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC ostial =", round(roc_ostial$auc, 2)), color = "darkred")

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/roc_ostial.png", plot_roc_ostial, width = 6, height = 5)

############################################################################################################
plot_ffr_ostial <- ggplot(baseline, aes(x = inv_ffrdobu, y = ccta_ostial_a)) +
  geom_smooth(method = "gam", formula = y ~ poly(x, 3), se = FALSE, linetype = "dashed", aes(color = "#d8d6d6"), alpha=0.5) +
  geom_point(aes(color = as.factor(funct_pos))) +
  scale_color_manual(values = c("0" = "darkblue", "1" = "orange")) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  annotate("text", x = 0.6, y = 7.7, label = "OLA >7.7 mm²", color = "darkgreen", vjust = -1) +
  geom_hline(yintercept = 7.7, linetype = "dashed", color = "darkgreen") +
  labs(x = "Fractional Flow Reserve", y = "Ostial lumen area [mm²]") +
  theme_classic()

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/plot_ffr_ostial.png", plot_ffr_ostial, width = 5, height = 4)

plot_ffr_ostial_w <- ggplot(baseline, aes(x = inv_ffrdobu, y = ccta_ostial_w)) +
  geom_smooth(method = "gam", formula = y ~ poly(x, 3), se = FALSE, linetype = "dashed", aes(color = "#d8d6d6"), alpha=0.5) +
  geom_jitter(aes(color = as.factor(funct_pos))) +
  scale_color_manual(values = c("0" = "darkblue", "1" = "orange"), labels = c("Negative", "Positive")) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
  annotate("text", x = 0.6, y = 1.75, label = "OLA minor axis ≥1.8 mm", color = "darkgreen", vjust = -1) +
  geom_hline(yintercept = 1.75, linetype = "dashed", color = "darkgreen") +
  labs(x = "Fractional Flow Reserve", y = "CCTA ostial minor axis [mm]", color = "Nuclear cardiac imaging") +
  theme_classic()

ggsave("C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/figures/plot_ffr_ostial_w.png", plot_ffr_ostial_w, width = 7, height = 5)

no_ffr <- baseline %>% filter(ffr_0.8 == 0)
ffr <- baseline %>% filter(ffr_0.8 == 1)

sum(ffr$funct_spect_mismatchcaa == 1)

baseline %>% filter(inv_cad_perc_rcaprox >=2 | inv_cad_perc_rcamid >=2 | inv_cad_perc_rcadist >=2) %>% nrow()
baseline %>% filter(inv_cad_perc_ladprox >=2 | inv_cad_perc_ladmid >=2 | inv_cad_perc_rcaprox >=2 | inv_cad_perc_rcamid >=2 | inv_cad_perc_rcadist >=2 | 
        inv_cad_perc_cxprox >=2 | inv_cad_perc_cxmid >=2 | inv_cad_perc_cxdist >=2 |  inv_cad_perc_diag1 >= 2 | inv_cad_perc_diag2 >=2) %>% nrow()

no_ffr %>% filter(inv_cad_perc_rcaprox >=2 | inv_cad_perc_rcamid >=2 | inv_cad_perc_rcadist >=2) %>% nrow()
no_ffr %>% filter(inv_cad_perc_ladprox >=2 | inv_cad_perc_ladmid >=2 | inv_cad_perc_rcaprox >=2 | inv_cad_perc_rcamid >=2 | inv_cad_perc_rcadist >=2 | 
        inv_cad_perc_cxprox >=2 | inv_cad_perc_cxmid >=2 | inv_cad_perc_cxdist >=2 |  inv_cad_perc_diag1 >= 2 | inv_cad_perc_diag2 >=2) %>% nrow()

no_ffr %>% select(ccta_ostial_a, ccta_ostial_w, inv_ffrdobu) %>% arrange(ccta_ostial_a, ascending=F)


chisq.test(a0.8$synp_treatment_type, b0.8$synp_treatment_type)