library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)
library(readxl)
library(pROC)
library(yardstick)
library(writexl)
library(caret)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# invasive <- read_excel("C:/WorkingData/Documents/3_Research/IVUS_data.xlsx")

# Load the data
baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))

baseline <- baseline %>%
  filter(record_id != 109) %>%
  filter(!is.na(inv_ffrdobu))

baseline <- baseline %>% mutate(
  ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0),
  funct_pos = case_when(
    funct_spect_mismatchcaa == "yes" | funct_pet_mismatchcaa == "yes" ~ 1,
    funct_spect_mismatch == "no" | funct_pet_mismatch == "no" | funct_spect_mismatchcaa == "no" | funct_pet_mismatchcaa == "no" ~ 0,
    TRUE ~ NA_real_
  ),
)

baseline %>% select(record_id, patient_id, synp_inv_yn, synp_funct_yn, synp_funct_type) %>% View()

perform_logistic_regression <- function(response_var, data, ivus_vars) {
  logistic_results <- list()
  
  for (i in 1:length(ivus_vars)) {
    ivus_var <- ivus_vars[i]
    
    p_value <- NA
    odds_ratio <- NA
    ci_lower <- NA
    ci_upper <- NA
    
    tryCatch({
      ivus_glm <- glm(as.formula(paste(response_var, "~", ivus_var)), family = "binomial", data = data)
      
      p_value <- round(summary(ivus_glm)$coefficients[2, 4], 3)
      odds_ratio <- round(exp(coef(ivus_glm)[2]), 2)
      ci_lower <- round(exp(confint(ivus_glm)[2, 1]), 2)
      ci_upper <- round(exp(confint(ivus_glm)[2, 2]), 2)
    }, error = function(e) {
      # If error occurs, leave p_value and odds_ratio as NA
    })
    
    logistic_results[[i]] <- data.frame(
      Variable = ivus_var,
      P_Value = p_value,
      Odds_Ratio = odds_ratio,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    )
  }
  
  logistic_results_df <- do.call(rbind, logistic_results)
  
  return(logistic_results_df)
}

perform_roc_analysis <- function(response_var, data, ivus_vars) {
  roc_results <- list()
  
  for (i in 1:length(ivus_vars)) {
    ivus_var <- ivus_vars[i]
    
    auc <- NA
    sensitivity <- NA
    specificity <- NA
    threshold <- NA
    
    tryCatch({
      roc_model <- glm(as.formula(paste(response_var, "~", ivus_var)), data = data, family = binomial)
      
      roc_pred <- predict(roc_model, data, type = "response")
      
      roc_curve <- roc(data[[response_var]], roc_pred)
      
      auc <- round(auc(roc_curve), 2)
      
      roc_threshold <- roc(data[[response_var]], data[[ivus_var]])
      opt_threshold <- coords(roc_threshold, "best", ret = "threshold")
      threshold <- round(opt_threshold, 2)
      
      sens_spec <- coords(roc_curve, "best", ret = c("sensitivity", "specificity"))
      sensitivity <- round(sens_spec["sensitivity"], 2)
      specificity <- round(sens_spec["specificity"], 2)
    }, error = function(e) {
      # If error occurs, leave variables as NA
    })    

    roc_results[[i]] <- data.frame(
      Variable = ivus_var,
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

# Define the variables for comparison
rest_vars <- c("inv_ffrado", 
                "inv_rest_hr", 
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean")

dobu_vars <- c("inv_ffrdobu", 
                "inv_dobu_hr",
               "inv_dobu_aosys", 
               "inv_dobu_aodia", 
               "inv_dobu_aomean")

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

# Optionally, save the results to a CSV file
write.csv(results_df, "statistical_analysis/data/invasive_change.csv", row.names = FALSE)

# baseline IVUS variables to predict FFR change
ivus_vars <- c("inv_ivusrest_mla", 
               "inv_ivusrest_mla_ellip",
               "inv_ivusrest_mla_ln", 
               "inv_ffrado")

logistic_results_ffr_0.8 <- perform_logistic_regression("ffr_0.8", baseline, ivus_vars)
logistic_results_risk_mla_any <- perform_logistic_regression("risk_mla_any", baseline, ivus_vars)

results_list <- list(
  "ffr_0.8" = logistic_results_ffr_0.8,
  "risk_mla_any" = logistic_results_risk_mla_any
)

# Write the list of data frames to an Excel file
write_xlsx(results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ivus_logistic_regression_results.xlsx")

######### ROC analysis for IVUS variables ###############################################################
# Perform ROC analysis for all the specified response variables
roc_results_ffr_0.8 <- perform_roc_analysis("ffr_0.8", baseline, ivus_vars)
roc_results_risk_mla_any <- perform_roc_analysis("risk_mla_any", baseline, ivus_vars)

# Create a list of data frames to be written to the Excel file
roc_results_list <- list(
  "ffr_0.8" = roc_results_ffr_0.8,
  "risk_mla_any" = roc_results_risk_mla_any
)

# Write the list of data frames to an Excel file
write_xlsx(roc_results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ivus_roc_analysis_results.xlsx")

############################################################################################################
# best cutoffs
# Fit models and get prediction data
mdl_mla_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla")[[1]]
mdl_mla_ellip_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ellip")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[1]]
mdl_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[1]]

prediction_data_mla <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla")[[2]]
prediction_data_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ellip")[[2]]
prediction_data_mla_ln <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[2]]
prediction_data_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[2]]

data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla) %>% drop_na()
roc_mla <- roc(data_ivus$ffr_0.8, mdl_mla_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_ellip) %>% drop_na()
roc_mla_ellip <- roc(data_ivus$ffr_0.8, mdl_mla_ellip_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ffrado) %>% drop_na()
roc_ffrado <- roc(data_ivus$ffr_0.8, mdl_ffrado$fitted.values)

# Define test sets
test_set <- list(
  coords(roc_mla), coords(roc_mla_ellip),
  coords(roc_mla_ln), coords(roc_ffrado)
)

prediction_data <- list(
  prediction_data_mla, prediction_data_ellip,
  prediction_data_mla_ln, prediction_data_ffrado
)

variables <- c(
  "IVUS rest MLA", "IVUS rest elliptic ratio",
  "IVUS rest MLA luminal narrowing", "FFR adenosine"
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
    
    if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla", "inv_ffrado")) {
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
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/"
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

ggsave(paste0(output_dir,"/roc_mla.png"), ggroc_mla, width = 6, height = 5)

shapiro.test(baseline$inv_rfr)
median(baseline$inv_rfr, na.rm = T)
quantile(baseline$inv_rfr, c(0.25, 0.75), na.rm = T)

# ANOVA paired for inv_rfr, inv_ffrado, inv_ffrdobu
anova_test <- aov(inv_rfr ~ inv_ffrado + inv_ffrdobu, data = baseline)
summary(anova_test)
# post hoc
wilcox.test(baseline$inv_rfr, baseline$inv_ffrado, paired = T)
wilcox.test(baseline$inv_rfr, baseline$inv_ffrdobu, paired = T)
wilcox.test(baseline$inv_ffrado, baseline$inv_ffrdobu, paired = T)

# same with inv_ivusrest_mla_ln, inv_ivusado_mla_ln, inv_ivusdobu_mla_ln_any
anova_test <- aov(inv_ivusrest_mla_ln ~ inv_ivusado_mla_ln + inv_ivusdobu_mla_ln_any, data = baseline)
summary(anova_test)
# post hoc
wilcox.test(baseline$inv_ivusrest_mla_ln, baseline$inv_ivusado_mla_ln, paired = T)
wilcox.test(baseline$inv_ivusrest_mla_ln, baseline$inv_ivusdobu_mla_ln_any, paired = T)
wilcox.test(baseline$inv_ivusado_mla_ln, baseline$inv_ivusdobu_mla_ln_any, paired = T)

median(baseline$inv_ivusrest_mla_ln, na.rm = T)
quantile(baseline$inv_ivusrest_mla_ln, c(0.25, 0.75), na.rm = T)
median(baseline$inv_ivusado_mla_ln, na.rm = T)
quantile(baseline$inv_ivusado_mla_ln, c(0.25, 0.75), na.rm = T)
median(baseline$inv_ivusdobu_mla_ln_any, na.rm = T)
quantile(baseline$inv_ivusdobu_mla_ln_any, c(0.25, 0.75), na.rm = T)

############################################################################################################
# ROC Cross-Validation for FFR 0.8
cross <- baseline %>% drop_na(ffr_0.8, all_of(ivus_vars))
cross$ffr_0.8 <- factor(cross$ffr_0.8, levels = c(0, 1), labels = c("Class0", "Class1"))

# Create a list of 50 seeds
set.seed(69)
seeds <- sample(1:1000, 100)

results_list <- list()

# Function to train and evaluate models
train_and_evaluate <- function(seed) {
  set.seed(seed)
  
  # Split the data into training (80%) and testing (20%) sets
  train_indices <- createDataPartition(cross$ffr_0.8, p = 0.78, list = FALSE)
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
  for (ivus_var in ivus_vars) {
    formula <- as.formula(paste("ffr_0.8 ~", ivus_var))
    
    # Model fitting with error handling
    model <- tryCatch({
      glm(formula, data = train_data, family = "binomial")
    }, error = function(e) {
      message(sprintf("Error fitting model for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(model)) next

    # Predictions with error handling
    predictions <- tryCatch({
      predict(model, test_data, type = "response")
    }, error = function(e) {
      message(sprintf("Error predicting for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(predictions)) next

    # ROC curve calculation with error handling
    roc_curve <- tryCatch({
      roc(test_data$ffr_0.8, predictions)
    }, error = function(e) {
      message(sprintf("Error calculating ROC for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(roc_curve)) next

    best_coords <- tryCatch({
      coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
    }, error = function(e) {
      message(sprintf("Error in coords for variable %s: %s", ivus_var, e))
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
      value <- test_data[[ivus_var]][closest_index]

      table <- if (ivus_var %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla", "inv_ffrado")) {
        test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) > value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      } else {
        test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) < value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      }

      # Confusion matrix with error handling
      if (nrow(table) != 2 || ncol(table) != 2) {
        message(sprintf("Confusion matrix not square for variable %s at seed %d", ivus_var, seed))
        next
      }

      confusion <- conf_mat(table)
      result <- summary(confusion, event_level = "second")

      sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
      spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
      ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
      npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
      acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()

      # Naming convention ivus_var_1, ivus_var_2, ivus_var_3
      seed_results[[paste(ivus_var, i, sep = "_")]] <- list(
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
  
  for (ivus_var in ivus_vars) {
    for (i in 1:3) {
      var_name <- paste(ivus_var, i, sep = "_")
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
write_xlsx(roc_summary_df_ffr0.8_rounded, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/roc_summary_ffr0.8_rounded_summary.xlsx")


############################################################################################################
# cut-off 0.8 for ffr_ado
# Fit models and get prediction data
baseline <- baseline %>% mutate(
  ffr_ado0.8 = ifelse(inv_ffrado <= 0.8, 1, 0)
)

table <- baseline %>%
  select(ffr_ado0.8, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

############################################################################################################

# roc curves
roc_invffrado <- roc(baseline$ffr_0.8, baseline$inv_ffrado)
ggroc(roc_invffrado, legacy.axes = TRUE, color = "darkblue") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("FFR adenosine") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC FFR adenosine =", round(roc_invffrado$auc, 2)))

# forest plot
# for hemodynamic and anatomic relevance as predictor
forest_df_bsa <- data.frame(
    var = c(
            "Minimal lumen area",
            "Minimal lumen elliptic ratio",
            "Minimal lumen narrowing [%]", 
            "FFR adenosine"
            ),
    odds_ratio = c(
            exp(mdl_mla_ffr$coefficients[2]),
            exp(mdl_mla_ellip_ffr$coefficients[2]),
            exp(mdl_mla_ln_ffr$coefficients[2]),
            exp(mdl_ffrado$coefficients[2])),
    lower_ci = c(
            exp(confint(mdl_mla_ffr))[2],
            exp(confint(mdl_mla_ellip_ffr))[2],
            exp(confint(mdl_mla_ln_ffr))[2], 
            exp(confint(mdl_ffrado))[2]), 
    upper_ci = c(
            exp(confint(mdl_mla_ffr))[4],
            exp(confint(mdl_mla_ellip_ffr))[4],
            exp(confint(mdl_mla_ln_ffr))[4],  
            exp(confint(mdl_ffrado))[4])
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

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_ffr_0.8.png", forest_ffr_0.8, width = 5, height = 2)

cross <- baseline %>% drop_na(ffr_0.8, all_of(ivus_vars))
cross <- cross %>% select(-inv_ffrado)
cross$ffr_0.8 <- factor(cross$ffr_0.8, levels = c(0, 1), labels = c("Class0", "Class1"))

set.seed(69)

n_seeds <- seq(10, 120, 10)

# Initialize dataframe to store convergence results
convergence_results <- data.frame(
  Seed = numeric(),
  Variable = character(),
  Mean_AUC = numeric()
)

# Function to train and evaluate models
train_and_evaluate <- function(seed) {
  set.seed(seed)
  
  # Split data into training (80%) and testing (20%) sets
  train_indices <- createDataPartition(cross$ffr_0.8, p = 0.78, list = FALSE)
  train_data <- cross[train_indices, ]
  test_data <- cross[-train_indices, ]
  
  # Check class balance in training and test sets
  if (length(unique(train_data$ffr_0.8)) < 2 || length(unique(test_data$ffr_0.8)) < 2) {
    message(sprintf("Seed %d: Training or Test set does not contain both classes.", seed))
    return(NULL)
  }

  # Store results for this seed
  seed_results <- list()
  
  # Train models and evaluate ROC for each IVUS variable
  for (ivus_var in ivus_vars) {
    formula <- as.formula(paste("ffr_0.8 ~", ivus_var))
    
    # Fit model
    model <- tryCatch({
      glm(formula, data = train_data, family = "binomial")
    }, error = function(e) {
      message(sprintf("Error fitting model for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(model)) next

    # Predict and evaluate ROC
    predictions <- tryCatch({
      predict(model, test_data, type = "response")
    }, error = function(e) {
      message(sprintf("Error predicting for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(predictions)) next

    # Calculate ROC and AUC
    roc_curve <- tryCatch({
      roc(test_data$ffr_0.8, predictions)
    }, error = function(e) {
      message(sprintf("Error calculating ROC for variable %s: %s", ivus_var, e))
      return(NULL)
    })

    if (is.null(roc_curve)) next

    seed_results[[ivus_var]] <- as.numeric(auc(roc_curve))
  }
  
  return(seed_results)
}

# Main loop for different seed counts
for (seed_count in n_seeds) {
  seeds <- sample(1:1000, seed_count)

  results_list <- vector("list", length(seeds))
  names(results_list) <- seeds
  
  # Run the training and evaluation for each seed
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

  # Aggregate results
  aggregate_results <- function(results_list) {
    aggregate_list <- list()
    
    for (ivus_var in ivus_vars) {
      aucs <- sapply(results_list, function(x) if (!is.null(x[[ivus_var]])) x[[ivus_var]] else NA)
      
      aggregate_list[[ivus_var]] <- data.frame(
        Seed = seed_count,
        Variable = ivus_var,
        Mean_AUC = mean(aucs, na.rm = TRUE)
      )
    }
    
    return(aggregate_list)
  }

  final_results <- aggregate_results(results_list)

  # Convert list to dataframe
  seed_results_df <- do.call(rbind, final_results)
  convergence_results <- rbind(convergence_results, seed_results_df)
}

# Save results to Excel
# write_xlsx(convergence_results, "roc_auc_convergence_results.xlsx")

ggplot(convergence_results, aes(x = Seed, y = Mean_AUC, color = Variable)) +
  # geom_point() +
  geom_smooth(method = "loess", se = FALSE) +
  scale_x_continuous(breaks = seq(0, 200, 10)) +
  scale_y_continuous(breaks = seq(0.6, 1.0, 0.01)) +
  labs(title = "AUC Convergence by Seed Count",
       x = "Number of Seeds",
       y = "Mean AUC",
       color = "IVUS Variable") +
  theme_classic()

# Save the plot
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/auc_convergence.png", width = 6, height = 4)

above0.8 <- baseline %>% filter(inv_ffrdobu > 0.8)
below0.8 <- baseline %>% filter(inv_ffrdobu <= 0.8)

median(above0.8$inv_ffrado, na.rm = T)
quantile(above0.8$inv_ffrado, c(0.25, 0.75), na.rm = T)

median(below0.8$inv_ffrado, na.rm = T)
quantile(below0.8$inv_ffrado, c(0.25, 0.75), na.rm = T)

median(baseline$inv_ivusrest_mla_ellip, na.rm = T)
quantile(baseline$inv_ivusrest_mla_ellip, c(0.25, 0.75), na.rm = T)

median(above0.8$inv_ivusrest_mla_ellip, na.rm = T)
quantile(above0.8$inv_ivusrest_mla_ellip, c(0.25, 0.75), na.rm = T)

median(below0.8$inv_ivusrest_mla_ellip, na.rm = T)
quantile(below0.8$inv_ivusrest_mla_ellip, c(0.25, 0.75), na.rm = T)

wilcox.test(above0.8$inv_ffrado, below0.8$inv_ffrado, paired = F)

no_treatment <- baseline %>% filter(synp_surgery == "no" & synp_stent == "no" & inv_ffrdobu <= 0.8)
treatment <- baseline %>% filter(synp_surgery == "yes" | synp_stent == "yes")

median(no_treatment$inv_ffrdobu, na.rm = T)
quantile(no_treatment$inv_ffrdobu, c(0.25, 0.75), na.rm = T)

median(treatment$inv_ffrdobu, na.rm = T)
quantile(treatment$inv_ffrdobu, c(0.25, 0.75), na.rm = T)

wilcox.test(no_treatment$inv_ffrdobu, treatment$inv_ffrdobu, paired = F)

test <- baseline %>% filter(!is.na(inv_ffrado)) %>% select(patient_id, inv_ffrado, inv_ffrdobu, inv_ivusrest_mla, inv_date)

test %>% filter(inv_ffrado <= 0.8)

test %>% filter(inv_ffrado > 0.8) %>% print(n = 50)

test %>% filter(inv_ffrado >0.8 & inv_ffrdobu <= 0.8)