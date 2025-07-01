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
  filter(caa_course___0 == "yes") %>%
  filter(!is.na(inv_ivusdobu_mla) & !is.na(inv_ffrdobu))

baseline <- baseline %>%
  mutate(percent_stenosis = cut(inv_ivusdobu_mla_ln_any, 
                        breaks = c(-Inf, 50, 70, 90, Inf), 
                        labels = c("<50", "50-70", "70-90", ">90")))
# Map the bins to specific sizes
size_values <- c("<50" = 2, "50-70" = 5, "70-90" = 10, ">90" = 20)

# ostium response variable
baseline <- baseline %>% mutate(
  risk_mla_any = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_mla_ln_any) ~ NA_real_,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_any < 70 ~ 0,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_any >= 70 ~ 1,
    inv_ffrdobu <= 0.8 ~ 1
  ),
)

baseline <- baseline %>% mutate(
  ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0)
)
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

# FFR data
n <- nrow(baseline)
inv_var1 <- pull(baseline, inv_ffrado)
inv_var2 <- pull(baseline, inv_ffrdobu)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ffrado <- ifelse(inv_ffrado <= 0.8, 1, 0)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
inv_ffrdobu <- ifelse(inv_ffrdobu <= 0.8, 1, 0)
method_1 <- rep("FFRadenosine", n)
method_2 <- rep("FFRdobutamine", n)
ivus_df <- data.frame(method = c(method_1, method_2), 
                      value = c(inv_var1, inv_var2),
                      ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
ivus_df$method <- factor(ivus_df$method, levels = c("FFRadenosine", "FFRdobutamine"), ordered = TRUE)
colnames(ivus_df) <- c("method", "value", "ffr_value")

plot_ffr_line <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  # stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(0, 1, 0.05))  + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + # Plot 1
  ylab("FFR")

# iFR plot
# new dataset ifr data with only non missing values for iFR_mean_rest and iFR_mean_dobu
ifr_data <- baseline %>% filter(!is.na(iFR_mean_rest) & !is.na(iFR_mean_dobu))
n <- nrow(ifr_data)
inv_var1 <- pull(ifr_data, iFR_mean_rest)
inv_var2 <- pull(ifr_data, iFR_mean_dobu)
iFR_mean_rest <- pull(ifr_data, iFR_mean_rest)
iFR_mean_rest <- ifelse(iFR_mean_rest <= 0.8, 1, 0)
iFR_mean_dobu <- pull(ifr_data, iFR_mean_dobu)
iFR_mean_dobu <- ifelse(iFR_mean_dobu <= 0.8, 1, 0)
method_1 <- rep("iFR rest", n)
method_2 <- rep("iFR dobu", n)
ivus_df <- data.frame(method = c(method_1, method_2), 
                      value = c(inv_var1, inv_var2),
                      ffr_value = c(iFR_mean_rest, iFR_mean_dobu))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("iFR > 0.8", "iFR <= 0.8"))
ivus_df$method <- factor(ivus_df$method, levels = c("iFR rest", "iFR dobu"), ordered = TRUE)
colnames(ivus_df) <- c("method", "value", "ffr_value")

plot_ifr_line <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(0, 1, 0.05))  + 
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") + # Plot 1
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "orange") + # Plot 1
  geom_hline(yintercept = 0.89, linetype = "dashed", color = "blue") + # Plot 1
  ylab("iFR")

# FFRdobu versus iFRdobu
# new dataset ifr data with only non missing values for iFR_mean_rest and iFR_mean_dobu
comp_data <- baseline %>% filter(!is.na(iFR_mean_dobu) & !is.na(inv_ffrdobu))
n <- nrow(comp_data)
inv_var1 <- pull(comp_data, inv_ffrdobu)
inv_var2 <- pull(comp_data, iFR_mean_dobu)
inv_ffrdobu <- pull(comp_data, inv_ffrdobu)
inv_ffrdobu <- ifelse(inv_ffrdobu <= 0.8, 1, 0)
iFR_mean_dobu <- pull(comp_data, iFR_mean_dobu)
iFR_mean_dobu <- ifelse(iFR_mean_dobu <= 0.8, 1, 0)
method_1 <- rep("FFRdobu", n)
method_2 <- rep("iFRdobu", n)
ivus_df <- data.frame(method = c(method_1, method_2), 
                      value = c(inv_var1, inv_var2),
                      ffr_value = c(inv_ffrdobu, iFR_mean_dobu))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("pdpa > 0.8", "pdpa <= 0.8"))
ivus_df$method <- factor(ivus_df$method, levels = c("FFRdobu", "iFRdobu"), ordered = TRUE)
colnames(ivus_df) <- c("method", "value", "ffr_value")

plot_ffr_ifr_line <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(0, 1, 0.05))  + 
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") + # Plot 1
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "orange") + # Plot 1
  ylab("pdpa")

# FFRado versus iFRado
# new dataset ifr data with only non missing values for iFR_mean_rest and iFR_mean_dobu
comp_data <- baseline %>% filter(!is.na(iFR_mean_ado) & !is.na(inv_ffrado))
n <- nrow(baseline)
inv_var1 <- pull(baseline, inv_ffrado)
inv_var2 <- pull(baseline, iFR_mean_ado)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ffrado <- ifelse(inv_ffrado <= 0.8, 1, 0)
iFR_mean_ado <- pull(baseline, iFR_mean_ado)
iFR_mean_ado <- ifelse(iFR_mean_ado <= 0.8, 1, 0)
method_1 <- rep("FFRado", n)
method_2 <- rep("iFRado", n)
ivus_df <- data.frame(method = c(method_1, method_2), 
                      value = c(inv_var1, inv_var2),
                      ffr_value = c(inv_ffrado, iFR_mean_ado))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("pdpa > 0.8", "pdpa <= 0.8"))
ivus_df$method <- factor(ivus_df$method, levels = c("FFRado", "iFRado"), ordered = TRUE)
colnames(ivus_df) <- c("method", "value", "ffr_value")

plot_ffr_ifr_line_ado <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(0, 1, 0.05))  + 
  geom_hline(yintercept = 0.75, linetype = "dashed", color = "red") + # Plot 1
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "orange") + # Plot 1
  ylab("pdpa")

pressure_rest_ado_dobu <- create_boxplot(baseline, c("inv_rfr", "inv_ffrado", "inv_ffrdobu"), seq(0, 1, 0.1))
ln_rest_ado_dobu <- create_boxplot(baseline, c("inv_ivusrest_mla_ln", "inv_ivusado_mla_ln", "inv_ivusdobu_mla_ln_any"), seq(0, 100, 10))

ggsave(paste0(output_dir,"/plot_ffr_line.png"), plot_ffr_line, width = 4, height = 5)
ggsave(paste0(output_dir,"/pressure_rest_ado_dobu.png"), pressure_rest_ado_dobu, width = 6, height = 5)
ggsave(paste0(output_dir,"/ln_rest_ado_dobu.png"), ln_rest_ado_dobu, width = 6, height = 5)
ggsave(paste0(output_dir,"/plot_ifr_line.png"), plot_ifr_line, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_ffr_ifr_line.png"), plot_ffr_ifr_line, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_ffr_ifr_line_ado.png"), plot_ffr_ifr_line_ado, width = 4, height = 5)

# linear relationships
mla_ffr <- ggplot(baseline, aes(x = inv_ivusrest_mla, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = percent_stenosis, alpha = 0.7, color = inv_ivusrest_mla_ellip)) +  # Map color to the points
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
  area_categories = ifelse(inv_ivusrest_mla < 2, "<2", 
                           ifelse(inv_ivusrest_mla < 4, "2-4", 
                                  ifelse(inv_ivusrest_mla < 6, "4-6", 
                                         ifelse(inv_ivusrest_mla < 8, "6-8", ">8"))))
  )

area_sizes <- c("<2" = 0.5, "2-4" = 1, "4-6" = 3, "6-8" = 6, ">8" = 10)

mln_ffr <- ggplot(baseline, aes(x = inv_ivusrest_mla_ln, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = area_categories, alpha = 0.7, color = inv_ivusrest_mla_ellip)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_gradient(low = "#070775", high = "#ff9100") +
  scale_size_manual(values = area_sizes) +
  labs(size = "MLA categories [mm²]") +
  xlab("Maximal lumen narrowing [%]") +
  ylab("FFR")

mla_ellip_ffr <- ggplot(baseline, aes(x = inv_ivusrest_mla_ellip, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = area_categories, alpha = 0.7, color = inv_ivusrest_mla_ln)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_x_continuous(breaks = seq(0, 7, 0.5)) +
  scale_color_gradient(low = "#070775", high = "#ff9100") +
  scale_size_manual(values = area_sizes) +
  labs(size = "MLA categories [mm²]") +
  xlab("Minimal lumen elliptic ratio") +
  ylab("FFR")

ggsave(paste0(output_dir,"/mla_ffr.png"), mla_ffr, width = 6, height = 5)
ggsave(paste0(output_dir,"/mln_ffr.png"), mln_ffr, width = 6, height = 5)
ggsave(paste0(output_dir,"/mla_ellip_ffr.png"), mla_ellip_ffr, width = 6, height = 5)

# Define the variables for comparison
rest_vars <- c("inv_ffrado",
                "pdpa_mean_ado",
                "iFR_mean_rest",
                "mid_systolic_ratio_mean_rest",
                "inv_rest_hr", 
                "inv_rest_aosys",
                "inv_rest_aodia",
               "inv_rest_aomean"
               )

dobu_vars <- c("inv_ffrdobu", 
                "pdpa_mean_dobu",
                "iFR_mean_dobu",
                "mid_systolic_ratio_mean_dobu",
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

baseline <- baseline %>% mutate(
  ffr_change = inv_ffrdobu - inv_ffrado,
  inv_hr_change = inv_dobu_hr - inv_rest_hr,
  inv_aosys_change = inv_dobu_aosys - inv_rest_aosys,
  inv_aodia_change = inv_dobu_aodia - inv_rest_aodia,
  inv_aomean_change = inv_dobu_aomean - inv_rest_aomean,
  ifr_change = iFR_mean_dobu - iFR_mean_rest,
  mid_systolic_ratio_change = mid_systolic_ratio_mean_dobu - mid_systolic_ratio_mean_rest,
)

relevant <- baseline %>% filter(inv_ffrdobu <= 0.8)
nonrelevant <- baseline %>% filter(inv_ffrdobu > 0.8)

inv_changes <- c( "ffr_change",
"inv_hr_change",
"inv_aosys_change",
"inv_aodia_change",
"inv_aomean_change",
"ifr_change",
"mid_systolic_ratio_change")

p_values <- data.frame(
  Variable = c("FFR", "Heart rate", "Aortic systolic pressure", "Aortic diastolic pressure", "Aortic mean pressure", "iFR", "Mid-systolic ratio"),
  Test = NA,
  P_Value = NA
)

for (i in 1:length(inv_changes)) {
  if (shapiro.test(relevant[[inv_changes[i]]])$p.value > 0.05 & shapiro.test(nonrelevant[[inv_changes[i]]])$p.value > 0.05) {
    p_value <- t.test(relevant[[inv_changes[i]]], nonrelevant[[inv_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "t-test"
    p_values$P_Value[i] <- p_value
  } else {
    p_value <- wilcox.test(relevant[[inv_changes[i]]], nonrelevant[[inv_changes[i]]], paired = FALSE)$p.value
    p_values$Test[i] <- "Wilcoxon"
    p_values$P_Value[i] <- p_value
  }
}

# adjust p-values
p_values$P_Value <- p.adjust(p_values$P_Value, method = "fdr")

results_list <- list(
  "Invasive Changes" = results_df,
  "P-Values" = p_values
)

# Optionally, save the results to a CSV file
write_xlsx(results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/invasive_change.xlsx")

# baseline IVUS variables to predict FFR change
ivus_vars <- c("inv_ivusrest_mla", 
                "inv_ivusrest_imla_w",
               "inv_ivusrest_mla_ellip",
               "inv_ivusrest_mla_ln",
              #  "inv_ivusrest_ostial_a",
                # "inv_ivusrest_ostial_w",
                # "iFR_mean_rest",
                # "mid_systolic_ratio_mean_rest",
                # "iFR_mean_ado",
                # "mid_systolic_ratio_mean_ado",
               "inv_ffrado")

logistic_results_ffr_0.8 <- perform_logistic_regression("ffr_0.8", baseline, ivus_vars)
logistic_results_risk_mla_any <- perform_logistic_regression("risk_mla_any", baseline, ivus_vars)

results_list <- list(
  "ffr_0.8" = logistic_results_ffr_0.8,
  "risk_mla_any" = logistic_results_risk_mla_any
)

# Write the list of data frames to an Excel file
write_xlsx(results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/ivus_logistic_regression_results.xlsx")

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
write_xlsx(roc_results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/ivus_roc_analysis_results.xlsx")

############################################################################################################
# best cutoffs
# Fit models and get prediction data
mdl_mla_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla")[[1]]
mdl_mla_ellip_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ellip")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[1]]
mdl_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[1]]
mdl_ifr_rest <- simple_logistic_regression(baseline, "ffr_0.8", "iFR_mean_rest")[[1]]
mdl_ifr_ado <- simple_logistic_regression(baseline, "ffr_0.8", "iFR_mean_ado")[[1]]
mdl_imla_w <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_imla_w")[[1]]
mdl_ostial_a <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_a")[[1]]
mdl_ostial_w <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_w")[[1]]

prediction_data_mla <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla")[[2]]
prediction_data_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ellip")[[2]]
prediction_data_mla_ln <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[2]]
prediction_data_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[2]]
prediction_data_ifr_rest <- simple_logistic_regression(baseline, "ffr_0.8", "iFR_mean_rest")[[2]]
prediction_data_ifr_ado <- simple_logistic_regression(baseline, "ffr_0.8", "iFR_mean_ado")[[2]]
prediction_data_imla_w <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_imla_w")[[2]]
prediction_data_ostial_a <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_a")[[2]]
prediction_data_ostial_w <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_w")[[2]]

data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla) %>% drop_na()
roc_mla <- roc(data_ivus$ffr_0.8, mdl_mla_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_ellip) %>% drop_na()
roc_mla_ellip <- roc(data_ivus$ffr_0.8, mdl_mla_ellip_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ffrado) %>% drop_na()
roc_ffrado <- roc(data_ivus$ffr_0.8, mdl_ffrado$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, iFR_mean_rest) %>% drop_na()
roc_ifr_rest <- roc(data_ivus$ffr_0.8, mdl_ifr_rest$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, iFR_mean_ado) %>% drop_na()
roc_ifr_ado <- roc(data_ivus$ffr_0.8, mdl_ifr_ado$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_imla_w) %>% drop_na()
roc_imla_w <- roc(data_ivus$ffr_0.8, mdl_imla_w$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_ostial_a) %>% drop_na()
roc_ostial_a <- roc(data_ivus$ffr_0.8, mdl_ostial_a$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_ostial_w) %>% drop_na()
roc_ostial_w <- roc(data_ivus$ffr_0.8, mdl_ostial_w$fitted.values)

# Define test sets
test_set <- list(
  coords(roc_mla), coords(roc_mla_ellip),
  coords(roc_mla_ln), coords(roc_ffrado)
  # ,
  # coords(roc_ifr_rest), coords(roc_ifr_ado),
  # coords(roc_imla_w), coords(roc_ostial_a),
  # coords(roc_ostial_w)
)

prediction_data <- list(
  prediction_data_mla, prediction_data_ellip,
  prediction_data_mla_ln, prediction_data_ffrado
  # ,
  # prediction_data_ifr_rest, prediction_data_ifr_ado,
  # prediction_data_imla_w, prediction_data_ostial_a,
  # prediction_data_ostial_w
)

variables <- c(
  "IVUS rest MLA", "IVUS rest elliptic ratio",
  "IVUS rest MLA luminal narrowing", "FFR adenosine"
  # ,
  # "iFR rest", "iFR adenosine",
  # "IVUS rest MLA minor axis", "IVUS rest ostial area",
  # "IVUS rest ostial minor axis"
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
    
    if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla", "inv_ffrado", "iFR_mean_rest", "iFR_mean_ado")) {
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
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/"
for (i in 1:length(cutoff_tables)) {
  name <- variables[i]
  table <- cutoff_tables[[i]]
  write_csv(table, paste0(path, name, ".csv"))
}

roc_list_mla <- list(roc_mla_ln, roc_mla_ellip)
ggroc_mla <- ggroc(roc_list_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkblue", "darkred")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("IVUS") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC MLN =", round(roc_mla_ln$auc, 2)), color = "darkblue") +
    annotate("text", x = 0.7, y = 0.1, label = paste("AUC elliptic ratio =", round(roc_mla_ellip$auc, 2)), color = "darkred")

roc_list_ffr <- list(roc_mla, roc_ffrado)
ggroc_ffr <- ggroc(roc_list_ffr, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC FFRAdenosine =", round(roc_ffrado$auc, 2)), color = "darkblue")

ggsave(paste0(output_dir,"/roc_ivus.png"), ggroc_mla, width = 6, height = 5)
ggsave(paste0(output_dir,"/roc_ffrado.png"), ggroc_ffr, width = 6, height = 5)

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

      table <- if (ivus_var %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla", "inv_ffrado", "iFR_mean_rest", "iFR_mean_ado")) {
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
write_xlsx(roc_summary_df_ffr0.8_rounded, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/roc_summary_ffr0.8_rounded_summary.xlsx")


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

baseline <- baseline %>% mutate(
  ifr_rest0.89 = ifelse(iFR_mean_rest <= 0.89, 1, 0),
  ifr_ado0.8 = ifelse(iFR_mean_ado <= 0.8, 1, 0),
  ifrado0.8_ffrado0.8 = ifelse(inv_ffrado <= 0.8 | inv_ffrado <= 0.8, 1, 0)
)

table <- baseline %>%
  select(ifr_rest0.89, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

table <- baseline %>%
  select(ifr_ado0.8, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

table <- baseline %>%
  select(ifrado0.8_ffrado0.8, ffr_0.8) %>%
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

roc_ifrrest <- roc(baseline$ffr_0.8, baseline$iFR_mean_rest)
ggroc(roc_ifrrest, legacy.axes = TRUE, color = "darkblue") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("iFR") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC iFR =", round(roc_ifrrest$auc, 2)))

roc_ifr_ado <- roc(baseline$ffr_0.8, baseline$iFR_mean_ado)
ggroc(roc_ifr_ado, legacy.axes = TRUE, color = "darkblue") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("iFR adenosine") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC iFR adenosine =", round(roc_ifr_ado$auc, 2)))

# forest plot
# for hemodynamic and anatomic relevance as predictor
forest_df_bsa <- data.frame(
    var = c(
        "Minimal lumen area",
        "Minimal lumen elliptic ratio",
        "Maximal lumen narrowing [%]",
        "FFR adenosine"
        # ,
        # "iFR rest",
        # "iFR adenosine"
    ),
    odds_ratio = c(
        exp(mdl_mla_ffr$coefficients[2]),
        exp(mdl_mla_ellip_ffr$coefficients[2]),
        exp(mdl_mla_ln_ffr$coefficients[2]),
        exp(mdl_ffrado$coefficients[2])
        # ,
        # exp(mdl_ifr_rest$coefficients[2]),
        # exp(mdl_ifr_ado$coefficients[2])
    ),
    lower_ci = c(
        exp(confint(mdl_mla_ffr))[2],
        exp(confint(mdl_mla_ellip_ffr))[2],
        exp(confint(mdl_mla_ln_ffr))[2], 
        exp(confint(mdl_ffrado))[2]
        # , 
        # exp(confint(mdl_ifr_rest))[2],
        # exp(confint(mdl_ifr_ado))[2]
    ),
    upper_ci = c(
        exp(confint(mdl_mla_ffr))[4],
        exp(confint(mdl_mla_ellip_ffr))[4],
        exp(confint(mdl_mla_ln_ffr))[4], 
        exp(confint(mdl_ffrado))[4]
        # ,
        # exp(confint(mdl_ifr_rest))[4],
        # exp(confint(mdl_ifr_ado))[4]
    )
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

cross <- baseline %>% drop_na(ffr_0.8, all_of(ivus_vars))
# cross <- cross %>% select(-inv_ffrado)
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
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/figures/auc_convergence.png", width = 6, height = 4)

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
min(no_treatment$inv_ffrdobu, na.rm = T)
max(no_treatment$inv_ffrdobu, na.rm = T)

treatment <- treatment %>% filter(inv_ffrdobu <= 0.8)
median(treatment$inv_ffrdobu, na.rm = T)
quantile(treatment$inv_ffrdobu, c(0.25, 0.75), na.rm = T)
min(treatment$inv_ffrdobu, na.rm = T)
max(treatment$inv_ffrdobu, na.rm = T)

wilcox.test(no_treatment$inv_ffrdobu, treatment$inv_ffrdobu, paired = F)

test <- baseline %>% filter(!is.na(inv_ffrado)) %>% select(patient_id, inv_ffrado, inv_ffrdobu, inv_ivusrest_mla)

test %>% filter(inv_ffrado <= 0.8)

test %>% filter(inv_ffrado > 0.8) %>% print(n = 50)

test %>% filter(inv_ffrado > 0.8 & inv_ffrdobu <= 0.8)

chisq.test(table(above0.8$inv_access_init___0), table(below0.8$inv_access_init___0))

# invasive changes above0.8
# Define the variables for comparison
rest_vars <- c("inv_ffrado", 
                "iFR_mean_rest",
                "inv_rest_hr", 
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean")

dobu_vars <- c("inv_ffrdobu", 
                "iFR_mean_dobu",
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
    rest_shapiro <- shapiro.test(as.numeric(above0.8[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(above0.8[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(above0.8[[rest_var]], above0.8[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(above0.8[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(above0.8[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(above0.8[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(above0.8[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(above0.8[[rest_var]], above0.8[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(above0.8[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(above0.8[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(above0.8[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(above0.8[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(above0.8[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(above0.8[[dobu_var]], 0.75, na.rm = TRUE)
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
write.csv(results_df, "statistical_analysis/data/invasive_change_above0.8.csv", row.names = FALSE)

# invasive changes below0.8
# Define the variables for comparison
rest_vars <- c("inv_ffrado", 
                "iFR_mean_rest",
                "inv_rest_hr", 
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean")

dobu_vars <- c("inv_ffrdobu", 
                "iFR_mean_dobu",
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
    rest_shapiro <- shapiro.test(as.numeric(below0.8[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(below0.8[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(below0.8[[rest_var]], below0.8[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(below0.8[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(below0.8[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(below0.8[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(below0.8[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(below0.8[[rest_var]], below0.8[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(below0.8[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(below0.8[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(below0.8[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(below0.8[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(below0.8[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(below0.8[[dobu_var]], 0.75, na.rm = TRUE)
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
write.csv(results_df, "statistical_analysis/data/invasive_change_below0.8.csv", row.names = FALSE)

graphical_mla <- ggroc(roc_mla, legacy.axes = TRUE, color = "black") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen area [mm²]") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC =", round(roc_mla$auc, 2)), color = "black")

graphical_ffr <- ggroc(roc_ffrado, legacy.axes = TRUE, color = "black") +
    scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("FFRadenosine") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC =", round(roc_ffrado$auc, 2)), color = "black")

ggsave(paste0(output_dir,"/graphical_mla.png"), graphical_mla, width = 2.5, height = 2.5)
ggsave(paste0(output_dir,"/graphical_ffr.png"), graphical_ffr, width = 2.5, height = 2.5)

test <- baseline %>% select(record_id, ffr_0.8, inv_ivusrest_mla, inv_ffrado)

write.csv(test, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/statistical_analysis/data/test.csv", row.names = FALSE)

# baseline IVUS variables to predict FFR change
resting_pressure_vars <- c("iFR_mean_rest", 
                "mid_systolic_ratio_mean_rest",
               "pdpa_mean_rest",
               "iFR_mean_ado",
               "mid_systolic_ratio_mean_ado",
                "pdpa_mean_ado",
                "integral_aortic_rest",
                "integral_distal_rest",
                "integral_diff_rest",
                "diastolic_integral_aortic_rest",
               "diastolic_integral_distal_rest",
               "diastolic_integral_diff_rest",
               "systolic_integral_aortic_rest",
               "systolic_integral_distal_rest",
               "systolic_integral_diff_rest",
               "integral_aortic_ado",
                "integral_distal_ado",
                "integral_diff_ado",
                "diastolic_integral_aortic_ado",
                "diastolic_integral_distal_ado",
                "diastolic_integral_diff_ado",
                "systolic_integral_aortic_ado",
                "systolic_integral_distal_ado",
                "systolic_integral_diff_ado")

logistic_results_ffr_0.8 <- perform_logistic_regression("ffr_0.8", baseline, resting_pressure_vars)

baseline %>% 
  select(inclusion_date, inv_date) %>% 
  arrange(inv_date) %>% 
  # filter(inv_date > "2020-01-01" & inv_date < "2025-01-01") %>%
  mutate(year = lubridate::year(inv_date)) %>%
  group_by(year) %>%
  summarise(count = n()) %>%
  View()

baseline <- baseline %>% mutate(inv_date = as.Date(inv_date))

ggplot(baseline, aes(x = inv_date)) +
  geom_histogram(binwidth = 120) +
  theme_minimal() +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", limits = as.Date(c("2019-01-01", "2025-01-01"))) +
  labs(title = "Distribution of Invasive Dates",
       x = "Invasive Date",
       y = "Count")
