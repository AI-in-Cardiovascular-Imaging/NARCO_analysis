library(tidyverse)
library(broom)
library(ggplot2)
library(yardstick)
library(pROC)
library(randomForest)
library(car)
library(readxl)
library(lmtest)
library(corrplot)
library(yardstick)

baseline <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")
baseline <- baseline %>%
  filter(record_id != 110, record_id != 111, record_id != 5, record_id != 75, record_id != 109) # left aaoca and <18 years

baseline <- baseline %>% mutate(ffr_0.8 = case_when(
    ffr_0.8 == "yes" ~ 1,
    ffr_0.8 == "no" ~ 0,
    TRUE ~ NA_real_
))

# functions needed
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

logistic_regression_categorical <- function(baseline_data, response, explanatory, family = "binomial") {
    data <- baseline_data %>% 
        select(!!sym(response), !!sym(explanatory)) %>% 
        drop_na()

    formula <- as.formula(paste(response, "~", explanatory))
    mdl <- glm(formula, data = data, family = family)
    odds_ratio <- exp(coef(mdl))[2]
    lower_ci <- exp(confint(mdl))[2]
    upper_ci <- exp(confint(mdl))[4]
    print(formula)
    tryCatch(
        {
            outcome <- table(baseline_data[[explanatory]], baseline_data[[response]])
            confusion <- conf_mat(outcome)
            plot_acc <- autoplot(confusion)
        },
        error = function(e) {
            print("Unfactor the explanatory variable first for confusion matrix and accuracy plot!")
        }
    )
    return(list(mdl, odds_ratio, lower_ci, upper_ci, confusion, plot_acc))
}

roc_table_creater <- function(list_roc, list_prediction_data, wgt) {
    roc_table <- data.frame()
    for (i in 1:length(list_roc)) {
        roc <- list_roc[[i]]
        prediction_data <- list_prediction_data[[i]]
        best_coords <- coords(roc, x = "best", best.method = "youden", ret = "all")
        wgt_coords <- coords(roc, x = "best", best.method = "youden", best.weights = wgt, ret = "all")
        dist_coords <- coords(roc, x = "best", best.method = "closest.topleft", ret = "all")

        value_cutoff <- prediction_data %>% filter(prediction_data[,2] > best_coords$threshold)
        value_cutoff_wgt <- prediction_data %>% filter(prediction_data[,2] > wgt_coords$threshold)
        value_cutoff_opt <- prediction_data %>% filter(prediction_data[,2] > dist_coords$threshold)
        
        roc_table <- rbind(roc_table, data.frame(
            "Predictor" = names(list_roc)[i],
            "AUC" = round(roc$auc, 3),
            "Best threshold" = round(best_coords$threshold, 3),
            "Accuracy" = round(best_coords$accuracy*100, 0),
            "Best sensitivity" = round(best_coords$sensitivity *100, 0),
            "Best specificity" = round(best_coords$specificity*100, 0),
            "Best PPV" = round(best_coords$ppv*100, 0),
            "Best NPV" = round(best_coords$npv*100, 0),
            "Cutoff value" = round(max(value_cutoff[,1]), 2),
            "Weighted threshold" = round(wgt_coords$threshold, 3),
            "Weighted accuracy" = round(wgt_coords$accuracy*100, 0),
            "Weighted sensitivity" = round(wgt_coords$sensitivity*100, 0),
            "Weighted specificity" = round(wgt_coords$specificity*100, 0),
            "Weighted PPV" = round(wgt_coords$ppv*100, 0),
            "Weighted NPV" = round(wgt_coords$npv*100, 0),
            "Weighted cutoff value" = round(max(value_cutoff_wgt[,1]), 2),
            "Distance threshold" = round(dist_coords$threshold, 3),
            "Distance accuracy" = round(dist_coords$accuracy*100, 0),
            "Distance sensitivity" = round(dist_coords$sensitivity*100, 0),
            "Distance specificity" = round(dist_coords$specificity*100, 0),
            "Distance PPV" = round(dist_coords$ppv*100, 0),
            "Distance NPV" = round(dist_coords$npv*100, 0),
            "Distance cutoff value" = round(max(value_cutoff_opt[,1]), 2)
        ))
    }
    return(roc_table)
}

# all continuous models
# ostium
all_ost_a_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a")
mdl_ost_a_bsa_ffr <- all_ost_a_bsa[[1]] # significant
all_ost_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")
mdl_ost_ellip_ffr <- all_ost_ellip[[1]] # significant
all_ost_h_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_h")
mdl_ost_h_bsa_ffr <- all_ost_h_bsa[[1]] # not significant!!!!
all_ost_w_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_w")
mdl_ost_w_bsa_ffr <- all_ost_w_bsa[[1]] # significant
all_ost_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")
mdl_ost_pn_ffr <- all_ost_pn[[1]] # significant

# mla
all_mla_a_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")
mdl_mla_a_bsa_ffr <- all_mla_a_bsa[[1]] # significant
all_mla_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")
mdl_mla_ellip_ffr <- all_mla_ellip[[1]] # not significant!!!
all_mla_h_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_h")
mdl_mla_h_bsa_ffr <- all_mla_h_bsa[[1]] # not significant!!!
all_mla_w_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_w")
mdl_mla_w_bsa_ffr <- all_mla_w_bsa[[1]] # significant
all_mla_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_pn_dist")
mdl_mla_pn_ffr <- all_mla_pn[[1]] # significant

# rest of the variables
all_imc_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length")
mdl_imc_bsa_ffr <- all_imc_bsa[[1]] # not significant!!!
all_aa_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree")
mdl_aa_bsa_ffr <- all_aa_bsa[[1]] # not significant!!!

all_ost_y_pos <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_stj_rca_bsa")
mdl_ost_y_pos_ffr <- all_ost_y_pos[[1]] # not significant!!!
all_ost_y_pos_scaled <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_stj_rca_scaled")
mdl_ost_y_pos_scaled_ffr <- all_ost_y_pos_scaled[[1]] # not significant!!!

all_ost_x_pos <- simple_logistic_regression(baseline, "ffr_0.8", "bogen_rca_bsa")
mdl_ost_x_pos_ffr <- all_ost_x_pos[[1]] # not significant!!!
all_ost_x_pos_scaled <- simple_logistic_regression(baseline, "ffr_0.8", "bogen_rca_bsa_scaled")
mdl_ost_x_pos_scaled_ffr <- all_ost_x_pos_scaled[[1]] # not significant!!!

# all categorical models
# ostium
mdl_slo_ffr <- glm(ffr_0.8 ~ caa_slo, data = baseline, family = "binomial") # perfect seperation
summary(mdl_slo_ffr) # not significant!!! AIC 49.59
lrtest(mdl_slo_ffr, mdl_ost_a_bsa_ffr) # mdl_ost_a_bsa_ffr is better
lrtest(mdl_slo_ffr, mdl_ost_ellip_ffr) # mdl_ost_ellip_ffr is better
# mdl_ellip_ffr <- glm(ffr_0.8 ~ caa_elliptic, data = baseline, family = "binomial")
# summary(mdl_ellip_ffr) # not significant!!! perfect seperation
mdl_pn_ffr <- glm(ffr_0.8 ~ caa_pn, data = baseline, family = "binomial")
summary(mdl_pn_ffr) # not significant!!! AIC 62.0
exp(coef(mdl_pn_ffr))[2]
exp(confint(mdl_pn_ffr))[2]
exp(confint(mdl_pn_ffr))[4]
summary(mdl_ost_pn_ffr) # significant!!! AIC 54.1
summary(mdl_mla_pn_ffr) # not significant!!! AIC 59.4
lrtest(mdl_pn_ffr, mdl_ost_pn_ffr) # mdl_ost_pn_ffr is better
lrtest(mdl_pn_ffr, mdl_mla_pn_ffr) # mdl_mla_pn_ffr is better

# forest plot
# bsa adjusted
forest_df_bsa <- data.frame(
    var = c(
            "Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial major axis (BSA) [mm]", "Ostial minor axis (BSA) [mm]", "Ostial lumen narrowing [%]", 
            "Minimal lumen area (BSA) [mm2]", "MLA elliptic ratio", 
            "IM major axis (BSA) [mm]", "IM minor axis (BSA) [mm]", 
            "Minimal lumen narrowing [%]", "Intramural length (BSA) [mm]", "Take-off angle (BSA) [deg]"),
    odds_ratio = c(
            exp(mdl_ost_a_bsa_ffr$coefficients[2]), 
            exp(mdl_ost_ellip_ffr$coefficients[2]), 
            exp(mdl_ost_h_bsa_ffr$coefficients[2]), 
            exp(mdl_ost_w_bsa_ffr$coefficients[2]), 
            exp(mdl_ost_pn_ffr$coefficients[2]), 
            exp(mdl_mla_a_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_ellip_ffr$coefficients[2]), 
            exp(mdl_mla_h_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_w_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_pn_ffr$coefficients[2]),
            exp(mdl_imc_bsa_ffr$coefficients[2]), 
            exp(mdl_aa_bsa_ffr$coefficients[2])),
    lower_ci = c(
            exp(confint(mdl_ost_a_bsa_ffr))[2], 
            exp(confint(mdl_ost_ellip_ffr))[2], 
            exp(confint(mdl_ost_h_bsa_ffr))[2], 
            exp(confint(mdl_ost_w_bsa_ffr))[2], 
            exp(confint(mdl_ost_pn_ffr))[2], 
            exp(confint(mdl_mla_a_bsa_ffr))[2], 
            exp(confint(mdl_mla_ellip_ffr))[2], 
            exp(confint(mdl_mla_h_bsa_ffr))[2], 
            exp(confint(mdl_mla_w_bsa_ffr))[2], 
            exp(confint(mdl_mla_pn_ffr))[2], 
            exp(confint(mdl_imc_bsa_ffr))[2], 
            exp(confint(mdl_aa_bsa_ffr))[2]),
    upper_ci = c(
            exp(confint(mdl_ost_a_bsa_ffr))[4], 
            exp(confint(mdl_ost_ellip_ffr))[4], 
            exp(confint(mdl_ost_h_bsa_ffr))[4], 
            exp(confint(mdl_ost_w_bsa_ffr))[4], 
            exp(confint(mdl_ost_pn_ffr))[4], 
            exp(confint(mdl_mla_a_bsa_ffr))[4], 
            exp(confint(mdl_mla_ellip_ffr))[4], 
            exp(confint(mdl_mla_h_bsa_ffr))[4], 
            exp(confint(mdl_mla_w_bsa_ffr))[4], 
            exp(confint(mdl_mla_pn_ffr))[4], 
            exp(confint(mdl_imc_bsa_ffr))[4], 
            exp(confint(mdl_aa_bsa_ffr))[4])
)

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_cut_ffr_bsa <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 7, 1), labels = parse(text = seq(0, 7, 1))) +
    geom_vline(xintercept = 6, linetype = "solid", color = "white", lwd = 2) +
    geom_segment(aes(x = 6, xend = 7, y = 3, yend = 3), linetype = "dashed") +        # adds dotted line for CI, Chris preferred arrow
    # geom_segment(aes(x = 6, xend = 7, y = 11, yend = 11), linetype = "dashed") +
    geom_segment(aes(x = 6, xend = 7, y = 8, yend = 8), linetype = "dashed") +
    # geom_segment(aes(x = 5.5, xend = 7, y = 3, yend = 3), arrow = arrow(length = unit(0.5, 'cm')), lwd = 0.3) +
    # geom_segment(aes(x = 5.5, xend = 7, y = 11, yend = 11), arrow = arrow(length = unit(0.5, 'cm')), lwd = 0.4) +
    # geom_segment(aes(x = 5.5, xend = 7, y = 8, yend = 8), arrow = arrow(length = unit(0.5, 'cm')), lwd = 0.4) +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_ccta_ffr.png", forest_cut_ffr_bsa, width = 4, height = 5)

# regression models #################################################################################################
# exploration of continuous variables that were significant
## ostium
### area
summary(mdl_ost_a_bsa_ffr) # significant
prediction_data_ost_a <- all_ost_a_bsa[[2]]
print(prediction_data_ost_a, n = 50)
plot_log_ost_a <- all_ost_a_bsa[[3]]
plot_odds_ost_a <- all_ost_a_bsa[[4]]
confusion_ost_a <- all_ost_a_bsa[[5]]
summary(confusion_ost_a, event_level = "second")
plot_acc_ost_a <- all_ost_a_bsa[[6]]

### elliptic
summary(mdl_ost_ellip_ffr) # significant
prediction_data_ost_ellip <- all_ost_ellip[[2]]
print(prediction_data_ost_ellip, n = 50)
plot_log_ost_ellip <- all_ost_ellip[[3]]
plot_odds_ost_ellip <- all_ost_ellip[[4]]

### width
summary(mdl_ost_w_bsa_ffr) # significant
prediction_data_ost_w <- all_ost_w_bsa[[2]]
print(prediction_data_ost_w, n = 50)
plot_log_ost_w <- all_ost_w_bsa[[3]]
plot_odds_ost_w <- all_ost_w_bsa[[4]]

### proximal narrowing
summary(mdl_ost_pn_ffr) # not significant
prediction_data_ost_pn <- all_ost_pn[[2]]
print(prediction_data_ost_pn, n = 50)
plot_log_ost_pn <- all_ost_pn[[3]]
plot_odds_ost_pn <- all_ost_pn[[4]]

## mla
### area
summary(mdl_mla_a_bsa_ffr) # significant
prediction_data_mla <- all_mla_a_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla <- all_mla_a_bsa[[3]]
plot_odds_mla <- all_mla_a_bsa[[4]]
confusion_mla <- all_mla_a_bsa[[5]]
summary(confusion_mla, event_level = "second")
plot_acc_mla <- all_mla_a_bsa[[6]]

### elliptic
summary(mdl_mla_ellip_ffr) # significant
prediction_data_mla_ellip <- all_mla_ellip[[2]]
print(prediction_data_mla_ellip, n = 50)
plot_log_mla_ellip <- all_mla_ellip[[3]]
plot_odds_mla_ellip <- all_mla_ellip[[4]]

### width
summary(mdl_mla_w_bsa_ffr) # significant
prediction_data_mla_w <- all_mla_w_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_w <- all_mla_w_bsa[[3]]
plot_odds_mla_w <- all_mla_w_bsa[[4]]

### proximal narrowing
summary(mdl_mla_pn_ffr) # not significant
prediction_data_mla_pn <- all_mla_pn[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_pn <- all_mla_pn[[3]]
plot_odds_mla_pn <- all_mla_pn[[4]]

# roc analysis ######################################################################################################
## ostium ###########################################################################################################
### area
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_a) %>% drop_na()
roc_ost_a <- roc(data_ostium$ffr_0.8, mdl_ost_a_bsa_ffr$fitted.values)

### elliptic
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_elliptic) %>% drop_na()
roc_ost_ellip <- roc(data_ostium$ffr_0.8, mdl_ost_ellip_ffr$fitted.values)

### width
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_w) %>% drop_na()
roc_ost_w <- roc(data_ostium$ffr_0.8, mdl_ost_w_bsa_ffr$fitted.values)

### proximal narrowing
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_pn) %>% drop_na()
roc_ost_pn <- roc(data_ostium$ffr_0.8, mdl_ost_pn_ffr$fitted.values)

## mla #############################################################################################################
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_a) %>% drop_na()
roc_mla <- roc(data_mla$ffr_0.8, mdl_mla_a_bsa_ffr$fitted.values)

### width
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_w) %>% drop_na()
roc_mla_w <- roc(data_mla$ffr_0.8, mdl_mla_w_bsa_ffr$fitted.values)

### proximal narrowing
data_mla <- baseline %>% select(ffr_0.8, ccta_pn_dist) %>% drop_na()
roc_mla_pn <- roc(data_mla$ffr_0.8, mdl_mla_pn_ffr$fitted.values)

# test different thresholds
list_prediction_data <- list("CCTA ostial area (BSA)" = all_ost_a_bsa[[2]], 
                            "CCTA ostial elliptic" = all_ost_ellip[[2]],
                            "CCTA ostial minor axis (BSA)" = all_ost_w_bsa[[2]],
                            "CCTA ostial proximal narrowing" = all_ost_pn[[2]],
                            "CCTA MLA (BSA)" = all_mla_a_bsa[[2]],
                            "CCTA MLA width (BSA)" = all_mla_w_bsa[[2]],
                            "CCTA MLA proximal narrowing" = all_mla_pn[[2]])
list_roc <- list("CCTA ostial area (BSA)" = roc_ost_a, 
                "CCTA ostial elliptic" = roc_ost_ellip,
                "CCTA ostial width (BSA)" = roc_ost_w,
                "CCTA ostial proximal narrowing" = roc_ost_pn,
                "CCTA MLA (BSA)" = roc_mla,
                "CCTA MLA width (BSA)" = roc_mla_w,
                "CCTA MLA proximal narrowing" = roc_mla_pn)


roc_table <- roc_table_creater(list_roc, list_prediction_data, c(1, 38/50))
write_csv(roc_table, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/roc_table.csv")

# combined ROC plots
list_roc_ost <- list("CCTA ostial elliptic" = roc_ost_ellip,
                    "CCTA ostial proximal narrowing" = roc_ost_pn,
                    "CCTA ostial area (BSA)" = roc_ost_a, 
                    "CCTA ostial minor axis (BSA)" = roc_ost_w)

roc_ost_combined <- ggroc(list_roc_ost, legacy.axes = TRUE) +
    scale_color_manual(values = c("#1f1fb8", "#1c771c", "#fa9600", "#a31800")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial measurements") +
    theme_classic() +
    # add all AUCs with color
    annotate("text", x = 0.5, y = 0.1, label = paste0("AUC elliptic ratio: ", round(auc(roc_ost_ellip), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.2, label = paste0("AUC prox. narrowing: ", round(auc(roc_ost_pn), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.3, label = paste0("AUC area: ", round(auc(roc_ost_a), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.4, label = paste0("AUC minor axis: ", round(auc(roc_ost_w), 3)), size = 3, hjust = 0)

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_combined.png", roc_ost_combined, width = 8, height = 6)



# combined ROC plots
list_roc_ost <- list(
                    "CCTA intramural area (BSA)" = roc_mla, 
                    "CCTA intramural proximal narrowing" = roc_mla_pn,
                    "CCTA intramural minor axis (BSA)" = roc_mla_w)

roc_mla_combined <- ggroc(list_roc_ost, legacy.axes = TRUE) +
    scale_color_manual(values = c("#fa9600", "#a31800", "#1c771c", "black")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    # geom_point(aes(x = x_weighted_mla, y = y_weighted_mla, color = "cut-off"), size = 2) +
    # geom_point(aes(x = x_weighted_mla_w, y = y_weighted_mla_w, color = "cut-off"), size = 2) +
    # geom_point(aes(x = x_weighted_mla_pn, y = y_weighted_mla_pn, color = "cut-off"), size = 2) +
    ggtitle("Intramural measurements") +
    theme_classic() +
    # add all AUCs with color
    annotate("text", x = 0.5, y = 0.2, label = paste0("AUC prox. narrowing: ", round(auc(roc_mla_pn), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.3, label = paste0("AUC area: ", round(auc(roc_mla), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.4, label = paste0("AUC minor axis: ", round(auc(roc_mla_w), 3)), size = 3, hjust = 0)

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_mla_combined.png", roc_mla_combined, width = 8, height = 6)

# correlation matrix of all high risk features
high_risk_features <- baseline %>% select(ccta_ostial_a, ccta_mla_a, ccta_ostial_elliptic, ccta_mla_elliptic, ccta_ostial_h, ccta_mla_h, 
                                        ccta_ostial_w, ccta_mla_w, ccta_ostial_pn, ccta_pn_dist, 
                                        ccta_imc_length, ccta_aa_degree)

colnames(high_risk_features) <- c("OLA (BSA)", "IMLA (BSA)", "OLA elliptic-ratio", "IMLA elliptic-ratio", "OLA major axis (BSA)", "IMLA major axis (BSA)", 
                                  "OLA minor axis (BSA)", "IMLA minor axis (BSA)", "Ostial lumen narrowing", "IM lumen narrowing", "IM length", "Take-off angle")

correlation_matrix <- cor(high_risk_features, use = "pairwise.complete.obs")
correlation_matrix <- round(correlation_matrix, 2)
# save as excel
write.csv(correlation_matrix, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/correlation_matrix.csv")

M <- cor(high_risk_features, use = "pairwise.complete.obs")

# http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram
cor.mtest <- function(mat, ...) {
    mat <- as.matrix(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
            tmp <- cor.test(mat[, i], mat[, j], ...)
            p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
        }
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat) 
  p.mat
}

N <- length(high_risk_features) -1
p.mat <- cor.mtest(high_risk_features)
head(p.mat[, 1:N])
ids <- seq(1,N) 

# plot correlation matrix
png(height = 800, width = 800, file = "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/correlation_matrix.png")
fig_correlation_matrix <- 
    corrplot(M, 
    type="upper", 
    order="hclust", 
    tl.pos=c("td"), 
    method="color",  
    tl.cex = 1.1, tl.col = 'black', 
    tl.srt = 45, 
    addCoef.col = "black", 
    number.cex = 1.1,
    number.digits = 2,
    diag = FALSE, 
    p.mat = p.mat, 
    sig.level = 0.05)

dev.off()

colnames(M) <- colnames(high_risk_features)
rownames(M) <- colnames(high_risk_features)
png(height = 800, width = 800, file = "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/correlation_matrix_full.png")
fig_correlation_matrix_full <- 
    corrplot(M, 
    type="full", 
    order="hclust", 
    method="color",  
    tl.cex = 1.1, tl.col = 'black', 
    tl.srt = 45, 
    addCoef.col = "black", 
    number.cex = 1.1,
    number.digits = 2,
    diag = TRUE, 
    p.mat = p.mat, 
    sig.level = 0.05)
dev.off()


# repeat whole forest plot as normalized
# normalize variables in baseline
normalized <- baseline %>% mutate(
    ccta_ostial_a = (ccta_ostial_a - mean(ccta_ostial_a, na.rm = TRUE)) / sd(ccta_ostial_a, na.rm = TRUE),
    ccta_ostial_elliptic = (ccta_ostial_elliptic - mean(ccta_ostial_elliptic, na.rm = TRUE)) / sd(ccta_ostial_elliptic, na.rm = TRUE),
    ccta_ostial_w = (ccta_ostial_w - mean(ccta_ostial_w, na.rm = TRUE)) / sd(ccta_ostial_w, na.rm = TRUE),
    ccta_ostial_h = (ccta_ostial_h - mean(ccta_ostial_h, na.rm = TRUE)) / sd(ccta_ostial_h, na.rm = TRUE),
    ccta_ostial_pn = (ccta_ostial_pn - mean(ccta_ostial_pn, na.rm = TRUE)) / sd(ccta_ostial_pn, na.rm = TRUE),
    ccta_mla_a = (ccta_mla_a - mean(ccta_mla_a, na.rm = TRUE)) / sd(ccta_mla_a, na.rm = TRUE),
    ccta_mla_elliptic = (ccta_mla_elliptic - mean(ccta_mla_elliptic, na.rm = TRUE)) / sd(ccta_mla_elliptic, na.rm = TRUE),
    ccta_mla_w = (ccta_mla_w - mean(ccta_mla_w, na.rm = TRUE)) / sd(ccta_mla_w, na.rm = TRUE),
    ccta_mla_h = (ccta_mla_h - mean(ccta_mla_h, na.rm = TRUE)) / sd(ccta_mla_h, na.rm = TRUE),
    ccta_pn_dist = (ccta_pn_dist - mean(ccta_pn_dist, na.rm = TRUE)) / sd(ccta_pn_dist, na.rm = TRUE),
    ccta_imc_length = (ccta_imc_length - mean(ccta_imc_length, na.rm = TRUE)) / sd(ccta_imc_length, na.rm = TRUE),
    ccta_aa_degree = (ccta_aa_degree - mean(ccta_aa_degree, na.rm = TRUE)) / sd(ccta_aa_degree, na.rm = TRUE)
)

# ostium
all_ost_a_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_a")
mdl_ost_a_bsa_ffr <- all_ost_a_bsa[[1]] # significant
all_ost_ellip <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_elliptic")
mdl_ost_ellip_ffr <- all_ost_ellip[[1]] # significant
all_ost_h_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_h")
mdl_ost_h_bsa_ffr <- all_ost_h_bsa[[1]] # not significant!!!!
all_ost_w_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_w")
mdl_ost_w_bsa_ffr <- all_ost_w_bsa[[1]] # significant
all_ost_pn <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_pn")
mdl_ost_pn_ffr <- all_ost_pn[[1]] # significant

# mla
all_mla_a_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_a")
mdl_mla_a_bsa_ffr <- all_mla_a_bsa[[1]] # significant
all_mla_ellip <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_elliptic")
mdl_mla_ellip_ffr <- all_mla_ellip[[1]] # not significant!!!
all_mla_h_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_h")
mdl_mla_h_bsa_ffr <- all_mla_h_bsa[[1]] # not significant!!!
all_mla_w_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_w")
mdl_mla_w_bsa_ffr <- all_mla_w_bsa[[1]] # significant
all_mla_pn <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_pn_dist")
mdl_mla_pn_ffr <- all_mla_pn[[1]] # significant

# rest of the variables
all_imc_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_imc_length")
mdl_imc_bsa_ffr <- all_imc_bsa[[1]] # not significant!!!
all_aa_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_aa_degree")
mdl_aa_bsa_ffr <- all_aa_bsa[[1]] # not significant!!!


forest_df_bsa <- data.frame(
    var = c("Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial major axis (BSA) [mm]", "Ostial minor axis (BSA) [mm]", 
            "Ostial lumen narrowing [%]", "IM area (BSA) [mm2]", "IM elliptic ratio", "IM major axis (BSA) [mm]", "IM minor axis (BSA) [mm]", 
            "IM lumen narrowing [%]", "Intramural length (BSA) [mm]", "Take-off angle (BSA) [deg]"),
    odds_ratio = c(exp(mdl_ost_a_bsa_ffr$coefficients[2]), exp(mdl_ost_ellip_ffr$coefficients[2]), 
            exp(mdl_ost_h_bsa_ffr$coefficients[2]), exp(mdl_ost_w_bsa_ffr$coefficients[2]), 
            exp(mdl_ost_pn_ffr$coefficients[2]), exp(mdl_mla_a_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_ellip_ffr$coefficients[2]), exp(mdl_mla_h_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_w_bsa_ffr$coefficients[2]), exp(mdl_mla_pn_ffr$coefficients[2]),
            exp(mdl_imc_bsa_ffr$coefficients[2]), exp(mdl_aa_bsa_ffr$coefficients[2])),
    lower_ci = c(exp(confint(mdl_ost_a_bsa_ffr))[2], exp(confint(mdl_ost_ellip_ffr))[2], 
            exp(confint(mdl_ost_h_bsa_ffr))[2], exp(confint(mdl_ost_w_bsa_ffr))[2], 
            exp(confint(mdl_ost_pn_ffr))[2], exp(confint(mdl_mla_a_bsa_ffr))[2], 
            exp(confint(mdl_mla_ellip_ffr))[2], exp(confint(mdl_mla_h_bsa_ffr))[2], 
            exp(confint(mdl_mla_w_bsa_ffr))[2], exp(confint(mdl_mla_pn_ffr))[2], 
            exp(confint(mdl_imc_bsa_ffr))[2], exp(confint(mdl_aa_bsa_ffr))[2]),
    upper_ci = c(exp(confint(mdl_ost_a_bsa_ffr))[4], exp(confint(mdl_ost_ellip_ffr))[4], 
            exp(confint(mdl_ost_h_bsa_ffr))[4], exp(confint(mdl_ost_w_bsa_ffr))[4], 
            exp(confint(mdl_ost_pn_ffr))[4], exp(confint(mdl_mla_a_bsa_ffr))[4], 
            exp(confint(mdl_mla_ellip_ffr))[4], exp(confint(mdl_mla_h_bsa_ffr))[4], 
            exp(confint(mdl_mla_w_bsa_ffr))[4], exp(confint(mdl_mla_pn_ffr))[4], 
            exp(confint(mdl_imc_bsa_ffr))[4], exp(confint(mdl_aa_bsa_ffr))[4])
)

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_cut_ffr_bsa <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 17, 1), labels = parse(text = seq(0, 17, 1))) +
    geom_vline(xintercept = 6, linetype = "solid", color = "white", lwd = 2) +
    geom_segment(aes(x = 6, xend = 7, y = 8, yend = 8), linetype = "dashed") +
    geom_segment(aes(x = 6, xend = 7, y = 3, yend = 3), linetype = "dashed") +
    # geom_segment(aes(x = 6, xend = 7, y = 11, yend = 11), linetype = "dashed") +
    theme_classic()

####################################################################################################################
# distribution with different cut-offs
test_set <- list(coords(roc_ost_a), coords(roc_ost_ellip), 
            coords(roc_ost_w), coords(roc_ost_pn), 
            coords(roc_mla), coords(roc_mla_w), coords(roc_mla_pn))

prediction_data <- list(prediction_data_ost_a, prediction_data_ost_ellip, 
                        prediction_data_ost_w, prediction_data_ost_pn, 
                        prediction_data_mla, prediction_data_mla_w, prediction_data_mla_pn)

variables <- c("CCTA ostial area (BSA)", "CCTA ostial elliptic", 
                "CCTA ostial minor axis (BSA)", "CCTA ostial proximal narrowing", 
                "CCTA MLA (BSA)", "CCTA MLA width (BSA)", "CCTA MLA proximal narrowing")

thresholds <- list()

for (i in 1:length(test_set)){
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

cutoff_tables <- list()

for (i in 1:length(thresholds)){
    thresholds_testing <- thresholds[[i]]
    prediction <- prediction_data[[i]]
    name <- names(prediction)[1]

    table_thresholds <- tibble()
    for (j in 1:length(thresholds_testing)){
        threshold <- thresholds_testing[j]
        closest_index <- which.min(abs(prediction$ffr_0.8 - threshold))
        value <- prediction[[name]][closest_index]
        print(paste0("Threshold: ", threshold, "Value: ", value, " iteration: ", j))
        if (name %in% c("ccta_ostial_a", "ccta_mla_a", "ccta_ostial_w", "ccta_mla_w")){
            table <- baseline %>% 
            mutate(pred_label = ifelse(!!sym(name) > value, 0, 1)) %>% # careful sometimes it's > and sometimes <
            select(pred_label, ffr_0.8) %>% 
            drop_na() %>% 
            table()
        } else {
            table <- baseline %>% 
            mutate(pred_label = ifelse(!!sym(name) < value, 0, 1)) %>% # careful sometimes it's > and sometimes <
            select(pred_label, ffr_0.8) %>% 
            drop_na() %>% 
            table()
        }
        print(table)

        confusion <- conf_mat(table)
        result <- summary(confusion, event_level = "second")

        sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
        spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
        ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
        npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
        acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()
        
        df <- tibble(
                value = as.numeric(value),
                sensitivity =as.numeric(sens),
                specificity = as.numeric(spec),
                ppv = as.numeric(ppv),
                npv = as.numeric(npv),
                accuracy = as.numeric(acc),
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

# save all as csv
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/"
for (i in 1:length(cutoff_tables)){
    name <- names(cutoff_tables)[i]
    table <- cutoff_tables[[i]]
    write_csv(table, paste0(path, name, ".csv"))
}

baseline %>% select(ccta_ostial_a, ffr_0.8) %>% mutate(cutoff = ifelse(ccta_ostial_a > 4.13, 0, 1)) %>% drop_na() %>% View()
baseline %>% select(inv_ffrdobu, record_id, patient_id) %>% drop_na() %>% mutate(patho = ifelse(inv_ffrdobu <=0.8, 1, 0)) %>% View()

mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_ln")[[1]]

roc_imc <- roc(baseline$ffr_0.8, baseline$ccta_imc_length)
data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)
roc_list <- list(roc_mla, roc_imc, roc_mla_ln)
graphical_abstract <- ggroc(roc_list, legacy.axes = TRUE) + 
scale_color_manual(values = c("#1f1fb8", "red", "darkgreen")) + 
theme(legend.position = "none") + 
geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
ggtitle("Minmal lumen area (BSA) vs. IM length (BSA)") + 
theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/graphical_abstract.png", graphical_abstract, width = 6, height = 5)








############################################################################################################
# best cutoffs
# Fit models and get prediction data
# baseline <- baseline %>% mutate(
#     ffr_0.8 = as.numeric(ffr_0.8) - 1,
# )

baseline <- baseline %>% mutate(
    ffr_0.8 = ffr_0.8 + 1
)

mdl_mla_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_ln")[[1]]

prediction_data_mla <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")[[2]]
prediction_data_mla_ln <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_ln")[[2]]


data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_a) %>% drop_na()
roc_mla <- roc(data_ivus$ffr_0.8, mdl_mla_bsa_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, ccta_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)

# Define test sets
test_set <- list(
    coords(roc_mla), coords(roc_mla_ln)
)

prediction_data <- list(
    prediction_data_mla, prediction_data_mla_ln
)

variables <- c(
    "CCTA MLA (BSA)", "CCTA MLA luminal narrowing"
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
    
    if (name %in% c("ccta_mla_a")) {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) < value, 0, 1)) %>%
        select(pred_label, ffr_0.8) %>%
        drop_na() %>%
        table()
    } else {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) > value, 0, 1)) %>%
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

# save all as csv
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/"
for (i in 1:length(cutoff_tables)){
    name <- names(cutoff_tables)[i]
    table <- cutoff_tables[[i]]
    write_csv(table, paste0(path, name, ".csv"))
}

baseline %>% filter(!is.na(ffr_0.8)) %>% select(record_id, patient_id, ccta_mla_a, ccta_mla_a, ccta_dist_a, ccta_mla_ln, ffr_0.8) %>% View()

ggplot(baseline, aes(ccta_mla_a, inv_ivusrest_mla)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    scale_y_continuous(limits = c(1, 10, 1)) +
    scale_x_continuous(limits = c(1, 10, 1)) +
    theme_classic()
