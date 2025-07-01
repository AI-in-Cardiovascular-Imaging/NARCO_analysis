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
heartflow <- read_excel("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/NARCO_FFRct.xlsx")

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

baseline <- baseline %>% mutate(
  ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0)
)

heartflow <- heartflow %>% 
  mutate(
    FFRct_wire = ifelse(FFRct_wire == "NA" | is.na(FFRct_wire), FFRct_prox, FFRct_wire)
  )

baseline <- baseline %>% select(patient_id, inv_ffrado, inv_ffrdobu)
baseline <- baseline %>% rename(record_id = patient_id)
baseline <- baseline %>% mutate(record_id = as.double(gsub("NARCO_", "", record_id)))
heartflow <- heartflow %>% select(`NARCO ID`, FFRct_wire, FFRct_prox, FFRct_mid, FFRct_dist)
heartflow <- heartflow %>% rename(record_id = `NARCO ID`)
heartflow <- heartflow %>% mutate(FFRct_wire = as.numeric(FFRct_wire),
                                  FFRct_prox = as.numeric(FFRct_prox),
                                  FFRct_mid = as.numeric(FFRct_mid),
                                  FFRct_dist = as.numeric(FFRct_dist))
# remove NARCO_ from record_id and make numeric
heartflow$record_id <- as.double(gsub("NARCO_", "", heartflow$record_id))

data <- left_join(baseline, heartflow, by = "record_id")
data <- data %>% mutate(inv_ffrado_0.8 = ifelse(inv_ffrado < 0.8, 1, 0),
                        inv_ffrdobu_0.8 = ifelse(inv_ffrdobu < 0.8, 1, 0),
                        ffrct_0.8 = ifelse(FFRct_wire < 0.8, 1, 0),
                        ffrct_prox_0.8 = ifelse(FFRct_prox < 0.8, 1, 0),
                        ffrct_mid_0.8 = ifelse(FFRct_mid < 0.8, 1, 0),
                        ffrct_dist_0.8 = ifelse(FFRct_dist < 0.8, 1, 0))

scatter <- ggplot(data, aes(x = inv_ffrado, y = FFRct_wire)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRadenosine", y = "FFRct") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRadenosine vs FFRct") +
  theme_minimal()

scatter_ado_ctffrprox <- ggplot(data, aes(x = inv_ffrado, y = FFRct_prox)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRadenosine", y = "FFRct_prox") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRadenosine vs FFRct_prox") +
  theme_minimal()

scatter_ado_ctffrmid <- ggplot(data, aes(x = inv_ffrado, y = FFRct_mid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRadenosine", y = "FFRct_mid") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRadenosine vs FFRct_mid") +
  theme_minimal()

scatter_dobu <- ggplot(data, aes(x = inv_ffrdobu, y = FFRct_wire)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRdobutamine", y = "FFRct") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRdobutamine vs FFRct") +
  theme_minimal()

scatter_dobu_ctffrprox <- ggplot(data, aes(x = inv_ffrdobu, y = FFRct_prox)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRdobutamine", y = "FFRct_prox") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRdobutamine vs FFRct_prox") +
  theme_minimal()

scatter_dobu_ctffrmid <- ggplot(data, aes(x = inv_ffrdobu, y = FFRct_mid)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(x = "FFRdobutamine", y = "FFRct_mid") +
  xlim(0.7, 1) +
  ylim(0.7, 1) +
  ggtitle("FFRdobutamine vs FFRct_mid") +
  theme_minimal()

ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrado.png", scatter, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrdobu.png", scatter_dobu, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrado_prox.png", scatter_ado_ctffrprox, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrado_mid.png", scatter_ado_ctffrmid, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrdobu_prox.png", scatter_dobu_ctffrprox, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/ffrct_vs_ffrdobu_mid.png", scatter_dobu_ctffrmid, width = 6, height = 6)

roc_ado_ffrct <- roc(data$inv_ffrado_0.8, data$FFRct_wire)
roc_ado_ffrctprox <- roc(data$inv_ffrado_0.8, data$FFRct_prox)
roc_ado_ffrctmid <- roc(data$inv_ffrado_0.8, data$FFRct_mid)
roc_dobu_ffrct <- roc(data$inv_ffrdobu_0.8, data$FFRct_wire)
roc_dobu_ffrctprox <- roc(data$inv_ffrdobu_0.8, data$FFRct_prox)
roc_dobu_ffrctmid <- roc(data$inv_ffrdobu_0.8, data$FFRct_mid)

# summarize the roc model
roc_ado_ffrct
roc_ado_ffrctprox
roc_ado_ffrctmid
roc_dobu_ffrct
roc_dobu_ffrctprox
roc_dobu_ffrctmid

# get auc value
auc_ffrct <- roc_ado_ffrct$auc
auc_ffrctprox <- roc_ado_ffrctprox$auc
auc_ffrctmid <- roc_ado_ffrctmid$auc
auc_ffrct_dobu <- roc_dobu_ffrct$auc
auc_ffrctprox_dobu <- roc_dobu_ffrctprox$auc
auc_ffrctmid_dobu <- roc_dobu_ffrctmid$auc

ggroc_ado_ffrct <- ggroc(roc_ado_ffrct) +
  ggtitle("ROC curve for FFRct") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrct, 2))) +
  theme_minimal()

ggroc_ado_ffrctprox <- ggroc(roc_ado_ffrctprox) +
  ggtitle("ROC curve for FFRct_prox") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrctprox, 2))) +
  theme_minimal()

ggroc_ado_ffrctmid <- ggroc(roc_ado_ffrctmid) +
  ggtitle("ROC curve for FFRct_mid") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrctmid, 2))) +
  theme_minimal()

ggroc_dobu_ffrct <- ggroc(roc_dobu_ffrct) +
  ggtitle("ROC curve for FFRct") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrct_dobu, 2))) +
  theme_minimal()

ggroc_dobu_ffrctprox <- ggroc(roc_dobu_ffrctprox) +
  ggtitle("ROC curve for FFRct_prox") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrctprox_dobu, 2))) +
  theme_minimal()

ggroc_dobu_ffrctmid <- ggroc(roc_dobu_ffrctmid) +
  ggtitle("ROC curve for FFRct_mid") +
  annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(auc_ffrctmid_dobu, 2))) +
  theme_minimal()

ggsave("C:/Users/ansel/Downloads/roc_ffrct.png", roc, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/roc_ffrct_prox.png", roc_ado_ffrctprox, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/roc_ffrct_mid.png", roc_ado_ffrctmid, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/roc_ffrct_dobu.png", roc_dobu_ffrct, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/roc_ffrct_prox_dobu.png", roc_dobu_ffrctprox, width = 6, height = 6)
ggsave("C:/Users/ansel/Downloads/roc_ffrct_mid_dobu.png", roc_dobu_ffrctmid, width = 6, height = 6)


# confusion matrix
table <- data %>%
  select(inv_ffrado_0.8, ffrct_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

table <- data %>%
  select(inv_ffrado_0.8, ffrct_prox_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result_prox <- summary(confusion, event_level = "second")

table <- data %>%
  select(inv_ffrado_0.8, ffrct_mid_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result_mid <- summary(confusion, event_level = "second")

table <- data %>%
  select(inv_ffrdobu_0.8, ffrct_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result_dobu <- summary(confusion, event_level = "second")

table <- data %>%
  select(inv_ffrdobu_0.8, ffrct_prox_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result_prox_dobu <- summary(confusion, event_level = "second")

table <- data %>%
  select(inv_ffrdobu_0.8, ffrct_mid_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result_mid_dobu <- summary(confusion, event_level = "second")

write_xlsx(list(
  "FFRado versus CT-FFR (wire)" = result,
  "Confusion Matrix FFRado versus CT-FFR (wire)" = as.data.frame(confusion$table),
  "FFRado versus CT-FFR (Prox)" = result_prox,
  "Confusion Matrix FFRado versus CT-FFR (Prox)" = as.data.frame(confusion$table),
  "FFRado versus CT-FFR (Mid)" = result_mid,
  "Confusion Matrix FFRado versus CT-FFR (Mid)" = as.data.frame(confusion$table),
  "FFRdobu versus CT-FFR (wire)" = result_dobu,
  "Confusion Matrix FFRdobu versus CT-FFR (wire)" = as.data.frame(confusion$table),
  "FFRdobu versus CT-FFR (Prox)" = result_prox_dobu,
  "Confusion Matrix FFRdobu versus CT-FFR (Prox)" = as.data.frame(confusion$table),
  "FFRdobu versus CT-FFR (Mid)" = result_mid_dobu,
  "Confusion Matrix FFRdobu versus CT-FFR (Mid)" = as.data.frame(confusion$table)
), "C:/Users/ansel/Downloads/confusion_matrix.xlsx")
