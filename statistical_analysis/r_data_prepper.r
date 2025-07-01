library(yaml)
library(readxl)
library(tidyverse)


## DATA LOADING AND PREPPING ################################################################################
yaml <- yaml.load_file(file.path(getwd(), "config.yaml")) # only works after activating venv

# baseline <- read_excel(paste0(yaml$first_analysis$output_dir, "/baseline.xlsx"))
baseline <- read_excel("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis_team/dataframes/baseline_marius.xlsx")
invasive <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/IVUS_data.xlsx")
pressure <- read_excel("C:/Users/ansel/OneDrive/Dokumente/3_Research/results.xlsx")

# remove the first row, which is test data
baseline <- baseline[-1, ]

## MUTATE ##########################################################################################################
baseline <- baseline %>% 
  mutate(
    sports_classification = case_when(
      sports_static_component == 0 & sports_dynamic_component == 0 ~ 1,
      (sports_static_component == 1 & sports_dynamic_component == 0) | (sports_static_component == 0 & sports_dynamic_component ==1) ~ 2,
      (sports_static_component == 2 & sports_dynamic_component == 0) | (sports_static_component == 1 & sports_dynamic_component == 1) | (sports_static_component == 0 & sports_dynamic_component == 2) ~ 3,
      (sports_static_component == 2 & sports_dynamic_component == 1) | (sports_static_component == 1 & sports_dynamic_component ==2) ~ 4,
      sports_static_component == 2 & sports_dynamic_component == 2 ~ 5,
    ),
    sports_classification_ch = case_when(
      sports_static_component_ch == 0 & sports_dynamic_comp_ch == 0 ~ 1,
      (sports_static_component_ch == 1 & sports_dynamic_comp_ch == 0) | (sports_static_component_ch == 0 & sports_dynamic_comp_ch ==1) ~ 2,
      (sports_static_component_ch == 2 & sports_dynamic_comp_ch == 0) | (sports_static_component_ch == 1 & sports_dynamic_comp_ch == 1) | (sports_static_component_ch == 0 & sports_dynamic_comp_ch == 2) ~ 3,
      (sports_static_component_ch == 2 & sports_dynamic_comp_ch == 1) | (sports_static_component_ch == 1 & sports_dynamic_comp_ch ==2) ~ 4,
      sports_static_component_ch == 2 & sports_dynamic_comp_ch == 2 ~ 5,
    ),
    ffr_0.8 = ifelse(inv_ffrdobu <= 0.8, 1, 0),
    slo_cheezum = ccta_ostial_w / (2 * sqrt(ccta_dist_a / pi)), 
    right_0_left_1 = ifelse(caa_origin___0 == 1 | caa_origin___4 == 1, 0, 1),
    cad_any_rca = ifelse(inv_cad_loc___9 == 1 | inv_cad_loc___10 == 1 | inv_cad_loc___11 == 1, 1, 0),
    cad_any_lcx = ifelse(inv_cad_loc___6 == 1 | inv_cad_loc___7 == 1 | inv_cad_loc___8 == 1, 1, 0),
    cad_any_lad = ifelse(inv_cad_loc___0 == 1 | inv_cad_loc___1 == 1 | inv_cad_loc___2 == 1 | inv_cad_loc___3 == 1 | inv_cad_loc___4 == 1 | inv_cad_loc___5 == 1, 1, 0),
    cad_relevant_rca = ifelse(inv_cad_perc_rcaprox > 1 | inv_cad_perc_rcamid > 1 | inv_cad_perc_rcadist > 1, 1, 0),
    cad_relevant_cx = ifelse(inv_cad_perc_cxprox > 1 | inv_cad_perc_cxmid > 1 | inv_cad_perc_cxdist > 1, 1, 0),
    cad_relevant_lad = ifelse(inv_cad_perc_lm > 1 | inv_cad_perc_ladprox > 1 | inv_cad_perc_ladmid > 1 | inv_cad_perc_laddist > 1 | inv_cad_perc_diag1 > 1 | inv_cad_perc_diag2 > 1, 1, 0),
    cad_any_anomalous = ifelse((right_0_left_1 == 0 & cad_any_rca == 1) | (right_0_left_1 == 1 & (cad_any_lcx == 1 | cad_any_lad == 1)), 1, 0),
    cad_relevant_anomalous = ifelse((right_0_left_1 == 0 & cad_relevant_rca == 1) | (right_0_left_1 == 1 & (cad_relevant_cx == 1 | cad_relevant_lad == 1)), 1, 0),
    ccta_ostial_pn = (1 - (ccta_ostial_a / ccta_dist_a)) * 100,
    ccta_pn_dist = (1 - (ccta_mla_a / ccta_dist_a)) * 100,
  )

additional_ccta <- read_excel("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/additional_ccta_data.xlsx")

# keep only columns record_id, length_sinus, bogen_rca, bogen_lca, scalor_width, scalor_height
additional_ccta <- additional_ccta %>% 
    select(record_id, length_sinus, bogen_rca, bogen_lca, scalor_width, scalor_height)

# merge dataframes by record_id
baseline <- baseline %>% 
    left_join(additional_ccta, by = "record_id")

# quick data prepping
baseline <- baseline %>% 
    mutate(
        # ccta_ostial_pn = (1 - (ccta_ostial_a / ccta_dist_a)) * 100,
        # ccta_pn_dist = (1 - (ccta_mla_a / ccta_dist_a)) * 100,
        ccta_quali = as.numeric(ccta_quali),
        ccta_quali = ifelse(ccta_photon == "yes", ccta_quali + 1, ifelse(is.na(ccta_quali), 2, ccta_quali)),
        ccta_ostial_a_bsa = ccta_ostial_a / bsa,
        ccta_mla_a_bsa = ccta_mla_a / bsa,
        ccta_ostial_h_bsa = ccta_ostial_h / bsa,
        ccta_mla_h_bsa = ccta_mla_h / bsa,
        ccta_ostial_w_bsa = ccta_ostial_w / bsa,
        ccta_mla_w_bsa = ccta_mla_w / bsa,
        ccta_stj_rca_bsa = ccta_stj_rca / bsa,
        bogen_rca_bsa = bogen_rca / bsa,
        ccta_stj_rca_scaled = ccta_stj_rca * scalor_height,
        bogen_rca_bsa_scaled = bogen_rca * scalor_height,
        ccta_pn_dist = ifelse(ccta_pn_dist < 0, 0, ccta_pn_dist),
        ccta_ostial_pn = ifelse(ccta_ostial_pn < 0, 0, ccta_ostial_pn)
    )

invasive <- invasive %>% rename(
    patient_id = `NARCO ID`,
    inv_ivusrest_mla = MLA_rest,
    inv_ivusrest_mla_c = MLA_circ_rest,
    inv_ivusrest_mla_h = MLA_h_rest,
    inv_ivusrest_mla_w = MLA_w_rest,
    inv_ivusrest_mla_ellip = MLA_elllip_rest,
    inv_ivusrest_ostial_a = ostial_a_rest,
    inv_ivusrest_ostial_c = ostial_circ_rest,
    inv_ivusrest_ostial_h = ostial_h_rest,
    inv_ivusrest_ostial_w = ostial_w_rest,
    inv_ivusrest_ostial_ellip = ostial_ellip_rest,
    inv_ivusrest_dist_a = reference_a_rest,
    inv_ivusrest_dist_c = reference_circ_rest,
    inv_ivusdobu_mla = MLA_stress,
    inv_ivusdobu_mla_c = MLA_circ_stress,
    inv_ivusdobu_mla_h = MLA_h_stress,
    inv_ivusdobu_mla_w = MLA_w_stress,
    inv_ivusdobu_mla_ellip = MLA_elllip_stress,
    inv_ivusdobu_ostial_a = ostial_a_stress,
    inv_ivusdobu_ostial_c = ostial_circ_stress,
    inv_ivusdobu_ostial_h = ostial_h_stress,
    inv_ivusdobu_ostial_w = ostial_w_stress,
    inv_ivusdobu_ostial_ellip = ostial_ellip_stress,
    inv_ivusdobu_dist_a = reference_a_stress,
    inv_ivusdobu_dist_c = reference_circ_stress,
    inv_ivusado_mla = MLA_adenosine,
    inv_ivusado_mla_c = MLA_circ_adenosine,
    inv_ivusado_mla_h = MLA_h_adenosine,
    inv_ivusado_mla_w = MLA_w_adenosine,
    inv_ivusado_mla_ellip = MLA_elllip_adenosine,
    inv_ivusado_ostial_a = ostial_a_adenosine,
    inv_ivusado_ostial_c = ostial_circ_adenosine,
    inv_ivusado_ostial_h = ostial_h_adenosine,
    inv_ivusado_ostial_w = ostial_w_adenosine,
    inv_ivusado_ostial_ellip = ostial_ellip_adenosine,
    inv_ivusado_dist_a = reference_a_adenosine,
    inv_ivusado_dist_c = reference_circ_adenosine, 
) %>% mutate(inv_ivusdobu_mla = as.numeric(inv_ivusdobu_mla))

invasive <- invasive %>% mutate(
    patient_id = paste0("NARCO_", patient_id)
)

# drop last 3 columns
invasive <- invasive %>% select(-c(38:40))

# merge with baseline by record_id
baseline <- left_join(as_tibble(baseline), as_tibble(invasive), by = "patient_id")
baseline <- baseline %>% mutate(
    inv_ivusrest_mla = ifelse(is.na(inv_ivusrest_mla.y), inv_ivusrest_mla.x, inv_ivusrest_mla.y),
    inv_ivusrest_mla_ellip = ifelse(is.na(inv_ivusrest_mla_ellip.y), inv_ivusrest_mla_ellip.x, inv_ivusrest_mla_ellip.y),
    inv_ivusrest_dist_a = ifelse(is.na(inv_ivusrest_dist_a.y), inv_ivusrest_dist_a.x, inv_ivusrest_dist_a.y),
    inv_ivusdobu_mla = ifelse(is.na(inv_ivusdobu_mla.y), inv_ivusdobu_mla.x, inv_ivusdobu_mla.y),
    inv_ivusdobu_mla_ellip = ifelse(is.na(inv_ivusdobu_mla_ellip.y), inv_ivusdobu_mla_ellip.x, inv_ivusdobu_mla_ellip.y),
    inv_ivusdobu_dist_a = ifelse(is.na(inv_ivusdobu_dist_a.y), inv_ivusdobu_dist_a.x, inv_ivusdobu_dist_a.y),
    inv_ivusado_mla = ifelse(is.na(inv_ivusado_mla.y), inv_ivusado_mla.x, inv_ivusado_mla.y),
    inv_ivusado_mla_ellip = ifelse(is.na(inv_ivusado_mla_ellip.y), inv_ivusado_mla_ellip.x, inv_ivusado_mla_ellip.y),
    inv_ivusado_dist_a = ifelse(is.na(inv_ivusado_dist_a.y), inv_ivusado_dist_a.x, inv_ivusado_dist_a.y),
) %>% 
select(-c("inv_ivusrest_mla.y", "inv_ivusrest_mla.x", "inv_ivusrest_dist_a.y", "inv_ivusrest_dist_a.x", 
          "inv_ivusdobu_mla.y", "inv_ivusdobu_mla.x", "inv_ivusdobu_mla_ellip.y", "inv_ivusdobu_mla_ellip.x", 
          "inv_ivusdobu_dist_a.y", "inv_ivusdobu_dist_a.x", "inv_ivusado_mla.y", "inv_ivusado_mla.x",))

# quick data prepping
baseline <- baseline %>% 
    mutate(
        inv_ivus_dist = pmax(inv_ivusrest_dist_a, inv_ivusdobu_dist_a),
        inv_ivus_dist_c = pmax(inv_ivusrest_dist_c, inv_ivusdobu_dist_c),
        inv_ivusrest_pn = inv_ivusrest_mla / inv_ivus_dist,
        inv_ivusdobu_pn = inv_ivusdobu_mla / inv_ivus_dist,
        inv_ivusrest_mla_bsa = inv_ivusrest_mla / bsa,
        inv_ivusdobu_mla_bsa = inv_ivusdobu_mla / bsa,
        inv_ivusrest_mla_h_bsa = inv_ivusrest_mla_h / bsa,
        inv_ivusdobu_mla_h_bsa = inv_ivusdobu_mla_h / bsa,
        inv_ivusrest_mla_w_bsa = inv_ivusrest_mla_w / bsa,
        inv_ivusdobu_mla_w_bsa = inv_ivusdobu_mla_w / bsa,
        inv_ivusrest_hypo = inv_ivusrest_mla_c / inv_ivus_dist_c,
        inv_ivusdobu_hypo = inv_ivusdobu_mla_c / inv_ivus_dist_c,
        inv_ivusrest_ostial_hypo = inv_ivusrest_ostial_c / inv_ivus_dist_c,
        inv_ivusdobu_ostial_hypo = inv_ivusdobu_ostial_c / inv_ivus_dist_c,
        inv_ivus_delta = inv_ivusdobu_mla - inv_ivusrest_mla,
        inv_ivusrest_ostial_a_bsa = inv_ivusrest_ostial_a / bsa,
        inv_ivusdobu_ostial_a_bsa = inv_ivusdobu_ostial_a / bsa,
        inv_ivusrest_ostial_h_bsa = inv_ivusrest_ostial_h / bsa,
        inv_ivusdobu_ostial_h_bsa = inv_ivusdobu_ostial_h / bsa,
        inv_ivusrest_ostial_w_bsa = inv_ivusrest_ostial_w / bsa,
        inv_ivusdobu_ostial_w_bsa = inv_ivusdobu_ostial_w / bsa,
        inv_ivusrest_ostial_pn = 100 * (1 - inv_ivusrest_ostial_a / inv_ivus_dist),
        inv_ivusado_ostial_pn = 100 * (1 - inv_ivusado_ostial_a / inv_ivus_dist),
        inv_ivusdobu_ostial_pn = 100 * (1 - inv_ivusdobu_ostial_a / inv_ivus_dist),
        inv_ivusrest_mla_pn = 100 * (1 - inv_ivusrest_mla / inv_ivus_dist),
        inv_ivusado_mla_pn = 100 * (1 - inv_ivusado_mla / inv_ivus_dist),
        inv_ivusdobu_mla_pn = 100 * (1 - inv_ivusdobu_mla / inv_ivus_dist),
        inv_ivusrest_ostial_pn = ifelse(inv_ivusrest_ostial_pn < 0, 0, inv_ivusrest_ostial_pn),
        inv_ivusado_ostial_pn = ifelse(inv_ivusado_ostial_pn < 0, 0, inv_ivusado_ostial_pn),
        inv_ivusdobu_ostial_pn = ifelse(inv_ivusdobu_ostial_pn < 0, 0, inv_ivusdobu_ostial_pn),
        inv_ivusrest_mla_pn = ifelse(inv_ivusrest_mla_pn < 0, 0, inv_ivusrest_mla_pn),
        inv_ivusado_mla_pn = ifelse(inv_ivusado_mla_pn < 0, 0, inv_ivusado_mla_pn),
        inv_ivusdobu_mla_pn = ifelse(inv_ivusdobu_mla_pn < 0, 0, inv_ivusdobu_mla_pn),
        delta_ivus_ostial_a = inv_ivusdobu_ostial_a - inv_ivusrest_ostial_a,
        delta_ivus_mla = inv_ivusdobu_mla - inv_ivusrest_mla,
        delta_ffr = inv_ffrdobu - inv_ffrado,
        percent_ffr = (1 - inv_ffrdobu/inv_ffrado) * 100,
        percent_pressure = (1 - inv_ffrdobu / inv_pdpa) * 100,
        percent_compression = (1 - inv_ivusdobu_ostial_a / inv_ivusrest_ostial_a) * 100,
        delta_percent_stenosis = 1 - inv_ivusrest_ostial_pn / inv_ivusdobu_ostial_pn,
        inv_aomean_change = inv_dobu_aomean - inv_rest_aomean,
        inv_aosys_change = inv_dobu_aosys - inv_rest_aosys,
        inv_aodia_change = inv_dobu_aodia - inv_rest_aodia,
        inv_hr_change = inv_dobu_hr - inv_rest_hr,
        inv_ffr_change = inv_ffrado - inv_ffrdobu,
        inv_ivus_change = inv_ivusrest_ostial_a - inv_ivusdobu_ostial_a,
        inv_ivusrest_ostial_50 = ifelse(inv_ivusrest_ostial_pn >= 50, 1, 0),
        inv_ivusrest_imla_50 = ifelse(inv_ivusrest_mla_pn >= 50, 1, 0),
        inv_ivusrest_ostial_70 = ifelse(inv_ivusrest_ostial_pn >= 70, 1, 0),
        inv_ivusrest_imla_70 = ifelse(inv_ivusrest_mla_pn >= 70, 1, 0),
        inv_ivusdobu_ostial_50 = ifelse(inv_ivusdobu_ostial_pn >= 50, 1, 0),
        inv_ivusdobu_imla_50 = ifelse(inv_ivusdobu_mla_pn >= 50, 1, 0),
        inv_ivusdobu_ostial_70 = ifelse(inv_ivusdobu_ostial_pn >= 70, 1, 0),
        inv_ivusdobu_imla_70 = ifelse(inv_ivusdobu_mla_pn >= 70, 1, 0),
    )

# rename every _mla with _imla
baseline <- baseline %>% rename_with(~str_replace(., "_mla", "_imla"), contains("mla"))

baseline <- baseline %>% mutate(
  inv_ivusrest_mla = pmin(inv_ivusrest_imla, inv_ivusrest_ostial_a),
  inv_ivusado_mla = pmin(inv_ivusado_imla, inv_ivusado_ostial_a),
  inv_ivusdobu_mla = pmin(inv_ivusdobu_imla, inv_ivusdobu_ostial_a),
  inv_ivusrest_mla_bsa = inv_ivusrest_mla / bsa,
  inv_ivusado_mla_bsa = inv_ivusado_mla / bsa,
  inv_ivusdobu_mla_bsa = inv_ivusdobu_mla / bsa,
  inv_ivusrest_mla_ellip = pmax(inv_ivusrest_imla_ellip, inv_ivusrest_ostial_ellip),
  inv_ivusdobu_mla_ellip = pmax(inv_ivusdobu_imla_ellip, inv_ivusdobu_ostial_ellip),
  inv_ivusrest_mla_ln = pmax(inv_ivusrest_imla_pn, inv_ivusrest_ostial_pn),
  inv_ivusado_mla_ln = pmax(inv_ivusado_imla_pn, inv_ivusado_ostial_pn),
  inv_ivusdobu_mla_ln_any = pmax(inv_ivusdobu_imla_pn, inv_ivusdobu_ostial_pn),
  inv_ivusdobu_mla_ln_loc = ifelse(inv_ivusrest_imla > inv_ivusrest_ostial_a, inv_ivusdobu_imla_pn, inv_ivusdobu_ostial_pn),
  inv_ivusrest_mla_50 = ifelse(inv_ivusrest_mla_ln >= 50, 1, 0),
  inv_ivusrest_mla_70 = ifelse(inv_ivusrest_mla_ln >= 70, 1, 0),
  inv_ivusdobu_mla_50 = ifelse(inv_ivusdobu_mla_ln_any >= 50, 1, 0),
  inv_ivusdobu_mla_70 = ifelse(inv_ivusdobu_mla_ln_any >= 70, 1, 0),
)

## FACTORIZE ##########################################################################################################
baseline$sex_0_male_1_female <- factor(baseline$sex_0_male_1_female, levels = c(0, 1), labels = c("male", "female"))
baseline$caa_modality <- factor(baseline$caa_modality, levels = c(0, 1, 2, 3, 4, 5), labels = c("ccta", "angiography", "echocardiography", "cmr", "surgery", "other"))
baseline$hx_sym_ap_ccs <- factor(baseline$hx_sym_ap_ccs, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$hx_sym_dysp_nyha <- factor(baseline$hx_sym_dysp_nyha, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$sports_level <- factor(baseline$sports_level, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_level_ch <- factor(baseline$sports_level_ch, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_static_component <- factor(baseline$sports_static_component, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_static_component_ch <- factor(baseline$sports_static_component_ch, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_dynamic_component <- factor(baseline$sports_dynamic_component, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_dynamic_comp_ch <- factor(baseline$sports_dynamic_comp_ch, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_classification <- factor(baseline$sports_classification, levels = 1:5, ordered = TRUE, labels = c("very low", "low", "moderate", "high", "very high"))
baseline$caa_malignancy <- factor(baseline$caa_malignancy, levels = c(0, 1), labels = c("benign", "malign"))
baseline$caa_ostia <- factor(baseline$caa_ostia, levels = c(0, 1), labels = c("one", "two"))
baseline$dgn_echo_diafun <- factor(baseline$dgn_echo_diafun, levels = c(1, 2, 3, 4, 5, 6), labels = c("normal", "grade 1", "grade 2", "grade 3", "undefined arrythmia", "undefined"))
baseline$dgn_ecg_classification <- factor(baseline$dgn_ecg_classification, levels = c(0, 1), labels = c("normal", "abnormal"))
baseline$ccta_quali <- factor(baseline$ccta_quali, levels = c(0, 1, 2, 3), ordered = TRUE, labels = c("1", "2", "3", "4"))
baseline$ccta_dominance <- factor(baseline$ccta_dominance, levels = c(0, 1, 2), labels = c("right", "left", "balanced"))
baseline$funct_ergo_findings <- factor(baseline$funct_ergo_findings, levels = c(1, 2, 3, 4), labels = c("clinical and electrical negative", "clinical negative, electrical positive", "clinical positive, electrical negative", "clinical and electrical positive"))
baseline$inv_dominance <- factor(baseline$inv_dominance, levels = c(0, 1, 2), labels = c("right", "left", "balanced"))
baseline$synp_symp <- factor(baseline$synp_symp, levels = c(0, 1, 2, 3), labels = c("incidental finding", "symptoms linked to caa", "symptoms linked to cad", "symptoms linked to other"))
baseline$hx_dm <- factor(baseline$hx_dm, levels = c(2, 1, 3, 0), labels = c("DM Typ1", "DM Typ2", "Prediabetes", "No DM"))
baseline$hx_hxtobacco <- factor(baseline$hx_hxtobacco, levels = c(1, 2), labels = c("Hx Tobacco", "No Hx Tobacco"))
binary_vars <- names(baseline)[
  sapply(names(baseline), function(var_name) {
    all(baseline[[var_name]] %in% c(0, 1, NA))
  })
]
for (var in binary_vars) {
  baseline[[var]] <- factor(baseline[[var]], levels = c(0, 1), labels = c("no", "yes"))
}

baseline <- baseline %>% mutate(
  ccta_mla_a = pmin(ccta_imla_a, ccta_ostial_a),
  ccta_mla_a_bsa = pmin(ccta_imla_a_bsa, ccta_ostial_a_bsa),
  ccta_mla_dist = ccta_imla_dist,
  ccta_mla_c = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_c, ccta_ostial_c),
  ccta_mla_h = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_h, ccta_ostial_h),
  ccta_mla_h_bsa = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_h_bsa, ccta_ostial_h_bsa),
  ccta_mla_w = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_w, ccta_ostial_w),
  ccta_mla_w_bsa = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_w_bsa, ccta_ostial_w_bsa),
  ccta_mla_elliptic = ifelse(ccta_imla_a < ccta_ostial_a, ccta_imla_elliptic, ccta_ostial_elliptic),
  ccta_mla_ln = pmax(ccta_pn_dist, ccta_ostial_pn),
)

pressure <- pressure %>% mutate(
  patient_id = toupper(patient_id),
)

baseline <- left_join(as_tibble(baseline), as_tibble(pressure), by = "patient_id")


# applying the shift from invasive data to the pressure curves
baseline <- baseline %>% mutate(
  rest_diff = ifelse(is.na(inv_pdpa), inv_rfr - pdpa_mean_rest, inv_pdpa - pdpa_mean_rest),
  ado_diff = inv_ffrado - pdpa_mean_ado,
  dobu_diff = inv_ffrdobu - pdpa_mean_dobu
)

baseline <- baseline %>% mutate(
  iFR_mean_rest = iFR_mean_rest + rest_diff,
  mid_systolic_ratio_mean_rest = mid_systolic_ratio_mean_rest + rest_diff,
  pdpa_mean_rest = pdpa_mean_rest + rest_diff,
  iFR_mean_ado = iFR_mean_ado + ado_diff,
  mid_systolic_ratio_mean_ado = mid_systolic_ratio_mean_ado + ado_diff,
  pdpa_mean_ado = pdpa_mean_ado + ado_diff,
  iFR_mean_dobu = iFR_mean_dobu + dobu_diff
)

# if pdpa_mean_rest is NA, then take inv_pdpa or inv_rfr
baseline <- baseline %>% mutate(
  pdpa_mean_rest = ifelse(is.na(pdpa_mean_rest), ifelse(is.na(inv_pdpa), inv_rfr, inv_pdpa), pdpa_mean_rest)
)

saveRDS(baseline, file = paste0(yaml$demographics$output_dir_data, "/baseline.rds"))
# saveRDS(baseline, file= "C:/Users/ansel/OneDrive/Dokumente/3_Research/7_NARCOfunctional/baseline.rds")
