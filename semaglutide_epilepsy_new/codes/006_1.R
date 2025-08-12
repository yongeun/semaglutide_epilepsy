# 006_1
# 필요한 라이브러리 설치 및 로드
# 필수 패키지 설치 (필요한 경우)
required_packages <- c(
  "tidyverse", "bigrquery", "dplyr", "glue", "MatchIt", "survival", "survminer",
  "lubridate", "ggplot2", "lme4", "lmerTest", "broom.mixed", 
  "emmeans", "survey", "tmle", "SuperLearner", "gridExtra", "writexl", "DHARMa"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}

# 라이브러리 로드
library(tidyverse)
library(bigrquery)
library(dplyr)
library(tidyr)
library(glue)
library(MatchIt)
library(survival)
library(survminer)
library(lubridate)
library(ggplot2)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)
library(survey)
library(tmle)
library(gridExtra)
library(writexl)
library(DHARMa)

# Install and load additional packages for flowchart
flowchart_packages <- c("DiagrammeR", "rsvg")
for (pkg in flowchart_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.rstudio.com/")
  }
}
library(DiagrammeR)
library(rsvg)

# =============================================================================
# MODIFICATION NOTE: This script has been modified to include calendar year 
# of the index date (index_year) as a covariate in both IPWT and TMLE analyses.
# This helps control for temporal trends and changes in clinical practice over time.
# The index_year variable is created as a factor in both analytic dataframes.
# =============================================================================

# 데이터 로드
name_of_file_in_bucket <- 'merged_df_2.csv'
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
system(paste0("gsutil cp ", my_bucket, "/data/", name_of_file_in_bucket, " ."), intern=TRUE)
merged_df_2 <- read_csv(name_of_file_in_bucket)

# BMI 데이터 로드
bmi_file_name <- 'bmi_panel.csv'
system(paste0("gsutil cp ", my_bucket, "/data/", bmi_file_name, " ."), intern=TRUE)
bmi_panel_df_loaded <- read_csv(bmi_file_name) %>%
  mutate(weight_date = as.Date(weight_date)) %>% # 날짜 형식 확인
  dplyr::select(person_id, weight_date, bmi) %>% # 필요한 열만 선택
  distinct(person_id, weight_date, .keep_all = TRUE) # 개인-날짜별 중복 항목 제거

# 범주형 변수 변환
cols_to_factor <- c(
 "epilepsy_or_seizure", "adrd", "mci","stroke",
  "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "DPP4", "Diuretic",
  "Ezetimibe", "GLP1", "MRA", "OtherHTN", "RAAS", "SGLT2", "SU", "Statin",
  "TZD", "Insulin",
  "mi", "chf", "pvd", "cvd", "dem", "cpd", "ctd", "pud", "mld",
  "dm", "dmc", "ph", "rd", "cancer", "msld", "mc", "hiv"
)
merged_df_2[cols_to_factor] <- lapply(
  merged_df_2[cols_to_factor],
  function(x) factor(x, levels = c(0,1), labels = c("0","1"))
)

# 환경 변수 설정
DATASET <- Sys.getenv("WORKSPACE_CDR")
GOOGLE_PROJECT <- Sys.getenv("GOOGLE_PROJECT")
cutoff_date <- as.Date("2018-01-01")  # Semaglutide FDA approval date
# Start exposure window from cutoff date to avoid unnecessary queries
# This ensures we only capture exposures after semaglutide became available
exposure_window <- list(start = as.character(cutoff_date), end = "2021-12-31")
data_cut_date <- as.Date("2023-10-01")

# 약물 클래스 및 비교 정의
drug_classes <- list(
  GLP1 = c("semaglutide"),
  SGLT2 = c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin"),
  SU_DPP4 = c("glimepiride", "glipizide", "glyburide", "alogliptin", "linagliptin", "sitagliptin", "saxagliptin")
)

comparisons <- list(
  list(name = "GLP1 vs SGLT2", exposure = "GLP1", comparator = "SGLT2"),
  list(name = "GLP1 vs SU_DPP4", exposure = "GLP1", comparator = "SU_DPP4"),
  list(name = "SGLT2 vs SU_DPP4", exposure = "SGLT2", comparator = "SU_DPP4")
)

outcomes <- list(
  list(var = "epilepsy_or_seizure_start_date", label = "Epilepsy/Seizure", late_onset = FALSE, early_onset = FALSE),
  list(var = "epilepsy_or_seizure_start_date", label = "Early-onset Epilepsy/Seizure", late_onset = FALSE, early_onset = TRUE),
  list(var = "epilepsy_or_seizure_start_date", label = "Late-onset Epilepsy/Seizure", late_onset = TRUE, early_onset = FALSE),
  list(var = "adrd_start_date", label = "ADRD", late_onset = FALSE, early_onset = FALSE),
  list(var = "stroke_start_date", label = "stroke", late_onset = FALSE, early_onset = FALSE)
)

# 공변량 변수 목록
all_ps_vars <- c(
  "age", "sex_cat", "raceethnicity_cat",
  "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
  "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
  "TZD", "Insulin", "SU", "DPP4", "SGLT2", "GLP1",
  "mi", "chf", "pvd", "cvd", "dem", "cpd", "ctd", "pud", "mld",
  "ph", "rd", "cancer", "msld", "mc", "hiv"
)

# 결과 저장 리스트 초기화
ipwt_results <- list()
tmle_results <- list()

# Initialize patient flow tracking
patient_flow_tracker <- list()

# Function to create study flowchart
create_study_flowchart <- function(flow_data, comparison_name, outcome_name) {
  
  # Create flowchart using DiagrammeR
  flowchart <- grViz("
    digraph flowchart {
      
      # Graph attributes
      graph [layout = dot, rankdir = TB]
      
      # Node attributes
      node [shape = box, style = filled, fillcolor = lightblue, 
            fontname = Arial, fontsize = 10]
      
      # Edge attributes  
      edge [color = black, arrowhead = vee]
      
      # Define nodes
      A [label = '@@1', fillcolor = lightgray]
      B [label = '@@2', fillcolor = lightcoral]
      C [label = '@@3', fillcolor = lightblue]
      D [label = '@@4', fillcolor = lightcoral]
      E [label = '@@5', fillcolor = lightblue]
      F [label = '@@6', fillcolor = lightcoral]
      G [label = '@@7', fillcolor = lightblue]
      H [label = '@@8', fillcolor = lightcoral]
      I [label = '@@9', fillcolor = lightgreen]
      
      # Define edges
      A -> C
      A -> B [style = dashed]
      C -> E  
      C -> D [style = dashed]
      E -> G
      E -> F [style = dashed]
      G -> I
      G -> H [style = dashed]
    }
    
    [1]: paste0('Initial Cohort with T2DM\\n', 'N = ', flow_data$initial_dm)
    [2]: paste0('Excluded: Non-T2DM patients\\n', 'N = ', flow_data$excluded_non_dm)
    [3]: paste0('After T2DM Filter\\n', 'N = ', flow_data$after_dm_filter)
    [4]: paste0('Excluded: Pre-existing ', outcome_name, '\\n', 'N = ', flow_data$excluded_preexisting)
    [5]: paste0('After Excluding Pre-existing\\n', 'N = ', flow_data$after_preexisting_exclusion)
    [6]: paste0('Excluded: Missing follow-up data\\n', 'N = ', flow_data$excluded_followup)
    [7]: paste0('After Follow-up Requirements\\n', 'N = ', flow_data$after_followup_filter)
    [8]: paste0('Excluded: Age < 50 (Late-onset only)\\n', 'N = ', flow_data$excluded_age)
    [9]: paste0('Final Analysis Cohort\\n', comparison_name, ' - ', outcome_name, '\\n', 'N = ', flow_data$final_cohort, '\\n', 'Events = ', flow_data$total_events)
  ")
  
  return(flowchart)
}

# Function to save flowchart
save_flowchart <- function(flowchart, comparison_name, outcome_name) {
  filename <- paste0("flowchart_", 
                     gsub(" ", "_", comparison_name), "_", 
                     gsub("[/ ]", "_", outcome_name), ".png")
  
  # Export to PNG
  flowchart %>%
    export_svg() %>%
    charToRaw() %>%
    rsvg::rsvg_png(filename, width = 1200, height = 800)
  
  cat("Saved flowchart:", filename, "\n")
}

# ---- BigQuery 관련 함수 ----
download_data <- function(sql) {
  bq_table_download(bq_project_query(GOOGLE_PROJECT, sql))
}

get_drug_concepts <- function(names_vec) {
  pattern <- paste0(
    'LOWER(c.concept_name) LIKE "%', tolower(names_vec), '%"',
    collapse = " OR "
  )
  sql <- glue("
    SELECT DISTINCT c2.concept_id
    FROM `{DATASET}.concept` c
    JOIN `{DATASET}.concept_ancestor` ca ON c.concept_id = ca.ancestor_concept_id
    JOIN `{DATASET}.concept` c2 ON c2.concept_id = ca.descendant_concept_id
    WHERE c.concept_class_id = 'Ingredient' AND ({pattern})
  ")
  download_data(sql)
}

build_exposure_idx <- function(concepts_df, idx_name, use_exposure_window = TRUE, start_date = NULL, end_date = NULL) {
  id_string <- paste(concepts_df$concept_id, collapse = ", ")
  
  # Build the date constraint based on parameters
  date_constraint <- if(use_exposure_window) {
    glue("AND drug_exposure_start_date BETWEEN '{exposure_window$start}' AND '{exposure_window$end}'")
  } else if(!is.null(start_date) && !is.null(end_date)) {
    glue("AND drug_exposure_start_date BETWEEN '{start_date}' AND '{end_date}'")
  } else if(!is.null(start_date)) {
    glue("AND drug_exposure_start_date >= '{start_date}'")
  } else if(!is.null(end_date)) {
    glue("AND drug_exposure_start_date <= '{end_date}'")
  } else {
    "" # No date constraint
  }
  
  sql <- glue("
    WITH first_use AS (
      SELECT person_id, MIN(drug_exposure_start_date) AS index_date
      FROM `{DATASET}.drug_exposure`
      WHERE drug_concept_id IN ({id_string})
        {date_constraint}
      GROUP BY person_id
    ),
    prior_use AS (
      SELECT f.person_id
      FROM first_use f
      JOIN `{DATASET}.drug_exposure` d
        ON d.person_id = f.person_id
       AND d.drug_concept_id IN ({id_string})
       AND d.drug_exposure_start_date BETWEEN DATE_SUB(f.index_date, INTERVAL 365 DAY)
                                           AND DATE_SUB(f.index_date, INTERVAL 1 DAY)
      GROUP BY f.person_id
    )
    SELECT f.person_id, f.index_date
    FROM first_use f
    LEFT JOIN prior_use p ON f.person_id = p.person_id
    WHERE p.person_id IS NULL
  ")
  download_data(sql) %>% rename(!!idx_name := index_date)
}

followup_and_event <- function(df, outcome_var, data_cut_date, late_onset = FALSE, early_onset = FALSE) {
  cat("\n--- DEBUG followup_and_event function ---\n")
  cat("Input parameters:\n")
  cat("- outcome_var:", outcome_var, "\n")
  cat("- data_cut_date:", as.character(data_cut_date), "\n")
  cat("- late_onset:", late_onset, "\n")
  cat("- early_onset:", early_onset, "\n")
  cat("- Input data rows:", nrow(df), "\n")
  
  # DEBUG: Check input data
  outcome_dates_input <- df[[outcome_var]]
  cat("- Patients with outcome dates:", sum(!is.na(outcome_dates_input)), "\n")
  if(sum(!is.na(outcome_dates_input)) > 0) {
    cat("- Outcome date range:", min(outcome_dates_input, na.rm = TRUE), "to", max(outcome_dates_input, na.rm = TRUE), "\n")
  }
  
  result <- df %>%
    mutate(
      crossover_date = case_when(
        treatment == 1 ~ .data[[attr(df, "comparator_idx")]],
        treatment == 0 ~ .data[[attr(df, "exposure_idx")]],
        TRUE ~ NA_Date_
      ),
      outcome_date = as.Date(.data[[outcome_var]]),
      ehr_end = as.Date(EHRmaxDT),
      censor_date = pmin(ehr_end, data_cut_date, crossover_date - 1, na.rm = TRUE),
      event_time = as.numeric(pmin(outcome_date, censor_date, na.rm = TRUE) - index_date),
      raw_event = as.integer(!is.na(outcome_date) & outcome_date <= censor_date),
      age_at_event = age + event_time / 365.25,
      event = if (late_onset) {
        as.integer(raw_event == 1 & age_at_event >= 50)
      } else if (early_onset) {
        as.integer(raw_event == 1 & age_at_event < 50)
      } else {
        raw_event
      }
    ) %>%
    filter(event_time >= 0)
  
  # DEBUG: Check intermediate calculations
  cat("\nDEBUG intermediate calculations:\n")
  cat("- Crossover dates available:", sum(!is.na(result$crossover_date)), "\n")
  cat("- EHR end dates available:", sum(!is.na(result$ehr_end)), "\n")
  cat("- Raw events (before late-onset filter):", sum(result$raw_event == 1, na.rm = TRUE), "\n")
  cat("- Events after time >= 0 filter:", sum(result$event == 1, na.rm = TRUE), "\n")
  cat("- Patients filtered out due to event_time < 0:", nrow(df) - nrow(result), "\n")
  
  if(late_onset) {
    cat("- Late-onset applied: events with age >= 50:", sum(result$raw_event == 1 & result$age_at_event >= 50, na.rm = TRUE), "\n")
  }
  if(early_onset) {
    cat("- Early-onset applied: events with age < 50:", sum(result$raw_event == 1 & result$age_at_event < 50, na.rm = TRUE), "\n")
  }
  
  # DEBUG: Check event_time distribution
  cat("\nEvent time distribution:\n")
  cat("- Min event_time:", min(result$event_time, na.rm = TRUE), "\n")
  cat("- Max event_time:", max(result$event_time, na.rm = TRUE), "\n")
  cat("- Mean event_time:", mean(result$event_time, na.rm = TRUE), "\n")
  
  # DEBUG: Check censoring reasons
  censoring_summary <- result %>%
    mutate(
      censor_reason = case_when(
        !is.na(outcome_date) & outcome_date <= censor_date ~ "Event occurred",
        is.na(outcome_date) ~ "No event - administrative censoring",
        outcome_date > censor_date ~ "Event after censoring",
        TRUE ~ "Other"
      )
    )
  cat("\nCensoring reasons:\n")
  print(table(censoring_summary$censor_reason))
  
  cat("--- END DEBUG followup_and_event ---\n\n")
  
  return(result)
}

# ---- 1. IPWT 함수 ----
run_ipwt_and_cox <- function(df, exclude_vars = NULL, trim_threshold = 0.01) {
  # 공변량 목록 구성 - all_ps_vars 사용
  # Note: income, education, baseline_bmi, baseline_hba1c, index_year are not in all_ps_vars
  # so we add them separately
  additional_vars <- c("income", "education", "baseline_bmi", "baseline_hba1c", "index_year")
  all_covariates <- c(all_ps_vars, additional_vars)
  
  # Create formula using all covariates
  ps_formula <- as.formula(paste("treatment ~", paste(all_covariates, collapse = " + ")))
  
  # 제외할 변수 처리
  rhs_vars <- all_covariates
  if (!is.null(exclude_vars)) {
    rhs_vars <- setdiff(rhs_vars, exclude_vars)
  }

  # 가변적인 변수만 유지
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x)
    length(unique(na.omit(x))) > 1)]

  # 성향 점수 공식 구성
  ps_form <- reformulate(keep_vars, response = "treatment")
  
  # Handle missing values - create complete cases subset for modeling
  model_vars <- c("treatment", keep_vars)
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]
  
  if (nrow(df_complete) < nrow(df)) {
    cat("Note: Dropping", nrow(df) - nrow(df_complete), "rows with missing values in covariates\n")
  }
  
  # 성향 점수 모델 - fit on complete cases
  ps_model <- glm(ps_form, data = df_complete, family = binomial())
  
  # 성향 점수 계산 - only for complete cases
  df_complete$ps <- predict(ps_model, type = "response")
  
  # IPWT 가중치 계산
  df_complete$ipw <- ifelse(df_complete$treatment == 1, 
                  1/df_complete$ps, 
                  1/(1-df_complete$ps))
  
  # 극단 가중치 처리
  if (!is.null(trim_threshold)) {
    lower <- quantile(df_complete$ipw, trim_threshold)
    upper <- quantile(df_complete$ipw, 1 - trim_threshold)
    df_complete$ipw_trimmed <- pmin(pmax(df_complete$ipw, lower), upper)
  } else {
    df_complete$ipw_trimmed <- df_complete$ipw
  }
  
  # 표준화된 가중치
  df_complete$ipw_std <- df_complete$ipw_trimmed * nrow(df_complete) / sum(df_complete$ipw_trimmed)
  
  # 균형 평가
  balance_before <- calculate_balance_metrics(df_complete, keep_vars, NULL)
  balance_after <- calculate_balance_metrics(df_complete, keep_vars, df_complete$ipw_std)
  
  # 가중치 적용 Cox 모델
  weighted_cox <- coxph(Surv(event_time, event) ~ treatment, 
                       data = df_complete, weights = ipw_std)
  
  # 결과 반환 - return the complete cases dataframe
  list(
    balance_before = balance_before,
    balance_after = balance_after,
    cox = summary(weighted_cox),
    cohort = df_complete,  # Return complete cases only
    ps_model = ps_model,
    dropped_covariates = setdiff(rhs_vars, keep_vars),
    n_dropped = nrow(df) - nrow(df_complete)
  )
}

# (Moved to consolidated section at end of script)

# ---- Create Before vs After Weighting Table Function ----
create_weighted_characteristics_table <- function(result, comparison_name, outcome_name) {
  
  df <- result$cohort
  vars <- result$balance_before$variable
  
  # Remove level suffixes to get base variable names
  base_vars <- unique(gsub("1$", "", vars))
  base_vars <- base_vars[base_vars %in% names(df)]
  
  characteristics_list <- list()
  
  for (var in base_vars) {
    if (var %in% names(df)) {
      
      if (is.factor(df[[var]])) {
        # Categorical variables - calculate proportions
        
        # Before weighting
        before_summary <- df %>%
          group_by(treatment) %>%
          summarise(
            n = n(),
            prop = sum(as.numeric(.data[[var]]) == 2, na.rm = TRUE) / n() * 100,  # Level "1" 
            .groups = "drop"
          ) %>%
          mutate(
            treatment_group = ifelse(treatment == 1, "Semaglutide", "Comparator"),
            weighting = "Before"
          ) %>%
          select(treatment_group, weighting, n, prop)
        
        # After weighting  
        after_summary <- df %>%
          group_by(treatment) %>%
          summarise(
            weighted_n = sum(ipw_std, na.rm = TRUE),
            weighted_prop = sum((as.numeric(.data[[var]]) == 2) * ipw_std, na.rm = TRUE) / 
                           sum(ipw_std, na.rm = TRUE) * 100,
            .groups = "drop"
          ) %>%
          mutate(
            treatment_group = ifelse(treatment == 1, "Semaglutide", "Comparator"),
            weighting = "After IPWT"
          ) %>%
          rename(n = weighted_n, prop = weighted_prop) %>%
          select(treatment_group, weighting, n, prop)
        
        var_summary <- bind_rows(before_summary, after_summary) %>%
          mutate(
            variable = var,
            comparison = comparison_name,
            outcome = outcome_name,
            value_type = "Proportion (%)"
          )
        
      } else {
        # Continuous variables - calculate means
        
        # Before weighting
        before_summary <- df %>%
          group_by(treatment) %>%
          summarise(
            n = n(),
            mean_val = mean(.data[[var]], na.rm = TRUE),
            sd_val = sd(.data[[var]], na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(
            treatment_group = ifelse(treatment == 1, "Semaglutide", "Comparator"),
            weighting = "Before",
            formatted_val = sprintf("%.2f (%.2f)", mean_val, sd_val)
          ) %>%
          select(treatment_group, weighting, n, formatted_val)
        
        # After weighting
        after_summary <- df %>%
          group_by(treatment) %>%
          summarise(
            weighted_n = sum(ipw_std, na.rm = TRUE),
            weighted_mean = weighted.mean(.data[[var]], ipw_std, na.rm = TRUE),
            # Weighted SD calculation
            weighted_var = sum(ipw_std * (.data[[var]] - weighted_mean)^2, na.rm = TRUE) / sum(ipw_std, na.rm = TRUE),
            weighted_sd = sqrt(weighted_var),
            .groups = "drop"
          ) %>%
          mutate(
            treatment_group = ifelse(treatment == 1, "Semaglutide", "Comparator"),
            weighting = "After IPWT",
            formatted_val = sprintf("%.2f (%.2f)", weighted_mean, weighted_sd)
          ) %>%
          rename(n = weighted_n) %>%
          select(treatment_group, weighting, n, formatted_val)
        
        var_summary <- bind_rows(before_summary, after_summary) %>%
          mutate(
            variable = var,
            comparison = comparison_name,
            outcome = outcome_name,
            value_type = "Mean (SD)",
            prop = NA  # Add prop column for consistency
          ) %>%
          rename(prop = formatted_val)
      }
      
      characteristics_list[[var]] <- var_summary
    }
  }
  
  # Combine all variables
  final_table <- bind_rows(characteristics_list)
  
  return(final_table)
}

# 균형 메트릭 계산 함수
calculate_balance_metrics <- function(df, vars, weights) {
  balance_df <- data.frame(
    variable = character(),
    std_diff = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Define multilevel variables that should get overall SMD (not individual level SMDs)
  multilevel_overall_vars <- c("index_year", "income", "education", "raceethnicity_cat", "sex_cat", "bmi_category")
  
  for (var in vars) {
    if (is.factor(df[[var]]) && var %in% multilevel_overall_vars) {
      # Special handling for index_year - group 2018/2019 and 2020/2021
      if (var == "index_year") {
        # Create grouped index_year variable
        df_temp <- df %>%
          mutate(index_year_grouped = case_when(
            as.character(.data[[var]]) %in% c("2018", "2019") ~ "2018/2019",
            as.character(.data[[var]]) %in% c("2020", "2021") ~ "2020/2021",
            TRUE ~ as.character(.data[[var]])
          ))
        
        # Calculate SMD for each group
        groups <- c("2018/2019", "2020/2021")
        level_smds <- numeric()
        
        for (group in groups) {
          var_binary <- as.numeric(df_temp$index_year_grouped == group)
          
          if (is.null(weights)) {
            # 가중치 없음
            mean_t1 <- mean(var_binary[df_temp$treatment == 1], na.rm = TRUE)
            mean_t0 <- mean(var_binary[df_temp$treatment == 0], na.rm = TRUE)
            var_t1 <- var(var_binary[df_temp$treatment == 1], na.rm = TRUE)
            var_t0 <- var(var_binary[df_temp$treatment == 0], na.rm = TRUE)
          } else {
            # 가중치 있음
            w_t1 <- weights[df_temp$treatment == 1]
            w_t0 <- weights[df_temp$treatment == 0]
            
            mean_t1 <- weighted.mean(var_binary[df_temp$treatment == 1], w_t1, na.rm = TRUE)
            mean_t0 <- weighted.mean(var_binary[df_temp$treatment == 0], w_t0, na.rm = TRUE)
            
            # 가중 분산 계산
            var_t1 <- sum(w_t1 * (var_binary[df_temp$treatment == 1] - mean_t1)^2, na.rm = TRUE) / 
                       sum(w_t1)
            var_t0 <- sum(w_t0 * (var_binary[df_temp$treatment == 0] - mean_t0)^2, na.rm = TRUE) / 
                       sum(w_t0)
          }
          
          # 표준화된 차이 계산
          pooled_sd <- sqrt((var_t1 + var_t0) / 2)
          level_smd <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
          level_smds <- c(level_smds, level_smd)
        }
        
        # Use maximum absolute SMD as overall SMD for grouped index_year
        overall_smd <- ifelse(length(level_smds) > 0, max(abs(level_smds), na.rm = TRUE), 0)
        
        # Preserve the sign of the SMD with maximum absolute value
        if(length(level_smds) > 0) {
          max_idx <- which.max(abs(level_smds))
          overall_smd <- overall_smd * sign(level_smds[max_idx])
        }
        
      } else {
        # Regular multilevel categorical variables - calculate overall SMD using max absolute SMD approach
        all_levels <- levels(df[[var]])
        level_smds <- numeric()
        
        # Calculate SMD for each level
        for (lvl in all_levels) {
          var_binary <- as.numeric(df[[var]] == lvl)
          
          if (is.null(weights)) {
            # 가중치 없음
            mean_t1 <- mean(var_binary[df$treatment == 1], na.rm = TRUE)
            mean_t0 <- mean(var_binary[df$treatment == 0], na.rm = TRUE)
            var_t1 <- var(var_binary[df$treatment == 1], na.rm = TRUE)
            var_t0 <- var(var_binary[df$treatment == 0], na.rm = TRUE)
          } else {
            # 가중치 있음
            w_t1 <- weights[df$treatment == 1]
            w_t0 <- weights[df$treatment == 0]
            
            mean_t1 <- weighted.mean(var_binary[df$treatment == 1], w_t1, na.rm = TRUE)
            mean_t0 <- weighted.mean(var_binary[df$treatment == 0], w_t0, na.rm = TRUE)
            
            # 가중 분산 계산
            var_t1 <- sum(w_t1 * (var_binary[df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / 
                       sum(w_t1)
            var_t0 <- sum(w_t0 * (var_binary[df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / 
                       sum(w_t0)
          }
          
          # 표준화된 차이 계산
          pooled_sd <- sqrt((var_t1 + var_t0) / 2)
          level_smd <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
          level_smds <- c(level_smds, level_smd)
        }
        
        # Use maximum absolute SMD as overall SMD
        overall_smd <- ifelse(length(level_smds) > 0, max(abs(level_smds), na.rm = TRUE), 0)
        
        # Preserve the sign of the SMD with maximum absolute value
        if(length(level_smds) > 0) {
          max_idx <- which.max(abs(level_smds))
          overall_smd <- overall_smd * sign(level_smds[max_idx])
        }
      }
      
      balance_df <- rbind(balance_df, 
                        data.frame(
                          variable = var,
                          std_diff = overall_smd
                        ))
      
    } else if (is.factor(df[[var]])) {
      # Regular categorical variables - individual level SMDs
      for (lvl in levels(df[[var]])[-1]) {  # 첫번째 레벨은 기준
        var_binary <- as.numeric(df[[var]] == lvl)
        
        if (is.null(weights)) {
          # 가중치 없음
          mean_t1 <- mean(var_binary[df$treatment == 1], na.rm = TRUE)
          mean_t0 <- mean(var_binary[df$treatment == 0], na.rm = TRUE)
          var_t1 <- var(var_binary[df$treatment == 1], na.rm = TRUE)
          var_t0 <- var(var_binary[df$treatment == 0], na.rm = TRUE)
        } else {
          # 가중치 있음
          w_t1 <- weights[df$treatment == 1]
          w_t0 <- weights[df$treatment == 0]
          
          mean_t1 <- weighted.mean(var_binary[df$treatment == 1], w_t1, na.rm = TRUE)
          mean_t0 <- weighted.mean(var_binary[df$treatment == 0], w_t0, na.rm = TRUE)
          
          # 가중 분산 계산
          var_t1 <- sum(w_t1 * (var_binary[df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / 
                     sum(w_t1)
          var_t0 <- sum(w_t0 * (var_binary[df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / 
                     sum(w_t0)
        }
        
        # 표준화된 차이 계산
        pooled_sd <- sqrt((var_t1 + var_t0) / 2)
        # 분모가 0인 경우 처리
        std_diff <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
        
        balance_df <- rbind(balance_df, 
                          data.frame(
                            variable = paste0(var, lvl),
                            std_diff = std_diff
                          ))
      }
    } else {
      # 수치형 변수
      if (is.null(weights)) {
        mean_t1 <- mean(df[[var]][df$treatment == 1], na.rm = TRUE)
        mean_t0 <- mean(df[[var]][df$treatment == 0], na.rm = TRUE)
        var_t1 <- var(df[[var]][df$treatment == 1], na.rm = TRUE)
        var_t0 <- var(df[[var]][df$treatment == 0], na.rm = TRUE)
      } else {
        w_t1 <- weights[df$treatment == 1]
        w_t0 <- weights[df$treatment == 0]
        
        mean_t1 <- weighted.mean(df[[var]][df$treatment == 1], w_t1, na.rm = TRUE)
        mean_t0 <- weighted.mean(df[[var]][df$treatment == 0], w_t0, na.rm = TRUE)
        
        var_t1 <- sum(w_t1 * (df[[var]][df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / 
                   sum(w_t1)
        var_t0 <- sum(w_t0 * (df[[var]][df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / 
                   sum(w_t0)
      }
      
      pooled_sd <- sqrt((var_t1 + var_t0) / 2)
      std_diff <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
      
      balance_df <- rbind(balance_df, 
                        data.frame(
                          variable = var,
                          std_diff = std_diff
                        ))
    }
  }
  
  return(balance_df)
}

# ---- 2. TMLE 함수 ----
run_tmle_analysis <- function(df, outcome_var = "event", exclude_vars = NULL) {
  # 기본 공변량 - all_ps_vars 사용
  additional_vars <- c("income", "education", "baseline_bmi", "index_year")
  base_covs <- c(all_ps_vars, additional_vars)
  
  # 제외 변수 처리
  covs <- setdiff(base_covs, exclude_vars)
  
  # 가변적인 변수만 유지
  keep_covs <- covs[sapply(df[covs], function(x) length(unique(na.omit(x))) > 1)]
  
  # TMLE 데이터 준비
  model_data <- df %>%
    select(treatment, !!sym(outcome_var), all_of(keep_covs)) %>%
    na.omit()
  
  # 범주형 변수 처리
  factor_vars <- names(model_data)[sapply(model_data, is.factor)]
  model_data <- model_data %>%
    mutate_at(vars(all_of(factor_vars)), as.numeric) %>%
    mutate_at(vars(all_of(factor_vars)), function(x) x - 1)
  
  # SuperLearner 라이브러리 정의
  sl_lib <- c("SL.glm", "SL.mean")
  
  # TMLE 실행
  tmle_result <- tmle(
    Y = model_data[[outcome_var]],
    A = model_data$treatment,
    W = model_data %>% select(-treatment, -!!sym(outcome_var)),
    Q.SL.library = sl_lib,
    g.SL.library = sl_lib,
    family = "binomial"
  )
  
  # NNT 계산
  if (!is.null(tmle_result$estimates$ATE$psi) && tmle_result$estimates$ATE$psi != 0) {
    nnt <- 1 / abs(tmle_result$estimates$ATE$psi)
  } else {
    nnt <- Inf
  }
  
  # 결과 반환
  list(
    tmle_result = tmle_result,
    data = df,
    nnt = nnt,
    covariates = keep_covs
  )
}
# 결과 표시 및 요약 테이블 생성 함수
create_results_table <- function(type, result, outcome_label, comparison_name) {
  if (type == "psp") {
    # PSP 결과 추출 (수정된 부분)
    cox_model <- result$cox$conf.int
    hr <- cox_model[1, "exp(coef)"]
    lower_ci <- cox_model[1, "lower .95"]
    upper_ci <- cox_model[1, "upper .95"]
    p_value <- result$cox$coefficients[1, "Pr(>|z|)"]
    
    data.frame(
      method = "PSP 1:1 Matching",
      comparison = comparison_name,
      outcome = outcome_label,
      effect_measure = "Hazard Ratio",
      estimate = hr,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      p_value = p_value,
      n_total = nrow(result$cohort),
      n_events = sum(result$cohort$event),
      stringsAsFactors = FALSE
    )
    
  } else if (type == "ipwt") {
    # IPWT 결과 추출 (수정된 부분)
    cox_model <- result$cox$conf.int
    hr <- cox_model[1, "exp(coef)"]
    lower_ci <- cox_model[1, "lower .95"]
    upper_ci <- cox_model[1, "upper .95"]
    p_value <- result$cox$coefficients[1, "Pr(>|z|)"]
    
    data.frame(
      method = "IPWT",
      comparison = comparison_name,
      outcome = outcome_label,
      effect_measure = "Hazard Ratio",
      estimate = hr,
      lower_ci = lower_ci,
      upper_ci = upper_ci,
      p_value = p_value,
      n_total = nrow(result$cohort),
      n_events = sum(result$cohort$event),
      stringsAsFactors = FALSE
    )
    
  } else if (type == "tmle") {
    # TMLE 결과 추출
    tmle_res <- result$tmle_result
    
    # ATE (Risk Difference)
    rd_row <- data.frame(
      method = "TMLE",
      comparison = comparison_name,
      outcome = outcome_label,
      effect_measure = "Risk Difference",
      estimate = tmle_res$estimates$ATE$psi,
      lower_ci = tmle_res$estimates$ATE$CI[1],
      upper_ci = tmle_res$estimates$ATE$CI[2],
      p_value = tmle_res$estimates$ATE$pvalue,
      n_total = nrow(result$data),
      n_events = sum(result$data$event),
      stringsAsFactors = FALSE
    )
    
    # OR
    or_row <- data.frame(
      method = "TMLE",
      comparison = comparison_name,
      outcome = outcome_label,
      effect_measure = "Odds Ratio",
      estimate = tmle_res$estimates$OR$psi,
      lower_ci = tmle_res$estimates$OR$CI[1],
      upper_ci = tmle_res$estimates$OR$CI[2],
      p_value = NA,
      n_total = nrow(result$data),
      n_events = sum(result$data$event),
      stringsAsFactors = FALSE
    )
    
    rbind(rd_row, or_row)
  }
}

# 생존 곡선 플롯 생성 함수
create_survival_plot <- function(df, method_name, comparison_name, outcome_name, weight_var = NULL) {
  if (is.null(weight_var)) {
    # 가중치 없는 경우 - PSP
    surv_fit <- survfit(Surv(event_time, event) ~ treatment, data = df)
    plot_title <- paste0(method_name, " Survival: ", comparison_name, " - ", outcome_name)
  } else {
    # 가중치 있는 경우 - IPWT
    surv_fit <- survfit(Surv(event_time, event) ~ treatment, data = df, weights = df[[weight_var]])
    plot_title <- paste0(method_name, " Weighted Survival: ", comparison_name, " - ", outcome_name)
  }
  
  # 생존 곡선 플롯
  survival_plot <- ggplot() +
    # 생존 곡선
    geom_step(
      data = fortify(surv_fit),
      aes(x = time, y = surv, color = factor(strata)),
      size = 1
    ) +
    # 신뢰 구간
    geom_ribbon(
      data = fortify(surv_fit),
      aes(x = time, ymin = lower, ymax = upper, fill = factor(strata)),
      alpha = 0.2,
      linetype = 0
    ) +
    # 각종 설정
    scale_color_manual(
      values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
      labels = c("Comparator", "Semaglutide"),
      name = "Treatment"
    ) +
    scale_fill_manual(
      values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
      labels = c("Comparator", "Semaglutide"),
      name = "Treatment"
    ) +
    labs(
      x = "Days from index date",
      y = "Event-free probability",
      title = plot_title
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(survival_plot)
}

analyse_hba1c_trajectory <- function(results,
                                       panel_path      = "a1c_panel.csv",
                                       years_window    = 2,
                                       drop_crossovers = TRUE,
                                       sample_max      = 20000,
                                       smoothing_k     = NULL,
                                       seed            = 1) {

  ## 0. packages -----------------------------------------------------
  pkgs <- c("tidyverse", "lubridate", "lme4", "broom.mixed",
            "emmeans", "DHARMa", "viridis")
  suppressPackageStartupMessages(lapply(pkgs, library, character.only = TRUE))

  ## 1. load panel ---------------------------------------------------
  panel <- read_csv(panel_path) |>
           mutate(date_of_measurement = lubridate::parse_date_time(
                    date_of_measurement,
                    orders = c("ymd HMS", "ymd HM", "ymd")))

  ## 2. *dplyr*::select to avoid stats::select -----------------------
  matched_ids <- results$cohort |>
                 dplyr::select(person_id,
                              index_date,
                              treatment,
                              baseline_hba1c,
                              dplyr::any_of(c("crossover_date", "ipw_std")))

  ## 3. join & filter ------------------------------------------------
  a1c_df <- panel |>
            inner_join(matched_ids, by = "person_id") |>
            mutate(days_since_index = as.numeric(difftime(date_of_measurement,
                                                          index_date,
                                                          units = "days")),
                   years_since      = days_since_index / 365) |>
            filter(between(days_since_index, -90, years_window * 365))
  
  # Add baseline HbA1c as a measurement at time 0 if not already present
  baseline_rows <- matched_ids |>
    filter(!is.na(baseline_hba1c)) |>
    mutate(
      date_of_measurement = index_date,
      A1c = baseline_hba1c,
      days_since_index = 0,
      years_since = 0
    )
  
  # Select only columns that exist in both datasets
  common_cols <- intersect(names(a1c_df), names(baseline_rows))
  baseline_rows <- baseline_rows |> dplyr::select(all_of(common_cols))
  
  # Combine with existing measurements, removing duplicates at day 0
  a1c_df <- a1c_df |>
    filter(days_since_index != 0) |>  # Remove any existing day 0 measurements
    bind_rows(baseline_rows) |>
    arrange(person_id, days_since_index)

  if (drop_crossovers && "crossover_date" %in% names(a1c_df)) {
    a1c_df <- a1c_df |>
              filter(is.na(crossover_date) |
                     date_of_measurement < crossover_date)
  }

  ## 4. sample & plot (unchanged) ------------------------------------
  set.seed(seed)
  a1c_df_plot <- a1c_df |>
                 group_by(treatment) |>
                 group_modify(~ slice_sample(.x,
                                             n = min(sample_max, nrow(.x)))) |>
                 ungroup()

  form_txt <- if (is.null(smoothing_k)) "y ~ s(x, bs = 'cs')" 
              else                       sprintf("y ~ s(x, bs = 'cs', k = %d)",
                                                 smoothing_k)

  # Option to show individual points - set to FALSE for cleaner plot
  show_points <- FALSE
  
  # Option to use GAM smoothing or LMM predictions
  use_lmm_predictions <- TRUE  # Set to TRUE to show LMM model predictions
  
  if (use_lmm_predictions) {
    # Create prediction data for LMM
    pred_data <- expand.grid(
      years_since = seq(0, years_window, by = 0.05),
      treatment = c(0, 1)
    )
    
    # Get predictions from the LMM (will be fitted later in the function)
    # For now, create empty plot that will be updated after model fitting
    hb_plot <- ggplot() +
               scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                                   labels  = c("Comparator", "GLP‑1 agonists"),
                                   name    = "Treatment") +
               labs(x = "Years since index",
                    y = "HbA1c (%)",
                    title = "HbA1c trajectories in the matched cohort") +
               theme_minimal() +
               theme(legend.position = "bottom",
                     panel.grid.minor = element_blank())
  } else {
    # Original GAM smoothing approach
    hb_plot <- ggplot(a1c_df_plot,
                      aes(years_since, A1c, colour = factor(treatment))) +
               {if(show_points) geom_point(alpha = 0.15, size = 0.6) else NULL} +
               stat_smooth(method   = "gam",
                           formula  = as.formula(form_txt),
                           se       = TRUE,
                           linewidth = 1.5,
                           alpha = 0.2) +
               scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                                   labels  = c("Comparator", "GLP‑1 agonists"),
                                   name    = "Treatment") +
               labs(x = "Years since index",
                    y = "HbA1c (%)",
                    title = "HbA1c trajectories in the matched cohort") +
               theme_minimal() +
               theme(legend.position = "bottom",
                     panel.grid.minor = element_blank())
  }

  ## 5. mixed‑effects model & emmeans (unweighted with baseline adjustment) --------------------
  # Create baseline HbA1c variable for all observations
  baseline_hba1c_df <- a1c_df |>
    dplyr::filter(years_since == 0) |>
    dplyr::select(person_id, A1c) |>
    dplyr::rename(baseline_a1c = A1c)
  
  a1c_df <- a1c_df |>
    dplyr::left_join(baseline_hba1c_df, by = "person_id")
  
  # Get baseline means by treatment for reference
  baseline_means <- a1c_df |>
    dplyr::filter(years_since == 0) |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(
      unweighted_baseline = mean(A1c, na.rm = TRUE),
      n_baseline = n(),
      .groups = 'drop'
    )
  
  # Print baseline means for reference
  cat("\nBaseline HbA1c means (unweighted):\n")
  print(baseline_means)
  cat("\n")
  
  # Primary model: Unweighted mixed model with baseline HbA1c as covariate
  lmm_hba1c <- lmer(
    A1c ~ baseline_a1c + years_since * treatment + (1 + years_since || person_id),
    data = a1c_df,
    REML = TRUE
  )
  
  cat("\nLinear Mixed Model Results (unweighted, adjusted for baseline HbA1c):\n")
  print(summary(lmm_hba1c)$coefficients)
  cat("\n")

  # Calculate mean baseline HbA1c for emmeans reference
  mean_baseline_a1c <- mean(a1c_df$baseline_a1c, na.rm = TRUE)
  
  emm <- emmeans(lmm_hba1c,
                 ~ treatment | years_since,
                 at = list(years_since = c(0, 1, 2),
                          baseline_a1c = mean_baseline_a1c),
                 re_formula = NA)

  dharma_res <- DHARMa::simulateResiduals(lmm_hba1c)

  # Update plot with LMM predictions if requested
  if (use_lmm_predictions) {
    # Create prediction data with baseline HbA1c
    pred_data <- expand.grid(
      years_since = seq(0, years_window, by = 0.05),
      treatment = c(0, 1),
      baseline_a1c = mean_baseline_a1c
    )
    
    # Get predictions from the fitted LMM
    pred_data$fit <- predict(lmm_hba1c, newdata = pred_data, re.form = NA)
    
    # Calculate confidence intervals using emmeans
    emm_smooth <- emmeans(lmm_hba1c,
                          ~ treatment | years_since,
                          at = list(years_since = seq(0, years_window, by = 0.05),
                                   baseline_a1c = mean_baseline_a1c),
                          re_formula = NA)
    emm_df <- as.data.frame(emm_smooth)
    
    # Merge predictions with confidence intervals
    pred_data <- pred_data |>
      dplyr::left_join(emm_df |> 
                  dplyr::select(treatment, years_since, emmean, asymp.LCL, asymp.UCL) |>
                  dplyr::rename(lower.CL = asymp.LCL, upper.CL = asymp.UCL),
                by = c("treatment", "years_since"))
    
    # Create the plot with LMM predictions
    hb_plot <- ggplot(pred_data, aes(x = years_since, colour = factor(treatment))) +
      geom_ribbon(aes(ymin = lower.CL, ymax = upper.CL, fill = factor(treatment)), 
                  alpha = 0.2, colour = NA) +
      geom_line(aes(y = emmean), linewidth = 1.5) +
      {if(show_points) geom_point(data = a1c_df_plot,
                                  aes(y = A1c),
                                  alpha = 0.15, size = 0.6) else NULL} +
      scale_colour_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                          labels = c("SGLT2i", "GLP-1 agonists"),
                          name = "Treatment") +
      scale_fill_manual(values = c("0" = "#1f77b4", "1" = "#ff7f0e"),
                        labels = c("SGLT2i", "GLP-1 agonists"),
                        name = "Treatment", guide = "none") +
      labs(x = "Years since index",
           y = "HbA1c (%)",
           title = "HbA1c trajectories in the IPTW-weighted cohort (Linear Mixed Model)") +
      theme_minimal() +
      theme(legend.position = "bottom",
            panel.grid.minor = element_blank())
  }

  list(
    plot        = hb_plot,
    model       = lmm_hba1c,
    diagnostic  = dharma_res,
    emmeans     = emm
  )
}

## ───────────────────────── 0.  LOAD EXTRA PKGS  ──────────────────────────
needs <- c("mediation", "regmedint", "data.table")
for (p in needs) if (!requireNamespace(p, quietly = TRUE))
    install.packages(p, repos = "https://cloud.r-project.org")
library(mediation)
library(regmedint)
library(data.table)

## ─────────────────────── 1.  GENERIC HELPERS  ───────────────────────────
get_hba1c_mediator <- function(cohort_df,
                               panel_path = "a1c_panel.csv",
                               win_days   = 365,
                               weight_var = NULL,
                               use_weighted_mean = FALSE) {

  ## 1. read HbA1c panel
  panel_dt <- data.table::fread(panel_path) |>
              dplyr::mutate(date_of_measurement =
                              lubridate::parse_date_time(
                                date_of_measurement,
                                orders = c("ymd HMS","ymd HM","ymd")))

  ## 2. Determine if we should use the new weighted mean approach
  if (use_weighted_mean && !is.null(weight_var)) {
    # Use the new weighted mean calculation
    cat("Using weighted mean HbA1c calculation with weights:", weight_var, "\n")
    result <- calculate_weighted_mean_hba1c(cohort_df, panel_path, win_days, weight_var)
    # Rename the variable to match expected name in mediation analysis
    result <- result |>
      dplyr::rename(hba1c_mean_6mo = hba1c_weighted_mean_1yr)
    return(result)
  }
  
  ## 3. Original unweighted calculation (backward compatibility)
  hba1c_med <- cohort_df |>
    dplyr::select(person_id, index_date) |>
    dplyr::left_join(panel_dt, by = "person_id", relationship = "many-to-many") |>
    dplyr::mutate(
      date_of_measurement = as.Date(date_of_measurement),
      index_date          = as.Date(index_date),
      days                = as.numeric(date_of_measurement - index_date)
    ) |>
    dplyr::filter(dplyr::between(days, 0, win_days)) |>
    dplyr::group_by(person_id) |>
    dplyr::summarise(hba1c_mean_6mo = mean(A1c, na.rm = TRUE),
                     .groups = "drop")

  ## 4. bring the mean back to the full cohort
  dplyr::right_join(hba1c_med, cohort_df, by = "person_id")
}

# New function to calculate 1-year weighted mean HbA1c using IPTW weights
calculate_weighted_mean_hba1c <- function(cohort_df,
                                         panel_path = "a1c_panel.csv",
                                         win_days   = 365,
                                         weight_var = "ipw_std") {
  
  cat("\n=== Calculating 1-Year Weighted Mean HbA1c ===\n")
  
  ## 1. read HbA1c panel
  panel_dt <- data.table::fread(panel_path) |>
              dplyr::mutate(date_of_measurement =
                              lubridate::parse_date_time(
                                date_of_measurement,
                                orders = c("ymd HMS","ymd HM","ymd")))
  
  ## 2. Check if weights exist in cohort
  cat("DEBUG: Available columns in cohort:", paste(names(cohort_df), collapse = ", "), "\n")
  cat("DEBUG: Looking for weight variable:", weight_var, "\n")
  
  if (!weight_var %in% names(cohort_df)) {
    cat("Warning: Weight variable", weight_var, "not found in cohort. Using unweighted mean.\n")
    cat("DEBUG: Available columns were:", paste(names(cohort_df), collapse = ", "), "\n")
    weight_var <- NULL
  } else {
    cat("DEBUG: Weight variable", weight_var, "found successfully.\n")
  }
  
  ## 3. join measurements with cohort data
  # Select columns conditionally based on whether weight_var exists
  if (!is.null(weight_var)) {
    hba1c_weighted <- cohort_df |>
      dplyr::select(person_id, index_date, treatment, all_of(weight_var))
  } else {
    hba1c_weighted <- cohort_df |>
      dplyr::select(person_id, index_date, treatment)
  }
  
  hba1c_weighted <- hba1c_weighted |>
    dplyr::left_join(panel_dt, by = "person_id", relationship = "many-to-many") |>
    dplyr::mutate(
      date_of_measurement = as.Date(date_of_measurement),
      index_date          = as.Date(index_date),
      days                = as.numeric(date_of_measurement - index_date)
    ) |>
    dplyr::filter(dplyr::between(days, 0, win_days))
  
  ## 4. Calculate weighted mean for each person
  if (!is.null(weight_var)) {
    # Calculate weighted mean using person-level weights
    hba1c_summary <- hba1c_weighted |>
      dplyr::group_by(person_id, treatment) |>
      dplyr::summarise(
        n_measurements = n(),
        hba1c_weighted_mean_1yr = weighted.mean(A1c, rep(first(.data[[weight_var]]), n()), na.rm = TRUE),
        hba1c_unweighted_mean = mean(A1c, na.rm = TRUE),
        weight_used = first(.data[[weight_var]]),
        .groups = "drop"
      )
  } else {
    # Fallback to unweighted mean
    hba1c_summary <- hba1c_weighted |>
      dplyr::group_by(person_id, treatment) |>
      dplyr::summarise(
        n_measurements = n(),
        hba1c_weighted_mean_1yr = mean(A1c, na.rm = TRUE),
        hba1c_unweighted_mean = mean(A1c, na.rm = TRUE),
        .groups = "drop"
      )
  }
  
  ## 5. Diagnostic output
  cat("\nHbA1c measurement summary:\n")
  measurement_summary <- hba1c_summary |>
    dplyr::group_by(treatment) |>
    dplyr::summarise(
      n_patients = n(),
      avg_measurements = mean(n_measurements),
      patients_with_1_plus = sum(n_measurements >= 1),
      patients_with_2_plus = sum(n_measurements >= 2),
      mean_weighted_hba1c = mean(hba1c_weighted_mean_1yr, na.rm = TRUE),
      sd_weighted_hba1c = sd(hba1c_weighted_mean_1yr, na.rm = TRUE),
      .groups = "drop"
    )
  print(measurement_summary)
  
  ## 6. Merge back to full cohort
  result <- dplyr::left_join(
    cohort_df,
    hba1c_summary |> dplyr::select(person_id, hba1c_weighted_mean_1yr, n_measurements),
    by = "person_id"
  )
  
  cat("\nPatients with weighted mean HbA1c:", sum(!is.na(result$hba1c_weighted_mean_1yr)), 
      "out of", nrow(result), "\n")
  
  return(result)
}

add_outcome_post_mediator <- function(df,
                                      outcome_var,
                                      data_cut,
                                      late_onset = FALSE,
                                      age_cut   = 50) {

  df %>%
    mutate(outcome_date = as.Date(.data[[outcome_var]]),
           start_fu     = index_date + 180,             # follow‑up starts after mediator
           censor_date  = pmin(as.Date(EHRmaxDT), data_cut, na.rm = TRUE),
           raw_event    = !is.na(outcome_date) &
                          outcome_date > start_fu &
                          outcome_date <= censor_date,
           time         = as.numeric(pmin(outcome_date, censor_date, na.rm = TRUE) -
                                     start_fu),
           event        = as.integer(raw_event)) %>%
    { if (late_onset) {
        mutate(.,
               age_at_event = age + time / 365.25,
               event        = as.integer(event == 1 & age_at_event >= age_cut))
      } else .
    } %>%
    filter(time >= 0)
}

# 수정된 중재 분석 함수
run_mediation <- function(df, covars, boot = 1000) {
  # 디버깅 정보 출력
  cat("\nMediation function input diagnostics:\n")
  cat("Total rows:", nrow(df), "\n")
  cat("Treatment distribution:", table(df$treatment), "\n")
  cat("Event distribution:", table(df$event), "\n")
  cat("HbA1c summary:", summary(df$hba1c_mean_6mo), "\n")
  
  # 필요한 열만 선택하고 누락된 행 제거
  keep_vars <- c("treatment", "hba1c_mean_6mo", "time", "event", covars)
  mod_df <- df[complete.cases(df[keep_vars]), keep_vars]
  
  cat("\nAfter complete case filtering:\n")
  cat("Remaining rows:", nrow(mod_df), "\n")
  cat("Treatment distribution:", table(mod_df$treatment), "\n")
  cat("Event distribution:", table(mod_df$event), "\n")

  if (length(unique(mod_df$treatment)) < 2 || sum(mod_df$event) == 0 || nrow(mod_df) < 10) {
    cat("[skip] 데이터 부족 - 치료군 변동 부족, 이벤트 없음, 또는 샘플 크기가 너무 작음\n")
    return(list(imai = NULL, reg = NULL, covs_used = covars))
  }
  
  # 1. 중재자 모델 - 변수 범위 문제 해결을 위해 직접 작성
  m_formula_str <- paste("hba1c_mean_6mo ~ treatment +", paste(covars, collapse = " + "))
  cat("중재자 모델 공식:", m_formula_str, "\n")
  
  m_fit <- tryCatch({
    lm(as.formula(m_formula_str), data = mod_df)
  }, error = function(e) {
    cat("중재자 모델 오류:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(m_fit)) {
    return(list(imai = NULL, reg = NULL, covs_used = covars))
  }
  
  cat("중재자 모델 성공!\n")
  cat("\n중재자 모델 요약:\n")
  print(summary(m_fit)$coefficients[1:2,])
  
  # 2. 결과 모델
  y_formula_str <- paste("Surv(time, event) ~ treatment + hba1c_mean_6mo +", 
                          paste(covars, collapse = " + "))
  cat("결과 모델 공식:", y_formula_str, "\n")
  
  y_fit <- tryCatch({
    survival::coxph(as.formula(y_formula_str), data = mod_df)
  }, error = function(e) {
    cat("결과 모델 오류:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(y_fit)) {
    return(list(imai = NULL, reg = NULL, covs_used = covars))
  }
  
  cat("결과 모델 성공!\n")
  cat("\n결과 모델 요약:\n")
  print(summary(y_fit)$coefficients[1:2,])
  
  # 3. mediate() 함수 직접 호출
  cat("\nmediate() 실행 중...\n")
  
  # 여기서 직접 mediate 함수 호출
  med_out <- NULL
  try({
    med_out <- mediation::mediate(
      model.m = m_fit,
      model.y = y_fit,
      treat = "treatment",
      mediator = "hba1c_mean_6mo",
      sims = boot,
      boot = TRUE
    )
    cat("mediate() 성공!\n")
  }, silent = FALSE)
  
  # 결과 반환
  list(imai = med_out, reg = NULL, covs_used = covars, m_fit = m_fit, y_fit = y_fit)
}

# 중재 효과 수동 계산 함수
calculate_mediation_effect <- function(med_result) {
  if (is.null(med_result$m_fit) || is.null(med_result$y_fit)) {
    cat("유효한 모델 결과가 없습니다.\n")
    return(NULL)
  }
  
  # 효과 추출
  a_path <- coef(med_result$m_fit)["treatment"]  # 치료 -> 중재자
  b_path <- coef(med_result$y_fit)["hba1c_mean_6mo"]  # 중재자 -> 결과
  c_path <- coef(med_result$y_fit)["treatment"]  # 직접 효과
  
  # 간접 효과 계산(근사치)
  indirect_effect <- a_path * b_path
  
  # 결과 표시
  cat("\n수동 계산된 중재 효과:\n")
  cat("A 경로 (치료 -> HbA1c):", a_path, "\n")
  cat("B 경로 (HbA1c -> 결과):", b_path, "\n")
  cat("C 경로 (직접 효과):", c_path, "\n")
  cat("간접 효과 (A*B, 근사치):", indirect_effect, "\n")
  cat("총 효과 (근사치):", c_path + indirect_effect, "\n")
  
  # 경로도 생성
  path_data <- data.frame(
    x = c(1, 3, 2),
    y = c(2, 2, 1),
    label = c("치료", "결과", "HbA1c")
  )
  
  path_plot <- ggplot() +
    # 노드
    geom_point(data = path_data, aes(x = x, y = y), size = 20, 
               color = c("#ff7f0e", "#1f77b4", "#2ca02c")) +
    # 노드 레이블
    geom_text(data = path_data, aes(x = x, y = y, label = label), 
              fontface = "bold", color = "white") +
    # 화살표
    geom_segment(aes(x = 1.3, y = 2, xend = 2.7, yend = 2),
                 arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    geom_segment(aes(x = 1, y = 1.7, xend = 1.7, yend = 1.3),
                 arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    geom_segment(aes(x = 2.3, y = 1.3, xend = 2.9, yend = 1.7),
                 arrow = arrow(length = unit(0.3, "cm")), size = 1) +
    # 효과 크기 레이블
    annotate("text", x = 2, y = 2.2, 
             label = sprintf("직접: %.3f", c_path),
             fontface = "bold") +
    annotate("text", x = 1.4, y = 1.4, 
             label = sprintf("%.3f", a_path),
             fontface = "bold") +
    annotate("text", x = 2.6, y = 1.4, 
             label = sprintf("%.3f", b_path),
             fontface = "bold") +
    # 제목 및 테마
    labs(title = "HbA1c를 통한 중재 효과 경로") +
    theme_void() +
    xlim(0, 4) +
    ylim(0, 3)
  
  print(path_plot)
  
  # 결과 반환
  list(
    a_path = a_path,
    b_path = b_path, 
    c_path = c_path,
    indirect_effect = indirect_effect,
    total_effect = c_path + indirect_effect,
    plot = path_plot
  )
}