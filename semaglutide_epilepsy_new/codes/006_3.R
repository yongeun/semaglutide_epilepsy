# 006_3
# ---- IPTW 분석 ----
late_cutoff <- 50

analyze_with_ipwt <- function() {
  results_table <- data.frame()
  mediation_results_list <- list()  # 중재 분석 결과 저장용 리스트
  
  for (cmp in comparisons) {
    # 약물 개념 ID 가져오기
    exposure_concepts <- get_drug_concepts(drug_classes[[cmp$exposure]])
    comparator_concepts <- get_drug_concepts(drug_classes[[cmp$comparator]])
    
    # 인덱스 날짜 구축
    idx_exposure <- build_exposure_idx(exposure_concepts, paste0(tolower(cmp$exposure), "_idx"))
    idx_comparator <- build_exposure_idx(comparator_concepts, paste0(tolower(cmp$comparator), "_idx"))
    
    # 코호트 구성
    exposure_wide <- full_join(idx_exposure, idx_comparator, by = "person_id")
    cohort_pair <- bind_rows(
      exposure_wide %>%
        filter(!is.na(.data[[paste0(tolower(cmp$exposure), "_idx")]]) & 
                 is.na(.data[[paste0(tolower(cmp$comparator), "_idx")]])) %>%
        mutate(index_date = .data[[paste0(tolower(cmp$exposure), "_idx")]], treatment = 1),
      exposure_wide %>%
        filter(!is.na(.data[[paste0(tolower(cmp$comparator), "_idx")]]) & 
                 is.na(.data[[paste0(tolower(cmp$exposure), "_idx")]])) %>%
        mutate(index_date = .data[[paste0(tolower(cmp$comparator), "_idx")]], treatment = 0)
    )
    
    # 속성 설정
    attr(cohort_pair, "exposure_idx") <- paste0(tolower(cmp$exposure), "_idx")
    attr(cohort_pair, "comparator_idx") <- paste0(tolower(cmp$comparator), "_idx")
    
    for (outc in outcomes) {
      cat("\n========= IPWT Analysis:", cmp$name, "-", outc$label, "=========\n")
      
      # Initialize flow tracking for this comparison-outcome pair
      flow_key <- paste(cmp$name, outc$label, sep = " - ")
      flow_data <- list()
      
      # 분석 데이터 준비 (age_at_end 포함)
      analytic_df_pre_bmi <- cohort_pair %>%
        left_join(merged_df_2, by = "person_id") %>%
        filter(dm == 1) %>%
        mutate(
          dob = as.Date(substr(date_of_birth, 1, 10), format = "%Y/%m/%d"),
          age_at_end = as.numeric(difftime(data_cut_date, dob, units = "days")) / 365.25
        )
      
      # DEBUG: Check initial cohort before any exclusions
      cat("\n=== INITIAL COHORT DEBUG ===\n")
      cat("Total patients after DM filter:", nrow(analytic_df_pre_bmi), "\n")
      cat("Treatment distribution:", table(analytic_df_pre_bmi$treatment), "\n")
      
      # TRACK: Initial numbers after DM filter
      flow_data$after_dm_filter <- nrow(analytic_df_pre_bmi)
      flow_data$excluded_non_dm <- "Not tracked in current pipeline"
      flow_data$initial_dm <- flow_data$after_dm_filter
      
      # DEBUG: Check raw outcome data before exclusions
      outcome_var_name <- outc$var
      cat("\nRaw outcome variable (", outcome_var_name, ") summary:\n")
      outcome_dates <- analytic_df_pre_bmi[[outcome_var_name]]
      cat("Total patients with outcome dates:", sum(!is.na(outcome_dates)), "\n")
      cat("Total patients without outcome dates:", sum(is.na(outcome_dates)), "\n")
      if(sum(!is.na(outcome_dates)) > 0) {
        cat("Outcome date range:", min(outcome_dates, na.rm = TRUE), "to", max(outcome_dates, na.rm = TRUE), "\n")
      }
      cat("===========================\n")
      
             # STEP 1: Apply exclusions for pre-existing conditions
      analytic_df_after_exclusion <- analytic_df_pre_bmi %>%
        # EXCLUDE PATIENTS WITH PRE-EXISTING CONDITIONS AT BASELINE FOR NEW ONSET ANALYSIS
        {
          before_exclusion <- nrow(.)
          cat("\nBefore excluding pre-existing", outc$label, "patients:", before_exclusion, "\n")
          
          # DEBUG: Check for pre-existing conditions using actual dates
          cat("\n=== DEBUGGING DATE-BASED EXCLUSIONS ===\n")
          if (outc$label %in% c("Epilepsy/Seizure", "Early-onset Epilepsy/Seizure", "Late-onset Epilepsy/Seizure")) {
            pre_existing <- !is.na(.$epilepsy_or_seizure_start_date) & 
                           as.Date(.$epilepsy_or_seizure_start_date) < as.Date(.$index_date)
            cat("Epilepsy/Seizure diagnosed BEFORE index date (excluding):", sum(pre_existing, na.rm = TRUE), "\n")
            cat("Patients eligible for analysis:", sum(!pre_existing, na.rm = TRUE), "\n")
          } else if (outc$label == "ADRD") {
            pre_existing <- !is.na(.$adrd_start_date) & 
                           as.Date(.$adrd_start_date) < as.Date(.$index_date)
            cat("ADRD diagnosed BEFORE index date (excluding):", sum(pre_existing, na.rm = TRUE), "\n")
            cat("Patients eligible for analysis:", sum(!pre_existing, na.rm = TRUE), "\n")
          } else if (outc$label == "stroke") {
            pre_existing <- !is.na(.$stroke_start_date) & 
                           as.Date(.$stroke_start_date) < as.Date(.$index_date)
            cat("Stroke diagnosed BEFORE index date (excluding):", sum(pre_existing, na.rm = TRUE), "\n")
            cat("Patients eligible for analysis:", sum(!pre_existing, na.rm = TRUE), "\n")
          }
          cat("=======================================\n")
          
          # PROPER PRE-EXISTING CONDITION EXCLUSION USING ACTUAL DATES
          # Compare diagnosis dates to index dates to identify true baseline conditions
          result <- filter(.,
            if (outc$label %in% c("Epilepsy/Seizure", "Early-onset Epilepsy/Seizure", "Late-onset Epilepsy/Seizure")) {
              # Exclude patients with epilepsy/seizure diagnosed BEFORE index date
              is.na(epilepsy_or_seizure_start_date) | 
              as.Date(epilepsy_or_seizure_start_date) >= as.Date(index_date)
            } else if (outc$label == "ADRD") {
              # Exclude patients with ADRD diagnosed BEFORE index date
              is.na(adrd_start_date) | 
              as.Date(adrd_start_date) >= as.Date(index_date)
            } else if (outc$label == "stroke") {
              # Exclude patients with stroke diagnosed BEFORE index date
              is.na(stroke_start_date) | 
              as.Date(stroke_start_date) >= as.Date(index_date)
            } else {
              # For other outcomes, keep everyone
              TRUE
            }
          )
          
          after_exclusion <- nrow(result)
          excluded_count <- before_exclusion - after_exclusion
          cat("After excluding pre-existing", outc$label, "patients:", after_exclusion, 
              "(excluded:", excluded_count, "patients)\n")
          
          # TRACK: Numbers after pre-existing exclusion
          flow_data$excluded_preexisting <- excluded_count
          flow_data$after_preexisting_exclusion <- after_exclusion
          
          # DEBUG: Check outcome data after exclusions but before followup_and_event
          cat("\n=== AFTER EXCLUSION DEBUG ===\n")
          outcome_dates_after <- result[[outcome_var_name]]
          cat("Patients with outcome dates after exclusion:", sum(!is.na(outcome_dates_after)), "\n")
          cat("Treatment distribution after exclusion:", table(result$treatment), "\n")
          if(sum(!is.na(outcome_dates_after)) > 0) {
            cat("Outcome date range after exclusion:", min(outcome_dates_after, na.rm = TRUE), "to", max(outcome_dates_after, na.rm = TRUE), "\n")
            # Check if outcome dates are after index dates
            valid_outcomes <- !is.na(outcome_dates_after) & !is.na(result$index_date) & 
                             as.Date(outcome_dates_after) > as.Date(result$index_date)
            cat("Outcomes occurring AFTER index date:", sum(valid_outcomes, na.rm = TRUE), "\n")
          }
          cat("=============================\n")
          
          result
        }

      # --- 기준 BMI 및 BMI 범주 추가 ---
      baseline_bmi_data_for_cohort <- analytic_df_after_exclusion %>%
        dplyr::select(person_id, index_date) %>%
        distinct() %>%
        left_join(bmi_panel_df_loaded, by = "person_id", relationship = "many-to-many") %>%
        filter(!is.na(weight_date) & !is.na(index_date)) %>%
        mutate(
          # Calculate absolute difference between BMI date and index date
          days_from_index = abs(as.numeric(difftime(weight_date, index_date, units = "days")))
        ) %>%
        group_by(person_id) %>%
        # Sort by smallest difference first, then by most recent date for ties
        arrange(days_from_index, desc(weight_date)) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(person_id, baseline_bmi = bmi, baseline_bmi_date = weight_date, days_from_index)
      
      # --- 기준 HbA1c 추가 ---
      # Load HbA1c panel data
      hba1c_panel_file <- "a1c_panel.csv"
      if (!file.exists(hba1c_panel_file)) {
        system(str_glue("gsutil cp gs://{WORKSPACE_BUCKET}/data/{hba1c_panel_file} ."), intern=T)
      }
      hba1c_panel_df_loaded <- read_csv(hba1c_panel_file, col_types = cols())
      
      baseline_hba1c_data_for_cohort <- analytic_df_after_exclusion %>%
        dplyr::select(person_id, index_date) %>%
        distinct() %>%
        left_join(hba1c_panel_df_loaded, by = "person_id", relationship = "many-to-many") %>%
        filter(!is.na(date_of_measurement) & !is.na(index_date)) %>%
        mutate(
          # Calculate absolute difference between HbA1c date and index date
          days_from_index = abs(as.numeric(difftime(date_of_measurement, index_date, units = "days")))
        ) %>%
        group_by(person_id) %>%
        # Sort by smallest difference first, then by most recent date for ties
        arrange(days_from_index, desc(date_of_measurement)) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(person_id, baseline_hba1c = A1c, baseline_hba1c_date = date_of_measurement, days_from_index_hba1c = days_from_index)

      # Add diagnostic output to check the baseline BMI selection
      cat("\nBaseline BMI Selection Summary:\n")
      bmi_timing_summary <- baseline_bmi_data_for_cohort %>%
        summarise(
          n_patients_with_bmi = n(),
          median_days_from_index = median(days_from_index, na.rm = TRUE),
          mean_days_from_index = mean(days_from_index, na.rm = TRUE),
          max_days_from_index = max(days_from_index, na.rm = TRUE),
          .groups = "drop"
        )
      print(bmi_timing_summary)
      
      # Add diagnostic output to check the baseline HbA1c selection
      cat("\nBaseline HbA1c Selection Summary:\n")
      hba1c_timing_summary <- baseline_hba1c_data_for_cohort %>%
        summarise(
          n_patients_with_hba1c = n(),
          median_days_from_index = median(days_from_index_hba1c, na.rm = TRUE),
          mean_days_from_index = mean(days_from_index_hba1c, na.rm = TRUE),
          max_days_from_index = max(days_from_index_hba1c, na.rm = TRUE),
          .groups = "drop"
        )
      print(hba1c_timing_summary)
      
      # STEP 2: Add BMI and HbA1c data to the cohort after exclusions
      analytic_df_with_bmi <- analytic_df_after_exclusion %>%
        left_join(baseline_bmi_data_for_cohort %>% dplyr::select(person_id, baseline_bmi), by = "person_id") %>%
        left_join(baseline_hba1c_data_for_cohort %>% dplyr::select(person_id, baseline_hba1c), by = "person_id") %>%
        mutate(
          # Add calendar year of index date as a covariate
          index_year = factor(year(as.Date(index_date))),
          
          bmi_category = case_when(
            baseline_bmi < 18.5 ~ "Underweight",
            baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
            baseline_bmi >= 25 & baseline_bmi < 30 ~ "Overweight",
            baseline_bmi >= 30 & baseline_bmi < 35 ~ "Obese I",
            baseline_bmi >= 35 & baseline_bmi < 40 ~ "Obese II",
            baseline_bmi >= 40 ~ "Obese III",
            TRUE ~ NA_character_
          ),
          bmi_category = factor(bmi_category, levels = c("Underweight", "Normal", "Overweight", "Obese I", "Obese II", "Obese III"))
        )
      
      # DEBUG: Check data before followup_and_event
      cat("\n=== BEFORE FOLLOWUP_AND_EVENT ===\n")
      cat("Total patients going into followup_and_event:", nrow(analytic_df_with_bmi), "\n")
      cat("Patients with outcome dates:", sum(!is.na(analytic_df_with_bmi[[outcome_var_name]])), "\n")
      
      # STEP 3: Apply followup_and_event function
      analytic_df <- analytic_df_with_bmi %>%
        followup_and_event(outcome_var = outc$var, data_cut_date = data_cut_date, late_onset = outc$late_onset, early_onset = outc$early_onset)
      
      # DEBUG: Check data after followup_and_event
      cat("\n=== AFTER FOLLOWUP_AND_EVENT ===\n")
      cat("Total patients after followup_and_event:", nrow(analytic_df), "\n")
      cat("Patients with events (event = 1):", sum(analytic_df$event == 1, na.rm = TRUE), "\n")
      cat("Patients without events (event = 0):", sum(analytic_df$event == 0, na.rm = TRUE), "\n")
      cat("Event distribution by treatment:\n")
      print(table(analytic_df$treatment, analytic_df$event, useNA = "ifany"))
      
      # Check event_time distribution
      cat("\nEvent time summary:\n")
      print(summary(analytic_df$event_time))
      cat("Patients with negative event_time (should be 0):", sum(analytic_df$event_time < 0, na.rm = TRUE), "\n")
      
      # Check raw_event vs final event
      if("raw_event" %in% names(analytic_df)) {
        cat("\nRaw events vs final events:\n")
        cat("Raw events:", sum(analytic_df$raw_event == 1, na.rm = TRUE), "\n")
        cat("Final events:", sum(analytic_df$event == 1, na.rm = TRUE), "\n")
        if(outc$late_onset) {
          cat("Late onset filter applied (age >= 50)\n")
          if("age_at_event" %in% names(analytic_df)) {
            cat("Events with age_at_event >= 50:", sum(analytic_df$raw_event == 1 & analytic_df$age_at_event >= 50, na.rm = TRUE), "\n")
          }
        }
        if(outc$early_onset) {
          cat("Early onset filter applied (age < 50)\n")
          if("age_at_event" %in% names(analytic_df)) {
            cat("Events with age_at_event < 50:", sum(analytic_df$raw_event == 1 & analytic_df$age_at_event < 50, na.rm = TRUE), "\n")
          }
        }
      }
      cat("=============================\n")
      
      # TRACK: Numbers after followup_and_event
      pre_followup_count <- nrow(analytic_df_with_bmi)
      post_followup_count <- nrow(analytic_df)
      flow_data$excluded_followup <- pre_followup_count - post_followup_count
      flow_data$after_followup_filter <- post_followup_count

      # --- 기준 BMI 비교 (가중치 적용 전) ---
      cat("\nBaseline BMI Comparison (Before Weighting):\n")
      summary_bmi_before <- analytic_df %>%
        group_by(treatment) %>%
        summarise(
          n = n(),
          mean_bmi = mean(baseline_bmi, na.rm = TRUE),
          sd_bmi = sd(baseline_bmi, na.rm = TRUE),
          median_bmi = median(baseline_bmi, na.rm = TRUE),
          missing_bmi = sum(is.na(baseline_bmi)),
          .groups = "drop"
        )
      print(summary_bmi_before)
      cat("\nBaseline BMI Category Distribution (Before Weighting):\n")
      print(table(analytic_df$treatment, analytic_df$bmi_category, useNA = "ifany"))
      
      # --- 기준 HbA1c 비교 (가중치 적용 전) ---
      cat("\nBaseline HbA1c Comparison (Before Weighting):\n")
      summary_hba1c_before <- analytic_df %>%
        group_by(treatment) %>%
        summarise(
          n = n(),
          mean_hba1c = mean(baseline_hba1c, na.rm = TRUE),
          sd_hba1c = sd(baseline_hba1c, na.rm = TRUE),
          median_hba1c = median(baseline_hba1c, na.rm = TRUE),
          missing_hba1c = sum(is.na(baseline_hba1c)),
          .groups = "drop"
        )
      print(summary_hba1c_before)
      
      # 전체 샘플 크기 저장
      total_sample_size <- nrow(analytic_df)
      total_events <- sum(analytic_df$event, na.rm = TRUE)
      cat("\nTotal sample size:", total_sample_size, "\n")
      cat("Total events:", total_events, "\n")
      
      # 인덱스 연도 분포 확인
      cat("\nIndex Year Distribution:\n")
      print(table(analytic_df$index_year, useNA = "ifany"))
      
      # 제외할 약물 클래스 변수
      split_comparator <- unlist(strsplit(cmp$comparator, "_"))
      exclude_vars <- unique(c(cmp$exposure, split_comparator))
      
      # Outcome-specific exclusions
      if (outc$label == "ADRD") {
        exclude_vars <- unique(c(exclude_vars, "dem"))
        cat("Excluding 'dem' covariate for ADRD outcome analysis\n")
      }
      if (outc$label == "stroke" || outc$label == "Stroke") {
        exclude_vars <- unique(c(exclude_vars, "cvd"))
        cat("Excluding 'cvd' covariate for Stroke outcome analysis\n")
      }
      
      # IPWT 분석
      set.seed(12345) 
      ipwt_result <- tryCatch({
        run_ipwt_and_cox(analytic_df, exclude_vars = exclude_vars, trim_threshold = 0.01)
      }, error = function(e) {
        cat("Error in IPWT analysis:", e$message, "\n")
        return(NULL)
      })
      
      if (is.null(ipwt_result)) {
        cat("Skipping this outcome due to error in IPWT analysis.\n")
        next
      }
      
      # --- late_onset일 때만, 분석 후 코호트에서 연령 필터 적용 ---
      pre_age_filter_count <- nrow(ipwt_result$cohort)
      if (isTRUE(outc$late_onset)) {
        ipwt_result$cohort <- ipwt_result$cohort %>%
          filter(age_at_end >= late_cutoff)
        cat("After age filtering (>=", late_cutoff, "): samples =", nrow(ipwt_result$cohort), 
            ", events =", sum(ipwt_result$cohort$event), "\n")
        flow_data$excluded_age <- pre_age_filter_count - nrow(ipwt_result$cohort)
      } else {
        flow_data$excluded_age <- 0
      }
      
      # TRACK: Final cohort numbers
      flow_data$final_cohort <- nrow(ipwt_result$cohort)
      flow_data$total_events <- sum(ipwt_result$cohort$event, na.rm = TRUE)
      
      # --- 기준 BMI 비교 (IPWT 가중치 적용 후) ---
      cat("\nBaseline BMI Comparison (After IPWT Weighting):\n")
      if (nrow(ipwt_result$cohort) > 0 && "ipw_std" %in% names(ipwt_result$cohort)) {
          summary_bmi_after <- ipwt_result$cohort %>%
            filter(!is.na(baseline_bmi)) %>% # NA 제거 후 가중 평균 계산
            group_by(treatment) %>%
            summarise(
              n_in_bmi_calc = n(), # 가중 평균 계산에 사용된 샘플 수
              mean_bmi_weighted = weighted.mean(baseline_bmi, ipw_std, na.rm = TRUE),
              # 가중 SD는 더 복잡하며, 필요 시 stats::weighted.sd 또는 유사 함수 사용
              .groups = "drop"
            )
          print(summary_bmi_after)
          
          cat("\nBaseline BMI Category Distribution (After IPWT Weighting - Weighted Proportions):\n")
          weighted_bmi_cat_summary <- ipwt_result$cohort %>%
            filter(!is.na(bmi_category)) %>% # NA 제거
            group_by(treatment, bmi_category) %>%
            summarise(weighted_n = sum(ipw_std, na.rm = TRUE), .groups = 'drop_last') %>%
            mutate(weighted_prop = weighted_n / sum(weighted_n, na.rm = TRUE)) %>%
            ungroup()
          print(weighted_bmi_cat_summary)
          
          # --- 기준 HbA1c 비교 (IPWT 가중치 적용 후) ---
          cat("\nBaseline HbA1c Comparison (After IPWT Weighting):\n")
          summary_hba1c_after <- ipwt_result$cohort %>%
            filter(!is.na(baseline_hba1c)) %>% # NA 제거 후 가중 평균 계산
            group_by(treatment) %>%
            summarise(
              n_in_hba1c_calc = n(), # 가중 평균 계산에 사용된 샘플 수
              mean_hba1c_weighted = weighted.mean(baseline_hba1c, ipw_std, na.rm = TRUE),
              .groups = "drop"
            )
          print(summary_hba1c_after)
      } else {
          cat("Skipping weighted BMI comparison: cohort is empty or weights are missing.\n")
      }
      
      # 결과 저장
      key <- paste0(cmp$name, " - ", outc$label)
      ipwt_results[[key]] <- ipwt_result
      
      # Create and save before vs after weighting characteristics table
      tryCatch({
        weighted_chars <- create_weighted_characteristics_table(ipwt_result, cmp$name, outc$label)
        
        # Save individual table
        write.csv(
          weighted_chars,
          paste0("ipwt_weighted_characteristics_", 
                 gsub(" ", "_", cmp$name), "_", 
                 gsub("[/ ]", "_", outc$label), ".csv"),
          row.names = FALSE
        )
        
        cat("Saved weighted characteristics comparison table\n")
        
      }, error = function(e) {
        cat("Error creating weighted characteristics table:", e$message, "\n")
      })
      
      # (Baseline characteristics tables will be generated in consolidated section at end)
      
      # Save cohort data using All of Us method
      # Replace df with THE NAME OF YOUR DATAFRAME
      my_dataframe <- ipwt_result$cohort
      
      # Replace 'test.csv' with THE NAME of the file you're going to store in the bucket
      destination_filename <- paste0("ipwt_cohort_", 
                                   gsub(" ", "_", cmp$name), "_", 
                                   gsub("[/ ]", "_", outc$label), ".csv")
      
      # store the dataframe in current workspace
      write_excel_csv(my_dataframe, destination_filename)
      
      # Get the bucket name
      my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
      
      # Copy the file from current workspace to the bucket
      system(paste0("gsutil cp ./", destination_filename, " ", my_bucket, "/data/"), intern=T)
      
      cat("Saved cohort data to bucket:", destination_filename, "\n")
      
      # 공변량 균형 확인
      cat("\nStandardized Mean Differences - Before Weighting (Top 10):\n")
      print(ipwt_result$balance_before %>% 
              arrange(desc(abs(std_diff))) %>% 
              head(10) %>%
              mutate(std_diff = round(std_diff, 4)))
      
      cat("\nStandardized Mean Differences - After Weighting (Top 10):\n")
      print(ipwt_result$balance_after %>% 
              arrange(desc(abs(std_diff))) %>% 
              head(10) %>%
              mutate(std_diff = round(std_diff, 4)))
      
      # 가중치 요약
      cat("\nIPWT Weight Summary:\n")
      print(summary(ipwt_result$cohort$ipw_std))
      
      # 가중 이벤트 수
      cat("\nWeighted Event Counts:\n")
      weighted_events <- ipwt_result$cohort %>%
        group_by(treatment) %>%
        summarise(
          n = n(),
          events = sum(event),
          weighted_n = sum(ipw_std),
          weighted_events = sum(event * ipw_std),
          raw_rate = events / n,
          weighted_rate = weighted_events / weighted_n
        )
      print(weighted_events)
      
      cat("\nWeighted Cox Results:\n")
      print(ipwt_result$cox)
      
      # Print SMD values in formatted table before visualization
      cat("\n", paste(rep("=", 80), collapse=""), "\n")
      cat("STANDARDIZED MEAN DIFFERENCES TABLE:", cmp$name, "-", outc$label, "\n")
      cat(paste(rep("=", 80), collapse=""), "\n")
      
      # Combine before and after SMD data for comparison
      smd_comparison_table <- merge(
        ipwt_result$balance_before[, c("variable", "std_diff")],
        ipwt_result$balance_after[, c("variable", "std_diff")],
        by = "variable",
        suffixes = c("_before", "_after")
      ) %>%
        mutate(
          improvement = abs(std_diff_before) - abs(std_diff_after),
          abs_smd_before = abs(std_diff_before),
          abs_smd_after = abs(std_diff_after)
        ) %>%
        arrange(desc(abs_smd_before))
      
      # Print detailed SMD table
      cat(sprintf("%-25s %10s %10s %10s %8s\n", 
                  "Variable", "Before", "After", "Abs_Before", "Abs_After"))
      cat(paste(rep("-", 75), collapse=""), "\n")
      
      for(i in 1:nrow(smd_comparison_table)) {
        row <- smd_comparison_table[i, ]
        cat(sprintf("%-25s %10.4f %10.4f %10.4f %10.4f\n",
                    substr(row$variable, 1, 25),  # Truncate long variable names
                    row$std_diff_before,
                    row$std_diff_after,
                    row$abs_smd_before,
                    row$abs_smd_after))
      }
      
      cat(paste(rep("-", 75), collapse=""), "\n")
      cat("Summary Statistics:\n")
      cat(sprintf("Variables with |SMD| > 0.1 before IPTW: %d\n", 
                  sum(abs(smd_comparison_table$std_diff_before) > 0.1)))
      cat(sprintf("Variables with |SMD| > 0.1 after IPTW: %d\n", 
                  sum(abs(smd_comparison_table$std_diff_after) > 0.1)))
      cat(sprintf("Variables improved (reduced |SMD|): %d\n", 
                  sum(smd_comparison_table$improvement > 0)))
      cat(sprintf("Mean |SMD| before: %.4f\n", 
                  mean(abs(smd_comparison_table$std_diff_before))))
      cat(sprintf("Mean |SMD| after: %.4f\n", 
                  mean(abs(smd_comparison_table$std_diff_after))))
      cat(sprintf("Maximum |SMD| before: %.4f (%s)\n", 
                  max(abs(smd_comparison_table$std_diff_before)),
                  smd_comparison_table$variable[which.max(abs(smd_comparison_table$std_diff_before))]))
      cat(sprintf("Maximum |SMD| after: %.4f (%s)\n", 
                  max(abs(smd_comparison_table$std_diff_after)),
                  smd_comparison_table$variable[which.max(abs(smd_comparison_table$std_diff_after))]))
      cat(paste(rep("=", 80), collapse=""), "\n")
      
      # Save SMD table to CSV
      smd_table_filename <- paste0("smd_comparison_table_", 
                                   gsub(" ", "_", cmp$name), "_", 
                                   gsub("[/ ]", "_", outc$label), ".csv")
      write.csv(smd_comparison_table, smd_table_filename, row.names = FALSE)
      cat("SMD comparison table saved to:", smd_table_filename, "\n\n")
      
      # 균형 평가 시각화 with SMD values displayed on plot
      balance_plot <- ggplot() +
        geom_point(
          data = ipwt_result$balance_before,
          aes(x = abs(std_diff), y = reorder(variable, abs(std_diff))),
          color = "red", size = 2, alpha = 0.7
        ) +
        geom_point(
          data = ipwt_result$balance_after,
          aes(x = abs(std_diff), y = reorder(variable, abs(std_diff))),
          color = "blue", size = 2, alpha = 0.7
        ) +
        # Add text labels showing SMD values
        geom_text(
          data = ipwt_result$balance_before,
          aes(x = abs(std_diff), y = reorder(variable, abs(std_diff)), 
              label = sprintf("%.3f", std_diff)),
          color = "darkred", size = 2.5, hjust = -0.1, vjust = -0.5
        ) +
        geom_text(
          data = ipwt_result$balance_after,
          aes(x = abs(std_diff), y = reorder(variable, abs(std_diff)), 
              label = sprintf("%.3f", std_diff)),
          color = "darkblue", size = 2.5, hjust = -0.1, vjust = 1.5
        ) +
        geom_vline(xintercept = 0.1, linetype = "dashed", color = "gray50") +
        labs(
          title = paste("Standardized Mean Differences:", cmp$name, "-", outc$label),
          subtitle = "Before (red) vs After (blue) weighting with SMD values displayed",
          x = "Absolute Standardized Mean Difference",
          y = "",
          caption = "Dashed line at |SMD| = 0.1 threshold"
        ) +
        theme_minimal() +
        theme(
          axis.text.y = element_text(size = 8),
          plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 10)
        ) +
        # Expand x-axis limits to accommodate text labels
        expand_limits(x = max(c(abs(ipwt_result$balance_before$std_diff), 
                               abs(ipwt_result$balance_after$std_diff))) * 1.3)
      
      print(balance_plot)
      
      # 가중치 적용 생존 곡선 플롯 (FI_ND_MASLD style)
      tryCatch({
        # 가중 생존 분석
        weighted_df <- ipwt_result$cohort
        
        # Debug: Check data structure
        cat("Debug: weighted_df columns:", names(weighted_df), "\n")
        cat("Debug: weighted_df rows:", nrow(weighted_df), "\n")
        
        # Convert event_time from days to months to match FI_ND_MASLD format
        weighted_df$event_time_months <- weighted_df$event_time / 30.44
        
        # Debug: Check converted time variable
        cat("Debug: event_time_months range:", range(weighted_df$event_time_months, na.rm = TRUE), "\n")
        
        # Fit the cumulative incidence model stratified by treatment with IPTW weights
        # (Create survival object directly in survfit call - matching FI_ND_MASLD approach)
        cum_incidence_fit <- survfit(Surv(time = event_time_months, event = event) ~ treatment, 
                                   data = weighted_df, 
                                   weights = weighted_df$ipw_std)
        
        # Debug: Check survival fit creation
        cat("Debug: cum_incidence_fit created with direct Surv() call\n")
        
        # Create comparison-specific treatment labels
        if (grepl("GLP1.*SGLT2", cmp$name)) {
          treatment_labels <- c("SGLT-2 inhibitors", "Semaglutide (GLP-1)")
        } else if (grepl("GLP1.*SU", cmp$name)) {
          treatment_labels <- c("SU/DPP-4i", "Semaglutide (GLP-1)")
        } else if (grepl("SGLT2.*SU", cmp$name)) {
          treatment_labels <- c("SU/DPP-4i", "SGLT-2 inhibitors")
        } else {
          treatment_labels <- c("Comparator", "Treatment")
        }
        
        
        # Create survival plot using exact FI_ND_MASLD format with confidence intervals and risk table
        surv_plot <- ggsurvplot(
          fit = cum_incidence_fit, 
          data = weighted_df,
          fun = "event", 
          conf.int = TRUE,
          risk.table = TRUE,
          risk.table.col = "strata",
          risk.table.fontsize = 3,
          tables.height = 0.3,
          legend.title = "Treatment Group", 
          legend.labs = treatment_labels,
          xlab = "Follow-up Time (months)", 
          ylab = "Cumulative Incidence", 
          title = paste("Cumulative Incidence of", outc$label, "Stratified by Treatment Group (IPTW-weighted)"), 
          break.time.by = 12,
          xlim = c(0, 60),
          lwd = 1.0,
          censor = FALSE,
          palette = "jama",
          ggtheme = theme_bw() +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank())
        )
        
        # Debug: Check if plot was created
        cat("Debug: surv_plot created successfully\n")
        
        print(surv_plot)
        
        # 그래프 저장 (simplified approach)
        tryCatch({
          plot_filename <- paste0("ipwt_survival_", 
                                gsub(" ", "_", cmp$name), "_", 
                                gsub("[/ ]", "_", outc$label), ".png")
          cat("Debug: Attempting to save plot as:", plot_filename, "\n")
          
          # Save the complete plot with risk table using png device
          png(filename = plot_filename, width = 12, height = 10, units = "in", res = 300)
          print(surv_plot)
          dev.off()
          cat("Debug: Plot saved successfully\n")
        }, error = function(save_error) {
          cat("Error saving plot:", save_error$message, "\n")
        })
        
      }, error = function(e) {
        cat("Error generating weighted survival plot:", e$message, "\n")
        cat("Skipping plot for this outcome.\n")
      })
      
      # 결과 테이블에 추가
      tryCatch({
        # Cox 모델에서 직접 HR과 CI 추출
        cox_model <- ipwt_result$cox$conf.int
        hr <- cox_model[1, "exp(coef)"]
        lower_ci <- cox_model[1, "lower .95"]
        upper_ci <- cox_model[1, "upper .95"]
        p_value <- ipwt_result$cox$coefficients[1, "Pr(>|z|)"]
        
        result_row <- data.frame(
          method = "IPWT",
          comparison = cmp$name,
          outcome = outc$label,
          effect_measure = "Hazard Ratio",
          estimate = hr,
          lower_ci = lower_ci,
          upper_ci = upper_ci,
          p_value = p_value,
          n_total = nrow(ipwt_result$cohort),
          n_events = sum(ipwt_result$cohort$event),
          stringsAsFactors = FALSE
        )
        
        results_table <- rbind(results_table, result_row)
        
        # 결과 저장 (각 분석별로)
        write.csv(
          result_row,
          paste0("ipwt_result_", 
                 gsub(" ", "_", cmp$name), "_", 
                 gsub("[/ ]", "_", outc$label), ".csv"),
          row.names = FALSE
        )
        
      }, error = function(e) {
        cat("Error creating results table:", e$message, "\n")
      })
      
      # ---- Generate Study Flowchart ----
      cat("\n=== Generating Study Flowchart ===\n")
      
      # Store flow data for this comparison-outcome
      patient_flow_tracker[[flow_key]] <- flow_data
      
      # Print flow summary
      cat("Study flow summary for", flow_key, ":\n")
      cat("- Initial cohort (T2DM):", flow_data$initial_dm, "\n")
      cat("- After pre-existing exclusion:", flow_data$after_preexisting_exclusion, 
          "(excluded:", flow_data$excluded_preexisting, ")\n")
      cat("- After follow-up requirements:", flow_data$after_followup_filter, 
          "(excluded:", flow_data$excluded_followup, ")\n")
      cat("- After age filtering:", flow_data$final_cohort, 
          "(excluded:", flow_data$excluded_age, ")\n")
      cat("- Final events:", flow_data$total_events, "\n")
      
      # Generate and save flowchart
      tryCatch({
        # Create flowchart
        flowchart <- create_study_flowchart(flow_data, cmp$name, outc$label)
        print(flowchart)
        
        # Save flowchart (requires rsvg package)
        if (requireNamespace("rsvg", quietly = TRUE)) {
          save_flowchart(flowchart, cmp$name, outc$label)
        } else {
          cat("rsvg package not available - flowchart displayed but not saved\n")
        }
        
        # Save flow data to CSV
        flow_df <- data.frame(
          comparison = cmp$name,
          outcome = outc$label,
          step = c("Initial T2DM", "After Pre-existing Exclusion", 
                   "After Follow-up Filter", "After Age Filter", "Final Cohort"),
          n_patients = c(flow_data$initial_dm, flow_data$after_preexisting_exclusion,
                        flow_data$after_followup_filter, 
                        flow_data$final_cohort, flow_data$final_cohort),
          n_events = c(NA, NA, NA, NA, flow_data$total_events),
          excluded = c(flow_data$excluded_non_dm, flow_data$excluded_preexisting,
                      flow_data$excluded_followup, flow_data$excluded_age, 0)
        )
        
        write.csv(flow_df, 
                  paste0("patient_flow_", 
                         gsub(" ", "_", cmp$name), "_", 
                         gsub("[/ ]", "_", outc$label), ".csv"),
                  row.names = FALSE)
        
      }, error = function(e) {
        cat("Error generating flowchart:", e$message, "\n")
      })
      
      # ---- HbA1c 중재 분석 (Mediation Analysis) - 모든 결과에서 수행 ----
      cat("\n--- HbA1c Mediation Analysis for", outc$label, "---\n")
      
      # Calculate weighted mean HbA1c using IPTW weights
      cat("\nCalculating 1-year weighted mean HbA1c for mediation analysis\n")
      cat("Total patients in IPWT cohort:", nrow(ipwt_result$cohort), "\n")
      
      # Get weighted HbA1c values
      med_df_temp <- get_hba1c_mediator(
        cohort_df = ipwt_result$cohort,
        panel_path = "a1c_panel.csv",
        win_days = 365,
        weight_var = "ipw_std",
        use_weighted_mean = TRUE
      )
      
      # Report HbA1c data availability
      cat("Patients with baseline HbA1c data:", sum(!is.na(med_df_temp$baseline_hba1c)), "\n")
      cat("Patients with weighted mean HbA1c data:", sum(!is.na(med_df_temp$hba1c_mean_6mo)), "\n")
      
      # Diagnostic: Compare baseline vs weighted mean HbA1c
      cat("\n--- HbA1c Distribution Comparison ---\n")
      cat("Baseline HbA1c - Mean:", round(mean(med_df_temp$baseline_hba1c, na.rm=TRUE), 2), 
          "SD:", round(sd(med_df_temp$baseline_hba1c, na.rm=TRUE), 2), "\n")
      cat("Weighted Mean HbA1c - Mean:", round(mean(med_df_temp$hba1c_mean_6mo, na.rm=TRUE), 2),
          "SD:", round(sd(med_df_temp$hba1c_mean_6mo, na.rm=TRUE), 2), "\n")
      cat("Correlation between baseline and weighted mean:", 
          round(cor(med_df_temp$baseline_hba1c, med_df_temp$hba1c_mean_6mo, use="complete.obs"), 3), "\n")
      
      # 수정된 landmark 기간으로 outcome 재정의
      add_outcome_post_mediator_modified <- function(df,
                                          outcome_var,
                                          data_cut,
                                          late_onset = FALSE,
                                          age_cut = 50,
                                          landmark_days = 90) { 
        df %>%
          mutate(outcome_date = as.Date(.data[[outcome_var]]),
                start_fu = index_date + landmark_days,  
                censor_date = pmin(as.Date(EHRmaxDT), data_cut, na.rm = TRUE),
                raw_event = !is.na(outcome_date) &
                            outcome_date > start_fu &
                            outcome_date <= censor_date,
                time = as.numeric(pmin(outcome_date, censor_date, na.rm = TRUE) - start_fu),
                event = as.integer(raw_event)) %>%
          { if (late_onset) {
              mutate(.,
                    age_at_event = age + time / 365.25,
                    event        = as.integer(event == 1 & age_at_event >= age_cut))
            } else .
          } %>%
          filter(time >= 0)
      }
      
      # 수정된 함수 적용 - IPWT cohort에서 추출한 변수 이용
      med_df_temp <- med_df_temp %>% 
                    add_outcome_post_mediator_modified(
                      outcome_var = outc$var,
                      data_cut = data_cut_date,
                      late_onset = outc$late_onset,
                      landmark_days = 90
                    )
      
      cat("\nDiagnostic after outcome redefinition but before HbA1c filtering:\n")
      cat("Total patients:", nrow(med_df_temp), "\n")
      cat("Total events:", sum(med_df_temp$event, na.rm=TRUE), "\n")
      
      # Filter to patients with baseline HbA1c - preserve IPWT weights
      med_df <- med_df_temp %>% 
                filter(!is.na(baseline_hba1c))
      
      cat("\nDiagnostic after HbA1c filtering:\n")
      cat("Remaining patients:", nrow(med_df), "\n")
      cat("Remaining events:", sum(med_df$event, na.rm=TRUE), "\n")
      cat("Events by treatment:", "\n")
      print(table(med_df$treatment, med_df$event))
      
      # 중재 분석 - 모든 결과에 대해 시도
      if(nrow(med_df) > 0 && sum(med_df$event, na.rm=TRUE) > 0 && length(unique(med_df$treatment)) >= 2) {
        base_covs <- c("age", "sex_cat", "raceethnicity_cat")
        
        cat("\nMediation analysis with IPWT approach:\n")
        
        # 1. 가중 선형 회귀 - 중재변수(HbA1c) 모델 구축
        # Use weighted mean HbA1c instead of baseline
        m_formula <- as.formula(paste("hba1c_mean_6mo ~ treatment +", paste(base_covs, collapse = " + ")))
        m_fit <- tryCatch({
          lm(m_formula, data = med_df, weights = ipw_std)
        }, error = function(e) {
          cat("Error in mediator model:", e$message, "\n")
          return(NULL)
        })
        
        if(is.null(m_fit)) {
          cat("Skipping mediation analysis due to error in mediator model.\n")
        } else {
          # HbA1c 모델 요약
          cat("\n중재자 모델 공식:", deparse(m_formula), "\n")
          cat("중재자 모델 성공! (가중치 적용됨)\n\n중재자 모델 요약:\n")
          m_summary <- summary(m_fit)
          print(m_summary$coefficients[c("(Intercept)", "treatment"), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")])
          
          # 가중 로지스틱 회귀 - 결과 모델 구축
          # Use weighted mean HbA1c instead of baseline
          y_formula <- as.formula(paste("event ~ treatment + hba1c_mean_6mo +", 
                          paste(base_covs, collapse = " + ")))
          
          # 이벤트가 너무 적을 경우 Firth's 로지스틱 회귀 시도
          y_fit <- tryCatch({
            glm(y_formula, family = binomial(), data = med_df, weights = ipw_std)
          }, error = function(e) {
            cat("Error in outcome model using glm:", e$message, "\n")
            cat("Trying with Firth's bias reduction...\n")
            
            # brglm 패키지가 로드되어 있는지 확인하고 없으면 설치
            if(!requireNamespace("brglm2", quietly = TRUE)) {
              cat("Installing brglm2 package for Firth's logistic regression...\n")
              install.packages("brglm2")
            }
            
            tryCatch({
              # Firth's method에는 weights 적용이 어려울 수 있음
              brglm2::brglm_fit(y_formula, family = binomial(), data = med_df)
            }, error = function(e2) {
              cat("Error in outcome model using Firth's method:", e2$message, "\n")
              return(NULL)
            })
          })
          
          if(is.null(y_fit)) {
            cat("Skipping mediation analysis due to error in outcome model.\n")
          } else {
            # 결과 모델 요약
            cat("\n결과 모델 공식:", deparse(y_formula), "\n")
            cat("결과 모델 성공! (가중치 적용됨)\n\n결과 모델 요약:\n")
            y_summary <- summary(y_fit)
            print(y_summary$coefficients[c("treatment", "hba1c_mean_6mo"), ])
            
            # 효과 추출
            a_path <- coef(m_fit)["treatment"]  # 치료 -> HbA1c
            b_path <- coef(y_fit)["hba1c_mean_6mo"]  # HbA1c -> 결과
            c_path <- coef(y_fit)["treatment"]  # 직접 효과
            
            # 간접 효과 계산
            indirect_effect <- a_path * b_path
            total_effect <- c_path + indirect_effect
            
            # 표준 오차 계산
            a_se <- m_summary$coefficients["treatment", "Std. Error"]
            b_se <- y_summary$coefficients["hba1c_mean_6mo", "Std. Error"]
            c_se <- y_summary$coefficients["treatment", "Std. Error"]
            
            # p-값 계산
            a_p <- m_summary$coefficients["treatment", "Pr(>|t|)"]
            
            # 표준 오차 및 p-값 열 이름이 다를 수 있음
            if("Pr(>|z|)" %in% colnames(y_summary$coefficients)) {
              b_p <- y_summary$coefficients["hba1c_mean_6mo", "Pr(>|z|)"]
              c_p <- y_summary$coefficients["treatment", "Pr(>|z|)"]
            } else if("Pr(>|t|)" %in% colnames(y_summary$coefficients)) {
              b_p <- y_summary$coefficients["hba1c_mean_6mo", "Pr(>|t|)"]
              c_p <- y_summary$coefficients["treatment", "Pr(>|t|)"]
            } else {
              # 마지막 열을 p-값으로 가정
              col_idx <- ncol(y_summary$coefficients)
              b_p <- y_summary$coefficients["hba1c_mean_6mo", col_idx]
              c_p <- y_summary$coefficients["treatment", col_idx]
            }
            
            # 소벨 테스트로 간접 효과의 p-값 계산
            sobel_z <- (a_path*b_path) / sqrt(b_path^2*a_se^2 + a_path^2*b_se^2)
            sobel_p <- 2 * (1 - pnorm(abs(sobel_z)))
            
            # 결과 출력
            cat("\n=== IPWT Mediation Analysis Summary for", outc$label, "===\n")
            cat("A path (Treatment -> HbA1c):", a_path, "p =", a_p, "\n")
            cat("B path (HbA1c -> Outcome):", b_path, "p =", b_p, "\n")
            cat("C path (Direct Effect):", c_path, "p =", c_p, "\n")
            cat("Indirect Effect (A*B):", indirect_effect, "p =", sobel_p, "(Sobel test)\n")
            cat("Total Effect (C + A*B):", total_effect, "\n")
            
            # Calculate standard error for indirect effect (Delta method)
            indirect_se <- sqrt(b_path^2 * a_se^2 + a_path^2 * b_se^2)
            
            # Calculate standard error for total effect
            # Assuming independence between direct and indirect effects
            total_se <- sqrt(c_se^2 + indirect_se^2)
            
            # Calculate hazard ratios and 95% CI (since outcome model is logistic, these are odds ratios)
            # For mediation effects, we report both Beta and OR
            hr_direct <- exp(c_path)
            hr_direct_lower <- exp(c_path - 1.96 * c_se)
            hr_direct_upper <- exp(c_path + 1.96 * c_se)
            
            hr_indirect <- exp(indirect_effect)
            hr_indirect_lower <- exp(indirect_effect - 1.96 * indirect_se)
            hr_indirect_upper <- exp(indirect_effect + 1.96 * indirect_se)
            
            hr_total <- exp(total_effect)
            hr_total_lower <- exp(total_effect - 1.96 * total_se)
            hr_total_upper <- exp(total_effect + 1.96 * total_se)
            
            # Calculate p-value for total effect
            z_total <- total_effect / total_se
            p_total <- 2 * (1 - pnorm(abs(z_total)))
            
            # 결과를 테이블에 저장 - 사용자가 요청한 형식으로
            # Calculate confidence intervals for A path
            a_lower <- a_path - 1.96 * a_se
            a_upper <- a_path + 1.96 * a_se
            
            med_result <- data.frame(
              method = "IPWT",
              comparison = cmp$name,
              outcome = outc$label,
              effect = c("Total Effect of Semaglutide for Prevention of Adult-Onset Epilepsy",
                        "Direct Effect of Semaglutide", 
                        "Indirect Effect mediated by HbA1c",
                        "A Path (Treatment → HbA1c)"),
              beta = sprintf("%.4f", c(total_effect, c_path, indirect_effect, a_path)),
              se = sprintf("%.4f", c(total_se, c_se, indirect_se, a_se)),
              beta_se = sprintf("%.4f (%.4f)", c(total_effect, c_path, indirect_effect, a_path), 
                                              c(total_se, c_se, indirect_se, a_se)),
              hr = c(sprintf("%.3f", c(hr_total, hr_direct, hr_indirect)), 
                     "N/A"),  # A path is not a hazard ratio
              hr_lower = c(sprintf("%.3f", c(hr_total_lower, hr_direct_lower, hr_indirect_lower)),
                          "N/A"),
              hr_upper = c(sprintf("%.3f", c(hr_total_upper, hr_direct_upper, hr_indirect_upper)),
                          "N/A"),
              hr_95ci = c(sprintf("%.3f (%.3f-%.3f)", 
                              c(hr_total, hr_direct, hr_indirect),
                              c(hr_total_lower, hr_direct_lower, hr_indirect_lower),
                              c(hr_total_upper, hr_direct_upper, hr_indirect_upper)),
                         sprintf("%.4f (%.4f-%.4f)", a_path, a_lower, a_upper)),  # Beta with 95% CI for A path
              p_value = sprintf("%.4f", c(p_total, c_p, sobel_p, a_p)),
              stringsAsFactors = FALSE
            )
            
            # 결과 리스트에 저장
            med_key <- paste(cmp$name, outc$label, sep = " - ")
            mediation_results_list[[med_key]] <- med_result
          }
        }
      } else {
        cat("\n[IPWT Mediation skipped for", outc$label, "] – No events remaining after HbA1c data filtering or insufficient treatment variation.\n")
        cat("This could be because:\n")
        cat("1. Few patients have HbA1c measurements\n")
        cat("2. Events occurred before the landmark period (90 days)\n")
        cat("3. Not enough events in both treatment groups\n")
      }
      cat("\n--- HbA1c Mediation Analysis End for", outc$label, "---\n")
      
      # ==========================================================
      # HbA1c trajectory analysis for all neurologic outcomes
      # ----------------------------------------------------------
      if (cmp$name %in% c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4") &&
          (outc$var %in% c("epilepsy_or_seizure_start_date", "adrd_start_date", "stroke_start_date"))) {
        
        cat("\n========= HbA1c Trajectory Analysis:", cmp$name, "-", outc$label, "=========\n")
        
        # HbA1c 데이터 확인 및 로드
        name_of_file_in_bucket <- "a1c_panel.csv"
        
        # 파일이 없으면 다운로드
        if(!file.exists(name_of_file_in_bucket)) {
          system(glue::glue("gsutil cp {my_bucket}/data/{name_of_file_in_bucket} ."),
                 intern = TRUE)
        }
        
        # HbA1c 궤적 분석 실행 - IPWT 결과 사용
        hba1c <- tryCatch({
          analyse_hba1c_trajectory(
            results        = ipwt_result,  # IPWT result object
            panel_path     = name_of_file_in_bucket,
            years_window   = 2,
            drop_crossovers = TRUE
          )
        }, error = function(e) {
          cat("Error in HbA1c trajectory analysis:", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(hba1c)) {
          # 분석 제목 설정
          analysis_title <- outc$label
          
          # 플롯 타이틀 업데이트
          hba1c$plot <- hba1c$plot + 
            labs(title = paste0("HbA1c Trajectories: ", cmp$name, " - ", analysis_title),
                 subtitle = "IPWT weighted cohort")
          
          print(hba1c$plot)
          
          # 파일 저장
          ggsave(paste0("ipwt_hba1c_trajectory_", gsub(" ", "_", cmp$name), "_", 
                        tolower(gsub("/", "_", analysis_title)), ".png"),
                 hba1c$plot, width = 8, height = 6)
          
          cat("\n==== HbA1c Trajectory LMM Results for", analysis_title, "(IPWT) ====\n")
          print(summary(hba1c$model)$coefficients[, c("Estimate", "Std. Error", "Pr(>|t|)")])
          
          # 상세 모델 결과
          model_results <- broom.mixed::tidy(hba1c$model, effects = "fixed", conf.int = TRUE)
          print(model_results)
          
          # 시간별 예측값
          cat("\nEstimated marginal means by time and treatment:\n")
          print(hba1c$emmeans)
          
          # 결과를 파일로 저장
          write.csv(model_results, 
                    paste0("ipwt_hba1c_lmm_", gsub(" ", "_", cmp$name), "_", 
                          tolower(gsub("/", "_", analysis_title)), ".csv"), 
                    row.names = FALSE)
          
          # 치료 효과와 상호작용 효과 요약
          cat("\nTreatment Effect (Baseline HbA1c difference):", 
              round(model_results$estimate[model_results$term == "treatment"], 3),
              "(p =", format.pval(model_results$p.value[model_results$term == "treatment"], digits = 3), ")\n")
          
          interaction_term <- model_results$term == "years_since:treatment"
          if(any(interaction_term)) {
            cat("Time×Treatment Interaction (HbA1c change over time):", 
                round(model_results$estimate[interaction_term], 3),
                "(p =", format.pval(model_results$p.value[interaction_term], digits = 3), ")\n")
          }
          
          # IPWT 가중치를 사용한 결과 비교를 위해 가중 평균 계산
          weighted_summary <- ipwt_result$cohort %>%
            group_by(treatment) %>%
            summarise(
              n = n(),
              weighted_n = sum(ipw_std),
              .groups = "drop"
            )
          
          cat("\nCohort sizes for HbA1c trajectory analysis:\n")
          print(weighted_summary)
        }
      }
    }
  }
  
  # 중재 분석 결과 요약 테이블 생성
  if (length(mediation_results_list) > 0) {
    mediation_summary_table <- do.call(rbind, mediation_results_list)
    cat("\n=== IPWT Mediation Analysis Summary Table ===\n")
    print(mediation_summary_table)
    
    # 중재 분석 결과 저장
    write.csv(mediation_summary_table, "ipwt_mediation_results.csv", row.names = FALSE)
    
    # 각 결과 유형별로 중재 분석 요약 출력
    for (outcome_type in unique(mediation_summary_table$outcome)) {
      outcome_mediation <- mediation_summary_table[mediation_summary_table$outcome == outcome_type, ]
      if (nrow(outcome_mediation) > 0) {
        cat("\n=== Mediation Results for", outcome_type, "===\n")
        # 경로별로 요약 - 새로운 형식의 출력
        print(outcome_mediation[, c("comparison", "effect", "beta_se", "hr_95ci", "p_value")])
      }
    }
  }
  
  # 전체 결과 저장
  write.csv(results_table, "ipwt_all_results.csv", row.names = FALSE)
  
  # ---- Generate Summary Flowchart for All Analyses ----
  cat("\n=== Generating Summary Patient Flow Tables ===\n")
  
  # Combine all flow data into a summary table
  if (length(patient_flow_tracker) > 0) {
    all_flow_data <- list()
    
    for (flow_name in names(patient_flow_tracker)) {
      flow_info <- patient_flow_tracker[[flow_name]]
      parts <- strsplit(flow_name, " - ")[[1]]
      comparison <- parts[1]
      outcome <- parts[2]
      
      # Check if all required fields exist and have values, otherwise set defaults
      safe_get <- function(field, default = 0) {
        if (is.null(flow_info[[field]]) || length(flow_info[[field]]) == 0) {
          return(default)
        } else {
          return(flow_info[[field]])
        }
      }
      
      initial_dm <- safe_get("initial_dm", 0)
      excluded_preexisting <- safe_get("excluded_preexisting", 0)
      after_preexisting <- safe_get("after_preexisting_exclusion", 0)
      excluded_followup <- safe_get("excluded_followup", 0)
      after_followup <- safe_get("after_followup_filter", 0)
      excluded_age <- safe_get("excluded_age", 0)
      final_cohort <- safe_get("final_cohort", 0)
      total_events <- safe_get("total_events", 0)
      
      # Calculate event rate safely
      event_rate <- if (final_cohort > 0) {
        round(total_events / final_cohort * 100, 2)
      } else {
        0
      }
      
      all_flow_data[[flow_name]] <- data.frame(
        comparison = comparison,
        outcome = outcome,
        initial_dm = initial_dm,
        excluded_preexisting = excluded_preexisting,
        after_preexisting = after_preexisting,
        excluded_followup = excluded_followup,
        after_followup = after_followup,
        excluded_age = excluded_age,
        final_cohort = final_cohort,
        total_events = total_events,
        event_rate = event_rate
      )
    }
    
    # Combine into one dataframe
    flow_summary <- do.call(rbind, all_flow_data)
    
    # Save summary table
    write.csv(flow_summary, "ipwt_all_patient_flows.csv", row.names = FALSE)
    
    # Print summary
    cat("\nPatient Flow Summary Table:\n")
    print(flow_summary)
  }
  
  # Combine all weighted characteristics tables
  tryCatch({
    weighted_char_files <- list.files(pattern = "^ipwt_weighted_characteristics_.*\\.csv$")
    if (length(weighted_char_files) > 0) {
      all_weighted_chars <- map_dfr(weighted_char_files, ~read_csv(.x, show_col_types = FALSE))
      write.csv(all_weighted_chars, "ipwt_all_weighted_characteristics.csv", row.names = FALSE)
      cat("\nSaved combined weighted characteristics table: ipwt_all_weighted_characteristics.csv\n")
    }
  }, error = function(e) {
    cat("Error combining weighted characteristics tables:", e$message, "\n")
  })
  
  return(list(results = results_table, mediation = mediation_results_list))
}

# IPWT 분석 실행 - 수정된 반환 구조
ipwt_output <- analyze_with_ipwt()
ipwt_table <- ipwt_output$results

# 결과 미리보기
cat("\n==== IPWT 분석 결과 요약 ====\n")
print(ipwt_table)

# Stroke 결과 확인
stroke_results <- ipwt_table %>% 
  filter(outcome == "Stroke") %>%
  arrange(comparison)

if(nrow(stroke_results) > 0) {
  cat("\n==== Stroke 분석 결과 ====\n")
  print(stroke_results)
  
  # Stroke 결과를 위한 별도 포레스트 플롯
  if(nrow(stroke_results) >= 2) {
    stroke_plot <- ggplot(stroke_results, 
                         aes(y = comparison, x = estimate)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
      geom_point(aes(size = n_events), color = "#ff7f0e") +
      geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci), height = 0.2, color = "#ff7f0e") +
      labs(
        title = "IPWT Analysis: Stroke Risk Across Comparisons",
        subtitle = paste("Total stroke events:", sum(stroke_results$n_events)),
        x = "Hazard Ratio",
        y = "",
        size = "Events"
      ) +
      coord_cartesian(xlim = c(0, 2)) +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    print(stroke_plot)
    ggsave("ipwt_stroke_results.png", stroke_plot, width = 9, height = 6)
  }
} else {
  cat("\n주의: Stroke 결과가 없습니다. outcomes 리스트에 Stroke가 포함되어 있는지 확인하세요.\n")
}


# 중재 분석 요약 시각화
if(length(ipwt_output$mediation) > 0) {
  # 모든 중재 분석 결과를 하나의 데이터프레임으로 병합
  all_mediations <- do.call(rbind, ipwt_output$mediation)
  
  # A, B 경로 효과 크기 비교 시각화
  paths_ab <- all_mediations %>%
    filter(effect %in% c("A Path (Treatment → HbA1c)", "Direct Effect of Semaglutide")) %>%
    mutate(
      beta_numeric = as.numeric(beta),
      p_value_numeric = as.numeric(p_value),
      significant = p_value_numeric < 0.05
    )
  
  if(nrow(paths_ab) > 0) {
    ab_plot <- ggplot(paths_ab, aes(x = outcome, y = beta_numeric, fill = effect, alpha = significant)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
      scale_fill_manual(values = c("A Path (Treatment → HbA1c)" = "#1f77b4", 
                                 "Direct Effect of Semaglutide" = "#ff7f0e")) +
      facet_wrap(~comparison) +
      labs(
        title = "IPWT Mediation Analysis: A & B Path Effects",
        subtitle = "Treatment -> HbA1c (A) and HbA1c -> Outcome (B) paths",
        x = "",
        y = "Effect Estimate",
        fill = "Path",
        alpha = "p < 0.05"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    print(ab_plot)
    ggsave("ipwt_mediation_ab_paths.png", ab_plot, width = 10, height = 7)
  }
}