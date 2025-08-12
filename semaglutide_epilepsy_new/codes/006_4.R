# 006_4
# ---- TMLE 분석 ----
late_cutoff <- 50

# Modified run_tmle_analysis_fixed to mirror IPTW covariate selection
run_tmle_analysis_fixed <- function(df, outcome_var, exclude_vars = NULL) {
  # 1. Define exactly the same base covariates as in IPTW:
  base_ps_vars    <- all_ps_vars
  additional_vars <- c("income", "education", "baseline_bmi", "baseline_hba1c", "index_year")
  
  # 2. Subset covariates exactly as IPTW did:
  requested_covs <- setdiff(c(base_ps_vars, additional_vars), exclude_vars)
  requested_covs <- requested_covs[requested_covs %in% colnames(df)]
  
  # 3. Filter out any covariate that has only a single level:
  keep_covs <- requested_covs[
    sapply(df[requested_covs], function(x) length(unique(na.omit(x))) > 1)
  ]
  
  # 4. Build a cleaned data frame for TMLE:
  clean_df <- df %>%
    dplyr::select(
      treatment,
      !!sym(outcome_var),
      all_of(keep_covs)
    ) %>%
    na.omit()
  
  # 5. Prepare W_vars; if empty, use intercept-only
  W_vars <- clean_df %>% dplyr::select(-treatment, -!!sym(outcome_var))
  if (ncol(W_vars) == 0) {
    W_vars <- data.frame(intercept = rep(1, nrow(clean_df)))
  }
  
  # 6. Run TMLE with SL.glm for both Q and g
  set.seed(12345)
  tmle_fit <- tmle(
    Y            = clean_df[[outcome_var]],
    A            = clean_df$treatment,
    W            = W_vars,
    family       = "binomial",
    Q.SL.library = c("SL.glm"),
    g.SL.library = c("SL.glm")
  )
  
  # 7. Compute NNT from the risk difference
  risk_diff <- tmle_fit$estimates$ATE$psi
  nnt <- ifelse(risk_diff != 0, 1 / abs(risk_diff), Inf)
  
  # 8. Return results
  list(
    tmle_result = tmle_fit,
    cohort      = df,
    nnt         = nnt,
    data        = clean_df
  )
}

analyze_with_tmle <- function() {
  results_table <- data.frame(
    method       = character(),
    comparison   = character(),
    outcome      = character(),
    effect_measure = character(),
    estimate     = numeric(),
    lower_ci     = numeric(),
    upper_ci     = numeric(),
    p_value      = numeric(),
    n_total      = integer(),
    n_events     = integer(),
    stringsAsFactors = FALSE
  )
  
  mediation_results_list <- list()
  
  for (cmp in comparisons) {
    # Get drug concept IDs
    exposure_concepts   <- get_drug_concepts(drug_classes[[cmp$exposure]])
    comparator_concepts <- get_drug_concepts(drug_classes[[cmp$comparator]])
    
    # Build index dates
    idx_exposure   <- build_exposure_idx(exposure_concepts, paste0(tolower(cmp$exposure), "_idx"))
    idx_comparator <- build_exposure_idx(comparator_concepts, paste0(tolower(cmp$comparator), "_idx"))
    
    # Construct cohort
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
    ) %>%
      filter(index_date >= cutoff_date)
    
    attr(cohort_pair, "exposure_idx")   <- paste0(tolower(cmp$exposure), "_idx")
    attr(cohort_pair, "comparator_idx") <- paste0(tolower(cmp$comparator), "_idx")
    
    for (outc in outcomes) {
      cat("\n========= TMLE Analysis:", cmp$name, "-", outc$label, "=========\n")
      
      # Prepare analytic dataframe
      analytic_df_pre_hba1c <- cohort_pair %>%
        left_join(merged_df_2, by = "person_id") %>%
        filter(dm == 1) %>%
        mutate(
          # Add calendar year of index date as a covariate
          index_year = factor(year(as.Date(index_date))),
          
          dob        = as.Date(substr(date_of_birth, 1, 10), format = "%Y/%m/%d"),
          age_at_end = as.numeric(difftime(data_cut_date, dob, units = "days")) / 365.25
        )
      
      # Add baseline HbA1c data - similar to IPTW section
      if (!exists("hba1c_panel_df_loaded")) {
        hba1c_panel_file <- "a1c_panel.csv"
        if (!file.exists(hba1c_panel_file)) {
          system(str_glue("gsutil cp gs://{WORKSPACE_BUCKET}/data/{hba1c_panel_file} ."), intern=T)
        }
        hba1c_panel_df_loaded <- read_csv(hba1c_panel_file, col_types = cols())
      }
      
      baseline_hba1c_data_for_tmle <- analytic_df_pre_hba1c %>%
        dplyr::select(person_id, index_date) %>%
        distinct() %>%
        left_join(hba1c_panel_df_loaded, by = "person_id", relationship = "many-to-many") %>%
        filter(!is.na(date_of_measurement) & !is.na(index_date)) %>%
        mutate(
          days_from_index = abs(as.numeric(difftime(date_of_measurement, index_date, units = "days")))
        ) %>%
        group_by(person_id) %>%
        arrange(days_from_index, desc(date_of_measurement)) %>%
        slice(1) %>%
        ungroup() %>%
        dplyr::select(person_id, baseline_hba1c = A1c)
      
      analytic_df <- analytic_df_pre_hba1c %>%
        left_join(baseline_hba1c_data_for_tmle, by = "person_id") %>%
        # EXCLUDE PATIENTS WITH PRE-EXISTING CONDITIONS AT BASELINE FOR NEW ONSET ANALYSIS
        {
          before_exclusion <- nrow(.)
          cat("\nBefore excluding pre-existing", outc$label, "patients:", before_exclusion, "\n")
          
          # DEBUG: Check for pre-existing conditions using actual dates
          cat("\n=== DEBUGGING DATE-BASED EXCLUSIONS (TMLE) ===\n")
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
          
          result
        } %>%
        followup_and_event(outcome_var = outc$var, data_cut_date = data_cut_date, late_onset = outc$late_onset, early_onset = outc$early_onset)
      
      total_sample_size <- nrow(analytic_df)
      total_events      <- sum(analytic_df$event, na.rm = TRUE)
      cat("\nTotal sample size before TMLE:", total_sample_size, "\n")
      cat("Total events:", total_events, "\n")
      
      # 인덱스 연도 분포 확인
      cat("\nIndex Year Distribution:\n")
      print(table(analytic_df$index_year, useNA = "ifany"))
      
      # Apply late-onset age filter if needed
      if (isTRUE(outc$late_onset)) {
        analytic_df <- analytic_df %>% filter(age_at_end >= late_cutoff)
        cat("After age filtering (≥", late_cutoff, "): samples =", nrow(analytic_df),
            ", events =", sum(analytic_df$event, na.rm = TRUE), "\n")
      }
      
      # Build exclude_vars exactly as in IPTW
      split_comparator <- unlist(strsplit(cmp$comparator, "_"))
      exclude_vars <- unique(c(cmp$exposure, split_comparator))
      if (outc$label == "ADRD") {
        exclude_vars <- unique(c(exclude_vars, "dem"))
        cat("Excluding 'dem' covariate for ADRD outcome\n")
      }
      if (outc$label == "stroke" || outc$label == "Stroke") {
        exclude_vars <- unique(c(exclude_vars, "cvd"))
        cat("Excluding 'cvd' covariate for Stroke outcome\n")
      }
      
      # Run TMLE
      set.seed(12345)
      tmle_result <- tryCatch({
        run_tmle_analysis_fixed(analytic_df, "event", exclude_vars)
      }, error = function(e) {
        cat("Error in TMLE analysis:", e$message, "\nTrying simplified model...\n")
        
        # Try a simplified model with only basic covariates
        simple_df <- analytic_df %>%
          dplyr::select(person_id, treatment, event, event_time, index_date, age, sex_cat, raceethnicity_cat)
        run_tmle_analysis_fixed(simple_df, "event", NULL)
      })
      
      if (is.null(tmle_result)) {
        cat("Skipping this outcome due to TMLE error.\n")
        next
      }
      
      key <- paste0(cmp$name, " - ", outc$label)
      tmle_results[[key]] <- tmle_result
      
      # Save cohort data using All of Us method
      # Replace df with THE NAME OF YOUR DATAFRAME
      my_dataframe <- tmle_result$cohort
      
      # Replace 'test.csv' with THE NAME of the file you're going to store in the bucket
      destination_filename <- paste0("tmle_cohort_", 
                                   gsub(" ", "_", cmp$name), "_", 
                                   gsub("[/ ]", "_", outc$label), ".csv")
      
      # store the dataframe in current workspace
      write_excel_csv(my_dataframe, destination_filename)
      
      # Get the bucket name (already defined earlier)
      # my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
      
      # Copy the file from current workspace to the bucket
      system(paste0("gsutil cp ./", destination_filename, " ", my_bucket, "/data/"), intern=T)
      
      cat("Saved cohort data to bucket:", destination_filename, "\n")
      
      # Print TMLE fit object
      cat("\nTMLE Results:\n")
      print(tmle_result$tmle_result)
      
      # Extract ATE (risk difference)
      ate     <- tmle_result$tmle_result$estimates$ATE
      ate_psi <- round(ate$psi, 4)
      ate_ci_lower <- round(ate$CI[1], 4)
      ate_ci_upper <- round(ate$CI[2], 4)
      ate_pval <- format.pval(ate$pvalue, digits = 3)
      
      cat("\nAverage Treatment Effect (Risk Difference):", ate_psi,
          "(95% CI:", ate_ci_lower, "to", ate_ci_upper, ")",
          "p =", ate_pval, "\n")
      
      # Extract OR
      or     <- tmle_result$tmle_result$estimates$OR
      or_psi <- round(or$psi, 2)
      or_ci_lower <- round(or$CI[1], 2)
      or_ci_upper <- round(or$CI[2], 2)
      or_pval <- or$pvalue
      
      cat("\nOdds Ratio:", or_psi,
          "(95% CI:", or_ci_lower, "to", or_ci_upper, ")",
          "p =", format.pval(or_pval, digits = 3), "\n")
      
      # Print NNT if finite
      if (!is.na(tmle_result$nnt) && is.finite(tmle_result$nnt)) {
        cat("\nNumber Needed to Treat (NNT):", round(tmle_result$nnt, 1), "\n")
      }
      
      # Observed event rates
      event_rates <- analytic_df %>%
        group_by(treatment) %>%
        summarise(
          n      = n(),
          events = sum(event),
          rate   = events / n
        )
      cat("\nObserved Event Rates:\n")
      print(event_rates)
      
      # Append to results table: Risk Difference
      results_table <- rbind(
        results_table,
        data.frame(
          method         = "TMLE",
          comparison     = cmp$name,
          outcome        = outc$label,
          effect_measure = "Risk Difference",
          estimate       = ate_psi,
          lower_ci       = ate_ci_lower,
          upper_ci       = ate_ci_upper,
          p_value        = ate$pvalue,
          n_total        = total_sample_size,
          n_events       = total_events,
          stringsAsFactors = FALSE
        )
      )
      
      # Append Odds Ratio row
      results_table <- rbind(
        results_table,
        data.frame(
          method         = "TMLE",
          comparison     = cmp$name,
          outcome        = outc$label,
          effect_measure = "Odds Ratio",
          estimate       = or_psi,
          lower_ci       = or_ci_lower,
          upper_ci       = or_ci_upper,
          p_value        = or_pval,
          n_total        = total_sample_size,
          n_events       = total_events,
          stringsAsFactors = FALSE
        )
      )
      
      # Save TMLE result to CSV
      write.csv(
        results_table %>%
          filter(comparison == cmp$name, outcome == outc$label),
        paste0("tmle_result_",
               gsub(" ", "_", cmp$name), "_",
               gsub("[/ ]", "_", outc$label), ".csv"),
        row.names = FALSE
      )
      
      # Create a bar plot of predicted event probabilities
      prob_df <- data.frame(
        Group       = c("Control", "Semaglutide"),
        Probability = c(
          tmle_result$tmle_result$estimates$EY0$psi,
          tmle_result$tmle_result$estimates$EY1$psi
        ),
        Lower = c(
          tmle_result$tmle_result$estimates$EY0$CI[1],
          tmle_result$tmle_result$estimates$EY1$CI[1]
        ),
        Upper = c(
          tmle_result$tmle_result$estimates$EY0$CI[2],
          tmle_result$tmle_result$estimates$EY1$CI[2]
        )
      )
      
      tmle_plot <- ggplot(prob_df, aes(x = Group, y = Probability, fill = Group)) +
        geom_bar(stat = "identity", alpha = 0.7) +
        geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2) +
        scale_fill_manual(values = c("Control" = "#1f77b4", "Semaglutide" = "#ff7f0e")) +
        labs(
          title    = paste("TMLE Estimated Event Probabilities for", outc$label),
          subtitle = paste(cmp$name, "- Risk Difference =", ate_psi, "(p =", ate_pval, ")"),
          y        = "Probability of Event",
          x        = ""
        ) +
        theme_minimal() +
        theme(legend.position = "none")
      
      print(tmle_plot)
      
      ggsave(
        filename = paste0("tmle_probs_",
                          gsub(" ", "_", cmp$name), "_",
                          gsub("[/ ]", "_", outc$label), ".png"),
        plot     = tmle_plot,
        width    = 8,
        height   = 6
      )
      
      # ---- HbA1c Mediation Analysis for TMLE ----
      cat("\n--- HbA1c Mediation Analysis for", outc$label, "---\n")
      
      # Calculate weighted mean HbA1c for TMLE mediation analysis
      cat("\nCalculating 1-year weighted mean HbA1c for mediation analysis\n")
      cat("Total patients in TMLE cohort:", nrow(analytic_df), "\n")
      
      # Get weighted HbA1c values (no weights for TMLE since it uses different approach)
      med_df_temp <- get_hba1c_mediator(
        cohort_df = analytic_df,
        panel_path = "a1c_panel.csv",
        win_days = 365,
        weight_var = NULL,  # No weights needed for TMLE
        use_weighted_mean = FALSE  # Use unweighted for TMLE
      )
      
      # Report HbA1c data availability
      cat("Patients with baseline HbA1c data:", sum(!is.na(med_df_temp$baseline_hba1c)), "\n")
      cat("Patients with mean HbA1c data:", sum(!is.na(med_df_temp$hba1c_mean_6mo)), "\n")
      
      # Diagnostic: Compare baseline vs mean HbA1c
      cat("\n--- HbA1c Distribution Comparison (TMLE) ---\n")
      cat("Baseline HbA1c - Mean:", round(mean(med_df_temp$baseline_hba1c, na.rm=TRUE), 2), 
          "SD:", round(sd(med_df_temp$baseline_hba1c, na.rm=TRUE), 2), "\n")
      cat("1-Year Mean HbA1c - Mean:", round(mean(med_df_temp$hba1c_mean_6mo, na.rm=TRUE), 2),
          "SD:", round(sd(med_df_temp$hba1c_mean_6mo, na.rm=TRUE), 2), "\n")
      cat("Correlation between baseline and 1-year mean:", 
          round(cor(med_df_temp$baseline_hba1c, med_df_temp$hba1c_mean_6mo, use="complete.obs"), 3), "\n")
      
      # Redefine outcome with 90-day landmark
      add_outcome_post_mediator_modified <- function(df,
                                                     outcome_var,
                                                     data_cut,
                                                     late_onset = FALSE,
                                                     age_cut = 50,
                                                     landmark_days = 90) {
        df %>%
          mutate(
            outcome_date = as.Date(.data[[outcome_var]]),
            start_fu     = index_date + landmark_days,
            censor_date  = pmin(as.Date(EHRmaxDT), data_cut, na.rm = TRUE),
            raw_event    = !is.na(outcome_date) &
                           outcome_date > start_fu &
                           outcome_date <= censor_date,
            time         = as.numeric(pmin(outcome_date, censor_date, na.rm = TRUE) - start_fu),
            event        = as.integer(raw_event)
          ) %>%
          { if (late_onset) {
              mutate(.,
                     age_at_event = age + time / 365.25,
                     event        = as.integer(event == 1 & age_at_event >= age_cut))
            } else .
          } %>%
          filter(time >= 0)
      }
      
      med_df_temp <- med_df_temp %>%
        add_outcome_post_mediator_modified(
          outcome_var   = outc$var,
          data_cut      = data_cut_date,
          late_onset    = outc$late_onset,
          landmark_days = 90
        )
      
      cat("\nDiagnostic after outcome redefinition:\n")
      cat("Total patients:", nrow(med_df_temp), "\n")
      cat("Total events:", sum(med_df_temp$event, na.rm = TRUE), "\n")
      
      # Filter to those with baseline HbA1c measured
      med_df <- med_df_temp %>% filter(!is.na(baseline_hba1c))
      
      cat("\nDiagnostic after HbA1c filtering:\n")
      cat("Remaining patients:", nrow(med_df), "\n")
      cat("Remaining events:", sum(med_df$event, na.rm = TRUE), "\n")
      cat("Events by treatment:\n")
      print(table(med_df$treatment, med_df$event))
      
      if (nrow(med_df) > 0 && sum(med_df$event, na.rm = TRUE) > 0 && length(unique(med_df$treatment)) >= 2) {
        base_covs <- c("age", "sex_cat", "raceethnicity_cat")
        cat("\nMediation with TMLE approach:\n")
        
        # Mediator model
        # Use mean HbA1c instead of baseline
        m_formula <- as.formula(paste("hba1c_mean_6mo ~ treatment +", paste(base_covs, collapse = " + ")))
        m_fit <- tryCatch({
          lm(m_formula, data = med_df)
        }, error = function(e) {
          cat("Error in mediator model:", e$message, "\n")
          return(NULL)
        })
        
        if (!is.null(m_fit)) {
          cat("\nMediator model formula:", deparse(m_formula), "\n")
          cat("Mediator model success!\n\nMediator summary:\n")
          m_summary <- summary(m_fit)
          print(m_summary$coefficients[c("(Intercept)", "treatment"), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")])
          
          # Outcome model (logistic)
          # Use mean HbA1c instead of baseline
          y_formula <- as.formula(paste("event ~ treatment + hba1c_mean_6mo +", paste(base_covs, collapse = " + ")))
          y_fit <- tryCatch({
            glm(y_formula, family = binomial(), data = med_df)
          }, error = function(e) {
            cat("Error in outcome model using glm:", e$message, "\nTrying Firth's method...\n")
            if (!requireNamespace("brglm2", quietly = TRUE)) {
              install.packages("brglm2")
            }
            tryCatch({
              brglm2::brglm_fit(y_formula, family = binomial(), data = med_df)
            }, error = function(e2) {
              cat("Error in Firth's method:", e2$message, "\n")
              return(NULL)
            })
          })
          
          if (!is.null(y_fit)) {
            cat("\nOutcome model formula:", deparse(y_formula), "\n")
            cat("Outcome model success!\n\nOutcome summary:\n")
            y_summary <- summary(y_fit)
            print(y_summary$coefficients[c("treatment", "hba1c_mean_6mo"), ])
            
            # Extract paths
            a_path <- coef(m_fit)["treatment"]
            b_path <- coef(y_fit)["hba1c_mean_6mo"]
            c_path <- coef(y_fit)["treatment"]
            indirect_effect <- a_path * b_path
            total_effect <- c_path + indirect_effect
            
            a_se <- m_summary$coefficients["treatment", "Std. Error"]
            b_se <- y_summary$coefficients["hba1c_mean_6mo", "Std. Error"]
            c_se <- y_summary$coefficients["treatment", "Std. Error"]
            a_p  <- m_summary$coefficients["treatment", "Pr(>|t|)"]
            
            if ("Pr(>|z|)" %in% colnames(y_summary$coefficients)) {
              b_p <- y_summary$coefficients["hba1c_mean_6mo", "Pr(>|z|)"]
              c_p <- y_summary$coefficients["treatment", "Pr(>|z|)"]
            } else if ("Pr(>|t|)" %in% colnames(y_summary$coefficients)) {
              b_p <- y_summary$coefficients["hba1c_mean_6mo", "Pr(>|t|)"]
              c_p <- y_summary$coefficients["treatment", "Pr(>|t|)"]
            } else {
              col_idx <- ncol(y_summary$coefficients)
              b_p <- y_summary$coefficients["hba1c_mean_6mo", col_idx]
              c_p <- y_summary$coefficients["treatment", col_idx]
            }
            
            sobel_z <- (a_path * b_path) / sqrt(b_path^2 * a_se^2 + a_path^2 * b_se^2)
            sobel_p <- 2 * (1 - pnorm(abs(sobel_z)))
            
            cat("\n=== TMLE Mediation Analysis Summary for", outc$label, "===\n")
            cat("A path (Treatment → HbA1c):", a_path, "p =", a_p, "\n")
            cat("B path (HbA1c → Outcome):", b_path, "p =", b_p, "\n")
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
              method = "TMLE",
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
            
            med_key <- paste(cmp$name, outc$label, sep = " - ")
            mediation_results_list[[med_key]] <- med_result
            
            write.csv(
              med_result,
              paste0("tmle_mediation_",
                     gsub(" ", "_", cmp$name), "_",
                     gsub("[/ ]", "_", outc$label), ".csv"),
              row.names = FALSE
            )
          } else {
            cat("Skipping mediation: outcome model failed.\n")
          }
        } else {
          cat("Skipping mediation: mediator model failed.\n")
        }
      } else {
        cat("\n[TMLE Mediation skipped for", outc$label, "] – insufficient data or no events after filtering.\n")
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
        
        # HbA1c 궤적 분석 실행 - TMLE 결과 사용
        hba1c <- tryCatch({
          analyse_hba1c_trajectory(
            results        = tmle_result,  # TMLE result object
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
                 subtitle = "TMLE adjusted cohort")
          
          print(hba1c$plot)
          
          # 파일 저장
          ggsave(paste0("tmle_hba1c_trajectory_", gsub(" ", "_", cmp$name), "_", 
                        tolower(gsub("/", "_", analysis_title)), ".png"),
                 hba1c$plot, width = 8, height = 6)
          
          cat("\n==== HbA1c Trajectory LMM Results for", analysis_title, "(TMLE) ====\n")
          print(summary(hba1c$model)$coefficients[, c("Estimate", "Std. Error", "Pr(>|t|)")])
          
          # 상세 모델 결과
          model_results <- broom.mixed::tidy(hba1c$model, effects = "fixed", conf.int = TRUE)
          print(model_results)
          
          # 시간별 예측값
          cat("\nEstimated marginal means by time and treatment:\n")
          print(hba1c$emmeans)
          
          # 결과를 파일로 저장
          write.csv(model_results, 
                    paste0("tmle_hba1c_lmm_", gsub(" ", "_", cmp$name), "_", 
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
          
          # TMLE 코호트 크기 정보
          cohort_summary <- tmle_result$cohort %>%
            group_by(treatment) %>%
            summarise(
              n = n(),
              events = sum(event),
              .groups = "drop"
            )
          
          cat("\nCohort sizes for HbA1c trajectory analysis:\n")
          print(cohort_summary)
        }
      }
      # End of HbA1c trajectory analysis section
    }
    # End of outcomes loop
  }
  # End of comparisons loop
  
  # Summarize mediation results if any
  if (length(mediation_results_list) > 0) {
    mediation_summary_table <- do.call(rbind, mediation_results_list)
    cat("\n=== TMLE Mediation Analysis Summary Table ===\n")
    print(mediation_summary_table)
    write.csv(mediation_summary_table, "tmle_all_mediation_results.csv", row.names = FALSE)
    
    for (outcome_type in unique(mediation_summary_table$outcome)) {
      outcome_mediation <- mediation_summary_table[mediation_summary_table$outcome == outcome_type, ]
      if (nrow(outcome_mediation) > 0) {
        cat("\n=== Mediation Results for", outcome_type, "===\n")
        # 경로별로 요약 - 새로운 형식의 출력
        print(outcome_mediation[, c("comparison", "effect", "beta_se", "hr_95ci", "p_value")])
      }
    }
  }
  
  # Save all TMLE results
  write.csv(results_table, "tmle_all_results.csv", row.names = FALSE)
  list(results = results_table, mediation = mediation_results_list)
}
# End of analyze_with_tmle function

# Execute TMLE analysis
tmle_output <- analyze_with_tmle()
tmle_table  <- tmle_output$results

# Preview TMLE results
cat("\n==== TMLE 분석 결과 요약 ====\n")
print(tmle_table)

# Note: Cohort data saving already happens within the analyze_with_tmle() function
# for each comparison and outcome combination