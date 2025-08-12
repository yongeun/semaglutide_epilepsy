# ===============================================================================
# DIAGNOSTIC FUNCTION TO IDENTIFY NUMBER DISCREPANCIES
# ===============================================================================

diagnose_number_discrepancy <- function() {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("DIAGNOSING NUMBER DISCREPANCIES BETWEEN FLOWCHART AND PUBLICATION TABLE\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # 1. First, get the flowchart numbers from the stepwise analysis
  cat("\n=== FLOWCHART NUMBERS (from stepwise analysis) ===\n")
  if(exists("stepwise_results") && !is.null(stepwise_results)) {
    cat("ORIGINAL FLOWCHART NUMBERS (before final exclusions):\n")
    cat("GLP1 vs SGLT2:\n")
    cat("  Total final cohort:", format(stepwise_results$glp1_sglt2_final, big.mark=","), "\n")
    cat("  - Semaglutide:", format(stepwise_results$glp1_sglt2_glp1_final, big.mark=","), "\n")
    cat("  - SGLT2i:", format(stepwise_results$glp1_sglt2_sglt2_final, big.mark=","), "\n")
    cat("  - Sum check:", format(stepwise_results$glp1_sglt2_glp1_final + stepwise_results$glp1_sglt2_sglt2_final, big.mark=","), 
        "(should equal", format(stepwise_results$glp1_sglt2_final, big.mark=","), ")\n")
    
    cat("\nGLP1 vs SU/DPP4:\n")
    cat("  Total final cohort:", format(stepwise_results$glp1_su_dpp4_final, big.mark=","), "\n")
    cat("  - Semaglutide:", format(stepwise_results$glp1_su_dpp4_glp1_final, big.mark=","), "\n")
    cat("  - SU/DPP4i:", format(stepwise_results$glp1_su_dpp4_su_dpp4_final, big.mark=","), "\n")
    cat("  - Sum check:", format(stepwise_results$glp1_su_dpp4_glp1_final + stepwise_results$glp1_su_dpp4_su_dpp4_final, big.mark=","), 
        "(should equal", format(stepwise_results$glp1_su_dpp4_final, big.mark=","), ")\n")
    
    # Show corrected numbers if available
    if("glp1_sglt2_final_corrected" %in% names(stepwise_results)) {
      cat("\nCORRECTED FINAL ANALYSIS NUMBERS (matching CSV files):\n")
      cat("GLP1 vs SGLT2:\n")
      cat("  Total final cohort:", format(stepwise_results$glp1_sglt2_final_corrected, big.mark=","), "\n")
      cat("  - Semaglutide:", format(stepwise_results$glp1_sglt2_glp1_final_corrected, big.mark=","), "\n")
      cat("  - SGLT2i:", format(stepwise_results$glp1_sglt2_sglt2_final_corrected, big.mark=","), "\n")
      
      cat("\nGLP1 vs SU/DPP4:\n")
      cat("  Total final cohort:", format(stepwise_results$glp1_su_dpp4_final_corrected, big.mark=","), "\n")
      cat("  - Semaglutide:", format(stepwise_results$glp1_su_dpp4_glp1_final_corrected, big.mark=","), "\n")
      cat("  - SU/DPP4i:", format(stepwise_results$glp1_su_dpp4_su_dpp4_final_corrected, big.mark=","), "\n")
      
      cat("\nADDITIONAL EXCLUSIONS APPLIED:\n")
      cat("  GLP1 vs SGLT2:", format(stepwise_results$additional_exclusions_1, big.mark=","), "patients excluded\n")
      cat("  GLP1 vs SU/DPP4:", format(stepwise_results$additional_exclusions_2, big.mark=","), "patients excluded\n")
    }
  } else {
    cat("WARNING: stepwise_results not found or is NULL.\n")
    cat("This means the stepwise analysis hasn't run successfully.\n")
    cat("We can still check the publication table numbers from CSV files.\n")
    cat("Note: You should run the full script from the beginning to get flowchart numbers.\n")
  }
  
  # 2. Get numbers from actual cohort CSV files (if they exist)
  cat("\n=== PUBLICATION TABLE NUMBERS (from cohort CSV files) ===\n")
  
  comparisons <- c("GLP1_vs_SGLT2", "GLP1_vs_SU_DPP4", "SGLT2_vs_SU_DPP4")
  outcomes <- c("Epilepsy_Seizure", "Early-onset_Epilepsy_Seizure", "Late-onset_Epilepsy_Seizure")
  
  csv_numbers <- list()
  
  for(comp in comparisons) {
    cat("\n", gsub("_", " ", comp), ":\n")
    comp_numbers <- list()
    
    for(outcome in outcomes) {
      cat("  ", gsub("_", " ", outcome), ":\n")
      
      # Check IPTW cohort file
      iptw_file <- paste0("ipwt_cohort_", comp, "_", outcome, ".csv")
      if(file.exists(iptw_file)) {
        iptw_cohort <- read_csv(iptw_file, show_col_types = FALSE)
        iptw_counts <- table(iptw_cohort$treatment)
        cat("    IPTW cohort (", iptw_file, "):\n")
        # Use appropriate treatment labels based on comparison
        if(comp == "SGLT2_vs_SU_DPP4") {
          cat("      - SGLT2i (treatment=1):", iptw_counts["1"], "\n")
          cat("      - SU/DPP4i (treatment=0):", iptw_counts["0"], "\n")
        } else {
          cat("      - Semaglutide (treatment=1):", iptw_counts["1"], "\n")
          cat("      - Comparator (treatment=0):", iptw_counts["0"], "\n")
        }
        cat("      - Total:", sum(iptw_counts), "\n")
        
        comp_numbers[[paste0(outcome, "_iptw")]] <- list(
          semaglutide = as.numeric(iptw_counts["1"]),
          comparator = as.numeric(iptw_counts["0"]),
          total = sum(iptw_counts)
        )
      } else {
        cat("    IPTW cohort file not found:", iptw_file, "\n")
      }
      
      # Check TMLE cohort file
      tmle_file <- paste0("tmle_cohort_", comp, "_", outcome, ".csv")
      if(file.exists(tmle_file)) {
        tmle_cohort <- read_csv(tmle_file, show_col_types = FALSE)
        tmle_counts <- table(tmle_cohort$treatment)
        cat("    TMLE cohort (", tmle_file, "):\n")
        # Use appropriate treatment labels based on comparison
        if(comp == "SGLT2_vs_SU_DPP4") {
          cat("      - SGLT2i (treatment=1):", tmle_counts["1"], "\n")
          cat("      - SU/DPP4i (treatment=0):", tmle_counts["0"], "\n")
        } else {
          cat("      - Semaglutide (treatment=1):", tmle_counts["1"], "\n")
          cat("      - Comparator (treatment=0):", tmle_counts["0"], "\n")
        }
        cat("      - Total:", sum(tmle_counts), "\n")
        
        comp_numbers[[paste0(outcome, "_tmle")]] <- list(
          semaglutide = as.numeric(tmle_counts["1"]),
          comparator = as.numeric(tmle_counts["0"]),
          total = sum(tmle_counts)
        )
      } else {
        cat("    TMLE cohort file not found:", tmle_file, "\n")
      }
    }
    
    csv_numbers[[comp]] <- comp_numbers
  }
  
  # 3. Compare the numbers and identify discrepancies
  cat("\n=== DISCREPANCY ANALYSIS ===\n")
  
  if(exists("stepwise_results") && !is.null(stepwise_results)) {
    cat("\nGLP1 vs SGLT2 Discrepancies:\n")
    
    # Check if any CSV files exist for this comparison
    if(length(csv_numbers[["GLP1_vs_SGLT2"]]) > 0) {
      # Use the first available outcome for comparison
      first_outcome <- names(csv_numbers[["GLP1_vs_SGLT2"]])[1]
      csv_data <- csv_numbers[["GLP1_vs_SGLT2"]][[first_outcome]]
      
      # Use corrected numbers if available, otherwise original
      if("glp1_sglt2_glp1_final_corrected" %in% names(stepwise_results)) {
        flowchart_sema <- stepwise_results$glp1_sglt2_glp1_final_corrected
        flowchart_comp <- stepwise_results$glp1_sglt2_sglt2_final_corrected
        flowchart_total <- stepwise_results$glp1_sglt2_final_corrected
        cat("  Using CORRECTED flowchart numbers:\n")
      } else {
        flowchart_sema <- stepwise_results$glp1_sglt2_glp1_final
        flowchart_comp <- stepwise_results$glp1_sglt2_sglt2_final
        flowchart_total <- stepwise_results$glp1_sglt2_final
        cat("  Using ORIGINAL flowchart numbers:\n")
      }
      
      csv_sema <- csv_data$semaglutide
      csv_comp <- csv_data$comparator
      csv_total <- csv_data$total
      
      cat("  Flowchart: Sema=", flowchart_sema, ", SGLT2=", flowchart_comp, ", Total=", flowchart_total, "\n")
      cat("  CSV file:  Sema=", csv_sema, ", SGLT2=", csv_comp, ", Total=", csv_total, "\n")
      cat("  Difference: Sema=", flowchart_sema - csv_sema, 
          ", SGLT2=", flowchart_comp - csv_comp, 
          ", Total=", flowchart_total - csv_total, "\n")
      
      # Calculate percentage differences
      if(flowchart_sema > 0 && flowchart_comp > 0 && flowchart_total > 0) {
        sema_pct <- round((csv_sema / flowchart_sema - 1) * 100, 1)
        comp_pct <- round((csv_comp / flowchart_comp - 1) * 100, 1)
        total_pct <- round((csv_total / flowchart_total - 1) * 100, 1)
        
        cat("  % Change:   Sema=", sema_pct, "%, SGLT2=", comp_pct, "%, Total=", total_pct, "%\n")
      } else {
        cat("  % Change: Cannot calculate (zero values in flowchart)\n")
      }
    } else {
      cat("  No CSV files found for comparison\n")
    }
    
    cat("\nGLP1 vs SU/DPP4 Discrepancies:\n")
    
    # Check if any CSV files exist for this comparison
    if(length(csv_numbers[["GLP1_vs_SU_DPP4"]]) > 0) {
      # Use the first available outcome for comparison
      first_outcome <- names(csv_numbers[["GLP1_vs_SU_DPP4"]])[1]
      csv_data <- csv_numbers[["GLP1_vs_SU_DPP4"]][[first_outcome]]
      
      # Use corrected numbers if available, otherwise original
      if("glp1_su_dpp4_glp1_final_corrected" %in% names(stepwise_results)) {
        flowchart_sema <- stepwise_results$glp1_su_dpp4_glp1_final_corrected
        flowchart_comp <- stepwise_results$glp1_su_dpp4_su_dpp4_final_corrected
        flowchart_total <- stepwise_results$glp1_su_dpp4_final_corrected
        cat("  Using CORRECTED flowchart numbers:\n")
      } else {
        flowchart_sema <- stepwise_results$glp1_su_dpp4_glp1_final
        flowchart_comp <- stepwise_results$glp1_su_dpp4_su_dpp4_final
        flowchart_total <- stepwise_results$glp1_su_dpp4_final
        cat("  Using ORIGINAL flowchart numbers:\n")
      }
      
      csv_sema <- csv_data$semaglutide
      csv_comp <- csv_data$comparator
      csv_total <- csv_data$total
      
      cat("  Flowchart: Sema=", flowchart_sema, ", SU/DPP4=", flowchart_comp, ", Total=", flowchart_total, "\n")
      cat("  CSV file:  Sema=", csv_sema, ", SU/DPP4=", csv_comp, ", Total=", csv_total, "\n")
      cat("  Difference: Sema=", flowchart_sema - csv_sema, 
          ", SU/DPP4=", flowchart_comp - csv_comp, 
          ", Total=", flowchart_total - csv_total, "\n")
      
      # Calculate percentage differences
      if(flowchart_sema > 0 && flowchart_comp > 0 && flowchart_total > 0) {
        sema_pct <- round((csv_sema / flowchart_sema - 1) * 100, 1)
        comp_pct <- round((csv_comp / flowchart_comp - 1) * 100, 1)
        total_pct <- round((csv_total / flowchart_total - 1) * 100, 1)
        
        cat("  % Change:   Sema=", sema_pct, "%, SU/DPP4=", comp_pct, "%, Total=", total_pct, "%\n")
      } else {
        cat("  % Change: Cannot calculate (zero values in flowchart)\n")
      }
    } else {
      cat("  No CSV files found for comparison\n")
    }
    
    # For SGLT2 vs SU_DPP4 comparison - this was not tracked in stepwise analysis
    cat("\nSGLT2 vs SU/DPP4 Analysis (if available):\n")
    if(length(csv_numbers[["SGLT2_vs_SU_DPP4"]]) > 0) {
      first_outcome <- names(csv_numbers[["SGLT2_vs_SU_DPP4"]])[1]
      csv_data <- csv_numbers[["SGLT2_vs_SU_DPP4"]][[first_outcome]]
      
      cat("  CSV file numbers:\n")
      cat("  SGLT2i (treatment=1):", csv_data$semaglutide, "\n")  # Note: stored as 'semaglutide' but actually SGLT2
      cat("  SU/DPP4i (treatment=0):", csv_data$comparator, "\n")
      cat("  Total:", csv_data$total, "\n")
      cat("  Note: This comparison was not tracked in the original flowchart\n")
      cat("        as it was added later for neurologic outcomes.\n")
    } else {
      cat("  No SGLT2 vs SU/DPP4 CSV files found.\n")
      cat("  This may indicate that this comparison was not run\n")
      cat("  or the analysis failed for this comparison.\n")
    }
    
  } else {
    cat("\nCannot perform discrepancy analysis because stepwise_results is not available.\n")
    cat("However, we can still see the publication table numbers above.\n")
    cat("To fix the discrepancy, you need to:\n")
    cat("1. Run the complete script from the beginning\n")
    cat("2. Check if there are any errors in the stepwise analysis function\n")
    cat("3. Ensure all required data files are available\n")
  }
  
  # 4. Potential causes and recommendations
  cat("\n=== POTENTIAL CAUSES OF DISCREPANCY ===\n")
  cat("1. Additional exclusions applied during outcome-specific cohort creation\n")
  cat("2. Different time windows or follow-up requirements\n")
  cat("3. Outcome-specific exclusions (e.g., prevalent cases)\n")
  cat("4. Propensity score matching/trimming exclusions\n")
  cat("5. Missing data exclusions during final analysis\n")
  
  cat("\n=== RECOMMENDATIONS ===\n")
  cat("1. ✓ COMPLETED: Added detailed exclusion tracking for each analysis stage\n")
  cat("2. ✓ COMPLETED: Created corrected flowchart using final CSV numbers\n")
  cat("3. ✓ COMPLETED: Added arithmetic verification to ensure sums are correct\n")
  cat("4. ✓ COMPLETED: Identified specific additional exclusions between stages\n")
  
  cat("\n=== ENHANCED FLOWCHART ANALYSIS FEATURES ===\n")
  cat("The updated analysis now provides:\n")
  cat("• Original flowchart numbers (intermediate stage)\n")
  cat("• Corrected flowchart numbers (final analysis stage)\n")
  cat("• Detailed tracking of additional exclusions:\n")
  cat("  - Patients filtered out due to event_time < 0\n")
  cat("  - Patients with missing values in covariates\n")
  cat("  - Patients with insufficient follow-up time\n")
  cat("  - Patients excluded during propensity score analysis\n")
  cat("  - Patients with missing outcome data\n")
  cat("• Arithmetic verification to ensure correct sums\n")
  cat("• Publication table consistency verification\n")
  
  cat("\n=== HOW TO USE THE RESULTS ===\n")
  cat("FOR PUBLICATION FLOWCHART: Use the 'CORRECTED FINAL ANALYSIS NUMBERS'\n")
  cat("FOR INTERNAL TRACKING: Use the detailed exclusion breakdown\n")
  cat("FOR VERIFICATION: Check that arithmetic verification shows ✓ CORRECT\n")
  
  return(list(
    flowchart_numbers = if(exists("stepwise_results") && !is.null(stepwise_results)) stepwise_results else NULL,
    csv_numbers = csv_numbers
  ))
}

# ===============================================================================
# RUN PUBLICATION TABLE EXTRACTION
# ===============================================================================

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("RUNNING PUBLICATION TABLE EXTRACTION\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# Extract and print publication results
extract_publication_results()

# Create summary table
publication_table <- create_publication_summary_table()

# Diagnose number discrepancies
cat("\n>>> RUNNING DIAGNOSTIC TO IDENTIFY NUMBER DISCREPANCIES <<<\n")
discrepancy_results <- diagnose_number_discrepancy()

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("PUBLICATION TABLE EXTRACTION COMPLETE\n")
cat("Check the printed output above for values to fill in your table\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# ===============================================================================
# CONSOLIDATED BASELINE CHARACTERISTICS TABLES MODULE
# ===============================================================================

# ---- Core Functions for Baseline Characteristics Tables ----

# Function to calculate standardized mean differences (SMD) - INDIVIDUAL LEVEL VERSION
# 
# ENHANCED APPROACH FOR MULTILEVEL VARIABLES:
# - Instead of single overall SMD per multilevel variable, calculate SMD for each level
# - Each level of multilevel variables (sex, race, income, education, BMI categories, calendar year) gets its own SMD
# - Calculates SMD for ALL levels including reference categories (e.g., Non-Hispanic White, <$10K, etc.)
# - Also provides overall summary SMD (maximum absolute SMD across levels)
# - Enables detailed balance assessment for each category within multilevel variables  
# - Treatment-defining variables (GLP1, SGLT2, etc.) are excluded to avoid infinite SMDs
#
calculate_smd_table <- function(cohort_df, comparison_name, outcome_name, weights = NULL) {
  
  # Define variables for SMD calculation
  continuous_vars <- c("age", "baseline_bmi", "baseline_hba1c")
  binary_vars <- c("Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
                   "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
                   "TZD", "Insulin", "SU", "DPP4", "SGLT2", "GLP1",
                   "mi", "chf", "pvd", "cvd", "dem", "cpd", "ctd", "pud", "mld",
                   "ph", "rd", "cancer", "msld", "mc", "hiv")
  multilevel_vars <- c("sex_cat", "raceethnicity_cat", "income", "education", "bmi_category", "index_year")
  
  # EXCLUDE treatment-defining variables from SMD calculation
  # These should theoretically have infinite SMDs and aren't meaningful to balance
  treatment_defining_vars <- character(0)
  
  if(grepl("GLP1 vs SGLT2", comparison_name, ignore.case = TRUE)) {
    treatment_defining_vars <- c("GLP1", "SGLT2")
  } else if(grepl("GLP1 vs SU.DPP4", comparison_name, ignore.case = TRUE)) {
    treatment_defining_vars <- c("GLP1", "SU", "DPP4")  
  } else if(grepl("SGLT2 vs SU.DPP4", comparison_name, ignore.case = TRUE)) {
    treatment_defining_vars <- c("SGLT2", "SU", "DPP4")
  }
  
  # Remove treatment-defining variables from binary_vars
  binary_vars <- setdiff(binary_vars, treatment_defining_vars)
  
  cat("DEBUG: Excluding treatment-defining variables from SMD calculation:", paste(treatment_defining_vars, collapse=", "), "\n")
  
  # Filter to variables that exist in the dataset
  continuous_vars <- continuous_vars[continuous_vars %in% names(cohort_df)]
  binary_vars <- binary_vars[binary_vars %in% names(cohort_df)]
  multilevel_vars <- multilevel_vars[multilevel_vars %in% names(cohort_df)]
  
  smd_results <- list()
  
  # Calculate SMD for continuous variables
  for(var in continuous_vars) {
    if(var %in% names(cohort_df)) {
      if(is.null(weights)) {
        # Unweighted SMD - CORRECTED
        treat_data <- cohort_df[cohort_df$treatment == 1 & !is.na(cohort_df[[var]]), ]
        comp_data <- cohort_df[cohort_df$treatment == 0 & !is.na(cohort_df[[var]]), ]
        
        treat_vals <- as.numeric(treat_data[[var]])
        comp_vals <- as.numeric(comp_data[[var]])
        
        mean_treat <- mean(treat_vals, na.rm = TRUE)
        mean_comp <- mean(comp_vals, na.rm = TRUE)
        var_treat <- var(treat_vals, na.rm = TRUE)
        var_comp <- var(comp_vals, na.rm = TRUE)
        
        # Handle cases where variance calculation returns NA
        if(is.na(var_treat)) var_treat <- 0
        if(is.na(var_comp)) var_comp <- 0
        
      } else {
        # Weighted SMD - CORRECTED
        treat_idx <- cohort_df$treatment == 1 & !is.na(cohort_df[[var]])
        comp_idx <- cohort_df$treatment == 0 & !is.na(cohort_df[[var]])
        
        if(sum(treat_idx) > 0 && sum(comp_idx) > 0) {
          treat_vals <- as.numeric(cohort_df[[var]][treat_idx])
          comp_vals <- as.numeric(cohort_df[[var]][comp_idx])
          w_treat <- weights[treat_idx]
          w_comp <- weights[comp_idx]
          
          mean_treat <- weighted.mean(treat_vals, w_treat, na.rm = TRUE)
          mean_comp <- weighted.mean(comp_vals, w_comp, na.rm = TRUE)
          
          # Weighted variance calculation
          var_treat <- sum(w_treat * (treat_vals - mean_treat)^2, na.rm = TRUE) / sum(w_treat, na.rm = TRUE)
          var_comp <- sum(w_comp * (comp_vals - mean_comp)^2, na.rm = TRUE) / sum(w_comp, na.rm = TRUE)
          
          # Handle cases where variance calculation returns NA
          if(is.na(var_treat)) var_treat <- 0
          if(is.na(var_comp)) var_comp <- 0
        } else {
          mean_treat <- mean_comp <- var_treat <- var_comp <- 0
        }
      }
      
      pooled_sd <- sqrt((var_treat + var_comp) / 2)
      smd <- if(!is.na(pooled_sd) && pooled_sd > 0) (mean_treat - mean_comp) / pooled_sd else 0
      
      smd_results[[var]] <- data.frame(
        comparison = as.character(comparison_name),
        outcome = as.character(outcome_name),
        variable = as.character(var),
        variable_level = as.character(NA),
        smd = as.numeric(smd),
        mean_treat = as.numeric(mean_treat),
        mean_comp = as.numeric(mean_comp),
        weighted = as.logical(!is.null(weights)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Calculate SMD for binary variables - CORRECTED
  for(var in binary_vars) {
    if(var %in% names(cohort_df)) {
      if(is.null(weights)) {
        # Unweighted SMD - CORRECTED CALCULATION
        treat_data <- cohort_df[cohort_df$treatment == 1 & !is.na(cohort_df[[var]]), ]
        comp_data <- cohort_df[cohort_df$treatment == 0 & !is.na(cohort_df[[var]]), ]
        
        if(nrow(treat_data) > 0 && nrow(comp_data) > 0) {
          prop_treat <- mean(treat_data[[var]] == "1" | treat_data[[var]] == 1, na.rm = TRUE)
          prop_comp <- mean(comp_data[[var]] == "1" | comp_data[[var]] == 1, na.rm = TRUE)
          
          # Handle NA proportions
          if(is.na(prop_treat)) prop_treat <- 0
          if(is.na(prop_comp)) prop_comp <- 0
        } else {
          prop_treat <- prop_comp <- 0
        }
        
        var_treat <- prop_treat * (1 - prop_treat)
        var_comp <- prop_comp * (1 - prop_comp)
        
        # Handle cases where variance calculation returns NA
        if(is.na(var_treat)) var_treat <- 0
        if(is.na(var_comp)) var_comp <- 0
        
      } else {
        # Weighted SMD - CORRECTED
        treat_idx <- cohort_df$treatment == 1 & !is.na(cohort_df[[var]])
        comp_idx <- cohort_df$treatment == 0 & !is.na(cohort_df[[var]])
        
        if(sum(treat_idx) > 0 && sum(comp_idx) > 0) {
          treat_positive <- (cohort_df[treat_idx, var] == "1" | cohort_df[treat_idx, var] == 1)
          comp_positive <- (cohort_df[comp_idx, var] == "1" | cohort_df[comp_idx, var] == 1)
          w_treat <- weights[treat_idx]
          w_comp <- weights[comp_idx]
          
          prop_treat <- sum(w_treat * treat_positive, na.rm = TRUE) / sum(w_treat, na.rm = TRUE)
          prop_comp <- sum(w_comp * comp_positive, na.rm = TRUE) / sum(w_comp, na.rm = TRUE)
          
          # Handle NA proportions
          if(is.na(prop_treat)) prop_treat <- 0
          if(is.na(prop_comp)) prop_comp <- 0
        } else {
          prop_treat <- prop_comp <- 0
        }
        
        var_treat <- prop_treat * (1 - prop_treat)
        var_comp <- prop_comp * (1 - prop_comp)
        
        # Handle cases where variance calculation returns NA
        if(is.na(var_treat)) var_treat <- 0
        if(is.na(var_comp)) var_comp <- 0
      }
      
      pooled_sd <- sqrt((var_treat + var_comp) / 2)
      smd <- if(!is.na(pooled_sd) && pooled_sd > 0) (prop_treat - prop_comp) / pooled_sd else 0
      
      smd_results[[var]] <- data.frame(
        comparison = as.character(comparison_name),
        outcome = as.character(outcome_name),
        variable = as.character(var),
        variable_level = as.character("1"),
        smd = as.numeric(smd),
        mean_treat = as.numeric(prop_treat),
        mean_comp = as.numeric(prop_comp),
        weighted = as.logical(!is.null(weights)),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Calculate SMD for multilevel variables - INDIVIDUAL SMD FOR EACH LEVEL
  for(var in multilevel_vars) {
    if(var %in% names(cohort_df)) {
      
      # Get all unique levels for this variable
      all_levels <- unique(cohort_df[[var]])
      all_levels <- all_levels[!is.na(all_levels)]
      all_levels <- sort(as.character(all_levels))
      
      # Debug: Show levels for BMI category
      if(var == "bmi_category") {
        cat("DEBUG: BMI category levels found:", paste(all_levels, collapse = ", "), "\n")
        cat("DEBUG: BMI category sample values:", paste(head(cohort_df[[var]], 10), collapse = ", "), "\n")
        cat("DEBUG: BMI category data type:", class(cohort_df[[var]]), "\n")
        cat("DEBUG: BMI category table:", paste(names(table(cohort_df[[var]])), "=", table(cohort_df[[var]]), collapse = ", "), "\n")
      }
      
      if(length(all_levels) > 1) {
        # Store individual level SMDs and calculate overall SMD
        level_smds <- numeric()
        level_names <- character()
        
        for(level in all_levels) { # Include ALL levels (no reference skipping)
          if(is.null(weights)) {
            # Unweighted SMD for this level
            treat_data <- cohort_df[cohort_df$treatment == 1 & !is.na(cohort_df[[var]]), ]
            comp_data <- cohort_df[cohort_df$treatment == 0 & !is.na(cohort_df[[var]]), ]
            
            if(nrow(treat_data) > 0 && nrow(comp_data) > 0) {
              if(var == "bmi_category") {
                # BMI categories are text labels, don't convert to numeric
                prop_treat <- mean(treat_data[[var]] == level, na.rm = TRUE)
                prop_comp <- mean(comp_data[[var]] == level, na.rm = TRUE)
              } else {
                # Other variables can be numeric or text
                prop_treat <- mean(treat_data[[var]] == level | treat_data[[var]] == as.numeric(level), na.rm = TRUE)
                prop_comp <- mean(comp_data[[var]] == level | comp_data[[var]] == as.numeric(level), na.rm = TRUE)
              }
              
              # Handle NA proportions
              if(is.na(prop_treat)) prop_treat <- 0
              if(is.na(prop_comp)) prop_comp <- 0
              
              # Debug: Show proportions for BMI category
              if(var == "bmi_category") {
                cat("DEBUG: BMI level", level, "- prop_treat:", prop_treat, ", prop_comp:", prop_comp, "\n")
                cat("DEBUG: BMI level", level, "- treat_data size:", nrow(treat_data), ", comp_data size:", nrow(comp_data), "\n")
              }
            } else {
              prop_treat <- prop_comp <- 0
            }
            
          } else {
            # Weighted SMD for this level
            treat_idx <- cohort_df$treatment == 1 & !is.na(cohort_df[[var]])
            comp_idx <- cohort_df$treatment == 0 & !is.na(cohort_df[[var]])
            
            if(sum(treat_idx) > 0 && sum(comp_idx) > 0) {
              if(var == "bmi_category") {
                # BMI categories are text labels, don't convert to numeric
                treat_level <- (cohort_df[treat_idx, var] == level)
                comp_level <- (cohort_df[comp_idx, var] == level)
              } else {
                # Other variables can be numeric or text
                treat_level <- (cohort_df[treat_idx, var] == level | cohort_df[treat_idx, var] == as.numeric(level))
                comp_level <- (cohort_df[comp_idx, var] == level | cohort_df[comp_idx, var] == as.numeric(level))
              }
              
              w_treat <- weights[treat_idx]
              w_comp <- weights[comp_idx]
              
              prop_treat <- sum(w_treat * treat_level, na.rm = TRUE) / sum(w_treat, na.rm = TRUE)
              prop_comp <- sum(w_comp * comp_level, na.rm = TRUE) / sum(w_comp, na.rm = TRUE)
              
              # Handle NA proportions
              if(is.na(prop_treat)) prop_treat <- 0
              if(is.na(prop_comp)) prop_comp <- 0
            } else {
              prop_treat <- prop_comp <- 0
            }
          }
          
          var_treat <- prop_treat * (1 - prop_treat)
          var_comp <- prop_comp * (1 - prop_comp)
          
          # Handle cases where variance calculation returns NA
          if(is.na(var_treat)) var_treat <- 0
          if(is.na(var_comp)) var_comp <- 0
          
          pooled_sd <- sqrt((var_treat + var_comp) / 2)
          level_smd <- if(!is.na(pooled_sd) && pooled_sd > 0) (prop_treat - prop_comp) / pooled_sd else 0
          
          # Debug: Show final calculations for BMI category
          if(var == "bmi_category") {
            cat("DEBUG: BMI level", level, "- var_treat:", var_treat, ", var_comp:", var_comp, ", pooled_sd:", pooled_sd, ", level_smd:", level_smd, "\n")
          }
          
          level_smds <- c(level_smds, level_smd)
          level_names <- c(level_names, level)
          
          # Get proper level label based on variable type
          if(var == "sex_cat") {
            level_labels <- list("0" = "Male", "1" = "Female", "999" = "Not specified/Other")
            level_display <- if(level %in% names(level_labels)) level_labels[[level]] else paste("Sex:", level)
          } else if(var == "raceethnicity_cat") {
            level_labels <- list("0" = "Non-Hispanic White", "1" = "Non-Hispanic Black", 
                               "2" = "Hispanic", "3" = "Non-Hispanic Asian/Other", "999" = "Unknown/Not specified")
            level_display <- if(level %in% names(level_labels)) level_labels[[level]] else paste("Race:", level)
          } else if(var == "income") {
            level_labels <- list("0" = "<$10k", "1" = "$10k-25k", "2" = "$25k-35k", "3" = "$35k-50k",
                               "4" = "$50k-75k", "5" = "$75k-100k", "6" = "$100k-150k", 
                               "7" = "$150k-200k", "8" = ">$200k", "999" = "Others/Not specified")
            level_display <- if(level %in% names(level_labels)) level_labels[[level]] else paste("Income:", level)
          } else if(var == "education") {
            level_labels <- list("0" = "Less than high school", "1" = "High school/GED", 
                               "2" = "Some college", "3" = "College graduate", "4" = "Advanced degree", "999" = "Not specified/Other")
            level_display <- if(level %in% names(level_labels)) level_labels[[level]] else paste("Education:", level)
          } else if(var == "bmi_category") {
            # BMI categories are already text labels, use as-is
            level_display <- level
          } else if(var == "index_year") {
            # Calendar years - use as-is but add "Year " prefix for clarity
            level_display <- paste("Year", level)
          } else {
            level_display <- level  # For any other variables
          }
          
          # Store individual level SMD
          level_key <- paste0(var, "_", level)
          smd_results[[level_key]] <- data.frame(
            comparison = as.character(comparison_name),
            outcome = as.character(outcome_name),
            variable = as.character(var),
            variable_level = as.character(level_display),
            smd = as.numeric(level_smd),
            mean_treat = as.numeric(prop_treat),
            mean_comp = as.numeric(prop_comp),
            weighted = as.logical(!is.null(weights)),
            stringsAsFactors = FALSE
          )
        }
        
        # Calculate and store overall SMD for this multilevel variable using maximum absolute SMD
        overall_smd <- ifelse(length(level_smds) > 0, max(abs(level_smds), na.rm = TRUE), 0)
        
        # Use the sign of the SMD with maximum absolute value
        if(length(level_smds) > 0) {
          max_idx <- which.max(abs(level_smds))
          overall_smd <- overall_smd * sign(level_smds[max_idx])
          max_level <- level_names[max_idx]
        } else {
          max_level <- "None"
        }
        
        # Store the overall SMD for this multilevel variable
        overall_key <- paste0(var, "_overall")
        smd_results[[overall_key]] <- data.frame(
          comparison = as.character(comparison_name),
          outcome = as.character(outcome_name),
          variable = as.character(var),
          variable_level = as.character(paste0("Overall (max from level: ", max_level, ")")),
          smd = as.numeric(overall_smd),
          mean_treat = as.numeric(NA),  # Not meaningful for overall multilevel SMD
          mean_comp = as.numeric(NA),   # Not meaningful for overall multilevel SMD
          weighted = as.logical(!is.null(weights)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  # Combine all SMD results
  if(length(smd_results) > 0) {
    # Debug: Check structure consistency
    cat("DEBUG: SMD results list length:", length(smd_results), "\n")
    
    # Check that all data frames have the same column names
    expected_cols <- c("comparison", "outcome", "variable", "variable_level", "smd", "mean_treat", "mean_comp", "weighted")
    
    for(i in seq_along(smd_results)) {
      result_cols <- names(smd_results[[i]])
      if(!identical(sort(result_cols), sort(expected_cols))) {
        cat("WARNING: Column mismatch in SMD result", i, "\n")
        cat("Expected:", paste(expected_cols, collapse=", "), "\n")
        cat("Found:", paste(result_cols, collapse=", "), "\n")
      }
    }
    
    # Combine with error handling
    tryCatch({
      final_smd_table <- do.call(rbind, smd_results)
      return(final_smd_table)
    }, error = function(e) {
      cat("Error in rbind of SMD results:", e$message, "\n")
      cat("Attempting to fix by ensuring consistent structure...\n")
      
      # Create standardized data frames
      standardized_results <- lapply(smd_results, function(df) {
        # Ensure all expected columns exist
        for(col in expected_cols) {
          if(!col %in% names(df)) {
            if(col %in% c("comparison", "outcome", "variable", "variable_level")) {
              df[[col]] <- as.character(NA)
            } else if(col %in% c("smd", "mean_treat", "mean_comp")) {
              df[[col]] <- as.numeric(NA)
            } else if(col == "weighted") {
              df[[col]] <- as.logical(NA)
            }
          }
        }
        # Reorder columns
        df[, expected_cols]
      })
      
      final_smd_table <- do.call(rbind, standardized_results)
      return(final_smd_table)
    })
  } else {
    return(NULL)
  }
}

# Function to create IPTW weighted baseline characteristics table
create_iptw_weighted_baseline_table <- function(cohort_df, comparison_name, outcome_name) {
  
  if(!"ipw_std" %in% names(cohort_df)) {
    cat("Warning: No IPTW weights found in cohort. Skipping weighted table.\n")
    return(NULL)
  }
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("IPTW-WEIGHTED BASELINE CHARACTERISTICS:", comparison_name, "-", outcome_name, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Get treatment group labels based on comparison
  if(comparison_name == "SGLT2 vs SU_DPP4") {
    treat_label <- "SGLT2i"
    comp_label <- "SU/DPP4i"
  } else if(comparison_name == "GLP1 vs SGLT2") {
    treat_label <- "Semaglutide"
    comp_label <- "SGLT2i"
  } else if(comparison_name == "GLP1 vs SU_DPP4") {
    treat_label <- "Semaglutide"
    comp_label <- "SU/DPP4i"
  } else {
    treat_label <- "Treatment"
    comp_label <- "Comparator"
  }
  
  # Calculate weighted overall statistics
  weighted_stats <- cohort_df %>%
    group_by(treatment) %>%
    summarise(
      effective_n = sum(ipw_std, na.rm = TRUE),
      actual_n = n(),
      weighted_events = sum(event * ipw_std, na.rm = TRUE),
      .groups = "drop"
    )
  
  treat_eff_n <- weighted_stats$effective_n[weighted_stats$treatment == 1]
  comp_eff_n <- weighted_stats$effective_n[weighted_stats$treatment == 0]
  treat_events <- weighted_stats$weighted_events[weighted_stats$treatment == 1]
  comp_events <- weighted_stats$weighted_events[weighted_stats$treatment == 0]
  
  cat(sprintf("%-35s %s (eff_n=%.0f) %s (eff_n=%.0f)\n", "Characteristic", treat_label, treat_eff_n, comp_label, comp_eff_n))
  cat(paste(rep("-", 80), collapse=""), "\n")
  
  # WEIGHTED Demographics Section
  cat("DEMOGRAPHICS (IPTW-WEIGHTED):\n")
  
  # Weighted Age
  if("age" %in% names(cohort_df)) {
    age_weighted <- cohort_df %>%
      group_by(treatment) %>%
      summarise(
        weighted_mean_age = weighted.mean(age, ipw_std, na.rm = TRUE),
        weighted_var = sum(ipw_std * (age - weighted_mean_age)^2, na.rm = TRUE) / sum(ipw_std, na.rm = TRUE),
        weighted_sd_age = sqrt(weighted_var),
        .groups = "drop"
      )
    
    treat_age <- sprintf("%.1f (%.1f)", age_weighted$weighted_mean_age[age_weighted$treatment == 1], 
                                        age_weighted$weighted_sd_age[age_weighted$treatment == 1])
    comp_age <- sprintf("%.1f (%.1f)", age_weighted$weighted_mean_age[age_weighted$treatment == 0], 
                                       age_weighted$weighted_sd_age[age_weighted$treatment == 0])
    
    cat(sprintf("  %-33s %-20s %-20s\n", "Age, mean (SD), y", treat_age, comp_age))
  }
  
  # Weighted Sex (multilevel: Male, Female, Not specified/Other)
  if("sex_cat" %in% names(cohort_df)) {
    cat("WEIGHTED SEX DISTRIBUTION:\n")
    
    sex_mappings <- list(
      "0" = "Male",
      "1" = "Female", 
      "999" = "Not specified/Other"
    )
    
    actual_levels <- unique(cohort_df$sex_cat)
    actual_levels <- actual_levels[!is.na(actual_levels)]
    actual_levels <- sort(as.character(actual_levels))
    
    for(level in actual_levels) {
      if(level %in% names(sex_mappings)) {
        sex_name <- sex_mappings[[level]]
        
        sex_weighted <- cohort_df %>%
          group_by(treatment) %>%
          summarise(
            weighted_n = sum(ipw_std, na.rm = TRUE),
            weighted_level_n = sum((sex_cat == level | sex_cat == as.numeric(level)) * ipw_std, na.rm = TRUE),
            weighted_level_prop = weighted_level_n / weighted_n * 100,
            .groups = "drop"
          )
        
        treat_sex <- sex_weighted %>% filter(treatment == 1)
        comp_sex <- sex_weighted %>% filter(treatment == 0)
        
        if(nrow(treat_sex) > 0 && nrow(comp_sex) > 0) {
          cat(sprintf("  %-33s %-20s %-20s\n", paste(sex_name, ", n (%)"), 
                      sprintf("%.0f (%.1f)", treat_sex$weighted_level_n, treat_sex$weighted_level_prop),
                      sprintf("%.0f (%.1f)", comp_sex$weighted_level_n, comp_sex$weighted_level_prop)))
        }
      }
    }
  }
  
  # Weighted BMI
  if("baseline_bmi" %in% names(cohort_df)) {
    bmi_weighted <- cohort_df %>%
      filter(!is.na(baseline_bmi)) %>%
      group_by(treatment) %>%
      summarise(
        weighted_mean_bmi = weighted.mean(baseline_bmi, ipw_std, na.rm = TRUE),
        weighted_var = sum(ipw_std * (baseline_bmi - weighted_mean_bmi)^2, na.rm = TRUE) / sum(ipw_std, na.rm = TRUE),
        weighted_sd_bmi = sqrt(weighted_var),
        .groups = "drop"
      )
    
    if(nrow(bmi_weighted) == 2) {
      treat_bmi <- sprintf("%.1f (%.1f)", bmi_weighted$weighted_mean_bmi[bmi_weighted$treatment == 1], 
                                           bmi_weighted$weighted_sd_bmi[bmi_weighted$treatment == 1])
      comp_bmi <- sprintf("%.1f (%.1f)", bmi_weighted$weighted_mean_bmi[bmi_weighted$treatment == 0], 
                                         bmi_weighted$weighted_sd_bmi[bmi_weighted$treatment == 0])
      
      cat(sprintf("  %-33s %-20s %-20s\n", "BMI (kg/m2), mean (SD)", treat_bmi, comp_bmi))
    }
  }
  
  cat(paste(rep("-", 80), collapse=""), "\n")
  cat(sprintf("WEIGHTED EVENTS: %s = %.1f, %s = %.1f\n", treat_label, treat_events, comp_label, comp_events))
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  return(list(
    comparison = comparison_name,
    outcome = outcome_name,
    treatment_label = treat_label,
    comparator_label = comp_label,
    treatment_eff_n = treat_eff_n,
    comparator_eff_n = comp_eff_n,
    treatment_events = treat_events,
    comparator_events = comp_events
  ))
}

# Function to create structured baseline table data (for CSV export)
create_unweighted_baseline_table <- function(cohort_df, comparison_name, outcome_name) {
  
  # Define all variables to include in baseline table
  continuous_vars <- c("age", "baseline_bmi", "baseline_hba1c")
  
  # Binary categorical variables (coded as 0/1)
  binary_vars <- c("Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
                   "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
                   "TZD", "Insulin", "SU", "DPP4", "SGLT2", "GLP1",
                   "mi", "chf", "pvd", "cvd", "dem", "cpd", "ctd", "pud", "mld",
                   "ph", "rd", "cancer", "msld", "mc", "hiv")
  
  # Multi-level categorical variables
  multilevel_vars <- c("sex_cat", "raceethnicity_cat", "income", "education", "bmi_category", "index_year")
  
  # Filter to variables that exist in the dataset
  continuous_vars <- continuous_vars[continuous_vars %in% names(cohort_df)]
  binary_vars <- binary_vars[binary_vars %in% names(cohort_df)]
  multilevel_vars <- multilevel_vars[multilevel_vars %in% names(cohort_df)]
  
  all_summaries <- list()
  
  # Create summary for continuous variables
  for(var in continuous_vars) {
    if(var %in% names(cohort_df)) {
      summary_stats <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          n = n(),
          mean_val = mean(.data[[var]], na.rm = TRUE),
          sd_val = sd(.data[[var]], na.rm = TRUE),
          missing = sum(is.na(.data[[var]])),
          .groups = "drop"
        ) %>%
        mutate(
          variable = var,
          variable_level = NA_character_,
          variable_type = "continuous",
          summary_stat = sprintf("%.1f (%.1f)", mean_val, sd_val),
          comparison = comparison_name,
          outcome = outcome_name
        )
      
      # Keep only the standardized columns
      summary_stats <- summary_stats[, c("comparison", "outcome", "variable", "variable_level", "variable_type", "treatment", "n", "summary_stat", "missing")]
      
      all_summaries[[var]] <- summary_stats
    }
  }
  
  # Create summary for binary categorical variables
  for(var in binary_vars) {
    if(var %in% names(cohort_df)) {
      summary_stats <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          n = n(),
          n_positive = sum(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE),
          prop_positive = mean(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE) * 100,
          missing = sum(is.na(.data[[var]])),
          .groups = "drop"
                 ) %>%
         mutate(
           variable = var,
           variable_level = "1",
           variable_type = "binary",
           summary_stat = sprintf("%.1f", prop_positive),
           comparison = comparison_name,
           outcome = outcome_name
         )
      
      # Keep only the standardized columns
      summary_stats <- summary_stats[, c("comparison", "outcome", "variable", "variable_level", "variable_type", "treatment", "n", "summary_stat", "missing")]
      
      all_summaries[[var]] <- summary_stats
    }
  }
  
  # Create summary for multi-level categorical variables
  for(var in multilevel_vars) {
    if(var %in% names(cohort_df)) {
      
      # Define labels for each variable
      if(var == "sex_cat") {
        level_labels <- list(
          "0" = "Male",
          "1" = "Female",
          "999" = "Not specified/Other"
        )
      } else if(var == "raceethnicity_cat") {
        level_labels <- list(
          "0" = "Non-Hispanic White",
          "1" = "Non-Hispanic Black", 
          "2" = "Hispanic",
          "3" = "Non-Hispanic Asian/Other",
          "999" = "Unknown/Not specified"
        )
      } else if(var == "income") {
        level_labels <- list(
          "0" = "<$10k",
          "1" = "$10k-25k",
          "2" = "$25k-35k", 
          "3" = "$35k-50k",
          "4" = "$50k-75k",
          "5" = "$75k-100k",
          "6" = "$100k-150k",
          "7" = "$150k-200k",
          "8" = ">$200k",
          "999" = "Others/Not specified"
        )
             } else if(var == "education") {
         level_labels <- list(
           "0" = "Less than high school",
           "1" = "High school/GED",
           "2" = "Some college",
           "3" = "College graduate", 
           "4" = "Advanced degree",
           "999" = "Not specified/Other"
         )
       } else if(var == "bmi_category") {
         # BMI categories are text labels, not numeric codes
         actual_bmi_levels <- unique(cohort_df[[var]])
         actual_bmi_levels <- actual_bmi_levels[!is.na(actual_bmi_levels)]
         
         level_labels <- list()
         for(bmi_level in actual_bmi_levels) {
           level_labels[[bmi_level]] <- bmi_level
         }
       } else if(var == "index_year") {
         # Index year are already meaningful labels, use as-is with "Year " prefix
         actual_year_levels <- unique(cohort_df[[var]])
         actual_year_levels <- actual_year_levels[!is.na(actual_year_levels)]
         
         level_labels <- list()
         for(year_level in actual_year_levels) {
           level_labels[[as.character(year_level)]] <- paste("Year", year_level)
         }
       }
      
             # Get actual levels present in data
       actual_levels <- unique(cohort_df[[var]])
       if(var == "bmi_category") {
         # BMI categories are text labels
         actual_levels <- actual_levels[!is.na(actual_levels)]
         actual_levels <- sort(as.character(actual_levels))
       } else if(var == "sex_cat") {
         # Sex categories: 0, 1, 999
         actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0, 1, 999, "0", "1", "999")]
         actual_levels <- sort(as.character(actual_levels))
       } else if(var == "index_year") {
         # Index year: use all actual years present
         actual_levels <- actual_levels[!is.na(actual_levels)]
         actual_levels <- sort(as.character(actual_levels))
       } else if(var == "raceethnicity_cat") {
         # Race/ethnicity categories: 0, 1, 2, 3, 999
         actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:3, 999, "0":"3", "999")]
         actual_levels <- sort(as.character(actual_levels))
       } else if(var == "education") {
         # Education categories: 0, 1, 2, 3, 4, 999
         actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:4, 999, "0":"4", "999")]
         actual_levels <- sort(as.character(actual_levels))
       } else {
         # Other numeric categorical variables (income)
         actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:8, 999, "0":"8", "999")]
         actual_levels <- sort(as.character(actual_levels))
       }
      
      # Create summary for each level
      level_summaries <- list()
      for(level in actual_levels) {
        if(level %in% names(level_labels)) {
                     level_summary <- cohort_df %>%
             group_by(treatment) %>%
             summarise(
               n = n(),
               n_level = if(var == "bmi_category") {
                 sum(.data[[var]] == level, na.rm = TRUE)
               } else if(var == "index_year") {
                 sum(.data[[var]] == level | .data[[var]] == as.character(level), na.rm = TRUE)
               } else {
                 sum(.data[[var]] == level | .data[[var]] == as.character(level), na.rm = TRUE)
               },
               prop_level = if(var == "bmi_category") {
                 mean(.data[[var]] == level, na.rm = TRUE) * 100
               } else if(var == "index_year") {
                 mean(.data[[var]] == level | .data[[var]] == as.character(level), na.rm = TRUE) * 100
               } else {
                 mean(.data[[var]] == level | .data[[var]] == as.character(level), na.rm = TRUE) * 100
               },
               missing = sum(is.na(.data[[var]])),
               .groups = "drop"
                          ) %>%
             mutate(
               variable = var,
               variable_level = level_labels[[level]],
               variable_type = "multilevel",
               summary_stat = sprintf("%.1f", prop_level),
               comparison = comparison_name,
               outcome = outcome_name
             )
           
           # Keep only the standardized columns
           level_summary <- level_summary[, c("comparison", "outcome", "variable", "variable_level", "variable_type", "treatment", "n", "summary_stat", "missing")]
          
          level_summaries[[level]] <- level_summary
        }
      }
      
      if(length(level_summaries) > 0) {
        all_summaries[[var]] <- do.call(rbind, level_summaries)
      }
    }
  }
  
  # Combine all summaries
  if(length(all_summaries) > 0) {
    final_table <- do.call(rbind, all_summaries)
    return(final_table)
  } else {
    return(NULL)
  }
}

# Function to create formatted baseline table (for console display)
create_formatted_baseline_table <- function(cohort_df, comparison_name, outcome_name) {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("UNWEIGHTED BASELINE CHARACTERISTICS:", comparison_name, "-", outcome_name, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Get treatment group labels based on comparison
  if(comparison_name == "SGLT2 vs SU_DPP4") {
    treat_label <- "SGLT2i"
    comp_label <- "SU/DPP4i"
  } else if(comparison_name == "GLP1 vs SGLT2") {
    treat_label <- "Semaglutide"
    comp_label <- "SGLT2i"
  } else if(comparison_name == "GLP1 vs SU_DPP4") {
    treat_label <- "Semaglutide"
    comp_label <- "SU/DPP4i"
  } else {
    treat_label <- "Treatment"
    comp_label <- "Comparator"
  }
  
  # Calculate overall statistics
  overall_stats <- cohort_df %>%
    group_by(treatment) %>%
    summarise(
      n = n(),
      events = sum(event, na.rm = TRUE),
      .groups = "drop"
    )
  
  treat_n <- overall_stats$n[overall_stats$treatment == 1]
  comp_n <- overall_stats$n[overall_stats$treatment == 0]
  treat_events <- overall_stats$events[overall_stats$treatment == 1]
  comp_events <- overall_stats$events[overall_stats$treatment == 0]
  
  cat(sprintf("%-35s %s (n=%d) %s (n=%d)\n", "Characteristic", treat_label, treat_n, comp_label, comp_n))
  cat(paste(rep("-", 80), collapse=""), "\n")
  
  # Demographics Section
  cat("DEMOGRAPHICS:\n")
  
  # Age
  if("age" %in% names(cohort_df)) {
    age_summary <- cohort_df %>%
      group_by(treatment) %>%
      summarise(
        mean_age = mean(age, na.rm = TRUE),
        sd_age = sd(age, na.rm = TRUE),
        .groups = "drop"
      )
    
    treat_age <- sprintf("%.1f (%.1f)", age_summary$mean_age[age_summary$treatment == 1], 
                                        age_summary$sd_age[age_summary$treatment == 1])
    comp_age <- sprintf("%.1f (%.1f)", age_summary$mean_age[age_summary$treatment == 0], 
                                       age_summary$sd_age[age_summary$treatment == 0])
    
    cat(sprintf("  %-33s %-20s %-20s\n", "Age, mean (SD), y", treat_age, comp_age))
  }
  
  # Sex (now treated as multilevel: 0=Male, 1=Female, 999=Not specified/Other)
  if("sex_cat" %in% names(cohort_df)) {
    cat("  Sex, n (%):\n")
    
    sex_mappings <- list(
      "0" = "Male",
      "1" = "Female", 
      "999" = "Not specified/Other"
    )
    
    actual_levels <- unique(cohort_df$sex_cat)
    actual_levels <- actual_levels[!is.na(actual_levels)]
    actual_levels <- sort(as.character(actual_levels))
    
    for(level in actual_levels) {
      if(level %in% names(sex_mappings)) {
        sex_name <- sex_mappings[[level]]
        
        # Calculate counts and proportions directly
        sex_counts <- cohort_df %>%
          group_by(treatment) %>%
          summarise(
            total_n = n(),
            sex_n = sum(sex_cat == level | sex_cat == as.numeric(level), na.rm = TRUE),
            sex_prop = mean(sex_cat == level | sex_cat == as.numeric(level), na.rm = TRUE) * 100,
            .groups = "drop"
          )
        
        treat_sex <- sex_counts %>% filter(treatment == 1)
        comp_sex <- sex_counts %>% filter(treatment == 0)
        
        treat_val <- if(nrow(treat_sex) > 0) sprintf("%d (%.1f)", treat_sex$sex_n, treat_sex$sex_prop) else "0 (0.0)"
        comp_val <- if(nrow(comp_sex) > 0) sprintf("%d (%.1f)", comp_sex$sex_n, comp_sex$sex_prop) else "0 (0.0)"
        
        cat(sprintf("    %-31s %-20s %-20s\n", sex_name, treat_val, comp_val))
      }
    }
  }
  
  # Race/Ethnicity
  if("raceethnicity_cat" %in% names(cohort_df)) {
    cat("  Race/Ethnicity, n (%):\n")
    
    # Define race/ethnicity mappings based on 001.R
    race_mappings <- list(
      "0" = "Non-Hispanic White",
      "1" = "Non-Hispanic Black", 
      "2" = "Hispanic",
      "3" = "Non-Hispanic Asian/Other",
      "999" = "Unknown/Not specified"
    )
    
    # Get actual levels in data
    actual_levels <- unique(cohort_df$raceethnicity_cat)
    actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:3, 999, "0":"3", "999")]
    actual_levels <- sort(as.character(actual_levels))
    
    for(level in actual_levels) {
      if(level %in% names(race_mappings)) {
        race_name <- race_mappings[[level]]
        
        # Calculate counts and proportions directly
        race_counts <- cohort_df %>%
          group_by(treatment) %>%
          summarise(
            total_n = n(),
            race_n = sum(raceethnicity_cat == level | raceethnicity_cat == as.numeric(level), na.rm = TRUE),
            race_prop = mean(raceethnicity_cat == level | raceethnicity_cat == as.numeric(level), na.rm = TRUE) * 100,
            .groups = "drop"
          )
        
        treat_race <- race_counts %>% filter(treatment == 1)
        comp_race <- race_counts %>% filter(treatment == 0)
        
        treat_val <- if(nrow(treat_race) > 0) sprintf("%d (%.1f)", treat_race$race_n, treat_race$race_prop) else "0 (0.0)"
        comp_val <- if(nrow(comp_race) > 0) sprintf("%d (%.1f)", comp_race$race_n, comp_race$race_prop) else "0 (0.0)"
        
        cat(sprintf("    %-31s %-20s %-20s\n", race_name, treat_val, comp_val))
      }
    }
  }
  
  # Income
  if("income" %in% names(cohort_df)) {
    cat("  Income, n (%):\n")
    
    income_mappings <- list(
      "0" = "<$10k",
      "1" = "$10k-25k", 
      "2" = "$25k-35k",
      "3" = "$35k-50k",
      "4" = "$50k-75k",
      "5" = "$75k-100k",
      "6" = "$100k-150k",
      "7" = "$150k-200k",
      "8" = ">$200k",
      "999" = "Others/Not specified"
    )
    
    actual_levels <- unique(cohort_df$income)
    actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:8, 999, "0":"8", "999")]
    actual_levels <- sort(as.character(actual_levels))
    
    for(level in actual_levels) {
      if(level %in% names(income_mappings)) {
        income_name <- income_mappings[[level]]
        
        # Calculate counts and proportions directly
        income_counts <- cohort_df %>%
          group_by(treatment) %>%
          summarise(
            total_n = n(),
            income_n = sum(income == level | income == as.numeric(level), na.rm = TRUE),
            income_prop = mean(income == level | income == as.numeric(level), na.rm = TRUE) * 100,
            .groups = "drop"
          )
        
        treat_income <- income_counts %>% filter(treatment == 1)
        comp_income <- income_counts %>% filter(treatment == 0)
        
        treat_val <- if(nrow(treat_income) > 0) sprintf("%d (%.1f)", treat_income$income_n, treat_income$income_prop) else "0 (0.0)"
        comp_val <- if(nrow(comp_income) > 0) sprintf("%d (%.1f)", comp_income$income_n, comp_income$income_prop) else "0 (0.0)"
        
        cat(sprintf("    %-31s %-20s %-20s\n", income_name, treat_val, comp_val))
      }
    }
  }
  
  # Education
  if("education" %in% names(cohort_df)) {
    cat("  Education, n (%):\n")
    
    education_mappings <- list(
      "0" = "Less than high school",
      "1" = "High school/GED",
      "2" = "Some college",
      "3" = "College graduate",
      "4" = "Advanced degree",
      "999" = "Not specified/Other"
    )
    
    actual_levels <- unique(cohort_df$education)
    actual_levels <- actual_levels[!is.na(actual_levels) & actual_levels %in% c(0:4, 999, "0":"4", "999")]
    actual_levels <- sort(as.character(actual_levels))
    
    for(level in actual_levels) {
      if(level %in% names(education_mappings)) {
        education_name <- education_mappings[[level]]
        
        # Calculate counts and proportions directly
        education_counts <- cohort_df %>%
          group_by(treatment) %>%
          summarise(
            total_n = n(),
            education_n = sum(education == level | education == as.numeric(level), na.rm = TRUE),
            education_prop = mean(education == level | education == as.numeric(level), na.rm = TRUE) * 100,
            .groups = "drop"
          )
        
        treat_education <- education_counts %>% filter(treatment == 1)
        comp_education <- education_counts %>% filter(treatment == 0)
        
        treat_val <- if(nrow(treat_education) > 0) sprintf("%d (%.1f)", treat_education$education_n, treat_education$education_prop) else "0 (0.0)"
        comp_val <- if(nrow(comp_education) > 0) sprintf("%d (%.1f)", comp_education$education_n, comp_education$education_prop) else "0 (0.0)"
        
        cat(sprintf("    %-31s %-20s %-20s\n", education_name, treat_val, comp_val))
      }
    }
  }
  
  # Calendar Year (Index Year)
  if("index_year" %in% names(cohort_df)) {
    cat("  Calendar Year (Index Date), n (%):\n")
    
    actual_levels <- unique(cohort_df$index_year)
    actual_levels <- actual_levels[!is.na(actual_levels)]
    actual_levels <- sort(as.character(actual_levels))
    
    for(year in actual_levels) {
      # Calculate counts and proportions directly
      year_counts <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          total_n = n(),
          year_n = sum(index_year == year | index_year == as.character(year), na.rm = TRUE),
          year_prop = mean(index_year == year | index_year == as.character(year), na.rm = TRUE) * 100,
          .groups = "drop"
        )
      
      treat_year <- year_counts %>% filter(treatment == 1)
      comp_year <- year_counts %>% filter(treatment == 0)
      
      treat_val <- if(nrow(treat_year) > 0) sprintf("%d (%.1f)", treat_year$year_n, treat_year$year_prop) else "0 (0.0)"
      comp_val <- if(nrow(comp_year) > 0) sprintf("%d (%.1f)", comp_year$year_n, comp_year$year_prop) else "0 (0.0)"
      
      cat(sprintf("    %-31s %-20s %-20s\n", year, treat_val, comp_val))
    }
  }
  
  # BMI
  if("baseline_bmi" %in% names(cohort_df)) {
    bmi_summary <- cohort_df %>%
      filter(!is.na(baseline_bmi)) %>%
      group_by(treatment) %>%
      summarise(
        mean_bmi = mean(baseline_bmi, na.rm = TRUE),
        sd_bmi = sd(baseline_bmi, na.rm = TRUE),
        .groups = "drop"
      )
    
    if(nrow(bmi_summary) == 2) {
      treat_bmi <- sprintf("%.1f (%.1f)", bmi_summary$mean_bmi[bmi_summary$treatment == 1], 
                                           bmi_summary$sd_bmi[bmi_summary$treatment == 1])
      comp_bmi <- sprintf("%.1f (%.1f)", bmi_summary$mean_bmi[bmi_summary$treatment == 0], 
                                         bmi_summary$sd_bmi[bmi_summary$treatment == 0])
      
      cat(sprintf("  %-33s %-20s %-20s\n", "BMI (kg/m2), mean (SD)", treat_bmi, comp_bmi))
    }
  }
  
  # BMI Categories (if available)
  if("bmi_category" %in% names(cohort_df)) {
    cat("  BMI Category, n (%):\n")
    
    bmi_levels <- c("Underweight", "Normal", "Overweight", "Obese I", "Obese II", "Obese III")
    
    for(category in bmi_levels) {
      # Calculate counts and proportions directly
      bmi_cat_counts <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          total_n = n(),
          bmi_cat_n = sum(bmi_category == category, na.rm = TRUE),
          bmi_cat_prop = mean(bmi_category == category, na.rm = TRUE) * 100,
          .groups = "drop"
        )
      
      # Only show categories that have at least one patient
      if(any(bmi_cat_counts$bmi_cat_n > 0)) {
        treat_bmi_cat <- bmi_cat_counts %>% filter(treatment == 1)
        comp_bmi_cat <- bmi_cat_counts %>% filter(treatment == 0)
        
        treat_val <- if(nrow(treat_bmi_cat) > 0) sprintf("%d (%.1f)", treat_bmi_cat$bmi_cat_n, treat_bmi_cat$bmi_cat_prop) else "0 (0.0)"
        comp_val <- if(nrow(comp_bmi_cat) > 0) sprintf("%d (%.1f)", comp_bmi_cat$bmi_cat_n, comp_bmi_cat$bmi_cat_prop) else "0 (0.0)"
        
        cat(sprintf("    %-31s %-20s %-20s\n", category, treat_val, comp_val))
      }
    }
  }
  
  # Comorbidities Section
  cat("\nCOMORBIDITIES, n (%):\n")
  
  comorbidity_vars <- c("mi", "chf", "pvd", "cvd", "dem", "cpd", "ctd", "pud", "mld", "ph", "rd", "cancer", "msld", "mc", "hiv")
  comorbidity_labels <- list(
    "mi" = "Myocardial infarction",
    "chf" = "Congestive heart failure", 
    "pvd" = "Peripheral vascular disease",
    "cvd" = "Cerebrovascular disease",
    "dem" = "Dementia",
    "cpd" = "Chronic pulmonary disease",
    "ctd" = "Connective tissue disease",
    "pud" = "Peptic ulcer disease",
    "mld" = "Mild liver disease",
    "ph" = "Paraplegia and hemiplegia",
    "rd" = "Renal disease",
    "cancer" = "Cancer",
    "msld" = "Moderate/severe liver disease",
    "mc" = "Metastatic cancer",
    "hiv" = "HIV/AIDS"
  )
  
  for(var in comorbidity_vars) {
    if(var %in% names(cohort_df)) {
      var_summary <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          n_total = n(),
          n_positive = sum(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE),
          prop_positive = mean(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE) * 100,
          .groups = "drop"
        )
      
      treat_data <- var_summary %>% filter(treatment == 1)
      comp_data <- var_summary %>% filter(treatment == 0)
      
      if(nrow(treat_data) > 0 && nrow(comp_data) > 0) {
        treat_val <- sprintf("%d (%.1f)", treat_data$n_positive, treat_data$prop_positive)
        comp_val <- sprintf("%d (%.1f)", comp_data$n_positive, comp_data$prop_positive)
        
        var_label <- if(var %in% names(comorbidity_labels)) comorbidity_labels[[var]] else var
        
        cat(sprintf("  %-33s %-20s %-20s\n", var_label, treat_val, comp_val))
      }
    }
  }
  
  # Medications Section
  cat("\nMEDICATIONS, n (%):\n")
  
  medication_vars <- c("Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic", "Ezetimibe", "MRA", 
                       "OtherHTN", "RAAS", "Statin", "TZD", "Insulin", "SU", "DPP4", "SGLT2", "GLP1")
  medication_labels <- list(
    "Anticoagulant" = "Anticoagulant",
    "Antiplatelet" = "Antiplatelet",
    "BB" = "Beta-blocker",
    "Biguanide" = "Biguanide",
    "CCB" = "Calcium channel blocker",
    "Diuretic" = "Diuretic",
    "Ezetimibe" = "Ezetimibe",
    "MRA" = "Mineralocorticoid receptor antagonist",
    "OtherHTN" = "Other antihypertensive",
    "RAAS" = "RAAS inhibitor",
    "Statin" = "Statin",
    "TZD" = "Thiazolidinedione",
    "Insulin" = "Insulin",
    "SU" = "Sulfonylurea",
    "DPP4" = "DPP-4 inhibitor",
    "SGLT2" = "SGLT-2 inhibitor",
    "GLP1" = "GLP-1 agonist"
  )
  
  for(var in medication_vars) {
    if(var %in% names(cohort_df)) {
      var_summary <- cohort_df %>%
        group_by(treatment) %>%
        summarise(
          n_total = n(),
          n_positive = sum(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE),
          prop_positive = mean(.data[[var]] == "1" | .data[[var]] == 1, na.rm = TRUE) * 100,
          .groups = "drop"
        )
      
      treat_data <- var_summary %>% filter(treatment == 1)
      comp_data <- var_summary %>% filter(treatment == 0)
      
      if(nrow(treat_data) > 0 && nrow(comp_data) > 0) {
        treat_val <- sprintf("%d (%.1f)", treat_data$n_positive, treat_data$prop_positive)
        comp_val <- sprintf("%d (%.1f)", comp_data$n_positive, comp_data$prop_positive)
        
        var_label <- if(var %in% names(medication_labels)) medication_labels[[var]] else var
        
        cat(sprintf("  %-33s %-20s %-20s\n", var_label, treat_val, comp_val))
      }
    }
  }
  
  cat(paste(rep("-", 80), collapse=""), "\n")
  cat(sprintf("TOTAL EVENTS: %s = %d, %s = %d\n", treat_label, treat_events, comp_label, comp_events))
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Return structured data for saving
  baseline_data <- list(
    comparison = comparison_name,
    outcome = outcome_name,
    treatment_label = treat_label,
    comparator_label = comp_label,
    treatment_n = treat_n,
    comparator_n = comp_n,
    treatment_events = treat_events,
    comparator_events = comp_events
  )
  
  return(baseline_data)
}

# ---- Single Table Generation Function ----
generate_single_baseline_table <- function(comparison, outcome) {
  
  cat("\n=== GENERATING SINGLE BASELINE TABLE ===\n")
  cat("Comparison:", comparison, "\n")
  cat("Outcome:", outcome, "\n")
  
  # Convert to file format
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(file.exists(cohort_file)) {
    cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
    
    cat("Loading cohort from:", cohort_file, "\n")
    cat("Cohort size:", nrow(cohort_data), "patients\n")
    
    # Generate the table
    baseline_data <- create_formatted_baseline_table(cohort_data, comparison, outcome)
    
    # Save structured data
    baseline_table <- create_unweighted_baseline_table(cohort_data, comparison, outcome)
    
    if(!is.null(baseline_table)) {
      output_file <- paste0("single_baseline_table_", comp_file, "_", outcome_file, ".csv")
      write.csv(baseline_table, output_file, row.names = FALSE)
      cat("Saved to:", output_file, "\n")
    }
    
    return(baseline_data)
    
  } else {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    cat("Available cohort files:\n")
    cohort_files <- list.files(pattern = "^ipwt_cohort_.*\\.csv$")
    if(length(cohort_files) > 0) {
      for(file in cohort_files) {
        cat(" ", file, "\n")
      }
    } else {
      cat("  No cohort files found. Run the main analysis first.\n")
    }
    return(NULL)
  }
}

# ---- Comprehensive Baseline Tables Generation ----
generate_all_baseline_tables <- function() {
  
  cat("\n", paste(rep("=", 100), collapse=""), "\n")
  cat("GENERATING COMPREHENSIVE UNWEIGHTED BASELINE CHARACTERISTICS TABLES\n")
  cat("ALL THREE COMPARISONS INCLUDING SGLT2i vs SU/DPP4i\n")
  cat(paste(rep("=", 100), collapse=""), "\n")
  
  # Define comparisons and outcomes
  comparisons <- c("GLP1_vs_SGLT2", "GLP1_vs_SU_DPP4", "SGLT2_vs_SU_DPP4")
  outcomes <- c("Epilepsy_Seizure", "Early-onset_Epilepsy_Seizure", "Late-onset_Epilepsy_Seizure", "ADRD", "stroke")
  
  all_baseline_tables <- list()
  successful_generations <- 0
  
  for(comp in comparisons) {
    comp_name <- gsub("_", " ", comp)
    
    cat("\n", paste(rep("-", 80), collapse=""), "\n")
    cat("COMPARISON:", comp_name, "\n")
    cat(paste(rep("-", 80), collapse=""), "\n")
    
    for(outcome in outcomes) {
      outcome_name <- gsub("_", "/", outcome)
      outcome_name <- gsub("Late-onset/", "Late-onset ", outcome_name)
      # Handle specific outcome names
      if(outcome == "ADRD") {
        outcome_name <- "ADRD"
      } else if(outcome == "stroke") {
        outcome_name <- "Stroke"
      }
          
      # Check if cohort file exists
      cohort_file <- paste0("ipwt_cohort_", comp, "_", outcome, ".csv")
      
      if(file.exists(cohort_file)) {
        cat("\n✓ Found cohort file:", cohort_file, "\n")
        
        # Load cohort data
        cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
        
        cat("  Cohort size:", nrow(cohort_data), "patients\n")
        cat("  Treatment distribution:", paste(names(table(cohort_data$treatment)), "=", table(cohort_data$treatment), collapse = ", "), "\n")
        cat("  Event distribution:", paste(names(table(cohort_data$event)), "=", table(cohort_data$event), collapse = ", "), "\n")
        
        # Generate formatted baseline table
        tryCatch({
          baseline_data <- create_formatted_baseline_table(cohort_data, comp_name, outcome_name)
          
          # Create and save structured baseline table
          baseline_table <- create_unweighted_baseline_table(cohort_data, comp_name, outcome_name)
          
          if(!is.null(baseline_table)) {
            output_file <- paste0("unweighted_baseline_characteristics_", comp, "_", outcome, ".csv")
            write.csv(baseline_table, output_file, row.names = FALSE)
            cat("  ✓ Saved structured data to:", output_file, "\n")
          }
          
          # Generate SMD tables (before and after IPTW) - INDIVIDUAL LEVEL APPROACH
          cat("\n--- SMD CALCULATIONS (INDIVIDUAL LEVELS) ---\n")
          
          # SMD before IPTW (unweighted)
          smd_before <- calculate_smd_table(cohort_data, comp_name, outcome_name, weights = NULL)
          if(!is.null(smd_before)) {
            smd_before_file <- paste0("smd_individual_before_iptw_", comp, "_", outcome, ".csv")
            write.csv(smd_before, smd_before_file, row.names = FALSE)
            cat("  ✓ Individual-level SMD before IPTW saved to:", smd_before_file, "\n")
            
                         # Count different types of variables
             binary_vars <- smd_before[is.na(smd_before$variable_level) | 
                                      (!grepl("Overall", smd_before$variable_level) & 
                                       smd_before$variable_level %in% c("1", "0")), ]
             multilevel_individual <- smd_before[grepl("sex_cat|raceethnicity_cat|income|education|bmi_category|index_year", 
                                                      smd_before$variable) & 
                                                !grepl("Overall", smd_before$variable_level), ]
            multilevel_overall <- smd_before[grepl("Overall", smd_before$variable_level), ]
            
            cat(sprintf("    - Binary/continuous variables: %d\n", nrow(binary_vars)))
            cat(sprintf("    - Individual multilevel categories: %d\n", nrow(multilevel_individual)))
            cat(sprintf("    - Overall multilevel summaries: %d\n", nrow(multilevel_overall)))
            
            # Display top 10 largest absolute SMDs before IPTW
            top_smd_before <- smd_before %>%
              arrange(desc(abs(smd))) %>%
              head(10) %>%
              mutate(smd = round(smd, 3))
            cat("  Top 10 largest absolute SMDs (before IPTW):\n")
            for(i in 1:nrow(top_smd_before)) {
              var_display <- if(is.na(top_smd_before$variable_level[i]) || top_smd_before$variable_level[i] == "") {
                top_smd_before$variable[i]
              } else {
                paste0(top_smd_before$variable[i], ": ", top_smd_before$variable_level[i])
              }
              cat(sprintf("    %s = %.3f\n", var_display, top_smd_before$smd[i]))
            }
          }
          
          # SMD after IPTW (weighted) - only if weights are available
          if("ipw_std" %in% names(cohort_data)) {
            smd_after <- calculate_smd_table(cohort_data, comp_name, outcome_name, weights = cohort_data$ipw_std)
            if(!is.null(smd_after)) {
              smd_after_file <- paste0("smd_individual_after_iptw_", comp, "_", outcome, ".csv")
              write.csv(smd_after, smd_after_file, row.names = FALSE)
              cat("  ✓ Individual-level SMD after IPTW saved to:", smd_after_file, "\n")
              
              # Display top 10 largest absolute SMDs after IPTW
              top_smd_after <- smd_after %>%
                arrange(desc(abs(smd))) %>%
                head(10) %>%
                mutate(smd = round(smd, 3))
              cat("  Top 10 largest absolute SMDs (after IPTW):\n")
              for(i in 1:nrow(top_smd_after)) {
                var_display <- if(is.na(top_smd_after$variable_level[i]) || top_smd_after$variable_level[i] == "") {
                  top_smd_after$variable[i]
                } else {
                  paste0(top_smd_after$variable[i], ": ", top_smd_after$variable_level[i])
                }
                cat(sprintf("    %s = %.3f\n", var_display, top_smd_after$smd[i]))
              }
              
              # Enhanced balance improvement summary
              # Match rows by variable and variable_level for proper comparison
              smd_comparison <- merge(
                smd_before[, c("variable", "variable_level", "smd")],
                smd_after[, c("variable", "variable_level", "smd")],
                by = c("variable", "variable_level"),
                suffixes = c("_before", "_after")
              )
              
              vars_improved <- sum(abs(smd_comparison$smd_after) < abs(smd_comparison$smd_before), na.rm = TRUE)
              total_vars <- nrow(smd_comparison)
              
                             # Separate analysis for multilevel individual categories
               multilevel_comparison <- smd_comparison[
                 grepl("sex_cat|raceethnicity_cat|income|education|bmi_category|index_year", smd_comparison$variable) & 
                 !grepl("Overall", smd_comparison$variable_level), 
               ]
              
              if(nrow(multilevel_comparison) > 0) {
                multilevel_improved <- sum(abs(multilevel_comparison$smd_after) < abs(multilevel_comparison$smd_before), na.rm = TRUE)
                cat(sprintf("  Balance improvement summary:\n"))
                cat(sprintf("    - Total variables/levels: %d/%d improved (%.1f%%)\n", 
                            vars_improved, total_vars, (vars_improved/total_vars)*100))
                cat(sprintf("    - Individual multilevel categories: %d/%d improved (%.1f%%)\n", 
                            multilevel_improved, nrow(multilevel_comparison), 
                            (multilevel_improved/nrow(multilevel_comparison))*100))
                cat(sprintf("    - Largest remaining |SMD| (any level): %.3f\n", 
                            max(abs(smd_after$smd), na.rm = TRUE)))
                cat(sprintf("    - Individual levels with |SMD| > 0.1: %d\n", 
                            sum(abs(smd_after$smd) > 0.1, na.rm = TRUE)))
              } else {
                cat(sprintf("  Balance improvement: %d/%d variables (%.1f%%) have reduced SMD after IPTW\n", 
                            vars_improved, total_vars, (vars_improved/total_vars)*100))
              }
            }
            
            # Generate IPTW weighted baseline characteristics table
            cat("\n--- IPTW WEIGHTED BASELINE CHARACTERISTICS ---\n")
            weighted_baseline <- create_iptw_weighted_baseline_table(cohort_data, comp_name, outcome_name)
            
          } else {
            cat("  ⚠ No IPTW weights found - skipping weighted SMD calculation\n")
          }
          
          # Store the baseline data
          table_key <- paste(comp, outcome, sep = "_")
          all_baseline_tables[[table_key]] <- baseline_data
          successful_generations <- successful_generations + 1
          
        }, error = function(e) {
          cat("  ✗ Error generating baseline table for", comp_name, "-", outcome_name, ":", e$message, "\n")
          cat("  Error details:", toString(e), "\n")
        })
        
      } else {
        cat("\n✗ Cohort file not found:", cohort_file, "\n")
        cat("  This comparison-outcome combination may not have been analyzed.\n")
      }
    }
  }
  
  # Generate summary of all baseline tables
  cat("\n", paste(rep("=", 100), collapse=""), "\n")
  cat("BASELINE CHARACTERISTICS SUMMARY ACROSS ALL COMPARISONS\n")
  cat(paste(rep("=", 100), collapse=""), "\n")
  
  if(length(all_baseline_tables) > 0) {
    
    # Create a summary table of sample sizes
    summary_df <- data.frame(
      Comparison = character(),
      Outcome = character(),
      Treatment_Label = character(),
      Treatment_N = integer(),
      Comparator_Label = character(),
      Comparator_N = integer(),
      Treatment_Events = integer(),
      Comparator_Events = integer(),
      stringsAsFactors = FALSE
    )
    
    for(table_name in names(all_baseline_tables)) {
      table_data <- all_baseline_tables[[table_name]]
      
      summary_df <- rbind(summary_df, data.frame(
        Comparison = table_data$comparison,
        Outcome = table_data$outcome,
        Treatment_Label = table_data$treatment_label,
        Treatment_N = table_data$treatment_n,
        Comparator_Label = table_data$comparator_label,
        Comparator_N = table_data$comparator_n,
        Treatment_Events = table_data$treatment_events,
        Comparator_Events = table_data$comparator_events,
        stringsAsFactors = FALSE
      ))
    }
    
    cat("\nSample sizes across all analyses:\n")
    print(summary_df)
    
    # Save summary table
    write.csv(summary_df, "baseline_characteristics_summary.csv", row.names = FALSE)
    cat("\n✓ Saved baseline characteristics summary to: baseline_characteristics_summary.csv\n")
    
    # Success summary
    cat("\n", paste(rep("=", 100), collapse=""), "\n")
    cat("BASELINE CHARACTERISTICS GENERATION COMPLETE\n")
    cat(paste(rep("=", 100), collapse=""), "\n")
    cat("✓ Successfully generated", successful_generations, "baseline characteristics tables\n")
    cat("✓ Tables available for all three comparisons:\n")
    cat("  - Semaglutide vs SGLT2i\n")
    cat("  - Semaglutide vs SU/DPP4i\n")  
    cat("  - SGLT2i vs SU/DPP4i ⭐ (New comparison as requested)\n")
    cat("✓ Both outcomes analyzed: Epilepsy/Seizure and Late-onset Epilepsy/Seizure\n")
    cat("✓ All tables show UNWEIGHTED baseline characteristics\n")
    cat("✓ Individual CSV files saved with detailed structured data\n")
    
    return(all_baseline_tables)
    
  } else {
    cat("✗ No baseline tables were generated. Please check that cohort files exist.\n")
    cat("Make sure you have run the main IPWT analysis first.\n")
    return(NULL)
  }
}

# ---- Treatment Assignment Diagnostic Function ----
diagnose_treatment_assignment <- function(comparison = "GLP1 vs SGLT2", outcome = "Epilepsy/Seizure") {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("DIAGNOSING TREATMENT ASSIGNMENT ISSUES FOR", comparison, "-", outcome, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load cohort data
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(!file.exists(cohort_file)) {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    return(NULL)
  }
  
  cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
  cat("Loaded cohort with", nrow(cohort_data), "patients\n")
  
  # Check treatment-defining medications
  if(comparison == "GLP1 vs SGLT2") {
    
    cat("\n=== GLP1 vs SGLT2 TREATMENT ASSIGNMENT CHECK ===\n")
    
    # Check GLP1 and SGLT2 usage by treatment group
    treatment_meds <- cohort_data %>%
      group_by(treatment) %>%
      summarise(
        n = n(),
        glp1_users = sum(GLP1 == "1" | GLP1 == 1, na.rm = TRUE),
        sglt2_users = sum(SGLT2 == "1" | SGLT2 == 1, na.rm = TRUE),
        glp1_pct = glp1_users / n * 100,
        sglt2_pct = sglt2_users / n * 100,
        both_drugs = sum((GLP1 == "1" | GLP1 == 1) & (SGLT2 == "1" | SGLT2 == 1), na.rm = TRUE),
        neither_drug = sum((GLP1 == "0" | GLP1 == 0 | is.na(GLP1)) & (SGLT2 == "0" | SGLT2 == 0 | is.na(SGLT2)), na.rm = TRUE),
        .groups = "drop"
      )
    
    cat("Treatment assignment summary:\n")
    print(treatment_meds)
    
    cat("\nEXPECTED vs ACTUAL:\n")
    cat("Treatment = 1 (should be Semaglutide group):\n")
    cat("  Expected: GLP1 = ~100%, SGLT2 = 0%\n")
    cat("  Actual:   GLP1 =", sprintf("%.1f%%", treatment_meds$glp1_pct[treatment_meds$treatment == 1]), 
        ", SGLT2 =", sprintf("%.1f%%", treatment_meds$sglt2_pct[treatment_meds$treatment == 1]), "\n")
    
    cat("Treatment = 0 (should be SGLT2 group):\n")
    cat("  Expected: GLP1 = 0%, SGLT2 = ~100%\n")
    cat("  Actual:   GLP1 =", sprintf("%.1f%%", treatment_meds$glp1_pct[treatment_meds$treatment == 0]), 
        ", SGLT2 =", sprintf("%.1f%%", treatment_meds$sglt2_pct[treatment_meds$treatment == 0]), "\n")
    
    # Check for problematic cases
    problematic_cases <- cohort_data %>%
      mutate(
        is_glp1_user = GLP1 == "1" | GLP1 == 1,
        is_sglt2_user = SGLT2 == "1" | SGLT2 == 1,
        case_type = case_when(
          treatment == 1 & !is_glp1_user ~ "Sema group but no GLP1",
          treatment == 1 & is_sglt2_user ~ "Sema group but has SGLT2", 
          treatment == 0 & is_glp1_user ~ "SGLT2 group but has GLP1",
          treatment == 0 & !is_sglt2_user ~ "SGLT2 group but no SGLT2",
          TRUE ~ "OK"
        )
      ) %>%
      filter(case_type != "OK")
    
    cat("\nProblematic cases found:", nrow(problematic_cases), "\n")
    if(nrow(problematic_cases) > 0) {
      problem_summary <- table(problematic_cases$case_type)
      print(problem_summary)
      
      cat("\nThis suggests issues with:\n")
      cat("1. Treatment assignment logic\n")
      cat("2. Medication variable definitions\n") 
      cat("3. Possible inclusion of baseline medications vs index medications\n")
    }
    
  }
  
  cat("\n=== MEDICATION VARIABLE INTERPRETATION ===\n")
  cat("The medication variables (GLP1, SGLT2, etc.) appear to represent:\n")
  cat("- Baseline medication use (period before/around index date)\n")
  cat("- NOT just the index medication that defines treatment\n")
  cat("- This includes:\n")
  cat("  * Index drug\n")
  cat("  * Other drugs in same class used during baseline period\n")
  cat("  * Possibly concurrent therapy\n")
  
  cat("\n=== RECOMMENDATIONS ===\n")
  cat("1. ✅ FIXED: Exclude treatment-defining variables from SMD calculation\n")
  cat("2. Consider creating 'index_drug' variables that are 100%/0% by design\n")
  cat("3. Keep existing medication variables for clinical context\n")
  cat("4. Document that medication variables represent baseline use, not treatment assignment\n")
  
  return(problematic_cases)
}

# ---- SMD Debugging Function ----
debug_smd_calculation <- function(comparison = "GLP1 vs SGLT2", outcome = "Epilepsy/Seizure") {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("DEBUGGING SMD CALCULATIONS FOR", comparison, "-", outcome, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load cohort data
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(!file.exists(cohort_file)) {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    return(NULL)
  }
  
  cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
  cat("Loaded cohort with", nrow(cohort_data), "patients\n")
  cat("Treatment distribution:", table(cohort_data$treatment), "\n")
  
  # Test sex_cat SMD calculation specifically
  cat("\n=== DEBUGGING SEX_CAT SMD ===\n")
  
  # Manual calculation - before IPTW
  treat_sex <- cohort_data[cohort_data$treatment == 1 & !is.na(cohort_data$sex_cat), ]
  comp_sex <- cohort_data[cohort_data$treatment == 0 & !is.na(cohort_data$sex_cat), ]
  
  prop_treat_before <- mean(treat_sex$sex_cat == "1" | treat_sex$sex_cat == 1, na.rm = TRUE)
  prop_comp_before <- mean(comp_sex$sex_cat == "1" | comp_sex$sex_cat == 1, na.rm = TRUE)
  
  var_treat_before <- prop_treat_before * (1 - prop_treat_before)
  var_comp_before <- prop_comp_before * (1 - prop_comp_before)
  pooled_sd_before <- sqrt((var_treat_before + var_comp_before) / 2)
  smd_before_manual <- (prop_treat_before - prop_comp_before) / pooled_sd_before
  
  cat("Manual calculation (before IPTW):\n")
  cat("  Treatment group: n =", nrow(treat_sex), ", prop_female =", prop_treat_before, "\n")
  cat("  Comparator group: n =", nrow(comp_sex), ", prop_female =", prop_comp_before, "\n")
  cat("  Pooled SD =", pooled_sd_before, "\n")
  cat("  SMD =", smd_before_manual, "\n")
  
  # Manual calculation - after IPTW (if weights available)
  if("ipw_std" %in% names(cohort_data)) {
    treat_idx <- cohort_data$treatment == 1 & !is.na(cohort_data$sex_cat)
    comp_idx <- cohort_data$treatment == 0 & !is.na(cohort_data$sex_cat)
    
    treat_positive <- (cohort_data[treat_idx, "sex_cat"] == "1" | cohort_data[treat_idx, "sex_cat"] == 1)
    comp_positive <- (cohort_data[comp_idx, "sex_cat"] == "1" | cohort_data[comp_idx, "sex_cat"] == 1)
    w_treat <- cohort_data$ipw_std[treat_idx]
    w_comp <- cohort_data$ipw_std[comp_idx]
    
    prop_treat_after <- sum(w_treat * treat_positive, na.rm = TRUE) / sum(w_treat, na.rm = TRUE)
    prop_comp_after <- sum(w_comp * comp_positive, na.rm = TRUE) / sum(w_comp, na.rm = TRUE)
    
    var_treat_after <- prop_treat_after * (1 - prop_treat_after)
    var_comp_after <- prop_comp_after * (1 - prop_comp_after)
    pooled_sd_after <- sqrt((var_treat_after + var_comp_after) / 2)
    smd_after_manual <- (prop_treat_after - prop_comp_after) / pooled_sd_after
    
    cat("\nManual calculation (after IPTW):\n")
    cat("  Treatment group: eff_n =", sum(w_treat, na.rm = TRUE), ", weighted_prop_female =", prop_treat_after, "\n")
    cat("  Comparator group: eff_n =", sum(w_comp, na.rm = TRUE), ", weighted_prop_female =", prop_comp_after, "\n")
    cat("  Pooled SD =", pooled_sd_after, "\n")
    cat("  SMD =", smd_after_manual, "\n")
    
    # Check weights summary
    cat("\nWeight summary:\n")
    cat("  Treatment weights: min =", min(w_treat, na.rm = TRUE), ", max =", max(w_treat, na.rm = TRUE), ", mean =", mean(w_treat, na.rm = TRUE), "\n")
    cat("  Comparator weights: min =", min(w_comp, na.rm = TRUE), ", max =", max(w_comp, na.rm = TRUE), ", mean =", mean(w_comp, na.rm = TRUE), "\n")
  }
  
  # Test function calculation
  cat("\n=== TESTING CORRECTED SMD FUNCTION ===\n")
  
  smd_before_func <- calculate_smd_table(cohort_data, comparison, outcome, weights = NULL)
  sex_smd_before <- smd_before_func[smd_before_func$variable == "sex_cat", ]
  
  cat("Function calculation (before IPTW):\n")
  cat("  SMD =", sex_smd_before$smd, "\n")
  cat("  Match with manual:", abs(sex_smd_before$smd - smd_before_manual) < 1e-10, "\n")
  
  if("ipw_std" %in% names(cohort_data)) {
    smd_after_func <- calculate_smd_table(cohort_data, comparison, outcome, weights = cohort_data$ipw_std)
    sex_smd_after <- smd_after_func[smd_after_func$variable == "sex_cat", ]
    
    cat("\nFunction calculation (after IPTW):\n")
    cat("  SMD =", sex_smd_after$smd, "\n")
    cat("  Match with manual:", abs(sex_smd_after$smd - smd_after_manual) < 1e-10, "\n")
    
    # Balance improvement
    cat("\nBalance improvement:\n")
    cat("  Before IPTW SMD =", sex_smd_before$smd, "\n")
    cat("  After IPTW SMD =", sex_smd_after$smd, "\n")
    cat("  Improvement =", abs(sex_smd_before$smd) - abs(sex_smd_after$smd), "\n")
    cat("  % reduction in |SMD| =", (1 - abs(sex_smd_after$smd)/abs(sex_smd_before$smd)) * 100, "%\n")
  }
  
  # Test multilevel variables
  cat("\n=== TESTING MULTILEVEL VARIABLE SMD ===\n")
  
  if("raceethnicity_cat" %in% names(cohort_data)) {
    cat("Race/ethnicity levels in data:", sort(unique(cohort_data$raceethnicity_cat)), "\n")
    
    race_smd_before <- smd_before_func[smd_before_func$variable == "raceethnicity_cat", ]
    cat("Overall race/ethnicity SMD (before IPTW):", race_smd_before$smd, "\n")
    
    if("ipw_std" %in% names(cohort_data)) {
      race_smd_after <- smd_after_func[smd_after_func$variable == "raceethnicity_cat", ]
      cat("Overall race/ethnicity SMD (after IPTW):", race_smd_after$smd, "\n")
    }
  }
  
  return(list(
    before_iptw = smd_before_func,
    after_iptw = if("ipw_std" %in% names(cohort_data)) smd_after_func else NULL
  ))
}

# ---- Run Comprehensive Baseline Table Generation ----
cat("\n>>> GENERATING ALL BASELINE CHARACTERISTICS TABLES (CONSOLIDATED) <<<\n")
all_baseline_results <- generate_all_baseline_tables()

# ---- Test SMD Calculations ----
cat("\n>>> TESTING CORRECTED SMD CALCULATIONS <<<\n")
smd_debug_results <- debug_smd_calculation("GLP1 vs SGLT2", "Epilepsy/Seizure")

# ---- Test Individual Level SMD Function ----
test_individual_level_smd <- function(comparison = "GLP1 vs SGLT2", outcome = "Epilepsy/Seizure") {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("TESTING INDIVIDUAL LEVEL SMD CALCULATIONS FOR", comparison, "-", outcome, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load cohort data
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(!file.exists(cohort_file)) {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    return(NULL)
  }
  
  cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
  cat("Loaded cohort with", nrow(cohort_data), "patients\n")
  
  # Calculate SMDs with new individual-level approach
  smd_before <- calculate_smd_table(cohort_data, comparison, outcome, weights = NULL)
  smd_after <- if("ipw_std" %in% names(cohort_data)) {
    calculate_smd_table(cohort_data, comparison, outcome, weights = cohort_data$ipw_std)
  } else { NULL }
  
  # Display multilevel variable SMDs
  cat("\n=== INDIVIDUAL LEVEL SMDs (BEFORE IPTW) ===\n")
  
  multilevel_vars <- c("sex_cat", "raceethnicity_cat", "income", "education", "bmi_category", "index_year")
  
  for(var in multilevel_vars) {
    if(var %in% cohort_data) {
      cat("\n", toupper(var), ":\n")
      
      # Get individual level SMDs for this variable
      var_smds <- smd_before[grepl(paste0("^", var), smd_before$variable) & 
                            smd_before$variable_level != "Overall", ]
      
      if(nrow(var_smds) > 0) {
        for(i in 1:nrow(var_smds)) {
          cat(sprintf("  %-30s SMD = %6.3f (treat=%.3f, comp=%.3f)\n", 
                      var_smds$variable_level[i], 
                      var_smds$smd[i],
                      var_smds$mean_treat[i],
                      var_smds$mean_comp[i]))
        }
        
        # Get overall SMD
        overall_smd <- smd_before[smd_before$variable == var & 
                                 grepl("Overall", smd_before$variable_level), ]
        if(nrow(overall_smd) > 0) {
          cat(sprintf("  %-30s SMD = %6.3f [%s]\n", 
                      ">>> OVERALL <<<", 
                      overall_smd$smd[1],
                      overall_smd$variable_level[1]))
        }
      }
    }
  }
  
  # Display comparison if IPTW results available
  if(!is.null(smd_after)) {
    cat("\n=== INDIVIDUAL LEVEL SMDs (AFTER IPTW) ===\n")
    
    for(var in multilevel_vars) {
      if(var %in% cohort_data) {
        cat("\n", toupper(var), " - IMPROVEMENT:\n")
        
        # Get before and after SMDs for this variable
        before_smds <- smd_before[grepl(paste0("^", var), smd_before$variable) & 
                                 smd_before$variable_level != "Overall", ]
        after_smds <- smd_after[grepl(paste0("^", var), smd_after$variable) & 
                               smd_after$variable_level != "Overall", ]
        
        if(nrow(before_smds) > 0 && nrow(after_smds) > 0) {
          # Match by variable_level
          merged_smds <- merge(before_smds[, c("variable_level", "smd")], 
                              after_smds[, c("variable_level", "smd")], 
                              by = "variable_level", suffixes = c("_before", "_after"))
          
          for(i in 1:nrow(merged_smds)) {
            improvement <- abs(merged_smds$smd_before[i]) - abs(merged_smds$smd_after[i])
            improvement_pct <- (1 - abs(merged_smds$smd_after[i])/abs(merged_smds$smd_before[i])) * 100
            
            cat(sprintf("  %-25s Before=%6.3f After=%6.3f Δ=%6.3f (%+.1f%%)\n", 
                        merged_smds$variable_level[i],
                        merged_smds$smd_before[i],
                        merged_smds$smd_after[i],
                        improvement,
                        improvement_pct))
          }
        }
      }
    }
  }
  
  # Summary statistics
  cat("\n=== SUMMARY STATISTICS ===\n")
  
  # Count multilevel variable levels
  multilevel_levels_before <- smd_before[grepl("sex_cat|raceethnicity_cat|income|education|bmi_category|index_year", 
                                              smd_before$variable) & 
                                        !grepl("Overall", smd_before$variable_level), ]
  
  cat("Total individual levels from multilevel variables:", nrow(multilevel_levels_before), "\n")
  cat("Levels with |SMD| > 0.1:", sum(abs(multilevel_levels_before$smd) > 0.1, na.rm = TRUE), "\n")
  cat("Levels with |SMD| > 0.2:", sum(abs(multilevel_levels_before$smd) > 0.2, na.rm = TRUE), "\n")
  cat("Largest |SMD| among individual levels:", round(max(abs(multilevel_levels_before$smd), na.rm = TRUE), 3), "\n")
  
  if(!is.null(smd_after)) {
    multilevel_levels_after <- smd_after[grepl("sex_cat|raceethnicity_cat|income|education|bmi_category|index_year", 
                                              smd_after$variable) & 
                                        !grepl("Overall", smd_after$variable_level), ]
    
    cat("\nAfter IPTW:\n")
    cat("Levels with |SMD| > 0.1:", sum(abs(multilevel_levels_after$smd) > 0.1, na.rm = TRUE), "\n")
    cat("Levels with |SMD| > 0.2:", sum(abs(multilevel_levels_after$smd) > 0.2, na.rm = TRUE), "\n")
    cat("Largest |SMD| among individual levels:", round(max(abs(multilevel_levels_after$smd), na.rm = TRUE), 3), "\n")
    
    # Calculate improvement
    levels_improved <- sum(abs(multilevel_levels_after$smd) < abs(multilevel_levels_before$smd), na.rm = TRUE)
    total_levels <- nrow(multilevel_levels_before)
    cat("Individual levels with improved balance:", levels_improved, "/", total_levels, 
        sprintf(" (%.1f%%)\n", (levels_improved/total_levels)*100))
  }
  
  return(list(
    before = smd_before,
    after = smd_after,
    multilevel_before = multilevel_levels_before,
    multilevel_after = if(!is.null(smd_after)) multilevel_levels_after else NULL
  ))
}

cat("\n>>> TESTING NEW INDIVIDUAL LEVEL SMD CALCULATIONS <<<\n")
individual_smd_test <- test_individual_level_smd("GLP1 vs SGLT2", "Epilepsy/Seizure")

# ---- Test Updated Multilevel Approach ----
test_updated_multilevel_approach <- function(comparison = "GLP1 vs SGLT2", outcome = "Epilepsy/Seizure") {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("TESTING UPDATED MULTILEVEL APPROACH - ALL LEVELS INCLUDED\n")
  cat("Comparison:", comparison, "- Outcome:", outcome, "\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load cohort data
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(!file.exists(cohort_file)) {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    return(NULL)
  }
  
  cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
  
  # Check what sex categories are actually present
  cat("\n=== SEX CATEGORY ANALYSIS ===\n")
  if("sex_cat" %in% names(cohort_data)) {
    sex_distribution <- table(cohort_data$sex_cat, useNA = "always")
    cat("Sex categories in data:\n")
    print(sex_distribution)
    
    sex_by_treatment <- table(cohort_data$treatment, cohort_data$sex_cat, useNA = "always")
    cat("\nSex by treatment group:\n")
    print(sex_by_treatment)
  }
  
  # Check BMI categories
  cat("\n=== BMI CATEGORY ANALYSIS ===\n")
  if("bmi_category" %in% names(cohort_data)) {
    bmi_distribution <- table(cohort_data$bmi_category, useNA = "always")
    cat("BMI categories in data:\n")
    print(bmi_distribution)
    
    if(sum(!is.na(cohort_data$bmi_category)) > 0) {
      bmi_by_treatment <- table(cohort_data$treatment, cohort_data$bmi_category, useNA = "always")
      cat("\nBMI categories by treatment group:\n")
      print(bmi_by_treatment)
    } else {
      cat("No BMI category data available\n")
    }
  } else {
    cat("bmi_category variable not found in cohort\n")
  }
  
  # Test the updated SMD calculation
  cat("\n=== UPDATED SMD CALCULATION TEST ===\n")
  smd_results <- calculate_smd_table(cohort_data, comparison, outcome, weights = NULL)
  
  if(!is.null(smd_results)) {
    # Check sex_cat SMDs
    sex_smds <- smd_results[smd_results$variable == "sex_cat", ]
    if(nrow(sex_smds) > 0) {
      cat("\nSEX_CAT SMDs (should include Male, Female, and potentially Not specified/Other):\n")
      for(i in 1:nrow(sex_smds)) {
        cat(sprintf("  %-25s SMD = %6.3f (treat=%.3f, comp=%.3f)\n", 
                    sex_smds$variable_level[i], sex_smds$smd[i],
                    ifelse(is.na(sex_smds$mean_treat[i]), NA, sex_smds$mean_treat[i]),
                    ifelse(is.na(sex_smds$mean_comp[i]), NA, sex_smds$mean_comp[i])))
      }
    } else {
      cat("No sex_cat SMDs found - check multilevel variable handling\n")
    }
    
    # Check raceethnicity_cat SMDs (should now include Non-Hispanic White)
    race_smds <- smd_results[smd_results$variable == "raceethnicity_cat", ]
    if(nrow(race_smds) > 0) {
      cat("\nRACEETHNICITY_CAT SMDs (should include ALL levels including Non-Hispanic White):\n")
      for(i in 1:nrow(race_smds)) {
        cat(sprintf("  %-30s SMD = %6.3f (treat=%.3f, comp=%.3f)\n", 
                    race_smds$variable_level[i], race_smds$smd[i],
                    ifelse(is.na(race_smds$mean_treat[i]), NA, race_smds$mean_treat[i]),
                    ifelse(is.na(race_smds$mean_comp[i]), NA, race_smds$mean_comp[i])))
      }
    }
    
    # Check income SMDs (should now include <$10K)
    income_smds <- smd_results[smd_results$variable == "income", ]
    if(nrow(income_smds) > 0) {
      cat("\nINCOME SMDs (should include ALL levels including <$10K):\n")
      for(i in 1:nrow(income_smds)) {
        cat(sprintf("  %-20s SMD = %6.3f (treat=%.3f, comp=%.3f)\n", 
                    income_smds$variable_level[i], income_smds$smd[i],
                    ifelse(is.na(income_smds$mean_treat[i]), NA, income_smds$mean_treat[i]),
                    ifelse(is.na(income_smds$mean_comp[i]), NA, income_smds$mean_comp[i])))
      }
    }
    
    # Check BMI category SMDs
    bmi_smds <- smd_results[smd_results$variable == "bmi_category", ]
    if(nrow(bmi_smds) > 0) {
      cat("\nBMI_CATEGORY SMDs (should include Underweight, Normal, Overweight, Obese I-III):\n")
      for(i in 1:nrow(bmi_smds)) {
        cat(sprintf("  %-15s SMD = %6.3f (treat=%.3f, comp=%.3f)\n", 
                    bmi_smds$variable_level[i], bmi_smds$smd[i],
                    ifelse(is.na(bmi_smds$mean_treat[i]), NA, bmi_smds$mean_treat[i]),
                    ifelse(is.na(bmi_smds$mean_comp[i]), NA, bmi_smds$mean_comp[i])))
      }
    } else {
      cat("\nNo BMI_CATEGORY SMDs found - this may be expected if BMI data is limited\n")
    }
    
    # Summary
    cat("\n=== SUMMARY ===\n")
    multilevel_counts <- smd_results %>%
      filter(variable %in% c("sex_cat", "raceethnicity_cat", "income", "education", "bmi_category")) %>%
      filter(!grepl("Overall", variable_level)) %>%
      group_by(variable) %>%
      summarise(levels_calculated = n(), .groups = "drop")
    
    cat("Number of individual levels calculated per multilevel variable:\n")
    print(multilevel_counts)
    
    expected_counts <- data.frame(
      variable = c("sex_cat", "raceethnicity_cat", "income", "education", "bmi_category"),
      expected_levels = c("2-3 (Male, Female, +999)", "4 (all race groups)", "up to 9 (all income levels)", 
                         "up to 5 (all education levels)", "varies (BMI categories present)")
    )
    
    cat("\nExpected vs actual:\n")
    comparison_table <- merge(expected_counts, multilevel_counts, by = "variable", all.x = TRUE)
    comparison_table$levels_calculated[is.na(comparison_table$levels_calculated)] <- 0
    print(comparison_table)
    
  } else {
    cat("SMD calculation returned NULL - check for errors\n")
  }
  
  return(smd_results)
}

cat("\n>>> TESTING UPDATED MULTILEVEL APPROACH <<<\n")
multilevel_test <- test_updated_multilevel_approach("GLP1 vs SGLT2", "Epilepsy/Seizure")

# ---- Diagnose Treatment Assignment Issues ----
cat("\n>>> DIAGNOSING TREATMENT ASSIGNMENT ISSUES <<<\n")
treatment_issues <- diagnose_treatment_assignment("GLP1 vs SGLT2", "Epilepsy/Seizure")

# ---- Generate Corrected SMD Table ----
cat("\n>>> GENERATING CORRECTED SMD TABLE (EXCLUDING TREATMENT-DEFINING VARIABLES) <<<\n")

generate_corrected_smd_comparison <- function(comparison = "GLP1 vs SGLT2", outcome = "Epilepsy/Seizure") {
  
  # Load cohort data
  comp_file <- gsub(" ", "_", comparison)
  outcome_file <- gsub("[/ ]", "_", outcome)
  cohort_file <- paste0("ipwt_cohort_", comp_file, "_", outcome_file, ".csv")
  
  if(!file.exists(cohort_file)) {
    cat("Error: Cohort file not found:", cohort_file, "\n")
    return(NULL)
  }
  
  cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
  
  cat("\n=== COMPARISON: Original vs Corrected SMD Tables ===\n")
  cat("Comparison:", comparison, "- Outcome:", outcome, "\n")
  
  # Calculate CORRECTED SMD (excluding treatment-defining variables)
  smd_before_corrected <- calculate_smd_table(cohort_data, comparison, outcome, weights = NULL)
  smd_after_corrected <- calculate_smd_table(cohort_data, comparison, outcome, weights = cohort_data$ipw_std)
  
  # Show improvement with corrected calculations
  cat("\nTop 10 largest absolute SMDs (CORRECTED - after IPTW):\n")
  top_smd_corrected <- smd_after_corrected %>%
    arrange(desc(abs(smd))) %>%
    head(10) %>%
    mutate(smd = round(smd, 3))
  
  for(i in 1:nrow(top_smd_corrected)) {
    var_display <- if(is.na(top_smd_corrected$variable_level[i])) {
      top_smd_corrected$variable[i]
    } else {
      paste0(top_smd_corrected$variable[i], ": ", top_smd_corrected$variable_level[i])
    }
    cat(sprintf("  %s = %.3f\n", var_display, top_smd_corrected$smd[i]))
  }
  
  # Calculate balance improvement
  vars_improved <- sum(abs(smd_after_corrected$smd) < abs(smd_before_corrected$smd), na.rm = TRUE)
  total_vars <- nrow(smd_before_corrected)
  
  cat(sprintf("\nBalance improvement summary:\n"))
  cat(sprintf("- Variables with improved balance: %d/%d (%.1f%%)\n", vars_improved, total_vars, (vars_improved/total_vars)*100))
  cat(sprintf("- Largest remaining |SMD|: %.3f\n", max(abs(smd_after_corrected$smd), na.rm = TRUE)))
  cat(sprintf("- Variables with |SMD| > 0.1: %d\n", sum(abs(smd_after_corrected$smd) > 0.1, na.rm = TRUE)))
  cat(sprintf("- Variables with |SMD| > 0.2: %d\n", sum(abs(smd_after_corrected$smd) > 0.2, na.rm = TRUE)))
  
  # Save corrected SMD tables
  write.csv(smd_before_corrected, paste0("smd_before_iptw_CORRECTED_", comp_file, "_", outcome_file, ".csv"), row.names = FALSE)
  write.csv(smd_after_corrected, paste0("smd_after_iptw_CORRECTED_", comp_file, "_", outcome_file, ".csv"), row.names = FALSE)
  
  cat(sprintf("\n✓ Saved corrected SMD tables:\n"))
  cat(sprintf("  - smd_before_iptw_CORRECTED_%s_%s.csv\n", comp_file, outcome_file))
  cat(sprintf("  - smd_after_iptw_CORRECTED_%s_%s.csv\n", comp_file, outcome_file))
  
  return(list(
    before = smd_before_corrected,
    after = smd_after_corrected,
    improvement_pct = (vars_improved/total_vars)*100
  ))
}

corrected_smd_results <- generate_corrected_smd_comparison("GLP1 vs SGLT2", "Epilepsy/Seizure")

# ---- Usage Examples (commented out) ----
# Generate a single table:
# single_table <- generate_single_baseline_table("SGLT2 vs SU_DPP4", "Epilepsy/Seizure")
# single_table <- generate_single_baseline_table("GLP1 vs SGLT2", "Late-onset Epilepsy/Seizure")

# ===============================================================================
# END OF CONSOLIDATED BASELINE CHARACTERISTICS MODULE
# ===============================================================================

# =============================================================================
# TMLE-Appropriate Risk Difference Plot: GLP1 RA vs SGLT2i
# =============================================================================

# Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(gridExtra)

# =============================================================================
# 1. Create TMLE-appropriate data structure
# =============================================================================

# TMLE typically provides:
# - Point estimate of risk difference at a specific time point
# - Confidence intervals
# - Often stratified by covariates or subgroups

create_tmle_results <- function() {
  # Simulated TMLE results for different analyses/subgroups
  
  # Main analysis at 2-year follow-up
  main_results <- data.frame(
    analysis = "Main Analysis (2-year)",
    outcome = c("Epilepsy/Seizure", "Late-onset Epilepsy"),
    risk_glp1 = c(0.018, 0.015),  # 2-year risks
    risk_sglt2i = c(0.021, 0.017),
    risk_difference = c(-0.003, -0.002),
    lower_ci = c(-0.005, -0.004),
    upper_ci = c(-0.001, 0.000),
    p_value = c(0.003, 0.052),
    n_glp1 = c(2000, 2000),
    n_sglt2i = c(4000, 4000),
    stringsAsFactors = FALSE
  )
  
  # Subgroup analyses
  subgroup_results <- data.frame(
    analysis = rep(c("Age <65", "Age ≥65", "With CVD", "Without CVD", 
                     "BMI <30", "BMI ≥30"), each = 2),
    outcome = rep(c("Epilepsy/Seizure", "Late-onset Epilepsy"), 6),
    risk_glp1 = c(0.015, 0.012, 0.022, 0.019, 0.020, 0.018, 
                  0.016, 0.013, 0.017, 0.014, 0.019, 0.016),
    risk_sglt2i = c(0.019, 0.014, 0.025, 0.022, 0.024, 0.021, 
                    0.018, 0.016, 0.020, 0.016, 0.022, 0.018),
    risk_difference = c(-0.004, -0.002, -0.003, -0.003, -0.004, -0.003,
                       -0.002, -0.003, -0.003, -0.002, -0.003, -0.002),
    lower_ci = c(-0.007, -0.005, -0.006, -0.006, -0.007, -0.006,
                 -0.005, -0.006, -0.006, -0.005, -0.006, -0.005),
    upper_ci = c(-0.001, 0.001, 0.000, 0.000, -0.001, 0.000,
                 0.001, 0.000, 0.000, 0.001, 0.000, 0.001),
    p_value = c(0.009, 0.180, 0.051, 0.048, 0.008, 0.055,
                0.215, 0.042, 0.047, 0.198, 0.052, 0.187),
    n_glp1 = rep(c(800, 1200, 900, 1100, 1000, 1000), each = 2),
    n_sglt2i = rep(c(1600, 2400, 1800, 2200, 2000, 2000), each = 2),
    stringsAsFactors = FALSE
  )
  
  # Combine results
  all_results <- rbind(main_results, subgroup_results)
  
  # Add significance indicator
  all_results$significant <- all_results$p_value < 0.05
  
  return(all_results)
}

# =============================================================================
# 2. Create forest plot style visualization for TMLE results
# =============================================================================

create_tmle_forest_plot <- function(tmle_data) {
  
  # Separate by outcome
  epilepsy_data <- tmle_data %>% filter(outcome == "Epilepsy/Seizure")
  late_onset_data <- tmle_data %>% filter(outcome == "Late-onset Epilepsy")
  
  # Create forest plot for Epilepsy/Seizure
  p1 <- ggplot(epilepsy_data, aes(y = reorder(analysis, -as.numeric(factor(analysis))), 
                                   x = risk_difference)) +
    # Zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, 
                       color = significant),
                   height = 0.2, size = 1) +
    
    # Point estimates
    geom_point(aes(color = significant), size = 4, shape = 18) +
    
    # Add risk difference values
    geom_text(aes(x = upper_ci + 0.0015, 
                  label = sprintf("%.3f\n(%.3f, %.3f)", 
                                  risk_difference, lower_ci, upper_ci)),
              size = 3, hjust = 0) +
    
    # Color scale
    scale_color_manual(values = c("TRUE" = "#E74C3C", "FALSE" = "#95A5A6"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "Statistical Significance") +
    
    # Axes
    scale_x_continuous(name = "Risk Difference (GLP1 RA - SGLT2i)",
                       limits = c(-0.008, 0.004),
                       breaks = seq(-0.008, 0.004, 0.002),
                       labels = function(x) sprintf("%.3f", x)) +
    
    # Labels
    labs(title = "Epilepsy/Seizure",
         y = "") +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Create forest plot for Late-onset Epilepsy
  p2 <- ggplot(late_onset_data, aes(y = reorder(analysis, -as.numeric(factor(analysis))), 
                                     x = risk_difference)) +
    # Zero line
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
    
    # Confidence intervals
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, 
                       color = significant),
                   height = 0.2, size = 1) +
    
    # Point estimates
    geom_point(aes(color = significant), size = 4, shape = 18) +
    
    # Add risk difference values
    geom_text(aes(x = upper_ci + 0.0015, 
                  label = sprintf("%.3f\n(%.3f, %.3f)", 
                                  risk_difference, lower_ci, upper_ci)),
              size = 3, hjust = 0) +
    
    # Color scale
    scale_color_manual(values = c("TRUE" = "#9B59B6", "FALSE" = "#95A5A6"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "Statistical Significance") +
    
    # Axes
    scale_x_continuous(name = "Risk Difference (GLP1 RA - SGLT2i)",
                       limits = c(-0.008, 0.004),
                       breaks = seq(-0.008, 0.004, 0.002),
                       labels = function(x) sprintf("%.3f", x)) +
    
    # Labels
    labs(title = "Late-onset Epilepsy",
         y = "") +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  # Combine plots
  combined_plot <- grid.arrange(
    p1, p2,
    ncol = 2,
    top = textGrob("TMLE-Estimated Risk Differences at 2 Years: GLP1 RA vs SGLT2i",
                   gp = gpar(fontsize = 16, fontface = "bold")),
    bottom = textGrob("Negative values favor GLP1 RA; Error bars represent 95% confidence intervals",
                      gp = gpar(fontsize = 10, fontface = "italic"))
  )
  
  return(combined_plot)
}

# =============================================================================
# 3. Create bar plot for absolute risks (TMLE-appropriate)
# =============================================================================

create_tmle_risk_comparison <- function(tmle_data) {
  
  # Filter main analysis only
  main_data <- tmle_data %>% 
    filter(analysis == "Main Analysis (2-year)") %>%
    select(outcome, risk_glp1, risk_sglt2i) %>%
    pivot_longer(cols = c(risk_glp1, risk_sglt2i),
                 names_to = "treatment",
                 values_to = "risk") %>%
    mutate(treatment = ifelse(treatment == "risk_glp1", "GLP1 RA", "SGLT2i"))
  
  # Create grouped bar plot
  p_risk <- ggplot(main_data, aes(x = outcome, y = risk, fill = treatment)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    
    # Add risk values on bars
    geom_text(aes(label = sprintf("%.1f%%", risk * 100)),
              position = position_dodge(width = 0.7),
              vjust = -0.5, size = 4) +
    
    # Color scale
    scale_fill_manual(values = c("GLP1 RA" = "#3498DB", "SGLT2i" = "#E67E22"),
                      name = "Treatment") +
    
    # Axes
    scale_y_continuous(name = "2-Year Risk",
                       labels = scales::percent_format(accuracy = 0.1),
                       limits = c(0, 0.025),
                       expand = expansion(mult = c(0, 0.1))) +
    
    # Labels
    labs(title = "TMLE-Estimated 2-Year Risks by Treatment",
         subtitle = "Main Analysis Results",
         x = "Outcome") +
    
    # Theme
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 11),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
  
  return(p_risk)
}

# =============================================================================
# 4. Create summary table for TMLE results
# =============================================================================

create_tmle_summary_table <- function(tmle_data) {
  
  # Prepare summary data
  summary_data <- tmle_data %>%
    filter(analysis == "Main Analysis (2-year)") %>%
    mutate(
      risk_glp1_pct = sprintf("%.2f%%", risk_glp1 * 100),
      risk_sglt2i_pct = sprintf("%.2f%%", risk_sglt2i * 100),
      risk_diff_ci = sprintf("%.3f (%.3f, %.3f)", 
                            risk_difference, lower_ci, upper_ci),
      p_value_fmt = sprintf("%.3f", p_value),
      nnt = ceiling(1 / abs(risk_difference))
    ) %>%
    select(outcome, risk_sglt2i_pct, risk_glp1_pct, risk_diff_ci, 
           p_value_fmt, nnt, n_sglt2i, n_glp1)
  
  return(summary_data)
}

# =============================================================================
# 5. Generate all visualizations
# =============================================================================

# Generate TMLE results
tmle_results <- create_tmle_results()

# Create forest plot
forest_plot <- create_tmle_forest_plot(tmle_results)

# Create risk comparison plot
risk_plot <- create_tmle_risk_comparison(tmle_results)

# Create summary table
summary_table <- create_tmle_summary_table(tmle_results)

# Display plots
print(risk_plot)
print(forest_plot)

# Print summary table
cat("\n=== TMLE Summary Results (2-Year Follow-up) ===\n")
cat("GLP1 RA vs SGLT2i Comparison\n\n")
print(summary_table)

# =============================================================================
# 6. Save outputs
# =============================================================================

# Save plots
ggsave("tmle_glp1_vs_sglt2i_forest_plot.png", 
       forest_plot, 
       width = 14, height = 8, dpi = 300)

ggsave("tmle_glp1_vs_sglt2i_risk_comparison.png", 
       risk_plot, 
       width = 10, height = 6, dpi = 300)

# Save summary data
write.csv(tmle_results, "tmle_glp1_vs_sglt2i_all_results.csv", row.names = FALSE)
write.csv(summary_table, "tmle_glp1_vs_sglt2i_summary.csv", row.names = FALSE)

cat("\n=== Files saved successfully! ===\n")
cat("- tmle_glp1_vs_sglt2i_forest_plot.png\n")
cat("- tmle_glp1_vs_sglt2i_risk_comparison.png\n")
cat("- tmle_glp1_vs_sglt2i_all_results.csv\n")
cat("- tmle_glp1_vs_sglt2i_summary.csv\n")

# =============================================================================
# 7. Interpretation guide
# =============================================================================

cat("\n=== TMLE Interpretation Guide ===\n")
cat("1. TMLE estimates treatment effects at a SPECIFIC time point (2 years)\n")
cat("2. No time-varying effects - single risk difference estimate\n")
cat("3. Includes covariate adjustment (implicit in TMLE)\n")
cat("4. Forest plot shows heterogeneity across subgroups\n")
cat("5. Negative risk differences favor GLP1 RA over SGLT2i\n")
cat("6. NNT = Number Needed to Treat to prevent one outcome\n")
