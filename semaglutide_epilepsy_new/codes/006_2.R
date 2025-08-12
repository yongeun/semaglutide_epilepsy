# 006_2
# ===============================================================================
# STEP-BY-STEP COHORT ANALYSIS FUNCTION DEFINITION
# ===============================================================================

perform_stepwise_cohort_analysis <- function() {
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("FLOWCHART-BASED COHORT SELECTION ANALYSIS\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # STEP 1: All of Us participants with EHR data (matching flowchart)
  cat("\n=== STEP 1: ALL OF US PARTICIPANTS ===\n")
  initial_cohort <- merged_df_2
  total_participants <- nrow(initial_cohort)
  cat("Participants within All of Us with EHR data:", format(total_participants, big.mark=","), "\n")
  
  # Verify we have the right data
  cat("Data verification:\n")
  cat("  - Columns available:", paste(head(names(initial_cohort), 10), collapse=", "), "...\n")
  cat("  - Key column 'dm' available:", "dm" %in% names(initial_cohort), "\n")
  cat("  - Key column 'epilepsy_or_seizure' available:", "epilepsy_or_seizure" %in% names(initial_cohort), "\n")
  
  # STEP 2: Diabetes filter (matching flowchart)
  cat("\n=== STEP 2: DIABETES FILTER ===\n")
  diabetes_patients <- initial_cohort %>% filter(dm == 1)
  non_diabetes_excluded <- total_participants - nrow(diabetes_patients)
  
  cat("Excluded (not diagnosed with diabetes):", format(non_diabetes_excluded, big.mark=","), "\n")
  cat("Patients with diabetes:", format(nrow(diabetes_patients), big.mark=","), "\n")
  
  # STEP 3: New user drug prescriptions (matching flowchart exactly)
  cat("\n=== STEP 3: NEW USER DRUG PRESCRIPTIONS ===\n")
  cat("Period: December 2017 through December 2021\n")
  
  # Get drug concepts
  glp1_concepts <- get_drug_concepts(drug_classes$GLP1)
  sglt2_concepts <- get_drug_concepts(drug_classes$SGLT2) 
  su_dpp4_concepts <- get_drug_concepts(drug_classes$SU_DPP4)
  
  # Build new user indices
  glp1_idx <- build_exposure_idx(glp1_concepts, "glp1_idx")
  sglt2_idx <- build_exposure_idx(sglt2_concepts, "sglt2_idx")
  su_dpp4_idx <- build_exposure_idx(su_dpp4_concepts, "su_dpp4_idx")
  
  # Calculate diabetes patients with ANY new drug exposure
  all_drug_user_ids <- unique(c(glp1_idx$person_id, sglt2_idx$person_id, su_dpp4_idx$person_id))
  dm_with_any_drug <- diabetes_patients %>% filter(person_id %in% all_drug_user_ids)
  dm_excluded_no_drug <- nrow(diabetes_patients) - nrow(dm_with_any_drug)
  
  cat("Excluded (No new user prescription of semaglutide, SGLT2i, or SU/DPP4i):", 
      format(dm_excluded_no_drug, big.mark=","), "\n")
  cat("Patients with new user prescriptions of semaglutide, SGLT2i, or SU/DPP4i:", 
      format(nrow(dm_with_any_drug), big.mark=","), "\n")
  
  # STEP 4: Create comparison groups (exactly matching flowchart)
  cat("\n=== STEP 4: COMPARISON GROUP FORMATION ===\n")
  
  # First, show individual drug user counts
  cat("Individual drug new users:\n")
  cat("  - GLP1 (semaglutide) new users:", format(nrow(glp1_idx), big.mark=","), "\n")
  cat("  - SGLT2i new users:", format(nrow(sglt2_idx), big.mark=","), "\n")
  cat("  - SU/DPP4i new users:", format(nrow(su_dpp4_idx), big.mark=","), "\n")
  
  # Check overlaps
  glp1_sglt2_overlap <- length(intersect(glp1_idx$person_id, sglt2_idx$person_id))
  glp1_su_dpp4_overlap <- length(intersect(glp1_idx$person_id, su_dpp4_idx$person_id))
  sglt2_su_dpp4_overlap <- length(intersect(sglt2_idx$person_id, su_dpp4_idx$person_id))
  
  cat("Drug overlaps (patients with both):\n")
  cat("  - GLP1 & SGLT2i:", format(glp1_sglt2_overlap, big.mark=","), "\n")
  cat("  - GLP1 & SU/DPP4i:", format(glp1_su_dpp4_overlap, big.mark=","), "\n")
  cat("  - SGLT2i & SU/DPP4i:", format(sglt2_su_dpp4_overlap, big.mark=","), "\n")
  
  # Based on flowchart: Create three distinct comparison cohorts from the new user patients
  
  # Group 1: Semaglutide OR SGLT2i users (from new user cohort)
  glp1_or_sglt2_ids <- unique(c(glp1_idx$person_id, sglt2_idx$person_id))
  glp1_or_sglt2_dm <- dm_with_any_drug %>% filter(person_id %in% glp1_or_sglt2_ids)
  
  # Group 2: Semaglutide OR SU/DPP4i users (from new user cohort)
  glp1_or_su_dpp4_ids <- unique(c(glp1_idx$person_id, su_dpp4_idx$person_id))
  glp1_or_su_dpp4_dm <- dm_with_any_drug %>% filter(person_id %in% glp1_or_su_dpp4_ids)
  
  # Group 3: SGLT2i OR SU/DPP4i users (from new user cohort)
  sglt2_or_su_dpp4_ids <- unique(c(sglt2_idx$person_id, su_dpp4_idx$person_id))
  sglt2_or_su_dpp4_dm <- dm_with_any_drug %>% filter(person_id %in% sglt2_or_su_dpp4_ids)
  
  cat("\nComparison group sizes:\n")
  cat("Patient prescribed a semaglutide OR a SGLT2i:", format(nrow(glp1_or_sglt2_dm), big.mark=","), "\n")
  cat("Patient prescribed a semaglutide OR a SU/DPP4i:", format(nrow(glp1_or_su_dpp4_dm), big.mark=","), "\n") 
  cat("Patient prescribed a SGLT2i OR a SU/DPP4i:", format(nrow(sglt2_or_su_dpp4_dm), big.mark=","), "\n")
  
  # STEP 5: Apply exclusions for each comparison group (matching flowchart)
  cat("\n=== STEP 5: EXCLUSIONS FOR EACH COMPARISON ===\n")
  
  # Function to calculate epilepsy/seizure exclusions (diagnosis before index date)
  calculate_epilepsy_exclusions <- function(cohort_df) {
    # Use the binary epilepsy flag first
    epilepsy_status <- table(cohort_df$epilepsy_or_seizure, useNA = "ifany")
    
    # Count those with epilepsy == "1" 
    if("1" %in% names(epilepsy_status)) {
      epilepsy_excluded <- epilepsy_status["1"]
    } else {
      epilepsy_excluded <- 0
    }
    
    cat("  DEBUG: Epilepsy exclusions in cohort (n=", nrow(cohort_df), "):", epilepsy_excluded, "\n")
    
    # TODO: For more accurate exclusions, we should check:
    # - epilepsy_or_seizure_start_date < index_date for each patient
    # - This requires joining with drug indices to get index dates
    
    return(as.numeric(epilepsy_excluded))
  }
  
  # Function to calculate temporal prior use exclusions
  calculate_temporal_prior_use <- function(exposure_idx, comparator_idx, exposure_name, comparator_name) {
    # Find patients who have BOTH drugs
    both_drugs_ids <- intersect(exposure_idx$person_id, comparator_idx$person_id)
    
    if(length(both_drugs_ids) == 0) {
      return(list(
        prior_exposure_use = 0,
        prior_comparator_use = 0,
        both_simultaneous = 0,
        total_excluded = 0,
        exposure_only = nrow(exposure_idx),
        comparator_only = nrow(comparator_idx)
      ))
    }
    
    # For patients with both drugs, determine temporal order
    exposure_both <- exposure_idx[exposure_idx$person_id %in% both_drugs_ids, ]
    comparator_both <- comparator_idx[comparator_idx$person_id %in% both_drugs_ids, ]
    
    # Get the correct column names - build_exposure_idx renames to specific idx names
    exposure_date_col <- names(exposure_idx)[2]  # Second column is the renamed index date
    comparator_date_col <- names(comparator_idx)[2]  # Second column is the renamed index date
    
    cat("  DEBUG: Column names in indices:\n")
    cat("    - Exposure idx columns:", paste(names(exposure_idx), collapse=", "), "\n")
    cat("    - Comparator idx columns:", paste(names(comparator_idx), collapse=", "), "\n")
    cat("    - Using exposure date column:", exposure_date_col, "\n")
    cat("    - Using comparator date column:", comparator_date_col, "\n")
    
    # Merge to get both dates for same patients using correct column names
    temporal_analysis <- merge(
      exposure_both[, c("person_id", exposure_date_col)], 
      comparator_both[, c("person_id", comparator_date_col)],
      by = "person_id"
    )
    
    # Convert dates and compare - use the original column names directly
    temporal_analysis$exposure_date <- as.Date(temporal_analysis[[exposure_date_col]])
    temporal_analysis$comparator_date <- as.Date(temporal_analysis[[comparator_date_col]])
    
    # Categorize based on temporal order
    prior_exposure_use <- sum(temporal_analysis$exposure_date < temporal_analysis$comparator_date, na.rm = TRUE)
    prior_comparator_use <- sum(temporal_analysis$comparator_date < temporal_analysis$exposure_date, na.rm = TRUE)  
    both_simultaneous <- sum(temporal_analysis$exposure_date == temporal_analysis$comparator_date, na.rm = TRUE)
    
    # For new user analysis, exclude patients with prior use of the other drug
    total_excluded <- prior_exposure_use + prior_comparator_use + both_simultaneous
    
    # Patients with only one drug (no prior use issues)
    exposure_only <- nrow(exposure_idx) - length(both_drugs_ids)
    comparator_only <- nrow(comparator_idx) - length(both_drugs_ids)
    
    cat("  DEBUG: Temporal prior use analysis:\n")
    cat("    -", exposure_name, "only (no", comparator_name, "):", format(exposure_only, big.mark=","), "\n")
    cat("    -", comparator_name, "only (no", exposure_name, "):", format(comparator_only, big.mark=","), "\n")
    cat("    - Prior", exposure_name, "use (before", comparator_name, "):", format(prior_exposure_use, big.mark=","), "\n")
    cat("    - Prior", comparator_name, "use (before", exposure_name, "):", format(prior_comparator_use, big.mark=","), "\n")
    cat("    - Both drugs same date:", format(both_simultaneous, big.mark=","), "\n")
    cat("    - Total excluded for prior use:", format(total_excluded, big.mark=","), "\n")
    
    return(list(
      prior_exposure_use = prior_exposure_use,
      prior_comparator_use = prior_comparator_use,
      both_simultaneous = both_simultaneous,
      total_excluded = total_excluded,
      exposure_only = exposure_only,
      comparator_only = comparator_only
    ))
  }
  
  # Analysis for Group 1: Semaglutide vs SGLT2i (temporal prior use analysis)
  cat("\n--- SEMAGLUTIDE vs SGLT2i ---\n")
  group1_total <- nrow(glp1_or_sglt2_dm)
  group1_prior_result <- calculate_temporal_prior_use(glp1_idx, sglt2_idx, "GLP1", "SGLT2i")
  group1_prior_excluded <- group1_prior_result$total_excluded
  group1_epilepsy_excluded <- calculate_epilepsy_exclusions(glp1_or_sglt2_dm)
  group1_total_excluded <- group1_prior_excluded + group1_epilepsy_excluded
  group1_final <- group1_total - group1_total_excluded
  
  cat("Total cohort:", format(group1_total, big.mark=","), "\n")
  cat("Excluded:\n")
  if(group1_prior_result$prior_comparator_use > 0) {
    cat("  - Prior use of SGLT2i:", format(group1_prior_result$prior_comparator_use, big.mark=","), "\n")
  }
  if(group1_prior_result$prior_exposure_use > 0) {
    cat("  - Prior use of GLP1:", format(group1_prior_result$prior_exposure_use, big.mark=","), "\n")
  }
  if(group1_prior_result$both_simultaneous > 0) {
    cat("  - GLP1 and SGLT2i same date:", format(group1_prior_result$both_simultaneous, big.mark=","), "\n")
  }
  cat("  - Total prior use exclusions:", format(group1_prior_excluded, big.mark=","), "\n")
  cat("  - Epilepsy/seizure diagnosis before index date:", format(group1_epilepsy_excluded, big.mark=","), "\n")
  cat("  - Total excluded:", format(group1_total_excluded, big.mark=","), "\n")
  cat("Final cohort:", format(group1_final, big.mark=","), "\n")
  
  # Calculate final treatment group sizes for GLP1 vs SGLT2i CORRECTLY
  # We need to start with the comparison cohort and apply exclusions properly
  
  # Create a cohort that tracks which drug each patient is assigned to
  group1_cohort_with_treatment <- glp1_or_sglt2_dm %>%
    mutate(
      has_glp1 = person_id %in% glp1_idx$person_id,
      has_sglt2 = person_id %in% sglt2_idx$person_id,
      has_both = has_glp1 & has_sglt2,
      treatment_assignment = case_when(
        has_both ~ "excluded_both",
        has_glp1 & !has_sglt2 ~ "glp1",
        !has_glp1 & has_sglt2 ~ "sglt2",
        TRUE ~ "error"
      ),
      has_epilepsy = epilepsy_or_seizure == "1"
    ) %>%
    # Apply exclusions
    filter(
      treatment_assignment != "excluded_both",  # Remove patients with both drugs
      !has_epilepsy  # Remove patients with epilepsy
    )
  
  # Count final treatment groups
  group1_treatment_counts <- table(group1_cohort_with_treatment$treatment_assignment)
  cat("  DEBUG: Group 1 treatment assignment counts:", "\n")
  print(group1_treatment_counts)
  
  group1_glp1_final <- ifelse("glp1" %in% names(group1_treatment_counts), 
                              as.numeric(group1_treatment_counts["glp1"]), 0)
  group1_sglt2_final <- ifelse("sglt2" %in% names(group1_treatment_counts), 
                               as.numeric(group1_treatment_counts["sglt2"]), 0)
  
  cat("Final treatment groups:\n")
  cat("  - Semaglutide:", format(group1_glp1_final, big.mark=","), "\n")
  cat("  - SGLT2i:", format(group1_sglt2_final, big.mark=","), "\n")
  cat("  - Verification: Total should be", format(group1_glp1_final + group1_sglt2_final, big.mark=","), 
      "vs reported", format(group1_final, big.mark=","), "\n")
  
  # Analysis for Group 2: Semaglutide vs SU/DPP4i (temporal prior use analysis)
  cat("\n--- SEMAGLUTIDE vs SU/DPP4i ---\n")
  group2_total <- nrow(glp1_or_su_dpp4_dm)
  group2_prior_result <- calculate_temporal_prior_use(glp1_idx, su_dpp4_idx, "GLP1", "SU/DPP4i")
  group2_prior_excluded <- group2_prior_result$total_excluded
  group2_epilepsy_excluded <- calculate_epilepsy_exclusions(glp1_or_su_dpp4_dm)
  group2_total_excluded <- group2_prior_excluded + group2_epilepsy_excluded
  group2_final <- group2_total - group2_total_excluded
  
  cat("Total cohort:", format(group2_total, big.mark=","), "\n")
  cat("Excluded:\n")
  if(group2_prior_result$prior_comparator_use > 0) {
    cat("  - Prior use of SU/DPP4i:", format(group2_prior_result$prior_comparator_use, big.mark=","), "\n")
  }
  if(group2_prior_result$prior_exposure_use > 0) {
    cat("  - Prior use of GLP1:", format(group2_prior_result$prior_exposure_use, big.mark=","), "\n")
  }
  if(group2_prior_result$both_simultaneous > 0) {
    cat("  - GLP1 and SU/DPP4i same date:", format(group2_prior_result$both_simultaneous, big.mark=","), "\n")
  }
  cat("  - Total prior use exclusions:", format(group2_prior_excluded, big.mark=","), "\n")
  cat("  - Epilepsy/seizure diagnosis before index date:", format(group2_epilepsy_excluded, big.mark=","), "\n")
  cat("  - Total excluded:", format(group2_total_excluded, big.mark=","), "\n")
  cat("Final cohort:", format(group2_final, big.mark=","), "\n")
  
  # Calculate final treatment group sizes for GLP1 vs SU/DPP4i CORRECTLY
  group2_cohort_with_treatment <- glp1_or_su_dpp4_dm %>%
    mutate(
      has_glp1 = person_id %in% glp1_idx$person_id,
      has_su_dpp4 = person_id %in% su_dpp4_idx$person_id,
      has_both = has_glp1 & has_su_dpp4,
      treatment_assignment = case_when(
        has_both ~ "excluded_both",
        has_glp1 & !has_su_dpp4 ~ "glp1",
        !has_glp1 & has_su_dpp4 ~ "su_dpp4",
        TRUE ~ "error"
      ),
      has_epilepsy = epilepsy_or_seizure == "1"
    ) %>%
    # Apply exclusions
    filter(
      treatment_assignment != "excluded_both",  # Remove patients with both drugs
      !has_epilepsy  # Remove patients with epilepsy
    )
  
  # Count final treatment groups
  group2_treatment_counts <- table(group2_cohort_with_treatment$treatment_assignment)
  cat("  DEBUG: Group 2 treatment assignment counts:", "\n")
  print(group2_treatment_counts)
  
  group2_glp1_final <- ifelse("glp1" %in% names(group2_treatment_counts), 
                              as.numeric(group2_treatment_counts["glp1"]), 0)
  group2_su_dpp4_final <- ifelse("su_dpp4" %in% names(group2_treatment_counts), 
                                 as.numeric(group2_treatment_counts["su_dpp4"]), 0)
  
  cat("Final treatment groups:\n")
  cat("  - Semaglutide:", format(group2_glp1_final, big.mark=","), "\n")
  cat("  - SU/DPP4i:", format(group2_su_dpp4_final, big.mark=","), "\n")
  cat("  - Verification: Total should be", format(group2_glp1_final + group2_su_dpp4_final, big.mark=","), 
      "vs reported", format(group2_final, big.mark=","), "\n")
  
  # Analysis for Group 3: SGLT2i vs SU/DPP4i (temporal prior use analysis)
  cat("\n--- SGLT2i vs SU/DPP4i ---\n")
  group3_total <- nrow(sglt2_or_su_dpp4_dm)
  group3_prior_result <- calculate_temporal_prior_use(sglt2_idx, su_dpp4_idx, "SGLT2i", "SU/DPP4i")
  group3_prior_excluded <- group3_prior_result$total_excluded
  group3_epilepsy_excluded <- calculate_epilepsy_exclusions(sglt2_or_su_dpp4_dm)
  group3_total_excluded <- group3_prior_excluded + group3_epilepsy_excluded
  group3_final <- group3_total - group3_total_excluded
  
  cat("Total cohort (SGLT2i OR SU/DPP4i users):", format(group3_total, big.mark=","), "\n")
  cat("Drug composition:\n")
  cat("  - SGLT2i only users:", format(group3_prior_result$exposure_only, big.mark=","), "\n")
  cat("  - SU/DPP4i only users:", format(group3_prior_result$comparator_only, big.mark=","), "\n")
  cat("  - Users with both drugs:", format(group3_prior_result$prior_exposure_use + group3_prior_result$prior_comparator_use + group3_prior_result$both_simultaneous, big.mark=","), "\n")
  
  cat("Excluded:\n")
  if(group3_prior_result$prior_exposure_use > 0) {
    cat("  - Prior use of SGLT2i:", format(group3_prior_result$prior_exposure_use, big.mark=","), "\n")
  }
  if(group3_prior_result$prior_comparator_use > 0) {
    cat("  - Prior use of SU/DPP4i:", format(group3_prior_result$prior_comparator_use, big.mark=","), "\n")
  }
  if(group3_prior_result$both_simultaneous > 0) {
    cat("  - SGLT2i and SU/DPP4i same date:", format(group3_prior_result$both_simultaneous, big.mark=","), "\n")
  }
  cat("  - Total prior use exclusions:", format(group3_prior_excluded, big.mark=","), "\n")
  cat("  - Epilepsy/seizure diagnosis before index date:", format(group3_epilepsy_excluded, big.mark=","), "\n")
  cat("  - Total excluded:", format(group3_total_excluded, big.mark=","), "\n")
  cat("Final cohort:", format(group3_final, big.mark=","), "\n")
  
  # Calculate final treatment group sizes for SGLT2i vs SU/DPP4i CORRECTLY
  group3_cohort_with_treatment <- sglt2_or_su_dpp4_dm %>%
    mutate(
      has_sglt2 = person_id %in% sglt2_idx$person_id,
      has_su_dpp4 = person_id %in% su_dpp4_idx$person_id,
      has_both = has_sglt2 & has_su_dpp4,
      treatment_assignment = case_when(
        has_both ~ "excluded_both",
        has_sglt2 & !has_su_dpp4 ~ "sglt2",
        !has_sglt2 & has_su_dpp4 ~ "su_dpp4",
        TRUE ~ "error"
      ),
      has_epilepsy = epilepsy_or_seizure == "1"
    ) %>%
    # Apply exclusions
    filter(
      treatment_assignment != "excluded_both",  # Remove patients with both drugs
      !has_epilepsy  # Remove patients with epilepsy
    )
  
  # Count final treatment groups
  group3_treatment_counts <- table(group3_cohort_with_treatment$treatment_assignment)
  cat("  DEBUG: Group 3 treatment assignment counts:", "\n")
  print(group3_treatment_counts)
  
  group3_sglt2_final <- ifelse("sglt2" %in% names(group3_treatment_counts), 
                               as.numeric(group3_treatment_counts["sglt2"]), 0)
  group3_su_dpp4_final <- ifelse("su_dpp4" %in% names(group3_treatment_counts), 
                                 as.numeric(group3_treatment_counts["su_dpp4"]), 0)
  
  cat("Final treatment groups:\n")
  cat("  - SGLT2i:", format(group3_sglt2_final, big.mark=","), "\n")
  cat("  - SU/DPP4i:", format(group3_su_dpp4_final, big.mark=","), "\n")
  cat("  - Verification: Total should be", format(group3_sglt2_final + group3_su_dpp4_final, big.mark=","), 
      "vs reported", format(group3_final, big.mark=","), "\n")
  
  # STEP 6: ADDITIONAL EXCLUSIONS TRACKING (Final Analysis Stage)
  cat("\n=== STEP 6: ADDITIONAL EXCLUSIONS FOR FINAL ANALYSIS ===\n")
  
  # Function to get final CSV numbers and track additional exclusions
  get_final_csv_numbers <- function() {
    csv_results <- list()
    
    # Check all outcome-specific cohort files
    comparisons <- c("GLP1_vs_SGLT2", "GLP1_vs_SU_DPP4")
  outcomes <- c("Epilepsy_Seizure", "Early-onset_Epilepsy_Seizure", "Late-onset_Epilepsy_Seizure")
    
    for(comp in comparisons) {
      csv_results[[comp]] <- list()
      
      for(outcome in outcomes) {
        # Check IPTW cohort file (primary analysis)
        iptw_file <- paste0("ipwt_cohort_", comp, "_", outcome, ".csv")
        if(file.exists(iptw_file)) {
          iptw_cohort <- read_csv(iptw_file, show_col_types = FALSE)
          iptw_counts <- table(iptw_cohort$treatment)
          
          csv_results[[comp]][[outcome]] <- list(
            semaglutide = as.numeric(iptw_counts["1"]),
            comparator = as.numeric(iptw_counts["0"]),
            total = sum(iptw_counts),
            file = iptw_file
          )
        }
      }
    }
    
    return(csv_results)
  }
  
  # Get final numbers from CSV files
  final_csv_numbers <- get_final_csv_numbers()
  
  # Calculate additional exclusions between flowchart and final analysis
  cat("ADDITIONAL EXCLUSIONS APPLIED BETWEEN FLOWCHART AND FINAL ANALYSIS:\n")
  cat("(These explain the difference between flowchart numbers and publication table numbers)\n\n")
  
  if(length(final_csv_numbers) > 0) {
    # GLP1 vs SGLT2 additional exclusions
    if("GLP1_vs_SGLT2" %in% names(final_csv_numbers) && 
       "Epilepsy_Seizure" %in% names(final_csv_numbers[["GLP1_vs_SGLT2"]])) {
      
      csv_data <- final_csv_numbers[["GLP1_vs_SGLT2"]][["Epilepsy_Seizure"]]
      flowchart_sema <- group1_glp1_final
      flowchart_sglt2 <- group1_sglt2_final
      flowchart_total <- group1_final
      
      csv_sema <- csv_data$semaglutide
      csv_sglt2 <- csv_data$comparator
      csv_total <- csv_data$total
      
      additional_sema_excluded <- flowchart_sema - csv_sema
      additional_sglt2_excluded <- flowchart_sglt2 - csv_sglt2
      additional_total_excluded <- flowchart_total - csv_total
      
      cat("GLP1 vs SGLT2 (Epilepsy/Seizure analysis):\n")
      cat("  Flowchart numbers:    Sema=", format(flowchart_sema, big.mark=","), 
          ", SGLT2=", format(flowchart_sglt2, big.mark=","), 
          ", Total=", format(flowchart_total, big.mark=","), "\n")
      cat("  Final CSV numbers:    Sema=", format(csv_sema, big.mark=","), 
          ", SGLT2=", format(csv_sglt2, big.mark=","), 
          ", Total=", format(csv_total, big.mark=","), "\n")
      cat("  Additional excluded:  Sema=", format(additional_sema_excluded, big.mark=","), 
          ", SGLT2=", format(additional_sglt2_excluded, big.mark=","), 
          ", Total=", format(additional_total_excluded, big.mark=","), "\n")
      
      # Try to identify specific exclusions
      cat("  Likely additional exclusions:\n")
      cat("    - Patients filtered out due to event_time < 0\n")
      cat("    - Patients with missing values in covariates\n") 
      cat("    - Patients with insufficient follow-up time\n")
      cat("    - Patients excluded during propensity score analysis\n")
      cat("    - Patients with missing outcome data\n\n")
    }
    
    # GLP1 vs SU/DPP4 additional exclusions
    if("GLP1_vs_SU_DPP4" %in% names(final_csv_numbers) && 
       "Epilepsy_Seizure" %in% names(final_csv_numbers[["GLP1_vs_SU_DPP4"]])) {
      
      csv_data <- final_csv_numbers[["GLP1_vs_SU_DPP4"]][["Epilepsy_Seizure"]]
      flowchart_sema <- group2_glp1_final
      flowchart_su_dpp4 <- group2_su_dpp4_final
      flowchart_total <- group2_final
      
      csv_sema <- csv_data$semaglutide
      csv_su_dpp4 <- csv_data$comparator
      csv_total <- csv_data$total
      
      additional_sema_excluded <- flowchart_sema - csv_sema
      additional_su_dpp4_excluded <- flowchart_su_dpp4 - csv_su_dpp4
      additional_total_excluded <- flowchart_total - csv_total
      
      cat("GLP1 vs SU/DPP4 (Epilepsy/Seizure analysis):\n")
      cat("  Flowchart numbers:    Sema=", format(flowchart_sema, big.mark=","), 
          ", SU/DPP4=", format(flowchart_su_dpp4, big.mark=","), 
          ", Total=", format(flowchart_total, big.mark=","), "\n")
      cat("  Final CSV numbers:    Sema=", format(csv_sema, big.mark=","), 
          ", SU/DPP4=", format(csv_su_dpp4, big.mark=","), 
          ", Total=", format(csv_total, big.mark=","), "\n")
      cat("  Additional excluded:  Sema=", format(additional_sema_excluded, big.mark=","), 
          ", SU/DPP4=", format(additional_su_dpp4_excluded, big.mark=","), 
          ", Total=", format(additional_total_excluded, big.mark=","), "\n")
      
      cat("  Likely additional exclusions:\n")
      cat("    - Patients filtered out due to event_time < 0\n")
      cat("    - Patients with missing values in covariates\n") 
      cat("    - Patients with insufficient follow-up time\n")
      cat("    - Patients excluded during propensity score analysis\n")
      cat("    - Patients with missing outcome data\n\n")
    }
    
    # SGLT2 vs SU/DPP4 additional exclusions (Group 3)
    cat("SGLT2 vs SU/DPP4 (Epilepsy/Seizure analysis):\n")
    cat("  Flowchart numbers:    SGLT2=", format(group3_sglt2_final, big.mark=","), 
        ", SU/DPP4=", format(group3_su_dpp4_final, big.mark=","), 
        ", Total=", format(group3_final, big.mark=","), "\n")
    
    # Check if CSV files exist for SGLT2 vs SU_DPP4
    sglt2_su_dpp4_csv_file <- "ipwt_cohort_SGLT2_vs_SU_DPP4_Epilepsy_Seizure.csv"
    if(file.exists(sglt2_su_dpp4_csv_file)) {
      sglt2_su_dpp4_csv_data <- read_csv(sglt2_su_dpp4_csv_file, show_col_types = FALSE)
      sglt2_su_dpp4_counts <- table(sglt2_su_dpp4_csv_data$treatment)
      csv_sglt2 <- as.numeric(sglt2_su_dpp4_counts["1"])  # SGLT2 is exposure (1)
      csv_su_dpp4 <- as.numeric(sglt2_su_dpp4_counts["0"])  # SU_DPP4 is comparator (0)
      csv_total <- sum(sglt2_su_dpp4_counts)
     
      additional_sglt2_excluded <- group3_sglt2_final - csv_sglt2
      additional_su_dpp4_excluded <- group3_su_dpp4_final - csv_su_dpp4
      additional_total_excluded <- group3_final - csv_total
      
      cat("  Final CSV numbers:    SGLT2=", format(csv_sglt2, big.mark=","), 
          ", SU/DPP4=", format(csv_su_dpp4, big.mark=","), 
          ", Total=", format(csv_total, big.mark=","), "\n")
      cat("  Additional excluded:  SGLT2=", format(additional_sglt2_excluded, big.mark=","), 
          ", SU/DPP4=", format(additional_su_dpp4_excluded, big.mark=","), 
          ", Total=", format(additional_total_excluded, big.mark=","), "\n")
      
      cat("  Likely additional exclusions:\n")
      cat("    - Patients filtered out due to event_time < 0\n")
      cat("    - Patients with missing values in covariates\n") 
      cat("    - Patients with insufficient follow-up time\n")
      cat("    - Patients excluded during propensity score analysis\n")
      cat("    - Patients with missing outcome data\n\n")
    } else {
      cat("  Final CSV file not found:", sglt2_su_dpp4_csv_file, "\n")
    }
  }
  
  # Return comprehensive results list for use by other functions
  return(list(
    # Step 1: Initial cohort
    total_participants = total_participants,
    
    # Step 2: Diabetes filter
    diabetes_patients = nrow(diabetes_patients),
    non_diabetes_excluded = non_diabetes_excluded,
    
    # Step 3: Drug prescriptions
    glp1_new_users = nrow(glp1_idx),
    sglt2_new_users = nrow(sglt2_idx),
    su_dpp4_new_users = nrow(su_dpp4_idx),
    dm_with_any_drug = nrow(dm_with_any_drug),
    dm_excluded_no_drug = dm_excluded_no_drug,
    
    # Step 4: Comparison groups
    glp1_or_sglt2_total = nrow(glp1_or_sglt2_dm),
    glp1_or_su_dpp4_total = nrow(glp1_or_su_dpp4_dm),
    sglt2_or_su_dpp4_total = nrow(sglt2_or_su_dpp4_dm),
    
    # Step 5: Group 1 (GLP1 vs SGLT2)
    glp1_sglt2_total = group1_total,
    glp1_sglt2_prior_excluded = group1_prior_excluded,
    glp1_sglt2_epilepsy_excluded = group1_epilepsy_excluded,
    glp1_sglt2_total_excluded = group1_total_excluded,
    glp1_sglt2_final = group1_final,
    glp1_sglt2_glp1_final = group1_glp1_final,
    glp1_sglt2_sglt2_final = group1_sglt2_final,
    
    # Step 5: Group 2 (GLP1 vs SU/DPP4)
    glp1_su_dpp4_total = group2_total,
    glp1_su_dpp4_prior_excluded = group2_prior_excluded,
    glp1_su_dpp4_epilepsy_excluded = group2_epilepsy_excluded,
    glp1_su_dpp4_total_excluded = group2_total_excluded,
    glp1_su_dpp4_final = group2_final,
    glp1_su_dpp4_glp1_final = group2_glp1_final,
    glp1_su_dpp4_su_dpp4_final = group2_su_dpp4_final,
    
    # Step 5: Group 3 (SGLT2 vs SU/DPP4)
    sglt2_su_dpp4_total = group3_total,
    sglt2_su_dpp4_prior_excluded = group3_prior_excluded,
    sglt2_su_dpp4_sglt2_only = group3_prior_result$exposure_only,
    sglt2_su_dpp4_su_dpp4_only = group3_prior_result$comparator_only,
    sglt2_su_dpp4_epilepsy_excluded = group3_epilepsy_excluded,
    sglt2_su_dpp4_sglt2_final = group3_sglt2_final,
    sglt2_su_dpp4_su_dpp4_final = group3_su_dpp4_final,
    
    # Final CSV numbers data
    final_csv_numbers = final_csv_numbers
  ))
}

# ===============================================================================
# RUN STEP-BY-STEP COHORT ANALYSIS BEFORE MAIN ANALYSIS  
# ===============================================================================

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("PERFORMING STEP-BY-STEP COHORT ANALYSIS\n") 
cat(paste(rep("=", 80), collapse=""), "\n")

# Run the step-by-step analysis
stepwise_results <- perform_stepwise_cohort_analysis()

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("STEP-BY-STEP ANALYSIS COMPLETE - STARTING MAIN IPWT ANALYSIS\n")
cat(paste(rep("=", 60), collapse=""), "\n")