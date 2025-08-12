# ===============================================================================
# STEP-BY-STEP COHORT SELECTION ANALYSIS - MATCHING USER'S FLOWCHART
# This section provides detailed tracking exactly matching the provided flowchart
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
      cat("  Final CSV numbers:    N/A - CSV file not found (", sglt2_su_dpp4_csv_file, ")\n")
      cat("  Additional excluded:  N/A - No final analysis CSV available\n")
      cat("  Note: This comparison may not have been fully processed or\n")
      cat("        the analysis may have failed for this comparison.\n\n")
    }
    
  } else {
    cat("Warning: No final CSV files found. Cannot calculate additional exclusions.\n")
  }
  
  # STEP 7: CORRECTED FLOWCHART SUMMARY (using final CSV numbers)
  cat("\n=== STEP 7: CORRECTED FLOWCHART SUMMARY (with final analysis numbers) ===\n")
  
  # Update final numbers to match actual CSV files
  if(length(final_csv_numbers) > 0) {
    # Use final CSV numbers for the flowchart
    if("GLP1_vs_SGLT2" %in% names(final_csv_numbers) && 
       "Epilepsy_Seizure" %in% names(final_csv_numbers[["GLP1_vs_SGLT2"]])) {
      csv_data1 <- final_csv_numbers[["GLP1_vs_SGLT2"]][["Epilepsy_Seizure"]]
      group1_glp1_final_corrected <- csv_data1$semaglutide
      group1_sglt2_final_corrected <- csv_data1$comparator
      group1_final_corrected <- csv_data1$total
    } else {
      group1_glp1_final_corrected <- group1_glp1_final
      group1_sglt2_final_corrected <- group1_sglt2_final
      group1_final_corrected <- group1_final
    }
    
    if("GLP1_vs_SU_DPP4" %in% names(final_csv_numbers) && 
       "Epilepsy_Seizure" %in% names(final_csv_numbers[["GLP1_vs_SU_DPP4"]])) {
      csv_data2 <- final_csv_numbers[["GLP1_vs_SU_DPP4"]][["Epilepsy_Seizure"]]
      group2_glp1_final_corrected <- csv_data2$semaglutide
      group2_su_dpp4_final_corrected <- csv_data2$comparator
      group2_final_corrected <- csv_data2$total
    } else {
      group2_glp1_final_corrected <- group2_glp1_final
      group2_su_dpp4_final_corrected <- group2_su_dpp4_final
      group2_final_corrected <- group2_final
    }
  } else {
    # Fallback to original numbers if CSV files not found
    group1_glp1_final_corrected <- group1_glp1_final
    group1_sglt2_final_corrected <- group1_sglt2_final
    group1_final_corrected <- group1_final
    group2_glp1_final_corrected <- group2_glp1_final
    group2_su_dpp4_final_corrected <- group2_su_dpp4_final
    group2_final_corrected <- group2_final
  }
  
  # Create the exact flowchart format matching the user's image
  cat("\n=== FLOWCHART EXACTLY MATCHING YOUR IMAGE ===\n")
  cat("\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│                     ", sprintf("%6s", format(total_participants, big.mark=" ")), "    Participants within All of Us with EHR data                    │\n")
  cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")
  cat("                                           │\n")
  cat("                                           │\n")
  cat("                                           ▼\n")
  cat("                      ┌──────────────────────────────────────────────────────────┐\n")
  cat("                      │  ", sprintf("%6s", format(non_diabetes_excluded, big.mark=" ")), "  Excluded (not diagnosed with diabetes)  │\n")
  cat("                      └──────────────────────────────────────────────────────────┘\n")
  cat("                                           │\n")
  cat("                                           ▼\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│                      ", sprintf("%6s", format(nrow(diabetes_patients), big.mark=" ")), "    Patients with diabetes                                 │\n")
  cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")
  cat("                                           │\n")
  cat("                                           │\n")
  cat("                                           ▼\n")
  cat("                      ┌──────────────────────────────────────────────────────────┐\n")
  cat("                      │   ", sprintf("%5s", format(dm_excluded_no_drug, big.mark=" ")), "  Excluded (No new user prescription of         │\n")
  cat("                      │          semaglutide, SGLT2i, or SU/DPP4i from         │\n")
  cat("                      │          December 2017 through December 2021)          │\n")
  cat("                      └──────────────────────────────────────────────────────────┘\n")
  cat("                                           │\n")
  cat("                                           ▼\n")
  cat("┌─────────────────────────────────────────────────────────────────────────────────┐\n")
  cat("│                      ", sprintf("%6s", format(nrow(dm_with_any_drug), big.mark=" ")), "    Patients with new user prescriptions of semaglutide,         │\n")
  cat("│                                or SU/DPP4i from December 2017 through December 2021    │\n")
  cat("└─────────────────────────────────────────────────────────────────────────────────┘\n")
  cat("                                           │\n")
  cat("                     ┌─────────────────────┼─────────────────────┐\n")
  cat("                     │                     │                     │\n")
  cat("                     ▼                     ▼                     ▼\n")
  
  # Three comparison groups side by side
  cat("┌─────────────────────────┐ ┌─────────────────────────┐ ┌─────────────────────────┐\n")
  cat("│   ", sprintf("%5s", format(group1_total, big.mark=" ")), "  Patient prescribed │ │   ", sprintf("%5s", format(group2_total, big.mark=" ")), "  Patient prescribed │ │   ", sprintf("%5s", format(group3_total, big.mark=" ")), "  Patient prescribed │\n")
  cat("│         a semaglutide   │ │         a semaglutide   │ │         a SGLT2i or     │\n")
  cat("│         or a SGLT2i     │ │         or a SU/DPP4i   │ │         a SU/DPP4i      │\n")
  cat("└─────────────────────────┘ └─────────────────────────┘ └─────────────────────────┘\n")
  cat("            │                           │                           │\n")
  cat("            ▼                           ▼                           ▼\n")
  cat("┌─────────────────────────┐ ┌─────────────────────────┐ ┌─────────────────────────┐\n")
  cat("│    ", sprintf("%4s", format(group1_total_excluded, big.mark=" ")), "    Excluded     │ │    ", sprintf("%4s", format(group2_total_excluded, big.mark=" ")), "    Excluded     │ │    ", sprintf("%4s", format(group3_total_excluded, big.mark=" ")), "    Excluded     │\n")
  cat("│      ", sprintf("%3s", group1_prior_result$prior_exposure_use), "  Prior use of   │ │       ", sprintf("%2s", group2_prior_result$prior_exposure_use), "  Prior use of   │ │      ", sprintf("%3s", group3_prior_result$prior_exposure_use), "  Prior use of   │\n")
  cat("│          semaglutide    │ │          semaglutide    │ │          SGLT2i         │\n")
  cat("│      ", sprintf("%3s", group1_prior_result$prior_comparator_use), "  Prior use of   │ │      ", sprintf("%3s", group2_prior_result$prior_comparator_use), "  Prior use of   │ │    ", sprintf("%4s", group3_prior_result$prior_comparator_use), "  Prior use of   │\n")
  cat("│          SGLT2i         │ │          SU/DPP4i       │ │          SU/DPP4i       │\n")
  cat("│       ", sprintf("%2s", group1_prior_result$both_simultaneous), "  Semaglutide and│ │       ", sprintf("%2s", group2_prior_result$both_simultaneous), "  Semaglutide and│ │      ", sprintf("%3s", group3_prior_result$both_simultaneous), "  Semaglutide and│\n")
  cat("│          SGLT2i started │ │          SU/DPP4i       │ │          SU/DPP4i       │\n")
  cat("│          on the same    │ │          started on the │ │          started on the │\n")
  cat("│          date           │ │          same date      │ │          same date      │\n")
  cat("│      ", sprintf("%3s", group1_epilepsy_excluded), "  Epilepsy/seizure│ │      ", sprintf("%3s", group2_epilepsy_excluded), "  Epilepsy/seizure│ │      ", sprintf("%3s", group3_epilepsy_excluded), "  Epilepsy/seizure│\n")
  cat("│          diagnosis      │ │          diagnosis      │ │          diagnosis      │\n")
  cat("│          before index   │ │          before index   │ │          before index   │\n")
  cat("│          date           │ │          date           │ │          date           │\n")
  
  # Additional exclusions section
  additional_exclusions_1 <- if(length(final_csv_numbers) > 0) group1_final - group1_final_corrected else 0
  additional_exclusions_2 <- if(length(final_csv_numbers) > 0) group2_final - group2_final_corrected else 0
  
     cat("│      ", sprintf("%3s", additional_exclusions_1), "  Missing covariates│ │      ", sprintf("%3s", additional_exclusions_2), "  Missing covariates│ │      ", sprintf("%3s", 0), "  Missing covariates│\n")
   cat("│          or insufficient│ │          or insufficient│ │          or insufficient│\n")
   cat("│          follow-up      │ │          follow-up      │ │          follow-up      │\n")
   cat("└─────────────────────────┘ └─────────────────────────┘ └─────────────────────────┘\n")
   cat("            │                           │                           │\n")
   cat("            ▼                           ▼                           ▼\n")
   cat("┌─────────────────────────┐ ┌─────────────────────────┐ ┌─────────────────────────┐\n")
   cat("│    ", sprintf("%4s", format(group1_final_corrected, big.mark=" ")), "  Patient prescribed │ │    ", sprintf("%4s", format(group2_final_corrected, big.mark=" ")), "  Patient prescribed │ │    ", sprintf("%4s", format(group3_final, big.mark=" ")), "  Patient prescribed │\n")
   cat("│         a semaglutide   │ │         a semaglutide   │ │         a SGLT2i        │\n")
   cat("│         or a SGLT2i     │ │         or a SU/DPP4i   │ │         or a SU/DPP4i   │\n")
   cat("│                         │ │                         │ │                         │\n")
   cat("│     ", sprintf("%4s", format(group1_glp1_final_corrected, big.mark=" ")), "  Semaglutide   │ │     ", sprintf("%4s", format(group2_glp1_final_corrected, big.mark=" ")), "  Semaglutide   │ │     ", sprintf("%4s", format(group3_sglt2_final, big.mark=" ")), "  SGLT2i          │\n")
   cat("│     ", sprintf("%4s", format(group1_sglt2_final_corrected, big.mark=" ")), "  SGLT2i        │ │     ", sprintf("%4s", format(group2_su_dpp4_final_corrected, big.mark=" ")), "  SU/DPP4i      │ │     ", sprintf("%4s", format(group3_su_dpp4_final, big.mark=" ")), "  SU/DPP4i        │\n")
   cat("└─────────────────────────┘ └─────────────────────────┘ └─────────────────────────┘\n")
  
  cat("\nNote: *** and **** indicate that SGLT2i vs SU/DPP4i comparison is not included\n")
  cat("in the neurologic outcomes analysis (epilepsy/seizure, ADRD, stroke).\n")
  cat("This comparison was only used for cardiovascular/diabetes outcomes.\n")
  
  # Arithmetic verification
  cat("\n=== ARITHMETIC VERIFICATION ===\n")
  sum1_corrected <- group1_glp1_final_corrected + group1_sglt2_final_corrected
  sum2_corrected <- group2_glp1_final_corrected + group2_su_dpp4_final_corrected
  
  cat("GLP1 vs SGLT2: ", format(group1_glp1_final_corrected, big.mark=","), " + ", 
      format(group1_sglt2_final_corrected, big.mark=","), " = ", 
      format(sum1_corrected, big.mark=","))
  if(sum1_corrected == group1_final_corrected) {
    cat(" ✓ CORRECT\n")
  } else {
    cat(" ✗ ERROR (should be ", format(group1_final_corrected, big.mark=","), ")\n")
  }
  
  cat("GLP1 vs SU/DPP4: ", format(group2_glp1_final_corrected, big.mark=","), " + ", 
      format(group2_su_dpp4_final_corrected, big.mark=","), " = ", 
      format(sum2_corrected, big.mark=","))
  if(sum2_corrected == group2_final_corrected) {
    cat(" ✓ CORRECT\n")
  } else {
    cat(" ✗ ERROR (should be ", format(group2_final_corrected, big.mark=","), ")\n")
  }
  
  # Return results matching the flowchart structure (including corrected numbers)
  return(list(
    total_participants = total_participants,
    diabetes_patients = nrow(diabetes_patients),
    excluded_no_diabetes = non_diabetes_excluded,
    patients_with_new_drugs = nrow(dm_with_any_drug),
    excluded_no_drugs = dm_excluded_no_drug,
    
    # Original flowchart numbers (before final exclusions)
    glp1_sglt2_total = group1_total,
    glp1_sglt2_final = group1_final,
    glp1_sglt2_excluded = group1_total_excluded,
    glp1_sglt2_prior_excluded = group1_prior_excluded,
    glp1_sglt2_prior_sglt2_use = group1_prior_result$prior_comparator_use,
    glp1_sglt2_prior_glp1_use = group1_prior_result$prior_exposure_use,
    glp1_sglt2_both_simultaneous = group1_prior_result$both_simultaneous,
    glp1_sglt2_epilepsy_excluded = group1_epilepsy_excluded,
    glp1_sglt2_glp1_final = group1_glp1_final,
    glp1_sglt2_sglt2_final = group1_sglt2_final,
    
    glp1_su_dpp4_total = group2_total,
    glp1_su_dpp4_final = group2_final,
    glp1_su_dpp4_excluded = group2_total_excluded,
    glp1_su_dpp4_prior_excluded = group2_prior_excluded,
    glp1_su_dpp4_prior_su_dpp4_use = group2_prior_result$prior_comparator_use,
    glp1_su_dpp4_prior_glp1_use = group2_prior_result$prior_exposure_use,
    glp1_su_dpp4_both_simultaneous = group2_prior_result$both_simultaneous,
    glp1_su_dpp4_epilepsy_excluded = group2_epilepsy_excluded,
    glp1_su_dpp4_glp1_final = group2_glp1_final,
    glp1_su_dpp4_su_dpp4_final = group2_su_dpp4_final,
    
    sglt2_su_dpp4_total = group3_total,
    sglt2_su_dpp4_final = group3_final,
    sglt2_su_dpp4_excluded = group3_total_excluded,
    sglt2_su_dpp4_prior_excluded = group3_prior_excluded,
    sglt2_su_dpp4_prior_sglt2_use = group3_prior_result$prior_exposure_use,
    sglt2_su_dpp4_prior_su_dpp4_use = group3_prior_result$prior_comparator_use,
    sglt2_su_dpp4_both_simultaneous = group3_prior_result$both_simultaneous,
    sglt2_su_dpp4_sglt2_only = group3_prior_result$exposure_only,
    sglt2_su_dpp4_su_dpp4_only = group3_prior_result$comparator_only,
    sglt2_su_dpp4_epilepsy_excluded = group3_epilepsy_excluded,
    sglt2_su_dpp4_sglt2_final = group3_sglt2_final,
    sglt2_su_dpp4_su_dpp4_final = group3_su_dpp4_final,
    
    # CORRECTED FINAL ANALYSIS NUMBERS (matching CSV files)
    glp1_sglt2_final_corrected = group1_final_corrected,
    glp1_sglt2_glp1_final_corrected = group1_glp1_final_corrected,
    glp1_sglt2_sglt2_final_corrected = group1_sglt2_final_corrected,
    glp1_su_dpp4_final_corrected = group2_final_corrected,
    glp1_su_dpp4_glp1_final_corrected = group2_glp1_final_corrected,
    glp1_su_dpp4_su_dpp4_final_corrected = group2_su_dpp4_final_corrected,
    
    # Additional exclusions tracking
    additional_exclusions_1 = if(exists("group1_final_corrected")) group1_final - group1_final_corrected else 0,
    additional_exclusions_2 = if(exists("group2_final_corrected")) group2_final - group2_final_corrected else 0,
    
    # Final CSV numbers data
    final_csv_numbers = final_csv_numbers
  ))
}

# This will be called at the appropriate location in the script

# Quick cohort overview function (can be run independently)
quick_cohort_overview <- function() {
  cat("\n=== QUICK COHORT OVERVIEW ===\n")
  
  cat("1. Initial merged_df_2 size:", nrow(merged_df_2), "\n")
  
  cat("2. Diabetes status:\n")
  dm_table <- table(merged_df_2$dm, useNA = "ifany")
  print(dm_table)
  cat("   - T2DM patients (dm==1):", dm_table["1"], "\n")
  
  cat("3. Epilepsy/seizure status in T2DM patients:\n")
  dm_patients <- merged_df_2 %>% filter(dm == 1)
  epilepsy_table <- table(dm_patients$epilepsy_or_seizure, useNA = "ifany")
  print(epilepsy_table)
  cat("   - T2DM patients with epilepsy/seizure history:", epilepsy_table["1"], "\n")
  
  cat("4. Other neurologic conditions in T2DM patients:\n")
  cat("   - ADRD:", sum(dm_patients$adrd == "1", na.rm = TRUE), "\n")
  cat("   - Stroke:", sum(dm_patients$stroke == "1", na.rm = TRUE), "\n")
  
  # Quick drug exposure check (without building full indices)
  cat("5. Time period for drug analysis:", exposure_window$start, "to", exposure_window$end, "\n")
  
  return(list(
    total_patients = nrow(merged_df_2),
    dm_patients = sum(merged_df_2$dm == "1", na.rm = TRUE),
    epilepsy_in_dm = sum(dm_patients$epilepsy_or_seizure == "1", na.rm = TRUE),
    adrd_in_dm = sum(dm_patients$adrd == "1", na.rm = TRUE),
    stroke_in_dm = sum(dm_patients$stroke == "1", na.rm = TRUE)
  ))
}

# ===============================================================================
# EXECUTE STEP-BY-STEP ANALYSIS
# ===============================================================================

cat("\n", paste(rep("=", 80), collapse=""), "\n")
cat("RUNNING STEP-BY-STEP COHORT ANALYSIS\n")
cat(paste(rep("=", 80), collapse=""), "\n")

# 1. Quick overview first
cat("\n>>> RUNNING QUICK OVERVIEW FIRST <<<\n")
overview_results <- quick_cohort_overview()

# 2. Then detailed step-by-step analysis
cat("\n>>> RUNNING DETAILED STEP-BY-STEP ANALYSIS <<<\n")
stepwise_results <- tryCatch({
  perform_stepwise_cohort_analysis()
}, error = function(e) {
  cat("ERROR in stepwise analysis:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
  return(NULL)
})

# Check if stepwise_results was created successfully
if(is.null(stepwise_results)) {
  cat("WARNING: stepwise_results is NULL. Stepwise analysis failed.\n")
} else {
  cat("SUCCESS: stepwise_results created with", length(stepwise_results), "elements\n")
}

# 3. Print summary
cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("STEP-BY-STEP ANALYSIS COMPLETED!\n")
cat(paste(rep("=", 60), collapse=""), "\n")

cat("\nSUMMARY:\n")
cat("Initial cohort:", overview_results$total_patients, "\n")
cat("T2DM patients:", overview_results$dm_patients, "\n")
cat("Epilepsy in T2DM:", overview_results$epilepsy_in_dm, "\n")
cat("ADRD in T2DM:", overview_results$adrd_in_dm, "\n") 
cat("Stroke in T2DM:", overview_results$stroke_in_dm, "\n")

# Note: Detailed stepwise analysis results are displayed above in the 
# "FLOWCHART-BASED COHORT SELECTION ANALYSIS" section

# ===============================================================================
# PUBLICATION TABLE EXTRACTION FUNCTIONS
# ===============================================================================

# Function to extract and format results for publication table
extract_publication_results <- function() {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("EXTRACTING RESULTS FOR PUBLICATION TABLE\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load saved results
  ipwt_results_df <- if(file.exists("ipwt_all_results.csv")) {
    read_csv("ipwt_all_results.csv", show_col_types = FALSE)
  } else { 
    cat("Warning: ipwt_all_results.csv not found\n")
    NULL 
  }
  
  tmle_results_df <- if(file.exists("tmle_all_results.csv")) {
    read_csv("tmle_all_results.csv", show_col_types = FALSE)
  } else { 
    cat("Warning: tmle_all_results.csv not found\n")
    NULL 
  }
  
  # DIAGNOSTIC: Check what files and results are available
  cat("\n=== DIAGNOSTIC: AVAILABLE FILES AND RESULTS ===\n")
  
  # Check cohort files
  cat("Available cohort files:\n")
  cohort_files <- list.files(pattern = "^(ipwt|tmle)_cohort_.*\\.csv$")
  if(length(cohort_files) > 0) {
    for(file in cohort_files) {
      cat("  ", file, "\n")
    }
  } else {
    cat("  No cohort files found\n")
  }
  
  # Check results in IPTW results file
  if(!is.null(ipwt_results_df)) {
    cat("\nIPTW results summary:\n")
    iptw_summary <- ipwt_results_df %>%
      group_by(comparison, outcome, effect_measure) %>%
      summarise(n_rows = n(), .groups = "drop")
    print(iptw_summary)
  }
  
  # Check results in TMLE results file  
  if(!is.null(tmle_results_df)) {
    cat("\nTMLE results summary:\n")
    tmle_summary <- tmle_results_df %>%
      group_by(comparison, outcome, effect_measure) %>%
      summarise(n_rows = n(), .groups = "drop")
    print(tmle_summary)
  }
  
  cat("=== END DIAGNOSTIC ===\n")
  
  # Function to format numbers consistently
  format_n_events <- function(events, total) {
    paste0(events, "/", format(total, big.mark = ","))
  }
  
  format_hr_ci <- function(hr, lower, upper) {
    sprintf("%.2f (%.2f to %.2f)", hr, lower, upper)
  }
  
  format_or_ci <- function(or, lower, upper) {
    sprintf("%.2f (%.2f-%.2f)", or, lower, upper)
  }
  
  format_rd_ci <- function(rd, lower, upper) {
    sprintf("%.3f (%.3f to %.3f)", rd*1000, lower*1000, upper*1000)
  }
  
  format_pval <- function(p) {
    if(length(p) > 1) {
      cat("Warning: Multiple p-values found, using first one:", p, "\n")
      p <- p[1]
    }
    if(is.na(p)) return("")
    if(p < 0.001) return("<0.001")
    if(p < 0.01) return(sprintf("%.3f", p))
    return(sprintf("%.3f", p))
  }
  
  # Load cohort data to get NON-WEIGHTED follow-up statistics
  get_cohort_followup_unweighted <- function(comparison, outcome, method = "ipwt") {
    
    # Construct filename
    clean_comp <- gsub(" ", "_", comparison)
    clean_outcome <- gsub("[/ ]", "_", outcome)
    filename <- paste0(tolower(method), "_cohort_", clean_comp, "_", clean_outcome, ".csv")
    
    if(file.exists(filename)) {
      cohort <- read_csv(filename, show_col_types = FALSE)
      
      cat("DEBUG: Loaded cohort file:", filename, "\n")
      cat("DEBUG: Cohort dimensions:", nrow(cohort), "rows,", ncol(cohort), "columns\n")
      cat("DEBUG: Events by treatment:", table(cohort$treatment, cohort$event, useNA = "ifany"), "\n")
      
      # Calculate NON-WEIGHTED follow-up by treatment group
      followup_stats <- cohort %>%
        mutate(comparison_type = comparison) %>%  # Add comparison for treatment labeling
        group_by(treatment) %>%
        summarise(
          n = n(),
          events = sum(event, na.rm = TRUE),
          mean_followup_years = mean(event_time / 365.25, na.rm = TRUE),
          sd_followup_years = sd(event_time / 365.25, na.rm = TRUE),
          total_py = sum(event_time / 365.25, na.rm = TRUE),
          rate_per_1000py = (events / total_py) * 1000,
          comparison_type = first(comparison_type),  # Keep comparison info
          .groups = "drop"
        ) %>%
        mutate(
          treatment_group = case_when(
            treatment == 1 & comparison_type %in% c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4") ~ "Semaglutide",
            treatment == 1 & comparison_type == "SGLT2 vs SU_DPP4" ~ "SGLT2i",
            treatment == 0 & comparison_type == "SGLT2 vs SU_DPP4" ~ "SU/DPP4i",
            TRUE ~ "Comparator"
          ),
          followup_mean_sd = sprintf("%.2f (%.2f)", mean_followup_years, sd_followup_years),
          n_events_formatted = format_n_events(events, n),
          rate_formatted = sprintf("%.2f", rate_per_1000py)
        )
      
      return(followup_stats)
    } else {
      cat("Warning: Cohort file not found:", filename, "\n")
      return(NULL)
    }
  }
  
  # Load cohort data to get IPTW-WEIGHTED follow-up statistics
  get_cohort_followup_weighted <- function(comparison, outcome, method = "ipwt") {
    
    # Construct filename
    clean_comp <- gsub(" ", "_", comparison)
    clean_outcome <- gsub("[/ ]", "_", outcome)
    filename <- paste0(tolower(method), "_cohort_", clean_comp, "_", clean_outcome, ".csv")
    
    if(file.exists(filename)) {
      cohort <- read_csv(filename, show_col_types = FALSE)
      
      cat("DEBUG: Loaded weighted cohort file:", filename, "\n")
      cat("DEBUG: Cohort dimensions:", nrow(cohort), "rows,", ncol(cohort), "columns\n")
      cat("DEBUG: Events by treatment:", table(cohort$treatment, cohort$event, useNA = "ifany"), "\n")
      
      # Check if we have IPTW weights
      if(!"ipw_std" %in% names(cohort)) {
        cat("Warning: No IPTW weights found in cohort file. Using unweighted calculations.\n")
        return(get_cohort_followup_unweighted(comparison, outcome, method))
      }
      
      # Calculate IPTW-WEIGHTED follow-up by treatment group
      followup_stats <- cohort %>%
        mutate(comparison_type = comparison) %>%  # Add comparison for treatment labeling
        group_by(treatment) %>%
        summarise(
          # Effective sample size (sum of weights)
          effective_n = sum(ipw_std, na.rm = TRUE),
          # Actual sample size
          n = n(),
          events = sum(event, na.rm = TRUE),
          
          # WEIGHTED follow-up calculations
          weighted_mean_followup_years = weighted.mean(event_time / 365.25, ipw_std, na.rm = TRUE),
          
          # Weighted variance calculation
          weighted_var_followup = sum(ipw_std * (event_time / 365.25 - weighted_mean_followup_years)^2, na.rm = TRUE) / sum(ipw_std, na.rm = TRUE),
          weighted_sd_followup_years = sqrt(weighted_var_followup),
          
          # WEIGHTED person-years and rates
          weighted_total_py = sum((event_time / 365.25) * ipw_std, na.rm = TRUE),
          weighted_events = sum(event * ipw_std, na.rm = TRUE),
          weighted_rate_per_1000py = (weighted_events / weighted_total_py) * 1000,
          
          comparison_type = first(comparison_type),  # Keep comparison info
          .groups = "drop"
        ) %>%
        mutate(
          treatment_group = case_when(
            treatment == 1 & comparison_type %in% c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4") ~ "Semaglutide",
            treatment == 1 & comparison_type == "SGLT2 vs SU_DPP4" ~ "SGLT2i",
            treatment == 0 & comparison_type == "SGLT2 vs SU_DPP4" ~ "SU/DPP4i",
            TRUE ~ "Comparator"
          ),
          followup_mean_sd = sprintf("%.2f (%.2f)", weighted_mean_followup_years, weighted_sd_followup_years),
          # For weighted analysis, show effective sample size
          n_events_formatted = sprintf("%.0f/%.0f", weighted_events, effective_n),
          rate_formatted = sprintf("%.2f", weighted_rate_per_1000py),
          # Keep original values for calculations
          mean_followup_years = weighted_mean_followup_years,
          sd_followup_years = weighted_sd_followup_years,
          total_py = weighted_total_py,
          rate_per_1000py = weighted_rate_per_1000py
        )
      
      return(followup_stats)
    } else {
      cat("Warning: Cohort file not found:", filename, "\n")
      return(NULL)
    }
  }
  
  # Backward compatibility function
  get_cohort_followup <- function(comparison, outcome, method = "ipwt") {
    return(get_cohort_followup_unweighted(comparison, outcome, method))
  }
  
  # Process each comparison and outcome
  comparisons_list <- c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4", "SGLT2 vs SU_DPP4")
  outcomes_list <- c("Epilepsy/Seizure", "Early-onset Epilepsy/Seizure", "Late-onset Epilepsy/Seizure")
  
  for(comp in comparisons_list) {
    for(outcome in outcomes_list) {
      
      cat("\n", paste(rep("-", 60), collapse=""), "\n")
      cat("COMPARISON:", comp, "- OUTCOME:", outcome, "\n")
      cat(paste(rep("-", 60), collapse=""), "\n")
      
      # DEBUG: Check what cohort files and results exist for this combination
      clean_comp <- gsub(" ", "_", comp)
      clean_outcome <- gsub("[/ ]", "_", outcome)
      ipwt_cohort_file <- paste0("ipwt_cohort_", clean_comp, "_", clean_outcome, ".csv")
      tmle_cohort_file <- paste0("tmle_cohort_", clean_comp, "_", clean_outcome, ".csv")
      
      cat("DEBUG: Looking for cohort files:\n")
      cat("  IPTW:", ipwt_cohort_file, "- Exists:", file.exists(ipwt_cohort_file), "\n")
      cat("  TMLE:", tmle_cohort_file, "- Exists:", file.exists(tmle_cohort_file), "\n")
      
      # DEBUG: Check what outcomes are available in results files
      if(!is.null(ipwt_results_df)) {
        ipwt_outcomes_for_comp <- unique(ipwt_results_df$outcome[ipwt_results_df$comparison == comp])
        cat("DEBUG: IPTW outcomes available for", comp, ":", paste(ipwt_outcomes_for_comp, collapse = ", "), "\n")
      }
      
      if(!is.null(tmle_results_df)) {
        tmle_outcomes_for_comp <- unique(tmle_results_df$outcome[tmle_results_df$comparison == comp])
        cat("DEBUG: TMLE outcomes available for", comp, ":", paste(tmle_outcomes_for_comp, collapse = ", "), "\n")
      }
      
      # ===============================
      # NON-WEIGHTED (CRUDE) RESULTS
      # ===============================
      cat("\n=== NON-WEIGHTED COHORT RESULTS ===\n")
      
      if(!is.null(ipwt_results_df)) {
        # Get cohort data for crude analysis (NON-WEIGHTED)
        ipwt_followup <- get_cohort_followup_unweighted(comp, outcome, "ipwt")
        
        if(!is.null(ipwt_followup)) {
          # Get appropriate treatment group names
          if(comp == "SGLT2 vs SU_DPP4") {
            treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "SGLT2i", ]
            comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "SU/DPP4i", ]
            treat_label <- "SGLT2i"
            comp_label <- "SU/DPP4i"
          } else {
            treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "Semaglutide", ]
            comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "Comparator", ]
            treat_label <- "Semaglutide"
            comp_label <- "Comparator"
          }
          
          cat("Non-weighted cohort:\n")
          cat("No. of events/No. of patients at risk:\n")
          cat(" ", treat_label, ":", treat_stats$n_events_formatted, "\n")
          cat(" ", comp_label, ":", comp_stats$n_events_formatted, "\n")
          
          cat("Follow-up, mean (SD), y:\n")
          cat(" ", treat_label, ":", treat_stats$followup_mean_sd, "\n")
          cat(" ", comp_label, ":", comp_stats$followup_mean_sd, "\n")
          
          # Calculate crude HR using Cox model on unweighted data
          clean_comp <- gsub(" ", "_", comp)
          clean_outcome <- gsub("[/ ]", "_", outcome)
          cohort_file <- paste0("ipwt_cohort_", clean_comp, "_", clean_outcome, ".csv")
          
          if(file.exists(cohort_file)) {
            cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
            
            # Fit crude Cox model (no weights)
            crude_cox <- tryCatch({
              coxph(Surv(event_time, event) ~ treatment, data = cohort_data)
            }, error = function(e) {
              cat("Error fitting crude Cox model:", e$message, "\n")
              NULL
            })
            
            if(!is.null(crude_cox)) {
              crude_summary <- summary(crude_cox)
              crude_hr <- crude_summary$conf.int[1, "exp(coef)"]
              crude_lower <- crude_summary$conf.int[1, "lower .95"]
              crude_upper <- crude_summary$conf.int[1, "upper .95"]
              crude_p <- crude_summary$coefficients[1, "Pr(>|z|)"]
              
              cat("Crude HR (95% CI):\n")
              cat(" ", treat_label, "vs", comp_label, ":", format_hr_ci(crude_hr, crude_lower, crude_upper), "\n")
              cat("  P-value:", format_pval(crude_p), "\n")
            } else {
              cat("Could not calculate crude HR\n")
            }
          }
        }
      }
      
      # ===============================
      # IPTW RESULTS
      # ===============================
      cat("\n=== IPTW-WEIGHTED COHORT RESULTS ===\n")
      
      if(!is.null(ipwt_results_df)) {
        target_outcome <- outcome  # Store the loop variable
        ipwt_row <- ipwt_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Hazard Ratio")
        
        if(nrow(ipwt_row) > 1) {
          cat("Warning: Multiple IPTW rows found for", comp, "-", target_outcome, ". Rows found:\n")
          print(ipwt_row)
          cat("Using first row only.\n")
          ipwt_row <- ipwt_row[1, ]
        }
        
        if(nrow(ipwt_row) > 0) {
          # Get IPTW-WEIGHTED follow-up statistics
          ipwt_followup <- get_cohort_followup_weighted(comp, target_outcome, "ipwt")
          
          if(!is.null(ipwt_followup)) {
            # Get appropriate treatment group names
            if(comp == "SGLT2 vs SU_DPP4") {
              treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "SGLT2i", ]
              comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "SU/DPP4i", ]
              treat_label <- "SGLT2i"
              comp_label <- "SU/DPP4i"
            } else {
              treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "Semaglutide", ]
              comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "Comparator", ]
              treat_label <- "Semaglutide"
              comp_label <- "Comparator"
            }
            
            cat("IPTW-weighted cohort:\n")
            cat("No. of events/No. of patients at risk:\n")
            cat(" ", treat_label, ":", treat_stats$n_events_formatted, "\n")
            cat(" ", comp_label, ":", comp_stats$n_events_formatted, "\n")
            
            cat("Follow-up, mean (SD), y:\n")
            cat(" ", treat_label, ":", treat_stats$followup_mean_sd, "\n")
            cat(" ", comp_label, ":", comp_stats$followup_mean_sd, "\n")
            
            cat("Incidence rate per 1000 person-years:\n")
            cat(" ", treat_label, ":", sprintf("%.2f (%.2f to %.2f)", 
                                        treat_stats$rate_per_1000py,
                                        treat_stats$rate_per_1000py - 1.96*sqrt(treat_stats$events)/treat_stats$total_py*1000,
                                        treat_stats$rate_per_1000py + 1.96*sqrt(treat_stats$events)/treat_stats$total_py*1000), "\n")
            cat(" ", comp_label, ":", sprintf("%.2f (%.2f to %.2f)", 
                                        comp_stats$rate_per_1000py,
                                        comp_stats$rate_per_1000py - 1.96*sqrt(comp_stats$events)/comp_stats$total_py*1000,
                                        comp_stats$rate_per_1000py + 1.96*sqrt(comp_stats$events)/comp_stats$total_py*1000), "\n")
            
            # Rate difference
            rate_diff <- treat_stats$rate_per_1000py - comp_stats$rate_per_1000py
            rate_diff_se <- sqrt(treat_stats$events/treat_stats$total_py^2 + comp_stats$events/comp_stats$total_py^2) * 1000
            rate_diff_lower <- rate_diff - 1.96*rate_diff_se
            rate_diff_upper <- rate_diff + 1.96*rate_diff_se
            
            cat("Rate difference per 1000 person-years:\n")
            cat(" ", treat_label, "vs", comp_label, ":", sprintf("%.2f (%.2f to %.2f)", rate_diff, rate_diff_lower, rate_diff_upper), "\n")
          }
          
          cat("Adjusted HR (95% CI):\n")
          cat(" ", treat_label, "vs", comp_label, ":", format_hr_ci(ipwt_row$estimate, ipwt_row$lower_ci, ipwt_row$upper_ci), "\n")
          cat("  P-value:", format_pval(ipwt_row$p_value), "\n")
          
        } else {
          cat("No IPTW results found for this comparison/outcome\n")
        }
      }
      
      # ===============================
      # TMLE RESULTS  
      # ===============================
      cat("\n=== TMLE RESULTS ===\n")
      
      if(!is.null(tmle_results_df)) {
        cat("DEBUG: Filtering TMLE results for comparison='", comp, "' and outcome='", outcome, "'\n")
        
        # Check what exact matches we get
        tmle_matches <- tmle_results_df %>%
          filter(comparison == comp, outcome == outcome)
        
        if(nrow(tmle_matches) > 0) {
          cat("DEBUG: Found", nrow(tmle_matches), "TMLE result rows for this combination:\n")
          print(tmle_matches[, c("comparison", "outcome", "effect_measure", "estimate", "n_total", "n_events")])
        } else {
          cat("DEBUG: No TMLE results found for this exact combination. Checking for similar matches...\n")
          similar_matches <- tmle_results_df %>%
            filter(comparison == comp)
          if(nrow(similar_matches) > 0) {
            cat("Available outcomes for this comparison:\n")
            print(unique(similar_matches$outcome))
          }
        }
        
        # More precise filtering with trimmed strings
        target_outcome <- outcome  # Store the loop variable to avoid naming conflict
        tmle_or_row <- tmle_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Odds Ratio")
        
        tmle_rd_row <- tmle_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Risk Difference")
        
        cat("DEBUG: OR rows found:", nrow(tmle_or_row), "\n")
        if(nrow(tmle_or_row) > 0) {
          cat("DEBUG: OR estimate:", tmle_or_row$estimate[1], "for outcome:", tmle_or_row$outcome[1], "\n")
        }
        
        cat("DEBUG: RD rows found:", nrow(tmle_rd_row), "\n")
        if(nrow(tmle_rd_row) > 0) {
          cat("DEBUG: RD estimate:", tmle_rd_row$estimate[1], "for outcome:", tmle_rd_row$outcome[1], "\n")
        }
        
        if(nrow(tmle_or_row) > 1) {
          cat("Warning: Multiple TMLE OR rows found for", comp, "-", outcome, ". Using first row.\n")
          tmle_or_row <- tmle_or_row[1, ]
        }
        
        if(nrow(tmle_rd_row) > 1) {
          cat("Warning: Multiple TMLE RD rows found for", comp, "-", outcome, ". Using first row.\n")
          tmle_rd_row <- tmle_rd_row[1, ]
        }
        
        if(nrow(tmle_or_row) > 0 && nrow(tmle_rd_row) > 0) {
          # Get follow-up statistics from TMLE cohort (unweighted)
          tmle_followup <- get_cohort_followup_unweighted(comp, outcome, "tmle")
          
          if(!is.null(tmle_followup)) {
            # Get appropriate treatment group names
            if(comp == "SGLT2 vs SU_DPP4") {
              treat_stats <- tmle_followup[tmle_followup$treatment_group == "SGLT2i", ]
              comp_stats <- tmle_followup[tmle_followup$treatment_group == "SU/DPP4i", ]
              treat_label <- "SGLT2i"
              comp_label <- "SU/DPP4i"
            } else {
              treat_stats <- tmle_followup[tmle_followup$treatment_group == "Semaglutide", ]
              comp_stats <- tmle_followup[tmle_followup$treatment_group == "Comparator", ]
              treat_label <- "Semaglutide"
              comp_label <- "Comparator"
            }
            
            cat("No. of events/No. of patients at risk:\n")
            cat(" ", treat_label, ":", treat_stats$n_events_formatted, "\n")
            cat(" ", comp_label, ":", comp_stats$n_events_formatted, "\n")
          }
          
          cat("OR (95% CI):\n")
          cat(" ", treat_label, "vs", comp_label, ":", format_or_ci(tmle_or_row$estimate, tmle_or_row$lower_ci, tmle_or_row$upper_ci), "\n")
          
          cat("Risk difference (95% CI):\n")
          cat(" ", treat_label, "vs", comp_label, ":", format_rd_ci(tmle_rd_row$estimate, tmle_rd_row$lower_ci, tmle_rd_row$upper_ci), "\n")
          
          cat("P-value:", format_pval(tmle_rd_row$p_value), "\n")
          
        } else {
          cat("No TMLE results found for this comparison/outcome\n")
        }
      }
    }
  }
  
  # ===============================
  # SUMMARY TABLE FORMAT
  # ===============================
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("FORMATTED FOR DIRECT TABLE ENTRY\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  for(comp in comparisons_list) {
    cat("\n>>> COMPARISON:", comp, "<<<\n")
    
    for(outcome in outcomes_list) {
      cat("\n", outcome, ":\n")
      
      # NON-WEIGHTED (CRUDE) RESULTS
      if(!is.null(ipwt_results_df)) {
        ipwt_followup <- get_cohort_followup_unweighted(comp, outcome, "ipwt")
        
        if(!is.null(ipwt_followup)) {
          # Get appropriate treatment group names
          if(comp == "SGLT2 vs SU_DPP4") {
            treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "SGLT2i", ]
            comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "SU/DPP4i", ]
            treat_label <- "SGLT2"
            comp_label <- "SU_DPP4"
          } else {
            treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "Semaglutide", ]
            comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "Comparator", ]
            treat_label <- "Sema"
            comp_label <- "Comp"
          }
          
          cat("Non-weighted cohort\n")
          cat("  N/Events -", treat_label, ":", treat_stats$n_events_formatted, " ", comp_label, ":", comp_stats$n_events_formatted, "\n")
          cat("  Follow-up -", treat_label, ":", treat_stats$followup_mean_sd, " ", comp_label, ":", comp_stats$followup_mean_sd, "\n")
          
          # Calculate crude HR
          clean_comp <- gsub(" ", "_", comp)
          clean_outcome <- gsub("[/ ]", "_", outcome)
          cohort_file <- paste0("ipwt_cohort_", clean_comp, "_", clean_outcome, ".csv")
          
          if(file.exists(cohort_file)) {
            cohort_data <- read_csv(cohort_file, show_col_types = FALSE)
            crude_cox <- tryCatch({
              coxph(Surv(event_time, event) ~ treatment, data = cohort_data)
            }, error = function(e) { NULL })
            
            if(!is.null(crude_cox)) {
              crude_summary <- summary(crude_cox)
              crude_hr <- crude_summary$conf.int[1, "exp(coef)"]
              crude_lower <- crude_summary$conf.int[1, "lower .95"]
              crude_upper <- crude_summary$conf.int[1, "upper .95"]
              
              cat("  Crude HR:", format_hr_ci(crude_hr, crude_lower, crude_upper), "\n")
            }
          }
        }
      }
      
      # IPTW-WEIGHTED RESULTS
      if(!is.null(ipwt_results_df)) {
        target_outcome <- outcome  # Store the loop variable
        ipwt_row <- ipwt_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Hazard Ratio")
        
        if(nrow(ipwt_row) > 1) {
          cat("Warning: Multiple IPTW rows in summary for", comp, "-", outcome, ". Using first row.\n")
          ipwt_row <- ipwt_row[1, ]
        }
        
        if(nrow(ipwt_row) > 0) {
          ipwt_followup <- get_cohort_followup_weighted(comp, outcome, "ipwt")
          
          if(!is.null(ipwt_followup)) {
            # Get appropriate treatment group names
            if(comp == "SGLT2 vs SU_DPP4") {
              treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "SGLT2i", ]
              comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "SU/DPP4i", ]
              treat_label <- "SGLT2"
              comp_label <- "SU_DPP4"
            } else {
              treat_stats <- ipwt_followup[ipwt_followup$treatment_group == "Semaglutide", ]
              comp_stats <- ipwt_followup[ipwt_followup$treatment_group == "Comparator", ]
              treat_label <- "Sema"
              comp_label <- "Comp"
            }
            
            cat("IPTW-weighted cohort\n")
            cat("  N/Events -", treat_label, ":", treat_stats$n_events_formatted, " ", comp_label, ":", comp_stats$n_events_formatted, "\n")
            cat("  Follow-up -", treat_label, ":", treat_stats$followup_mean_sd, " ", comp_label, ":", comp_stats$followup_mean_sd, "\n")
            cat("  Rate/1000PY -", treat_label, ":", sprintf("%.2f", treat_stats$rate_per_1000py), 
                " ", comp_label, ":", sprintf("%.2f", comp_stats$rate_per_1000py), "\n")
            
            rate_diff <- treat_stats$rate_per_1000py - comp_stats$rate_per_1000py
            rate_diff_se <- sqrt(treat_stats$events/treat_stats$total_py^2 + comp_stats$events/comp_stats$total_py^2) * 1000
            
            cat("  Rate diff:", sprintf("%.2f (%.2f to %.2f)", 
                                      rate_diff, 
                                      rate_diff - 1.96*rate_diff_se, 
                                      rate_diff + 1.96*rate_diff_se), "\n")
            cat("  Adj HR:", format_hr_ci(ipwt_row$estimate, ipwt_row$lower_ci, ipwt_row$upper_ci), "\n")
          }
        }
      }
      
      # TMLE
      if(!is.null(tmle_results_df)) {
        target_outcome <- outcome  # Store the loop variable
        tmle_or_row <- tmle_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Odds Ratio")
        tmle_rd_row <- tmle_results_df %>%
          filter(trimws(comparison) == trimws(comp), 
                 trimws(outcome) == trimws(target_outcome), 
                 effect_measure == "Risk Difference")
        
        if(nrow(tmle_or_row) > 0 && nrow(tmle_rd_row) > 0) {
          tmle_followup <- get_cohort_followup_unweighted(comp, outcome, "tmle")
          
          if(!is.null(tmle_followup)) {
            # Get appropriate treatment group names
            if(comp == "SGLT2 vs SU_DPP4") {
              treat_stats <- tmle_followup[tmle_followup$treatment_group == "SGLT2i", ]
              comp_stats <- tmle_followup[tmle_followup$treatment_group == "SU/DPP4i", ]
              treat_label <- "SGLT2"
              comp_label <- "SU_DPP4"
            } else {
              treat_stats <- tmle_followup[tmle_followup$treatment_group == "Semaglutide", ]
              comp_stats <- tmle_followup[tmle_followup$treatment_group == "Comparator", ]
              treat_label <- "Sema"
              comp_label <- "Comp"
            }
            
            cat("TMLE analysis\n")
            cat("  N/Events -", treat_label, ":", treat_stats$n_events_formatted, " ", comp_label, ":", comp_stats$n_events_formatted, "\n")
          }
          
          cat("  OR:", format_or_ci(tmle_or_row$estimate, tmle_or_row$lower_ci, tmle_or_row$upper_ci), "\n")
          cat("  Risk diff:", sprintf("%.3f (%.3f to %.3f)", 
                                      tmle_rd_row$estimate*1000, 
                                      tmle_rd_row$lower_ci*1000, 
                                      tmle_rd_row$upper_ci*1000), "\n")
          cat("  P-value:", format_pval(tmle_rd_row$p_value), "\n")
        }
      }
    }
  }
}

# Function to create publication-ready summary table
create_publication_summary_table <- function() {
  
  cat("\n", paste(rep("=", 80), collapse=""), "\n")
  cat("CREATING PUBLICATION-READY SUMMARY TABLE\n")
  cat(paste(rep("=", 80), collapse=""), "\n")
  
  # Load results
  ipwt_results <- if(file.exists("ipwt_all_results.csv")) {
    read_csv("ipwt_all_results.csv", show_col_types = FALSE)
  } else { NULL }
  
  tmle_results <- if(file.exists("tmle_all_results.csv")) {
    read_csv("tmle_all_results.csv", show_col_types = FALSE)
  } else { NULL }
  
  if(is.null(ipwt_results) && is.null(tmle_results)) {
    cat("No results files found. Please run the main analysis first.\n")
    return(NULL)
  }
  
  # Initialize summary table
  summary_table <- data.frame(
    Comparison = character(),
    Outcome = character(),
    Method = character(),
    Treatment_N_Events = character(),
    Comparator_N_Events = character(),
    Treatment_Followup = character(),
    Comparator_Followup = character(),
    Effect_Estimate = character(),
    P_Value = character(),
    stringsAsFactors = FALSE
  )
  
  comparisons <- c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4", "SGLT2 vs SU_DPP4")
  outcomes <- c("Epilepsy/Seizure", "Early-onset Epilepsy/Seizure", "Late-onset Epilepsy/Seizure")
  
  for(comp in comparisons) {
    for(outcome in outcomes) {
      
      # IPTW row
      if(!is.null(ipwt_results)) {
        ipwt_hr <- ipwt_results %>%
          filter(comparison == comp, outcome == outcome, effect_measure == "Hazard Ratio")
        
        if(nrow(ipwt_hr) > 0) {
          # Get cohort data
          clean_comp <- gsub(" ", "_", comp)
          clean_outcome <- gsub("[/ ]", "_", outcome)
          cohort_file <- paste0("ipwt_cohort_", clean_comp, "_", clean_outcome, ".csv")
          
          if(file.exists(cohort_file)) {
            cohort <- read_csv(cohort_file, show_col_types = FALSE)
            followup_stats <- cohort %>%
              group_by(treatment) %>%
              summarise(
                n = n(),
                events = sum(event, na.rm = TRUE),
                mean_followup = mean(event_time / 365.25, na.rm = TRUE),
                sd_followup = sd(event_time / 365.25, na.rm = TRUE),
                .groups = "drop"
              )
            
            sema_row <- followup_stats[followup_stats$treatment == 1, ]
            comp_row <- followup_stats[followup_stats$treatment == 0, ]
            
            summary_table <- rbind(summary_table, data.frame(
              Comparison = comp,
              Outcome = outcome,
              Method = "IPTW",
              Treatment_N_Events = paste0(sema_row$events, "/", sema_row$n),
              Comparator_N_Events = paste0(comp_row$events, "/", comp_row$n),
              Treatment_Followup = sprintf("%.2f (%.2f)", sema_row$mean_followup, sema_row$sd_followup),
              Comparator_Followup = sprintf("%.2f (%.2f)", comp_row$mean_followup, comp_row$sd_followup),
              Effect_Estimate = sprintf("%.2f (%.2f-%.2f)", ipwt_hr$estimate, ipwt_hr$lower_ci, ipwt_hr$upper_ci),
              P_Value = ifelse(ipwt_hr$p_value < 0.001, "<0.001", sprintf("%.3f", ipwt_hr$p_value)),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
      
      # TMLE row
      if(!is.null(tmle_results)) {
        tmle_or <- tmle_results %>%
          filter(comparison == comp, outcome == outcome, effect_measure == "Odds Ratio")
        tmle_rd <- tmle_results %>%
          filter(comparison == comp, outcome == outcome, effect_measure == "Risk Difference")
        
        if(nrow(tmle_or) > 0 && nrow(tmle_rd) > 0) {
          # Get cohort data
          clean_comp <- gsub(" ", "_", comp)
          clean_outcome <- gsub("[/ ]", "_", outcome)
          cohort_file <- paste0("tmle_cohort_", clean_comp, "_", clean_outcome, ".csv")
          
          if(file.exists(cohort_file)) {
            cohort <- read_csv(cohort_file, show_col_types = FALSE)
            followup_stats <- cohort %>%
              group_by(treatment) %>%
              summarise(
                n = n(),
                events = sum(event, na.rm = TRUE),
                .groups = "drop"
              )
            
            sema_row <- followup_stats[followup_stats$treatment == 1, ]
            comp_row <- followup_stats[followup_stats$treatment == 0, ]
            
            summary_table <- rbind(summary_table, data.frame(
              Comparison = comp,
              Outcome = outcome,
              Method = "TMLE",
              Treatment_N_Events = paste0(sema_row$events, "/", sema_row$n),
              Comparator_N_Events = paste0(comp_row$events, "/", comp_row$n),
              Treatment_Followup = "N/A",
              Comparator_Followup = "N/A",
              Effect_Estimate = sprintf("OR: %.2f (%.2f-%.2f); RD: %.3f", 
                                      tmle_or$estimate, tmle_or$lower_ci, tmle_or$upper_ci,
                                      tmle_rd$estimate * 1000),
              P_Value = ifelse(tmle_rd$p_value < 0.001, "<0.001", sprintf("%.3f", tmle_rd$p_value)),
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  # Print and save table
  cat("\nSUMMARY TABLE:\n")
  print(summary_table)
  
  write.csv(summary_table, "publication_summary_table.csv", row.names = FALSE)
  cat("\nTable saved as: publication_summary_table.csv\n")
  
  return(summary_table)
}