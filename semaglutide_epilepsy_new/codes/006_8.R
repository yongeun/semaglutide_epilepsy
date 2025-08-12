# =============================================================================
# INDIVIDUAL MEDICATION BREAKDOWN ANALYSIS FOR IPTW COHORTS
# =============================================================================

cat("\n================================================================================ \n")
cat("INDIVIDUAL MEDICATION BREAKDOWN ANALYSIS FOR EPILEPSY/SEIZURE OUTCOME\n")
cat("================================================================================ \n")

# Function to analyze individual medications within drug classes
analyze_individual_medications <- function(cohort_data, comparison_name, drug_class_name) {
  
  cat(paste("\n=== Individual", drug_class_name, "Medications in", comparison_name, "===\n"))
  
  # Determine treatment assignment based on comparison and drug class
  if (comparison_name == "GLP1 vs SGLT2" && drug_class_name == "SGLT2") {
    # In GLP1 vs SGLT2, SGLT2 patients have treatment = 0 (comparator)
    drug_class_patients <- cohort_data %>% filter(treatment == 0)
  } else if (comparison_name == "SGLT2 vs SU_DPP4" && drug_class_name == "SGLT2") {
    # In SGLT2 vs SU_DPP4, SGLT2 patients have treatment = 1 (exposure)
    drug_class_patients <- cohort_data %>% filter(treatment == 1)
  } else if (comparison_name == "SGLT2 vs SU_DPP4" && drug_class_name == "SU_DPP4") {
    # In SGLT2 vs SU_DPP4, SU_DPP4 patients have treatment = 0 (comparator)
    drug_class_patients <- cohort_data %>% filter(treatment == 0)
  } else if (comparison_name == "GLP1 vs SU_DPP4" && drug_class_name == "SU_DPP4") {
    # In GLP1 vs SU_DPP4, SU_DPP4 patients have treatment = 0 (comparator)
    drug_class_patients <- cohort_data %>% filter(treatment == 0)
  } else {
    cat("Invalid combination:", comparison_name, "-", drug_class_name, "\n")
    return(NULL)
  }
  
  if (nrow(drug_class_patients) == 0) {
    cat("No patients found in", drug_class_name, "group for", comparison_name, "\n")
    return(NULL)
  }
  
  cat(paste("Found", nrow(drug_class_patients), drug_class_name, "patients in", comparison_name, "comparison\n"))
  
  # Define individual medications for each class
  if (drug_class_name == "SGLT2") {
    # SGLT2 inhibitors
    individual_meds <- c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin")
  } else if (drug_class_name == "SU_DPP4") {
    # SU and DPP4 inhibitors
    individual_meds <- c("glimepiride", "glipizide", "glyburide", "alogliptin", "linagliptin", "sitagliptin", "saxagliptin")
  } else {
    cat("Unknown drug class:", drug_class_name, "\n")
    return(NULL)
  }
  
  # Query BigQuery to get individual medication usage
  medication_counts <- list()
  
  for (med in individual_meds) {
    cat(paste("Analyzing", med, "usage...\n"))
    
    # Get concept IDs for this medication
    med_concepts_sql <- glue("
      SELECT DISTINCT c2.concept_id
      FROM `{DATASET}.concept` c
      JOIN `{DATASET}.concept_ancestor` ca ON c.concept_id = ca.ancestor_concept_id
      JOIN `{DATASET}.concept` c2 ON c2.concept_id = ca.descendant_concept_id
      WHERE c.concept_class_id = 'Ingredient' 
        AND LOWER(c.concept_name) LIKE '%{tolower(med)}%'
    ")
    
    med_concepts <- tryCatch({
      download_data(med_concepts_sql)
    }, error = function(e) {
      cat(paste("Error getting concepts for", med, ":", e$message, "\n"))
      return(data.frame(concept_id = integer(0)))
    })
    
    if (nrow(med_concepts) > 0) {
      concept_ids <- paste(med_concepts$concept_id, collapse = ", ")
      
      # Get patients using this specific medication
      med_usage_sql <- glue("
        SELECT DISTINCT de.person_id
        FROM `{DATASET}.drug_exposure` de
        WHERE de.person_id IN ({paste(drug_class_patients$person_id, collapse = ', ')})
          AND de.drug_concept_id IN ({concept_ids})
          AND de.drug_exposure_start_date BETWEEN '{exposure_window$start}' AND '{exposure_window$end}'
      ")
      
      med_users <- tryCatch({
        download_data(med_usage_sql)
      }, error = function(e) {
        cat(paste("Error getting users for", med, ":", e$message, "\n"))
        return(data.frame(person_id = integer(0)))
      })
      
      medication_counts[[med]] <- nrow(med_users)
      cat(paste("-", med, ":", nrow(med_users), "patients\n"))
    } else {
      medication_counts[[med]] <- 0
      cat(paste("-", med, ": 0 patients (no concepts found)\n"))
    }
  }
  
  # Create summary table
  med_summary <- data.frame(
    Medication = names(medication_counts),
    N_Patients = unlist(medication_counts),
    Percentage = round(unlist(medication_counts) / nrow(drug_class_patients) * 100, 1)
  )
  
  cat(paste("\nTotal", drug_class_name, "patients:", nrow(drug_class_patients), "\n"))
  cat("\nIndividual medication breakdown:\n")
  print(med_summary)
  
  return(med_summary)
}

# Load IPTW cohort data for epilepsy/seizure outcome
cat("\nLoading IPTW cohort data for medication breakdown analysis...\n")

# GLP1 vs SGLT2 comparison
cat("\n--- Loading GLP1 vs SGLT2 cohort ---\n")
glp1_sglt2_file <- "ipwt_cohort_GLP1_vs_SGLT2_Epilepsy_Seizure.csv"
if (!file.exists(glp1_sglt2_file)) {
  system(paste0("gsutil cp ", my_bucket, "/data/", glp1_sglt2_file, " ."), intern = TRUE)
}

glp1_sglt2_cohort <- tryCatch({
  read_csv(glp1_sglt2_file, show_col_types = FALSE)
}, error = function(e) {
  cat("Could not load GLP1 vs SGLT2 cohort:", e$message, "\n")
  NULL
})

if (!is.null(glp1_sglt2_cohort)) {
  # Analyze SGLT2 medications in GLP1 vs SGLT2 comparison
  sglt2_in_glp1_comparison <- analyze_individual_medications(
    glp1_sglt2_cohort, 
    "GLP1 vs SGLT2", 
    "SGLT2"
  )
}

# SGLT2 vs SU_DPP4 comparison
cat("\n--- Loading SGLT2 vs SU_DPP4 cohort ---\n")
sglt2_su_file <- "ipwt_cohort_SGLT2_vs_SU_DPP4_Epilepsy_Seizure.csv"
if (!file.exists(sglt2_su_file)) {
  system(paste0("gsutil cp ", my_bucket, "/data/", sglt2_su_file, " ."), intern = TRUE)
}

sglt2_su_cohort <- tryCatch({
  read_csv(sglt2_su_file, show_col_types = FALSE)
}, error = function(e) {
  cat("Could not load SGLT2 vs SU_DPP4 cohort:", e$message, "\n")
  NULL
})

if (!is.null(sglt2_su_cohort)) {
  # Analyze SGLT2 medications in SGLT2 vs SU_DPP4 comparison
  sglt2_in_sglt2_comparison <- analyze_individual_medications(
    sglt2_su_cohort, 
    "SGLT2 vs SU_DPP4", 
    "SGLT2"
  )
  
  # Analyze SU_DPP4 medications in SGLT2 vs SU_DPP4 comparison
  su_dpp4_in_sglt2_comparison <- analyze_individual_medications(
    sglt2_su_cohort, 
    "SGLT2 vs SU_DPP4", 
    "SU_DPP4"
  )
}

# GLP1 vs SU_DPP4 comparison
cat("\n--- Loading GLP1 vs SU_DPP4 cohort ---\n")
glp1_su_file <- "ipwt_cohort_GLP1_vs_SU_DPP4_Epilepsy_Seizure.csv"
if (!file.exists(glp1_su_file)) {
  system(paste0("gsutil cp ", my_bucket, "/data/", glp1_su_file, " ."), intern = TRUE)
}

glp1_su_cohort <- tryCatch({
  read_csv(glp1_su_file, show_col_types = FALSE)
}, error = function(e) {
  cat("Could not load GLP1 vs SU_DPP4 cohort:", e$message, "\n")
  NULL
})

if (!is.null(glp1_su_cohort)) {
  # Analyze SU_DPP4 medications in GLP1 vs SU_DPP4 comparison
  su_dpp4_in_glp1_comparison <- analyze_individual_medications(
    glp1_su_cohort, 
    "GLP1 vs SU_DPP4", 
    "SU_DPP4"
  )
}

# Create comprehensive summary table
cat("\n================================================================================ \n")
cat("COMPREHENSIVE MEDICATION BREAKDOWN SUMMARY\n")
cat("================================================================================ \n")

# Combine all results into a comprehensive table
all_med_results <- list()

if (exists("sglt2_in_glp1_comparison") && !is.null(sglt2_in_glp1_comparison)) {
  sglt2_in_glp1_comparison$Comparison <- "GLP1 vs SGLT2"
  sglt2_in_glp1_comparison$Drug_Class <- "SGLT2"
  all_med_results[["sglt2_glp1"]] <- sglt2_in_glp1_comparison
}

if (exists("sglt2_in_sglt2_comparison") && !is.null(sglt2_in_sglt2_comparison)) {
  sglt2_in_sglt2_comparison$Comparison <- "SGLT2 vs SU_DPP4"
  sglt2_in_sglt2_comparison$Drug_Class <- "SGLT2"
  all_med_results[["sglt2_sglt2"]] <- sglt2_in_sglt2_comparison
}

if (exists("su_dpp4_in_sglt2_comparison") && !is.null(su_dpp4_in_sglt2_comparison)) {
  su_dpp4_in_sglt2_comparison$Comparison <- "SGLT2 vs SU_DPP4"
  su_dpp4_in_sglt2_comparison$Drug_Class <- "SU_DPP4"
  all_med_results[["su_dpp4_sglt2"]] <- su_dpp4_in_sglt2_comparison
}

if (exists("su_dpp4_in_glp1_comparison") && !is.null(su_dpp4_in_glp1_comparison)) {
  su_dpp4_in_glp1_comparison$Comparison <- "GLP1 vs SU_DPP4"
  su_dpp4_in_glp1_comparison$Drug_Class <- "SU_DPP4"
  all_med_results[["su_dpp4_glp1"]] <- su_dpp4_in_glp1_comparison
}

if (length(all_med_results) > 0) {
  # Combine all results
  comprehensive_med_breakdown <- do.call(rbind, all_med_results)
  rownames(comprehensive_med_breakdown) <- NULL
  
  # Reorder columns
  comprehensive_med_breakdown <- comprehensive_med_breakdown[, c("Comparison", "Drug_Class", "Medication", "N_Patients", "Percentage")]
  
  cat("\nComprehensive Individual Medication Breakdown:\n")
  print(comprehensive_med_breakdown)
  
  # Save to CSV
  tryCatch({
    write_csv(comprehensive_med_breakdown, "individual_medication_breakdown_epilepsy_seizure.csv")
    
    # Copy to bucket
    system(paste0("gsutil cp individual_medication_breakdown_epilepsy_seizure.csv ", 
                  Sys.getenv('WORKSPACE_BUCKET'), "/data/"), intern = TRUE)
    
    cat("\nSaved comprehensive medication breakdown to: individual_medication_breakdown_epilepsy_seizure.csv\n")
  }, error = function(e) {
    cat("Error saving medication breakdown file:", e$message, "\n")
  })
  
  # Create summary by drug class
  cat("\n=== SUMMARY BY DRUG CLASS ===\n")
  
  drug_class_summary <- comprehensive_med_breakdown %>%
    group_by(Comparison, Drug_Class) %>%
    summarise(
      Total_Patients = sum(N_Patients),
      N_Medications = n(),
      Top_Medication = Medication[which.max(N_Patients)],
      Top_Med_N = max(N_Patients),
      Top_Med_Pct = round(max(N_Patients) / sum(N_Patients) * 100, 1),
      .groups = "drop"
    )
  
  print(drug_class_summary)
  
} else {
  cat("No medication breakdown data available.\n")
}

cat("\n================================================================================ \n")
cat("INDIVIDUAL MEDICATION BREAKDOWN ANALYSIS COMPLETE\n")
cat("================================================================================ \n")

# End of script
cat("\n\nAll analyses completed successfully.\n")