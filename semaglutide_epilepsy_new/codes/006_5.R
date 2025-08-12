# 006_5
# ---- Libraries ----
# Make sure to load all required libraries at the top of your script
library(parallel)
library(doParallel)
library(foreach)
install.packages("doSNOW")
library(doSNOW) # ADDED: For the progress bar to work with foreach

# ---- Bootstrap TMLE Function (PARALLEL VERSION) ----
perform_tmle_bootstrap <- function(df, outcome_var = "event", exclude_vars = NULL,
                                   n_bootstrap = 1000, 
                                   seed = 12345, parallel = TRUE, n_cores = NULL) {

  cat("\nPerforming TMLE bootstrap with", n_bootstrap, "replicates")

  if (parallel) {
    # PARALLEL VERSION
    cat(" using parallel processing...\n")

    # Detect cores if not specified
    if(is.null(n_cores)) {
      n_cores <- min(detectCores() - 1, n_bootstrap)
      cat("Using", n_cores, "cores\n")
    }

    # Set up parallel backend
    cl <- makeCluster(n_cores)
    
    # CHANGED: Register with doSNOW for progress bar compatibility
    registerDoSNOW(cl)

    # Export necessary objects to workers
    clusterExport(cl, c("run_tmle_analysis_fixed", "all_ps_vars"), 
                  envir = .GlobalEnv)

    # Set up progress bar for parallel
    pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)

    # Parallel bootstrap
    bootstrap_results <- foreach(i = 1:n_bootstrap, 
                                 .combine = c,
                                 .packages = c("tmle", "dplyr", "SuperLearner", "tidyr"),
                                 .options.snow = opts, # This now works correctly
                                 .errorhandling = 'pass') %dopar% {

      # Set unique seed for each iteration
      set.seed(seed + i)

      # Create bootstrap sample with replacement
      boot_indices <- sample(1:nrow(df), nrow(df), replace = TRUE)
      boot_df <- df[boot_indices, ]

      # Run TMLE on bootstrap sample
      result <- tryCatch({
        boot_result <- run_tmle_analysis_fixed(boot_df, outcome_var, exclude_vars)
        
        # Extract risk difference (ATE)
        if (!is.null(boot_result$tmle_result$estimates$ATE$psi)) {
          boot_result$tmle_result$estimates$ATE$psi
        } else {
          NA
        }
      }, error = function(e) {
        NA # Return NA if error
      })
      
      result
    }

    close(pb)
    stopCluster(cl)

    # Handle errors in results
    bootstrap_samples <- unlist(bootstrap_results)
    bootstrap_samples <- bootstrap_samples[!is.na(bootstrap_samples)]

  } else {
    # SEQUENTIAL VERSION (Unchanged)
    cat(" using sequential processing...\n")
    set.seed(seed)

    bootstrap_samples <- numeric(n_bootstrap)
    pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)

    for (i in 1:n_bootstrap) {
      boot_indices <- sample(1:nrow(df), nrow(df), replace = TRUE)
      boot_df <- df[boot_indices, ]
      
      tryCatch({
        boot_result <- run_tmle_analysis_fixed(boot_df, outcome_var, exclude_vars)
        
        if (!is.null(boot_result$tmle_result$estimates$ATE$psi)) {
          bootstrap_samples[i] <- boot_result$tmle_result$estimates$ATE$psi
        } else {
          bootstrap_samples[i] <- NA
        }
      }, error = function(e) {
        bootstrap_samples[i] <- NA
        if (i <= 5) cat("\nBootstrap", i, "error:", e$message)
      })
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    bootstrap_samples <- bootstrap_samples[!is.na(bootstrap_samples)]
  }

  cat("\nSuccessful bootstrap replicates:", length(bootstrap_samples), "out of", n_bootstrap, "\n")

  # Calculate bootstrap statistics
  bootstrap_stats <- list(
    samples = bootstrap_samples,
    mean = mean(bootstrap_samples),
    median = median(bootstrap_samples),
    sd = sd(bootstrap_samples),
    ci_lower_2.5 = quantile(bootstrap_samples, 0.025),
    ci_upper_97.5 = quantile(bootstrap_samples, 0.975),
    n_successful = length(bootstrap_samples),
    n_attempted = n_bootstrap
  )

  return(bootstrap_stats)
}


# ---- Function to save bootstrap results ----
save_bootstrap_results <- function(bootstrap_stats, comparison, outcome, bucket = TRUE) {
  
  # Create filename
  filename <- paste0("tmle_bootstrap_",
                     gsub(" ", "_", comparison), "_",
                     gsub("[/ ]", "_", outcome), ".rds")
  
  # Save bootstrap object
  saveRDS(bootstrap_stats, filename)
  cat("\nSaved bootstrap results to:", filename, "\n")
  
  # Also save as CSV for the actual samples
  csv_filename <- paste0("tmle_bootstrap_samples_",
                         gsub(" ", "_", comparison), "_",
                         gsub("[/ ]", "_", outcome), ".csv")
  
  bootstrap_df <- data.frame(
    sample_id = 1:length(bootstrap_stats$samples),
    risk_difference = bootstrap_stats$samples
  )
  
  write_csv(bootstrap_df, csv_filename)
  cat("Saved bootstrap samples to:", csv_filename, "\n")
  
  # Upload to bucket if requested
  if (bucket) {
    my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
    system(paste0("gsutil cp ", filename, " ", my_bucket, "/data/"), intern = TRUE)
    system(paste0("gsutil cp ", csv_filename, " ", my_bucket, "/data/"), intern = TRUE)
    cat("Uploaded files to bucket\n")
  }
  
  return(list(rds_file = filename, csv_file = csv_filename))
}

# ---- TMLE BOOTSTRAP-FOCUSED ANALYSIS ----
late_cutoff <- 50

# --- MODIFICATION (3): Added 'baseline_hba1c' to the list of covariates.
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

analyze_bootstrap_with_tmle <- function() {
    # --- MODIFICATION (0): Filter to only "Late-onset Epilepsy/Seizure" outcome.
    filtered_outcomes <- list()
    for (outc in outcomes) {
      if (outc$label == "Epilepsy/Seizure") {
        filtered_outcomes[[length(filtered_outcomes) + 1]] <- outc
      }
    }

    # --- MODIFICATION (0): Filter comparisons to the two specified.
    filtered_comparisons <- list()
    for (cmp in comparisons) {
      if (cmp$name %in% c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4")) {
        filtered_comparisons[[length(filtered_comparisons) + 1]] <- cmp
      }
    }

    cat("\n=== Running TMLE Bootstrap-Focused Analysis ===\n")
    cat("Filtered to:\n")
    cat("- Outcomes:", paste(sapply(filtered_outcomes, function(x) x$label), collapse = ", "), "\n")
    cat("- Comparisons:", paste(sapply(filtered_comparisons, function(x) x$name), collapse = ", "), "\n\n")
                    
    results_table <- data.frame(
      method         = character(),
      comparison     = character(),
      outcome        = character(),
      effect_measure = character(),
      estimate       = numeric(),
      lower_ci       = numeric(),
      upper_ci       = numeric(),
      p_value        = numeric(),
      n_total        = integer(),
      n_events       = integer(),
      stringsAsFactors = FALSE
    )
  
    for (cmp in filtered_comparisons) {
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
        
        for (outc in filtered_outcomes) {
            cat("\n========= TMLE Bootstrap Analysis:", cmp$name, "-", outc$label, "=========\n")
            
            # --- MODIFICATION (3): Added data prep steps to include baseline HbA1c.
            # Prepare analytic dataframe (pre-HbA1c join)
            analytic_df_pre_hba1c <- cohort_pair %>%
                left_join(merged_df_2, by = "person_id") %>%
                filter(dm == 1) %>%
                mutate(
                    index_year = factor(year(as.Date(index_date))),
                    dob        = as.Date(substr(date_of_birth, 1, 10), format = "%Y/%m/%d"),
                    age_at_end = as.numeric(difftime(data_cut_date, dob, units = "days")) / 365.25
                )

            # Add baseline HbA1c data
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
                # END of baseline HbA1c data preparation
                # EXCLUDE PATIENTS WITH PRE-EXISTING CONDITIONS
                {
                    before_exclusion <- nrow(.)
                    cat("\nBefore excluding pre-existing", outc$label, "patients:", before_exclusion, "\n")
                    
                    result <- filter(.,
                        is.na(epilepsy_or_seizure_start_date) | 
                        as.Date(epilepsy_or_seizure_start_date) >= as.Date(index_date)
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

            # Build exclude_vars
            split_comparator <- unlist(strsplit(cmp$comparator, "_"))
            exclude_vars <- unique(c(cmp$exposure, split_comparator))

            # Run TMLE to get point estimate
            set.seed(12345)
            tmle_result <- tryCatch({
                run_tmle_analysis_fixed(analytic_df, "event", exclude_vars)
            }, error = function(e) {
                cat("Error in TMLE analysis:", e$message, "\n")
                return(NULL)
            })
            
            if (is.null(tmle_result)) {
                cat("Skipping this outcome due to TMLE error.\n")
                next
            }
            
            ate_psi <- round(tmle_result$tmle_result$estimates$ATE$psi, 4)

            # ---- Perform Bootstrap Analysis ----
            if (total_events >= 10) { 
                cat("\n=== Starting Bootstrap Analysis ===\n")

                bootstrap_stats <- perform_tmle_bootstrap(
                    df = analytic_df,
                    outcome_var = "event",
                    exclude_vars = exclude_vars,
                    n_bootstrap = 1000, # You can adjust this
                    seed = 12345
                )

                # Save bootstrap results
                saved_files <- save_bootstrap_results(
                    bootstrap_stats = bootstrap_stats,
                    comparison = cmp$name,
                    outcome = outc$label,
                    bucket = TRUE
                )

                # Print bootstrap summary
                cat("\nBootstrap Summary:\n")
                cat("Mean risk difference:", round(bootstrap_stats$mean, 4), "\n")
                cat("Bootstrap 95% CI: [", round(bootstrap_stats$ci_lower_2.5, 4), 
                    ",", round(bootstrap_stats$ci_upper_97.5, 4), "]\n")
                cat("Bootstrap SD:", round(bootstrap_stats$sd, 4), "\n")

                # Save plot
                png_filename <- paste0("bootstrap_hist_", gsub(" ", "_", cmp$name), "_", 
                                        gsub("[/ ]", "_", outc$label), ".png")
                png(png_filename)
                hist(bootstrap_stats$samples, 
                    main = paste("Bootstrap Distribution:", outc$label),
                    xlab = "Risk Difference",
                    col = "lightblue",
                    breaks = 30)
                abline(v = ate_psi, col = "red", lwd = 2, lty = 2)
                abline(v = 0, col = "black", lwd = 1, lty = 2)
                legend("topright", 
                        legend = c("TMLE Point Estimate", "Null Effect"),
                        col = c("red", "black"),
                        lty = c(2, 2),
                        lwd = c(2, 1))
                dev.off()
                cat("Saved bootstrap histogram to:", png_filename, "\n")

            } else {
                cat("\nSkipping bootstrap - insufficient events (", total_events, ")\n")
            }

            # --- MODIFICATION (1 & 2): ALL MEDIATION AND TRAJECTORY CODE HAS BEEN REMOVED.
        } # End of outcomes loop
    } # End of comparisons loop
  
    cat("\n\nBootstrap-focused analysis complete.\n")
    cat("Check saved files for bootstrap statistics and plots.\n")

    # The main purpose is saving bootstrap files, but we can return the point estimates for inspection.
    return(list(results = results_table))
}
  
# Execute the bootstrap-focused TMLE analysis
bootstrap_output <- analyze_bootstrap_with_tmle()

# ---- TMLE BOOTSTRAP-FOCUSED ANALYSIS ----
late_cutoff <- 50

# --- MODIFICATION (3): Added 'baseline_hba1c' to the list of covariates.
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

analyze_bootstrap_with_tmle <- function() {
    # --- MODIFICATION (0): Filter to only "Late-onset Epilepsy/Seizure" outcome.
    filtered_outcomes <- list()
    for (outc in outcomes) {
      if (outc$label == "Epilepsy/Seizure") {
        filtered_outcomes[[length(filtered_outcomes) + 1]] <- outc
      }
    }

    # --- MODIFICATION (0): Filter comparisons to the two specified.
    filtered_comparisons <- list()
    for (cmp in comparisons) {
      if (cmp$name %in% c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4")) {
        filtered_comparisons[[length(filtered_comparisons) + 1]] <- cmp
      }
    }

    cat("\n=== Running TMLE Bootstrap-Focused Analysis ===\n")
    cat("Filtered to:\n")
    cat("- Outcomes:", paste(sapply(filtered_outcomes, function(x) x$label), collapse = ", "), "\n")
    cat("- Comparisons:", paste(sapply(filtered_comparisons, function(x) x$name), collapse = ", "), "\n\n")
                    
    results_table <- data.frame(
      method         = character(),
      comparison     = character(),
      outcome        = character(),
      effect_measure = character(),
      estimate       = numeric(),
      lower_ci       = numeric(),
      upper_ci       = numeric(),
      p_value        = numeric(),
      n_total        = integer(),
      n_events       = integer(),
      stringsAsFactors = FALSE
    )
  
    for (cmp in filtered_comparisons) {
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
        
        for (outc in filtered_outcomes) {
            cat("\n========= TMLE Bootstrap Analysis:", cmp$name, "-", outc$label, "=========\n")
            
            # --- MODIFICATION (3): Added data prep steps to include baseline HbA1c.
            # Prepare analytic dataframe (pre-HbA1c join)
            analytic_df_pre_hba1c <- cohort_pair %>%
                left_join(merged_df_2, by = "person_id") %>%
                filter(dm == 1) %>%
                mutate(
                    index_year = factor(year(as.Date(index_date))),
                    dob        = as.Date(substr(date_of_birth, 1, 10), format = "%Y/%m/%d"),
                    age_at_end = as.numeric(difftime(data_cut_date, dob, units = "days")) / 365.25
                )

            # Add baseline HbA1c data
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
                # END of baseline HbA1c data preparation
                # EXCLUDE PATIENTS WITH PRE-EXISTING CONDITIONS
                {
                    before_exclusion <- nrow(.)
                    cat("\nBefore excluding pre-existing", outc$label, "patients:", before_exclusion, "\n")
                    
                    result <- filter(.,
                        is.na(epilepsy_or_seizure_start_date) | 
                        as.Date(epilepsy_or_seizure_start_date) >= as.Date(index_date)
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

            # Build exclude_vars
            split_comparator <- unlist(strsplit(cmp$comparator, "_"))
            exclude_vars <- unique(c(cmp$exposure, split_comparator))

            # Run TMLE to get point estimate
            set.seed(12345)
            tmle_result <- tryCatch({
                run_tmle_analysis_fixed(analytic_df, "event", exclude_vars)
            }, error = function(e) {
                cat("Error in TMLE analysis:", e$message, "\n")
                return(NULL)
            })
            
            if (is.null(tmle_result)) {
                cat("Skipping this outcome due to TMLE error.\n")
                next
            }
            
            ate_psi <- round(tmle_result$tmle_result$estimates$ATE$psi, 4)

            # ---- Perform Bootstrap Analysis ----
            if (total_events >= 10) { 
                cat("\n=== Starting Bootstrap Analysis ===\n")

                bootstrap_stats <- perform_tmle_bootstrap(
                    df = analytic_df,
                    outcome_var = "event",
                    exclude_vars = exclude_vars,
                    n_bootstrap = 1000, # You can adjust this
                    seed = 12345
                )

                # Save bootstrap results
                saved_files <- save_bootstrap_results(
                    bootstrap_stats = bootstrap_stats,
                    comparison = cmp$name,
                    outcome = outc$label,
                    bucket = TRUE
                )

                # Print bootstrap summary
                cat("\nBootstrap Summary:\n")
                cat("Mean risk difference:", round(bootstrap_stats$mean, 4), "\n")
                cat("Bootstrap 95% CI: [", round(bootstrap_stats$ci_lower_2.5, 4), 
                    ",", round(bootstrap_stats$ci_upper_97.5, 4), "]\n")
                cat("Bootstrap SD:", round(bootstrap_stats$sd, 4), "\n")

                # Save plot
                png_filename <- paste0("bootstrap_hist_", gsub(" ", "_", cmp$name), "_", 
                                        gsub("[/ ]", "_", outc$label), ".png")
                png(png_filename)
                hist(bootstrap_stats$samples, 
                    main = paste("Bootstrap Distribution:", outc$label),
                    xlab = "Risk Difference",
                    col = "lightblue",
                    breaks = 30)
                abline(v = ate_psi, col = "red", lwd = 2, lty = 2)
                abline(v = 0, col = "black", lwd = 1, lty = 2)
                legend("topright", 
                        legend = c("TMLE Point Estimate", "Null Effect"),
                        col = c("red", "black"),
                        lty = c(2, 2),
                        lwd = c(2, 1))
                dev.off()
                cat("Saved bootstrap histogram to:", png_filename, "\n")

            } else {
                cat("\nSkipping bootstrap - insufficient events (", total_events, ")\n")
            }

            # --- MODIFICATION (1 & 2): ALL MEDIATION AND TRAJECTORY CODE HAS BEEN REMOVED.
        } # End of outcomes loop
    } # End of comparisons loop
  
    cat("\n\nBootstrap-focused analysis complete.\n")
    cat("Check saved files for bootstrap statistics and plots.\n")

    # The main purpose is saving bootstrap files, but we can return the point estimates for inspection.
    return(list(results = results_table))
}
  
# Execute the bootstrap-focused TMLE analysis
bootstrap_output <- analyze_bootstrap_with_tmle()

# ---- Combined Visualization for Multiple Bootstrap Comparisons ----
visualize_combined_bootstrap <- function(comparisons, 
                                        outcome = "Late-onset Epilepsy/Seizure",
                                        colors = NULL,
                                        use_jama_style = TRUE,
                                        plot_type = "multi") {  # "single", "multi", or "both"
  
  # JAMA-style color palette (professional, muted colors with good contrast)
  jama_colors <- c(
    "#00274C",  # Navy blue (primary)
    "#C4622D",  # Burnt orange (secondary)
    "#4B7AA7",  # Steel blue (tertiary)
    "#969696",  # Medium gray
    "#2F7F7E"   # Teal
  )
  
  # Set default colors based on style choice
  if (is.null(colors)) {
    if (use_jama_style) {
      colors <- jama_colors
    } else {
      colors <- c("darkblue", "darkred", "darkgreen", "purple", "orange")
    }
  }
  
  # Set JAMA-style plot parameters
  if (use_jama_style) {
    # Save original parameters
    old_par <- par(no.readonly = TRUE)
    
    # Set JAMA-style parameters
    par(
      family = "sans",         # Clean sans-serif font
      cex.main = 1.1,         # Slightly smaller title
      cex.lab = 1.0,          # Axis labels
      cex.axis = 0.9,         # Axis tick labels
      font.main = 1,          # Plain (not bold) title
      las = 1,                # Horizontal axis labels
      mgp = c(2.5, 0.7, 0),   # Tighter margins
      tcl = -0.3              # Smaller tick marks
    )
  }
  
  # Storage for all bootstrap data
  all_bootstrap_data <- list()
  
  # Load data for each comparison
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    
    # Construct filenames
    rds_filename <- paste0("tmle_bootstrap_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".rds")
    
    csv_filename <- paste0("tmle_bootstrap_samples_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".csv")
    
    # Check if files exist
    if (!file.exists(rds_filename) && !file.exists(csv_filename)) {
      cat("Bootstrap files not found for", comparison, ". Attempting to download from bucket...\n")
      
      # Try to download from bucket
      my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
      system(paste0("gsutil cp ", my_bucket, "/data/", rds_filename, " ."), intern = TRUE)
      system(paste0("gsutil cp ", my_bucket, "/data/", csv_filename, " ."), intern = TRUE)
    }
    
    # Load the bootstrap results
    if (file.exists(rds_filename)) {
      bootstrap_stats <- readRDS(rds_filename)
      bootstrap_samples <- bootstrap_stats$samples
      cat("Loaded bootstrap results from RDS file for", comparison, "\n")
    } else if (file.exists(csv_filename)) {
      bootstrap_df <- read.csv(csv_filename)
      bootstrap_samples <- bootstrap_df$risk_difference
      cat("Loaded bootstrap results from CSV file for", comparison, "\n")
      
      # Recreate bootstrap_stats object
      bootstrap_stats <- list(
        samples = bootstrap_samples,
        mean = mean(bootstrap_samples),
        median = median(bootstrap_samples),
        sd = sd(bootstrap_samples),
        ci_lower_2.5 = quantile(bootstrap_samples, 0.025),
        ci_upper_97.5 = quantile(bootstrap_samples, 0.975),
        n_successful = length(bootstrap_samples),
        n_attempted = 20
      )
    } else {
      warning(paste("Could not find bootstrap results files for", comparison))
      next
    }
    
    # Store the data
    all_bootstrap_data[[comparison]] <- bootstrap_stats
  }
  
  # Check if we have any data
  if (length(all_bootstrap_data) == 0) {
    stop("No bootstrap data could be loaded")
  }
  
  # ---- Prepare plot data ----
  
  # Calculate plot limits
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(all_samples) * 1.1
  
  # Calculate densities for all comparisons
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1
  
  # ---- Create visualizations based on plot_type ----
  
  if (plot_type == "single" || plot_type == "both") {
    # Create single density plot
    plot(NULL, 
         xlim = x_range,
         ylim = c(0, y_max),
         main = if(use_jama_style) paste("Bootstrap Distributions:", outcome) else paste("Combined Bootstrap Densities:", outcome),
         xlab = "Risk Difference",
         ylab = "Density",
         type = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    # Add custom axes with JAMA style
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    
    # Add subtle grid lines (JAMA style)
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    # Add box around plot
    box(lwd = 0.5)
    
    # Add vertical line at zero (null effect)
    abline(v = 0, col = if(use_jama_style) "black" else "black", lwd = 1.5, lty = 1)
    
    # Plot density curves for each comparison
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      dens <- densities[[i]]
      color <- colors[i]
      
      # Plot density line with JAMA-appropriate line weight
      lines(dens, col = color, lwd = if(use_jama_style) 2 else 3)
      
      # Add shaded area under curve (with reduced transparency for JAMA)
      polygon(c(dens$x, rev(dens$x)), 
              c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(color, alpha.f = if(use_jama_style) 0.15 else 0.2),
              border = NA)
      
      # Add vertical line for mean (thinner for JAMA style)
      abline(v = all_bootstrap_data[[comparison_name]]$mean, 
             col = color, lwd = if(use_jama_style) 1 else 2, lty = 2)
    }
    
        # Add legend with JAMA styling - positioned outside plot area
        legend("topright",
               inset = c(-0.35, 0),  # Move outside plot to the right
               xpd = TRUE,           # Allow drawing outside plot region
               legend = c(paste(names(all_bootstrap_data), 
                               "\n  Mean:", 
                               sprintf("%.4f", sapply(all_bootstrap_data, function(x) x$mean))),
                         "Null effect"),
               col = c(colors[1:length(all_bootstrap_data)], "black"),
               lwd = c(rep(if(use_jama_style) 2 else 3, length(all_bootstrap_data)), 1.5),
               lty = c(rep(1, length(all_bootstrap_data)), 1),
               cex = if(use_jama_style) 0.75 else 0.8,
               bg = "white",
               box.lwd = if(use_jama_style) 0.5 else 1,
               box.col = if(use_jama_style) "gray50" else "black")
  }
  
  if (plot_type == "multi" || plot_type == "both") {
    # Create multi-panel visualization
    par(mfrow = c(2, 2))
  
    # Panel 1: Combined density plot
    plot(NULL, 
         xlim = x_range,
         ylim = c(0, y_max),
         main = if(use_jama_style) "A. Density Distributions" else "Combined Bootstrap Densities",
         xlab = "Risk Difference",
         ylab = "Density",
         type = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(v = 0, col = "black", lwd = 1.5, lty = 1)
    
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      dens <- densities[[i]]
      color <- colors[i]
      lines(dens, col = color, lwd = if(use_jama_style) 2 else 3)
      polygon(c(dens$x, rev(dens$x)), 
              c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(color, alpha.f = if(use_jama_style) 0.15 else 0.2),
              border = NA)
    }
    
    # Panel 2: Boxplots comparison
    boxplot_data <- lapply(all_bootstrap_data, function(x) x$samples)
    boxplot(boxplot_data,
            names = names(all_bootstrap_data),
            col = adjustcolor(colors[1:length(all_bootstrap_data)], alpha.f = if(use_jama_style) 0.3 else 0.5),
            border = colors[1:length(all_bootstrap_data)],
            main = if(use_jama_style) "B. Distribution Comparison" else "Bootstrap Distributions Comparison",
            ylab = "Risk Difference",
            las = 2,
            axes = FALSE,
            frame.plot = FALSE,
            lwd = if(use_jama_style) 1 else 1.5,
            boxwex = if(use_jama_style) 0.6 else 0.8,
            staplewex = 0.5,
            outwex = 0.5)
    
    axis(1, at = 1:length(all_bootstrap_data), 
         labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    
    # Panel 3: Confidence intervals comparison
    plot(NULL,
         xlim = c(0.5, length(all_bootstrap_data) + 0.5),
         ylim = x_range,
         main = if(use_jama_style) "C. 95% Confidence Intervals" else "95% Confidence Intervals",
         xlab = "",
         ylab = "Risk Difference",
         xaxt = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    axis(1, at = 1:length(all_bootstrap_data), 
         labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    
    for (i in seq_along(all_bootstrap_data)) {
      data <- all_bootstrap_data[[i]]
      color <- colors[i]
      
      # Draw CI with JAMA-appropriate line weight
      segments(i, data$ci_lower_2.5, i, data$ci_upper_97.5, 
               col = color, lwd = if(use_jama_style) 3 else 4)
      
      # Draw mean point
      points(i, data$mean, pch = 19, col = color, cex = if(use_jama_style) 1.5 else 2)
      
      # Draw median point (optional for JAMA style)
      if (!use_jama_style) {
        points(i, data$median, pch = 4, col = color, cex = 1.5, lwd = 2)
      }
    }
    
    # Panel 4: Statistical summary table as plot
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    title(if(use_jama_style) "D. Statistical Summary" else "Statistical Summary")
    
    # Create summary text with JAMA formatting
    y_pos <- 0.9
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      data <- all_bootstrap_data[[i]]
      color <- colors[i]
      
      text(0.05, y_pos, comparison_name, 
           col = color, font = if(use_jama_style) 1 else 2, adj = 0, cex = if(use_jama_style) 0.85 else 0.9)
      
      y_pos <- y_pos - 0.08
      text(0.1, y_pos, 
           paste("Mean:", sprintf("%.5f", data$mean)), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, 
           paste("95% CI: (", sprintf("%.5f", data$ci_lower_2.5), 
                 ", ", sprintf("%.5f", data$ci_upper_97.5), ")"), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, 
           paste("P(RD<0):", sprintf("%.3f", mean(data$samples < 0))), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.1
    }
    
    # Reset plot parameters
    par(mfrow = c(1, 1))
  }
  
  # Restore original parameters if JAMA style was used
  if (use_jama_style && exists("old_par")) {
    par(old_par)
  }
  
  # Print detailed comparison statistics
  cat("\n=== Bootstrap Distribution Comparison ===\n")
  cat("Outcome:", outcome, "\n\n")
  
  for (comparison_name in names(all_bootstrap_data)) {
    data <- all_bootstrap_data[[comparison_name]]
    cat("--- ", comparison_name, " ---\n")
    cat("  Sample size:", length(data$samples), "\n")
    cat("  Mean (SD):", round(data$mean, 5), "(", round(data$sd, 5), ")\n")
    cat("  Median:", round(data$median, 5), "\n")
    cat("  95% CI: [", round(data$ci_lower_2.5, 5), ",", 
        round(data$ci_upper_97.5, 5), "]\n")
    cat("  P(RD < 0):", round(mean(data$samples < 0), 4), "\n")
    cat("  Min/Max:", round(min(data$samples), 5), "/", 
        round(max(data$samples), 5), "\n\n")
  }
  
  # Perform pairwise comparisons if more than one group
  if (length(all_bootstrap_data) > 1) {
    cat("=== Pairwise Comparisons ===\n")
    comparison_names <- names(all_bootstrap_data)
    
    for (i in 1:(length(comparison_names) - 1)) {
      for (j in (i + 1):length(comparison_names)) {
        name1 <- comparison_names[i]
        name2 <- comparison_names[j]
        samples1 <- all_bootstrap_data[[name1]]$samples
        samples2 <- all_bootstrap_data[[name2]]$samples
        
        # Calculate difference
        mean_diff <- mean(samples1) - mean(samples2)
        
        # Approximate test (treating as independent - note this is an approximation)
        cat(name1, "vs", name2, ":\n")
        cat("  Mean difference:", round(mean_diff, 5), "\n")
        
        # Bootstrap comparison
        if (length(samples1) == length(samples2)) {
          paired_diff <- samples1 - samples2
          cat("  P(diff < 0):", round(mean(paired_diff < 0), 4), "\n")
        }
        cat("\n")
      }
    }
  }
  
  # Return all bootstrap data for further use
  invisible(all_bootstrap_data)
}

# ---- Example Usage ----

# Create only the single combined density plot (no redundancy)
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Late-onset Epilepsy/Seizure",
  use_jama_style = TRUE,
  plot_type = "single"  # Just the density plot
)

# Create only the 4-panel visualization
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Late-onset Epilepsy/Seizure",
  use_jama_style = TRUE,
  plot_type = "multi"  # Default - 4-panel figure
)

# You can also use custom JAMA-appropriate colors if needed
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Late-onset Epilepsy/Seizure",
  colors = c("#00274C", "#C4622D"),  # Navy blue and burnt orange
  use_jama_style = TRUE,
  plot_type = "single"
)

# Compare all three groups with just the density plot
# comparison_results <- visualize_combined_bootstrap(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4", "SGLT2 vs SU_DPP4"),
#   outcome = "Late-onset Epilepsy/Seizure",
#   use_jama_style = TRUE,
#   plot_type = "single"
# )

# For non-JAMA style (original colorful version)
# comparison_results <- visualize_combined_bootstrap(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
#   outcome = "Late-onset Epilepsy/Seizure",
#   use_jama_style = FALSE,
#   plot_type = "multi"
# )

# ---- Individual visualization function (kept for single comparison use) ----
# You can still use the original function for detailed single comparison analysis
visualize_bootstrap_results <- function(comparison = "GLP1 vs SGLT2", 
                                       outcome = "Late-onset Epilepsy/Seizure",
                                       tmle_estimate = NULL) {
  
  # [Original function code remains unchanged]
  # ... [rest of the original function as provided]
}

# ---- Create Publication-Quality JAMA-Style Figure ----
create_jama_figure <- function(comparisons, 
                              outcome = "Late-onset Epilepsy/Seizure",
                              figure_width = 7,
                              figure_height = 5,
                              save_file = NULL) {
  
  # JAMA color palette
  jama_colors <- c("#00274C", "#C4622D", "#4B7AA7")
  
  # Set up high-quality graphics device if saving
  if (!is.null(save_file)) {
    if (grepl("\\.pdf$", save_file)) {
      pdf(save_file, width = figure_width, height = figure_height, family = "Helvetica")
    } else if (grepl("\\.png$", save_file)) {
      png(save_file, width = figure_width * 300, height = figure_height * 300, res = 300)
    }
  }
  
  # Set JAMA-style parameters
  old_par <- par(no.readonly = TRUE)
  par(
    family = "sans",
    mar = c(4, 4, 2, 2),
    mgp = c(2.5, 0.7, 0),
    tcl = -0.3,
    cex = 1.0,
    lwd = 1
  )
  
  # Load all bootstrap data (simplified from main function)
  all_bootstrap_data <- list()
  
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    
    rds_filename <- paste0("tmle_bootstrap_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".rds")
    
    csv_filename <- paste0("tmle_bootstrap_samples_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".csv")
    
    if (file.exists(rds_filename)) {
      bootstrap_stats <- readRDS(rds_filename)
    } else if (file.exists(csv_filename)) {
      bootstrap_df <- read.csv(csv_filename)
      bootstrap_stats <- list(
        samples = bootstrap_df$risk_difference,
        mean = mean(bootstrap_df$risk_difference),
        ci_lower_2.5 = quantile(bootstrap_df$risk_difference, 0.025),
        ci_upper_97.5 = quantile(bootstrap_df$risk_difference, 0.975)
      )
    }
    
    all_bootstrap_data[[comparison]] <- bootstrap_stats
  }
  
  # Calculate plot limits
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(c(all_samples, 0)) * 1.1
  
  # Calculate densities
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1
  
  # Create the plot
  plot(NULL, 
       xlim = x_range,
       ylim = c(0, y_max),
       xlab = "Risk Difference",
       ylab = "Density",
       type = "n",
       axes = FALSE,
       frame.plot = FALSE)
  
  # Add axes
  axis(1, col = "black", col.axis = "black", lwd = 0.75)
  axis(2, col = "black", col.axis = "black", lwd = 0.75)
  
  # Add subtle grid
  abline(h = axTicks(2), col = "gray92", lty = 1, lwd = 0.5)
  abline(v = axTicks(1), col = "gray92", lty = 1, lwd = 0.5)
  
  # Add box
  box(lwd = 0.75)
  
  # Add null effect line
  abline(v = 0, col = "black", lwd = 1.5, lty = 1)
  
  # Plot densities
  for (i in seq_along(all_bootstrap_data)) {
    comparison_name <- names(all_bootstrap_data)[i]
    dens <- densities[[i]]
    color <- jama_colors[i]
    
    # Main density line
    lines(dens, col = color, lwd = 2.5)
    
    # Shaded area
    polygon(c(dens$x, rev(dens$x)), 
            c(dens$y, rep(0, length(dens$y))),
            col = adjustcolor(color, alpha.f = 0.15),
            border = NA)
    
    # Mean line (subtle)
    abline(v = all_bootstrap_data[[comparison_name]]$mean, 
           col = color, lwd = 1, lty = 2)
  }
  
  # Add legend
  legend_text <- paste0(
    names(all_bootstrap_data),
    " (", 
    sprintf("%.3f", sapply(all_bootstrap_data, function(x) x$mean)),
    ")"
  )
  
  legend(if(mean(all_samples) > 0) "topleft" else "topright",
         legend = legend_text,
         col = jama_colors[1:length(all_bootstrap_data)],
         lwd = 2.5,
         lty = 1,
         cex = 0.8,
         bg = "white",
         box.lwd = 0.5,
         box.col = "gray60",
         title = "Comparison (Mean RD)")
  
  # Add annotation for null effect
  text(0, y_max * 0.05, "Null", pos = 3, cex = 0.7, col = "gray40")
  
  # Restore parameters
  par(old_par)
  
  # Close device if saving
  if (!is.null(save_file)) {
    dev.off()
    cat("Figure saved to:", save_file, "\n")
  }
  
  invisible(all_bootstrap_data)
}

# Example: Create a publication-quality figure
# create_jama_figure(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
#   outcome = "Late-onset Epilepsy/Seizure",
#   save_file = "jama_bootstrap_comparison.pdf"
# )

# ---- Combined Visualization for Multiple Bootstrap Comparisons ----
visualize_combined_bootstrap <- function(comparisons, 
                                        outcome = "Epilepsy/Seizure",
                                        colors = NULL,
                                        use_jama_style = TRUE,
                                        plot_type = "multi") {  # "single", "multi", or "both"
  
  # JAMA-style color palette (professional, muted colors with good contrast)
  jama_colors <- c(
    "#00274C",  # Navy blue (primary)
    "#C4622D",  # Burnt orange (secondary)
    "#4B7AA7",  # Steel blue (tertiary)
    "#969696",  # Medium gray
    "#2F7F7E"   # Teal
  )
  
  # Set default colors based on style choice
  if (is.null(colors)) {
    if (use_jama_style) {
      colors <- jama_colors
    } else {
      colors <- c("darkblue", "darkred", "darkgreen", "purple", "orange")
    }
  }
  
  # Set JAMA-style plot parameters
  if (use_jama_style) {
    # Save original parameters
    old_par <- par(no.readonly = TRUE)
    
    # Set JAMA-style parameters
    par(
      family = "sans",         # Clean sans-serif font
      cex.main = 1.1,         # Slightly smaller title
      cex.lab = 1.0,          # Axis labels
      cex.axis = 0.9,         # Axis tick labels
      font.main = 1,          # Plain (not bold) title
      las = 1,                # Horizontal axis labels
      mgp = c(2.5, 0.7, 0),   # Tighter margins
      tcl = -0.3              # Smaller tick marks
    )
  }
  
  # Storage for all bootstrap data
  all_bootstrap_data <- list()
  
  # Load data for each comparison
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    
    # Construct filenames
    rds_filename <- paste0("tmle_bootstrap_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".rds")
    
    csv_filename <- paste0("tmle_bootstrap_samples_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".csv")
    
    # Check if files exist
    if (!file.exists(rds_filename) && !file.exists(csv_filename)) {
      cat("Bootstrap files not found for", comparison, ". Attempting to download from bucket...\n")
      
      # Try to download from bucket
      my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
      system(paste0("gsutil cp ", my_bucket, "/data/", rds_filename, " ."), intern = TRUE)
      system(paste0("gsutil cp ", my_bucket, "/data/", csv_filename, " ."), intern = TRUE)
    }
    
    # Load the bootstrap results
    if (file.exists(rds_filename)) {
      bootstrap_stats <- readRDS(rds_filename)
      bootstrap_samples <- bootstrap_stats$samples
      cat("Loaded bootstrap results from RDS file for", comparison, "\n")
    } else if (file.exists(csv_filename)) {
      bootstrap_df <- read.csv(csv_filename)
      bootstrap_samples <- bootstrap_df$risk_difference
      cat("Loaded bootstrap results from CSV file for", comparison, "\n")
      
      # Recreate bootstrap_stats object
      bootstrap_stats <- list(
        samples = bootstrap_samples,
        mean = mean(bootstrap_samples),
        median = median(bootstrap_samples),
        sd = sd(bootstrap_samples),
        ci_lower_2.5 = quantile(bootstrap_samples, 0.025),
        ci_upper_97.5 = quantile(bootstrap_samples, 0.975),
        n_successful = length(bootstrap_samples),
        n_attempted = 20
      )
    } else {
      warning(paste("Could not find bootstrap results files for", comparison))
      next
    }
    
    # Store the data
    all_bootstrap_data[[comparison]] <- bootstrap_stats
  }
  
  # Check if we have any data
  if (length(all_bootstrap_data) == 0) {
    stop("No bootstrap data could be loaded")
  }
  
  # ---- Prepare plot data ----
  
  # Calculate plot limits
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(all_samples) * 1.1
  
  # Calculate densities for all comparisons
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1
  
  # ---- Create visualizations based on plot_type ----
  
  if (plot_type == "single" || plot_type == "both") {
    # Create single density plot
    plot(NULL, 
         xlim = x_range,
         ylim = c(0, y_max),
         main = if(use_jama_style) paste("Bootstrap Distributions:", outcome) else paste("Combined Bootstrap Densities:", outcome),
         xlab = "Risk Difference",
         ylab = "Density",
         type = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    # Add custom axes with JAMA style
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    
    # Add subtle grid lines (JAMA style)
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    # Add box around plot
    box(lwd = 0.5)
    
    # Add vertical line at zero (null effect)
    abline(v = 0, col = if(use_jama_style) "black" else "black", lwd = 1.5, lty = 1)
    
    # Plot density curves for each comparison
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      dens <- densities[[i]]
      color <- colors[i]
      
      # Plot density line with JAMA-appropriate line weight
      lines(dens, col = color, lwd = if(use_jama_style) 2 else 3)
      
      # Add shaded area under curve (with reduced transparency for JAMA)
      polygon(c(dens$x, rev(dens$x)), 
              c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(color, alpha.f = if(use_jama_style) 0.15 else 0.2),
              border = NA)
      
      # Add vertical line for mean (thinner for JAMA style)
      abline(v = all_bootstrap_data[[comparison_name]]$mean, 
             col = color, lwd = if(use_jama_style) 1 else 2, lty = 2)
    }
    
    # Add legend with JAMA styling
    legend(if(mean(all_samples) > 0) "topleft" else "topright",
           legend = c(paste(names(all_bootstrap_data), 
                           "\n  Mean:", 
                           sprintf("%.4f", sapply(all_bootstrap_data, function(x) x$mean))),
                     "Null effect"),
           col = c(colors[1:length(all_bootstrap_data)], "black"),
           lwd = c(rep(if(use_jama_style) 2 else 3, length(all_bootstrap_data)), 1.5),
           lty = c(rep(1, length(all_bootstrap_data)), 1),
           cex = if(use_jama_style) 0.75 else 0.8,
           bg = "white",
           box.lwd = if(use_jama_style) 0.5 else 1,
           box.col = if(use_jama_style) "gray50" else "black")
  }
  
  if (plot_type == "multi" || plot_type == "both") {
    # Create multi-panel visualization
    par(mfrow = c(2, 2))
  
    # Panel 1: Combined density plot
    plot(NULL, 
         xlim = x_range,
         ylim = c(0, y_max),
         main = if(use_jama_style) "A. Density Distributions" else "Combined Bootstrap Densities",
         xlab = "Risk Difference",
         ylab = "Density",
         type = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(v = 0, col = "black", lwd = 1.5, lty = 1)
    
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      dens <- densities[[i]]
      color <- colors[i]
      lines(dens, col = color, lwd = if(use_jama_style) 2 else 3)
      polygon(c(dens$x, rev(dens$x)), 
              c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(color, alpha.f = if(use_jama_style) 0.15 else 0.2),
              border = NA)
    }
    
    # Panel 2: Boxplots comparison
    boxplot_data <- lapply(all_bootstrap_data, function(x) x$samples)
    boxplot(boxplot_data,
            names = names(all_bootstrap_data),
            col = adjustcolor(colors[1:length(all_bootstrap_data)], alpha.f = if(use_jama_style) 0.3 else 0.5),
            border = colors[1:length(all_bootstrap_data)],
            main = if(use_jama_style) "B. Distribution Comparison" else "Bootstrap Distributions Comparison",
            ylab = "Risk Difference",
            las = 2,
            axes = FALSE,
            frame.plot = FALSE,
            lwd = if(use_jama_style) 1 else 1.5,
            boxwex = if(use_jama_style) 0.6 else 0.8,
            staplewex = 0.5,
            outwex = 0.5)
    
    axis(1, at = 1:length(all_bootstrap_data), 
         labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    
    # Panel 3: Confidence intervals comparison
    plot(NULL,
         xlim = c(0.5, length(all_bootstrap_data) + 0.5),
         ylim = x_range,
         main = if(use_jama_style) "C. 95% Confidence Intervals" else "95% Confidence Intervals",
         xlab = "",
         ylab = "Risk Difference",
         xaxt = "n",
         axes = FALSE,
         frame.plot = FALSE)
    
    axis(1, at = 1:length(all_bootstrap_data), 
         labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    }
    
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    
    for (i in seq_along(all_bootstrap_data)) {
      data <- all_bootstrap_data[[i]]
      color <- colors[i]
      
      # Draw CI with JAMA-appropriate line weight
      segments(i, data$ci_lower_2.5, i, data$ci_upper_97.5, 
               col = color, lwd = if(use_jama_style) 3 else 4)
      
      # Draw mean point
      points(i, data$mean, pch = 19, col = color, cex = if(use_jama_style) 1.5 else 2)
      
      # Draw median point (optional for JAMA style)
      if (!use_jama_style) {
        points(i, data$median, pch = 4, col = color, cex = 1.5, lwd = 2)
      }
    }
    
    # Panel 4: Statistical summary table as plot
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))
    title(if(use_jama_style) "D. Statistical Summary" else "Statistical Summary")
    
    # Create summary text with JAMA formatting
    y_pos <- 0.9
    for (i in seq_along(all_bootstrap_data)) {
      comparison_name <- names(all_bootstrap_data)[i]
      data <- all_bootstrap_data[[i]]
      color <- colors[i]
      
      text(0.05, y_pos, comparison_name, 
           col = color, font = if(use_jama_style) 1 else 2, adj = 0, cex = if(use_jama_style) 0.85 else 0.9)
      
      y_pos <- y_pos - 0.08
      text(0.1, y_pos, 
           paste("Mean:", sprintf("%.5f", data$mean)), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, 
           paste("95% CI: (", sprintf("%.5f", data$ci_lower_2.5), 
                 ", ", sprintf("%.5f", data$ci_upper_97.5), ")"), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, 
           paste("P(RD<0):", sprintf("%.3f", mean(data$samples < 0))), 
           adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      
      y_pos <- y_pos - 0.1
    }
    
    # Reset plot parameters
    par(mfrow = c(1, 1))
  }
  
  # Restore original parameters if JAMA style was used
  if (use_jama_style && exists("old_par")) {
    par(old_par)
  }
  
  # Print detailed comparison statistics
  cat("\n=== Bootstrap Distribution Comparison ===\n")
  cat("Outcome:", outcome, "\n\n")
  
  for (comparison_name in names(all_bootstrap_data)) {
    data <- all_bootstrap_data[[comparison_name]]
    cat("--- ", comparison_name, " ---\n")
    cat("  Sample size:", length(data$samples), "\n")
    cat("  Mean (SD):", round(data$mean, 5), "(", round(data$sd, 5), ")\n")
    cat("  Median:", round(data$median, 5), "\n")
    cat("  95% CI: [", round(data$ci_lower_2.5, 5), ",", 
        round(data$ci_upper_97.5, 5), "]\n")
    cat("  P(RD < 0):", round(mean(data$samples < 0), 4), "\n")
    cat("  Min/Max:", round(min(data$samples), 5), "/", 
        round(max(data$samples), 5), "\n\n")
  }
  
  # Perform pairwise comparisons if more than one group
  if (length(all_bootstrap_data) > 1) {
    cat("=== Pairwise Comparisons ===\n")
    comparison_names <- names(all_bootstrap_data)
    
    for (i in 1:(length(comparison_names) - 1)) {
      for (j in (i + 1):length(comparison_names)) {
        name1 <- comparison_names[i]
        name2 <- comparison_names[j]
        samples1 <- all_bootstrap_data[[name1]]$samples
        samples2 <- all_bootstrap_data[[name2]]$samples
        
        # Calculate difference
        mean_diff <- mean(samples1) - mean(samples2)
        
        # Approximate test (treating as independent - note this is an approximation)
        cat(name1, "vs", name2, ":\n")
        cat("  Mean difference:", round(mean_diff, 5), "\n")
        
        # Bootstrap comparison
        if (length(samples1) == length(samples2)) {
          paired_diff <- samples1 - samples2
          cat("  P(diff < 0):", round(mean(paired_diff < 0), 4), "\n")
        }
        cat("\n")
      }
    }
  }
  
  # Return all bootstrap data for further use
  invisible(all_bootstrap_data)
}

# ---- Example Usage ----

# Create only the single combined density plot (no redundancy)
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  use_jama_style = TRUE,
  plot_type = "single"  # Just the density plot
)

# Create only the 4-panel visualization
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  use_jama_style = TRUE,
  plot_type = "multi"  # Default - 4-panel figure
)

# You can also use custom JAMA-appropriate colors if needed
comparison_results <- visualize_combined_bootstrap(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  colors = c("#00274C", "#C4622D"),  # Navy blue and burnt orange
  use_jama_style = TRUE,
  plot_type = "single"
)

# Compare all three groups with just the density plot
# comparison_results <- visualize_combined_bootstrap(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4", "SGLT2 vs SU_DPP4"),
#   outcome = "Late-onset Epilepsy/Seizure",
#   use_jama_style = TRUE,
#   plot_type = "single"
# )

# For non-JAMA style (original colorful version)
# comparison_results <- visualize_combined_bootstrap(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
#   outcome = "Late-onset Epilepsy/Seizure",
#   use_jama_style = FALSE,
#   plot_type = "multi"
# )

# ---- Individual visualization function (kept for single comparison use) ----
# You can still use the original function for detailed single comparison analysis
visualize_bootstrap_results <- function(comparison = "GLP1 vs SGLT2", 
                                       outcome = "Late-onset Epilepsy/Seizure",
                                       tmle_estimate = NULL) {
  
  # [Original function code remains unchanged]
  # ... [rest of the original function as provided]
}

# ---- Create Publication-Quality JAMA-Style Figure ----
create_jama_figure <- function(comparisons, 
                              outcome = "Late-onset Epilepsy/Seizure",
                              figure_width = 7,
                              figure_height = 5,
                              save_file = NULL) {
  
  # JAMA color palette
  jama_colors <- c("#00274C", "#C4622D", "#4B7AA7")
  
  # Set up high-quality graphics device if saving
  if (!is.null(save_file)) {
    if (grepl("\\.pdf$", save_file)) {
      pdf(save_file, width = figure_width, height = figure_height, family = "Helvetica")
    } else if (grepl("\\.png$", save_file)) {
      png(save_file, width = figure_width * 300, height = figure_height * 300, res = 300)
    }
  }
  
  # Set JAMA-style parameters
  old_par <- par(no.readonly = TRUE)
  par(
    family = "sans",
    mar = c(4, 4, 2, 2),
    mgp = c(2.5, 0.7, 0),
    tcl = -0.3,
    cex = 1.0,
    lwd = 1
  )
  
  # Load all bootstrap data (simplified from main function)
  all_bootstrap_data <- list()
  
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    
    rds_filename <- paste0("tmle_bootstrap_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".rds")
    
    csv_filename <- paste0("tmle_bootstrap_samples_",
                          gsub(" ", "_", comparison), "_",
                          gsub("[/ ]", "_", outcome), ".csv")
    
    if (file.exists(rds_filename)) {
      bootstrap_stats <- readRDS(rds_filename)
    } else if (file.exists(csv_filename)) {
      bootstrap_df <- read.csv(csv_filename)
      bootstrap_stats <- list(
        samples = bootstrap_df$risk_difference,
        mean = mean(bootstrap_df$risk_difference),
        ci_lower_2.5 = quantile(bootstrap_df$risk_difference, 0.025),
        ci_upper_97.5 = quantile(bootstrap_df$risk_difference, 0.975)
      )
    }
    
    all_bootstrap_data[[comparison]] <- bootstrap_stats
  }
  
  # Calculate plot limits
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(c(all_samples, 0)) * 1.1
  
  # Calculate densities
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1
  
  # Create the plot
  plot(NULL, 
       xlim = x_range,
       ylim = c(0, y_max),
       xlab = "Risk Difference",
       ylab = "Density",
       type = "n",
       axes = FALSE,
       frame.plot = FALSE)
  
  # Add axes
  axis(1, col = "black", col.axis = "black", lwd = 0.75)
  axis(2, col = "black", col.axis = "black", lwd = 0.75)
  
  # Add subtle grid
  abline(h = axTicks(2), col = "gray92", lty = 1, lwd = 0.5)
  abline(v = axTicks(1), col = "gray92", lty = 1, lwd = 0.5)
  
  # Add box
  box(lwd = 0.75)
  
  # Add null effect line
  abline(v = 0, col = "black", lwd = 1.5, lty = 1)
  
  # Plot densities
  for (i in seq_along(all_bootstrap_data)) {
    comparison_name <- names(all_bootstrap_data)[i]
    dens <- densities[[i]]
    color <- jama_colors[i]
    
    # Main density line
    lines(dens, col = color, lwd = 2.5)
    
    # Shaded area
    polygon(c(dens$x, rev(dens$x)), 
            c(dens$y, rep(0, length(dens$y))),
            col = adjustcolor(color, alpha.f = 0.15),
            border = NA)
    
    # Mean line (subtle)
    abline(v = all_bootstrap_data[[comparison_name]]$mean, 
           col = color, lwd = 1, lty = 2)
  }
  
  # Add legend
  legend_text <- paste0(
    names(all_bootstrap_data),
    " (", 
    sprintf("%.3f", sapply(all_bootstrap_data, function(x) x$mean)),
    ")"
  )
  
  legend(if(mean(all_samples) > 0) "topleft" else "topright",
         legend = legend_text,
         col = jama_colors[1:length(all_bootstrap_data)],
         lwd = 2.5,
         lty = 1,
         cex = 0.8,
         bg = "white",
         box.lwd = 0.5,
         box.col = "gray60",
         title = "Comparison (Mean RD)")
  
  # Add annotation for null effect
  text(0, y_max * 0.05, "Null", pos = 3, cex = 0.7, col = "gray40")
  
  # Restore parameters
  par(old_par)
  
  # Close device if saving
  if (!is.null(save_file)) {
    dev.off()
    cat("Figure saved to:", save_file, "\n")
  }
  
  invisible(all_bootstrap_data)
}

# Example: Create a publication-quality figure
# create_jama_figure(
#   comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
#   outcome = "Epilepsy/Seizure",
#   save_file = "jama_bootstrap_comparison.pdf"
# )

# ---- Save Bootstrap Density Plot as PNG ----

# Function to save the density plot as high-resolution PNG
save_density_plot <- function(comparisons, 
                             outcome = "Epilepsy/Seizure",
                             filename = NULL,
                             width = 10,
                             height = 6,
                             dpi = 300,
                             use_jama_style = TRUE,
                             plot_type = "single") {
  
  # Generate default filename if not provided
  if (is.null(filename)) {
    # Create timestamp
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    
    # Clean comparison names for filename
    comparison_str <- gsub(" ", "_", paste(comparisons, collapse="_vs_"))
    outcome_str <- gsub("[/ ]", "_", outcome)
    
    # Construct filename
    filename <- paste0("bootstrap_density_",
                      comparison_str, "_",
                      outcome_str, "_",
                      timestamp, ".png")
  }
  
  # Start PNG device with high resolution
  png(filename = filename, 
      width = width * dpi, 
      height = height * dpi, 
      res = dpi,
      type = "cairo-png")  # Use cairo for better quality if available
  
  # Generate the plot
  comparison_results <- visualize_combined_bootstrap(
    comparisons = comparisons,
    outcome = outcome,
    use_jama_style = use_jama_style,
    plot_type = plot_type
  )
  
  # Close the device
  dev.off()
  
  # Print confirmation
  cat("\n=== Plot saved successfully ===\n")
  cat("Filename:", filename, "\n")
  cat("Dimensions:", width, "x", height, "inches\n")
  cat("Resolution:", dpi, "DPI\n")
  cat("File size:", round(file.size(filename)/1024), "KB\n\n")
  
  # Return the filename for reference
  invisible(filename)
}

# ---- Quick save function (simpler version) ----
quick_save_plot <- function(comparisons, 
                           outcome = "Epilepsy/Seizure",
                           prefix = "figure") {
  
  # Simple filename with prefix
  filename <- paste0(prefix, "_bootstrap_density.png")
  
  # Save with default high-quality settings
  png(filename = filename, 
      width = 3000, 
      height = 1800, 
      res = 300)
  
  visualize_combined_bootstrap(
    comparisons = comparisons,
    outcome = outcome,
    use_jama_style = TRUE,
    plot_type = "single"
  )
  
  dev.off()
  
  cat("Plot saved as:", filename, "\n")
  invisible(filename)
}

# ---- USAGE EXAMPLES ----

# Method 1: Save with detailed control
saved_file <- save_density_plot(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  filename = "Figure_2_bootstrap_comparison.png",
  width = 10,      # Width in inches
  height = 6,      # Height in inches
  dpi = 300,       # Resolution for publication
  plot_type = "single"  # Just the density plot
)

# Method 2: Save with auto-generated filename
saved_file <- save_density_plot(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  plot_type = "single"
)

# Method 3: Quick save with simple filename
quick_save_plot(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  prefix = "manuscript_fig2"
)

# Method 4: Save the 4-panel figure
saved_file <- save_density_plot(
  comparisons = c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
  outcome = "Epilepsy/Seizure",
  filename = "bootstrap_4panel_analysis_epilepsy_seizure.png",
  width = 12,      # Wider for 4 panels
  height = 10,     # Taller for 4 panels
  plot_type = "multi"
)

# ---- Alternative: Save current plot after it's already displayed ----
# If you've already run the visualization and want to save what's currently displayed:

save_current_plot <- function(filename = "current_plot.png", 
                             width = 10, 
                             height = 6, 
                             dpi = 300) {
  
  # Save the current device content
  dev.copy(png, 
           filename = filename,
           width = width * dpi,
           height = height * dpi,
           res = dpi)
  dev.off()
  
  cat("Current plot saved as:", filename, "\n")
  invisible(filename)
}

# Usage: After running your visualization, just call:
# save_current_plot("my_density_plot.png")

# ---- Batch save multiple comparisons ----
save_all_comparisons <- function(outcome = "Epilepsy/Seizure") {
  
  # Define all comparison pairs
  all_comparisons <- list(
    c("GLP1 vs SGLT2", "GLP1 vs SU_DPP4"),
    c("GLP1 vs SGLT2", "SGLT2 vs SU_DPP4"),
    c("GLP1 vs SU_DPP4", "SGLT2 vs SU_DPP4")
  )
  
  # Save each comparison
  for (i in seq_along(all_comparisons)) {
    comparisons <- all_comparisons[[i]]
    
    filename <- paste0("bootstrap_comparison_", i, ".png")
    
    save_density_plot(
      comparisons = comparisons,
      outcome = outcome,
      filename = filename,
      plot_type = "single"
    )
  }
  
  cat("\nAll comparison plots saved successfully!\n")
}

# Usage:
# save_all_comparisons()