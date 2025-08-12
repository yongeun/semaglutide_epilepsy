library(bigrquery)
library(tidyverse)
library(lubridate)
library(glue)

#── Configuration ────────────────────────────────────────────────────────────
# get configuration info
DATASET <- Sys.getenv('WORKSPACE_CDR')    # e.g. "all_of_us_ehr_registered_2025q1"
Sys.setenv(GOOGLE_PROJECT = Sys.getenv('GOOGLE_PROJECT'))  # ensure billing project is set

# a function to extract data from BigQuery
download_data <- function(sql) {
  tb <- bq_project_query(Sys.getenv('GOOGLE_PROJECT'), sql)
  bq_table_download(tb, bigint = "integer64")
}

cdr <- DATASET  # alias for your dataset in glue()

#── Panel‐HbA1c SQL (Modified from BMI to extract HbA1c data) ────────────────
# OMOP measurement IDs for HbA1c: 3004410, 3005673, 3003309, 3007263

panel_sql <- str_glue("
WITH hba1c_data AS (
  SELECT
    m.person_id,
    m.measurement_datetime,
    m.value_as_number AS hba1c_value,
    m.unit_concept_id,
    m.measurement_concept_id,
    m.measurement_source_concept_id,
    m_ext.src_id AS data_source,
    -- Convert to standard percentage if needed
    CASE 
      WHEN m.unit_concept_id = 8554 THEN m.value_as_number  -- % (already percentage)
      WHEN m.unit_concept_id = 8840 THEN m.value_as_number / 10.929  -- mmol/mol to %
      ELSE m.value_as_number  -- assume percentage if unit unknown
    END AS hba1c_percent
  FROM `{cdr}.measurement` m
  LEFT JOIN `{cdr}.measurement_ext` m_ext
    ON m.measurement_id = m_ext.measurement_id
  WHERE m.measurement_concept_id IN (3004410, 3005673, 3003309, 3007263)  -- HbA1c concepts
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0  -- Ensure HbA1c is positive
)

-- Select and structure the data
SELECT
  person_id,
  measurement_datetime AS date_of_measurement,
  hba1c_value AS original_value,
  unit_concept_id,
  measurement_concept_id,
  hba1c_percent,
  CASE
    WHEN data_source IN ('Staff Portal: HealthPro', 'Participant Portal: PTSC') THEN 'PPI'
    ELSE 'EHR'
  END AS source_type,
  data_source AS detailed_source
FROM hba1c_data

-- Include participants who have at least one HbA1c measurement
WHERE person_id IN (
  SELECT DISTINCT person_id
  FROM `{cdr}.cb_search_all_events`
  WHERE is_standard = 1
    AND concept_id IN (3004410, 3005673, 3003309, 3007263)  -- HbA1c concepts
)

ORDER BY person_id, date_of_measurement
")

#── Download & Tidy ──────────────────────────────────────────────────────────
cat("Downloading HbA1c panel data (including both PPI and EHR sources)...\n")
hba1c_panel <- download_data(panel_sql) %>%
  mutate(
    date_of_measurement = as_date(date_of_measurement)
  ) %>%
  select(
    person_id,
    date_of_measurement,
    hba1c_value = hba1c_percent,
    original_value,
    unit_concept_id,
    measurement_concept_id,
    source_type,
    detailed_source
  ) %>%
  arrange(person_id, date_of_measurement)

#── Quick Summary by Source ──────────────────────────────────────────────────
cat("\n── Data Source Summary ──\n")
source_summary <- hba1c_panel %>%
  group_by(source_type) %>%
  summarise(
    n_measurements = n(),
    n_persons = n_distinct(person_id),
    mean_hba1c = round(mean(hba1c_value, na.rm = TRUE), 2),
    median_hba1c = round(median(hba1c_value, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(source_summary)

cat("\n── Detailed Source Summary ──\n")
detailed_source_summary <- hba1c_panel %>%
  filter(!is.na(detailed_source)) %>%
  group_by(source_type, detailed_source) %>%
  summarise(
    n_measurements = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_measurements)) %>%
  head(10)
print(detailed_source_summary)

#── Measurement Concept Summary ──────────────────────────────────────────────
cat("\n── HbA1c Measurement Concepts Used ──\n")
concept_summary <- hba1c_panel %>%
  group_by(measurement_concept_id) %>%
  summarise(
    n_measurements = n(),
    mean_value = round(mean(hba1c_value, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  mutate(
    concept_name = case_when(
      measurement_concept_id == 3004410 ~ "Hemoglobin A1c/Hemoglobin.total in Blood",
      measurement_concept_id == 3005673 ~ "Hemoglobin A1c [Mass/volume] in Blood",  
      measurement_concept_id == 3003309 ~ "Hemoglobin A1c/Hemoglobin.total in Blood by IFCC protocol",
      measurement_concept_id == 3007263 ~ "Hemoglobin A1c/Hemoglobin.total in Blood by Electrophoresis",
      TRUE ~ "Unknown"
    )
  )
print(concept_summary)

#── Quick Sanity Check ───────────────────────────────────────────────────────
hba1c_panel %>%
  group_by(person_id) %>%
  summarize(
    n_obs     = n(),
    first_obs = min(date_of_measurement),
    last_obs  = max(date_of_measurement),
    has_ppi   = any(source_type == "PPI"),
    has_ehr   = any(source_type == "EHR"),
    mean_hba1c = round(mean(hba1c_value, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  summarise(
    total_persons = n(),
    persons_with_ppi_only = sum(has_ppi & !has_ehr),
    persons_with_ehr_only = sum(!has_ppi & has_ehr),
    persons_with_both = sum(has_ppi & has_ehr),
    median_obs_per_person = median(n_obs),
    mean_hba1c_overall = round(mean(mean_hba1c, na.rm = TRUE), 2)
  ) %>%
  print()

# define plausible ranges for HbA1c (percentage)
min_hba1c <- 3.0   # Very low but possible
max_hba1c <- 20.0  # Very high but possible in extreme cases

# filter
hba1c_panel_clean <- hba1c_panel %>%
  filter(
    hba1c_value >= min_hba1c & hba1c_value <= max_hba1c
  )

n_before <- nrow(hba1c_panel)
n_after  <- nrow(hba1c_panel_clean)

cat("\nRows before filtering:", n_before, "\n")
cat("Rows after filtering:",  n_after,  "\n")
cat("Percent kept:", round(100 * n_after / n_before, 1), "%\n")

cat("\n── Filtered Data by Source ──\n")
filtered_source_summary <- hba1c_panel_clean %>%
  group_by(source_type) %>%
  summarise(
    n_measurements = n(),
    n_persons = n_distinct(person_id),
    mean_hba1c = round(mean(hba1c_value, na.rm = TRUE), 2),
    median_hba1c = round(median(hba1c_value, na.rm = TRUE), 2),
    sd_hba1c = round(sd(hba1c_value, na.rm = TRUE), 2),
    .groups = "drop"
  )
print(filtered_source_summary)

# Distribution check
cat("\n── HbA1c Distribution ──\n")
hba1c_dist <- hba1c_panel_clean %>%
  summarise(
    min = round(min(hba1c_value, na.rm = TRUE), 2),
    q1 = round(quantile(hba1c_value, 0.25, na.rm = TRUE), 2),
    median = round(median(hba1c_value, na.rm = TRUE), 2),
    q3 = round(quantile(hba1c_value, 0.75, na.rm = TRUE), 2),
    max = round(max(hba1c_value, na.rm = TRUE), 2),
    n_below_5.7 = sum(hba1c_value < 5.7),  # Normal
    n_5.7_to_6.4 = sum(hba1c_value >= 5.7 & hba1c_value < 6.5),  # Prediabetes
    n_6.5_and_above = sum(hba1c_value >= 6.5)  # Diabetes
  )
print(hba1c_dist)

# This snippet assumes that you run setup first

# This code saves your dataframe into a csv file in a "data" folder in Google Bucket

# Select and rename columns for final output to match expected structure
final_a1c_panel <- hba1c_panel_clean %>%
  select(
    person_id,
    date_of_measurement,
    A1c = hba1c_value
  )

# Replace df with THE NAME OF YOUR DATAFRAME
my_dataframe <- final_a1c_panel

# Replace 'test.csv' with THE NAME of the file you're going to store in the bucket (don't delete the quotation marks)
destination_filename <- 'a1c_panel.csv'

########################################################################
##
################# DON'T CHANGE FROM HERE ###############################
##
########################################################################

# store the dataframe in current workspace
write_excel_csv(my_dataframe, destination_filename)

# Get the bucket name
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

# Copy the file from current workspace to the bucket
system(paste0("gsutil cp ./", destination_filename, " ", my_bucket, "/data/"), intern=T)

# Check if file is in the bucket
system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern=T)