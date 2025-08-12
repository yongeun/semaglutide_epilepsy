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

#── Panel‐BMI SQL (Modified to include both EHR and PPI data) ────────────────
# Standard concept IDs:
# - Height: 3036277
# - Weight: 3025315  
# - BMI: 3038553

panel_sql <- str_glue("
WITH height_data AS (
  SELECT
    m.person_id,
    m.measurement_datetime,
    m.value_as_number AS height_cm,
    m.measurement_source_concept_id,
    m_ext.src_id AS data_source,
    ROW_NUMBER() OVER (
      PARTITION BY m.person_id
      ORDER BY m.measurement_datetime DESC
    ) AS rn
  FROM `{cdr}.measurement` m
  LEFT JOIN `{cdr}.measurement_ext` m_ext
    ON m.measurement_id = m_ext.measurement_id
  WHERE m.measurement_concept_id = 3036277  -- Standard concept for height
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0  -- Ensure height is positive
    AND m.unit_concept_id = 8582  -- centimeter
),
most_recent_height AS (
  SELECT person_id, height_cm, data_source AS height_source
  FROM height_data
  WHERE rn = 1
),
weight_data AS (
  SELECT
    m.person_id,
    m.measurement_datetime AS weight_datetime,
    m.value_as_number AS weight_kg,
    m.measurement_source_concept_id,
    m_ext.src_id AS data_source,
    -- Get the most recent height for this person
    h.height_cm,
    h.height_source,
    -- Calculate BMI with safety check
    CASE 
      WHEN h.height_cm IS NOT NULL AND h.height_cm > 0 THEN
        ROUND(m.value_as_number / POWER(h.height_cm/100, 2), 2)
      ELSE NULL
    END AS calculated_bmi
  FROM `{cdr}.measurement` m
  LEFT JOIN `{cdr}.measurement_ext` m_ext
    ON m.measurement_id = m_ext.measurement_id
  INNER JOIN most_recent_height h
    ON m.person_id = h.person_id
  WHERE m.measurement_concept_id = 3025315  -- Standard concept for weight
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0  -- Ensure weight is positive
    AND m.unit_concept_id = 9529  -- kilogram
),
direct_bmi AS (
  SELECT
    m.person_id,
    m.measurement_datetime AS bmi_datetime,
    m.value_as_number AS direct_bmi,
    m.measurement_source_concept_id,
    m_ext.src_id AS data_source
  FROM `{cdr}.measurement` m
  LEFT JOIN `{cdr}.measurement_ext` m_ext
    ON m.measurement_id = m_ext.measurement_id
  WHERE m.measurement_concept_id = 3038553  -- Standard concept for BMI
    AND m.value_as_number IS NOT NULL
    AND m.value_as_number > 0  -- Ensure BMI is positive
)

-- Combine calculated BMI from weight/height and direct BMI measurements
SELECT
  COALESCE(w.person_id, b.person_id) AS person_id,
  COALESCE(w.weight_datetime, b.bmi_datetime) AS measurement_date,
  w.weight_kg,
  w.height_cm,
  COALESCE(w.calculated_bmi, b.direct_bmi) AS bmi,
  CASE
    WHEN w.data_source IN ('Staff Portal: HealthPro', 'Participant Portal: PTSC') THEN 'PPI'
    WHEN b.data_source IN ('Staff Portal: HealthPro', 'Participant Portal: PTSC') THEN 'PPI'
    ELSE 'EHR'
  END AS source_type,
  COALESCE(w.data_source, b.data_source) AS detailed_source
FROM weight_data w
FULL OUTER JOIN direct_bmi b
  ON w.person_id = b.person_id
  AND DATE(w.weight_datetime) = DATE(b.bmi_datetime)

-- Include participants who have at least height data (for panel analysis)
WHERE COALESCE(w.person_id, b.person_id) IN (
  SELECT DISTINCT person_id
  FROM `{cdr}.cb_search_all_events`
  WHERE is_standard = 1
    AND concept_id IN (3036277, 3025315, 3038553)  -- Height, Weight, or BMI
)
  AND COALESCE(w.calculated_bmi, b.direct_bmi) IS NOT NULL  -- Ensure we have a valid BMI

ORDER BY person_id, measurement_date
")

#── Download & Tidy ──────────────────────────────────────────────────────────
cat("Downloading BMI panel data (including both PPI and EHR sources)...\n")
bmi_panel <- download_data(panel_sql) %>%
  mutate(
    weight_date = as_date(measurement_date)
  ) %>%
  select(
    person_id,
    weight_date,
    weight_kg,
    height_cm,
    bmi,
    source_type,
    detailed_source
  ) %>%
  arrange(person_id, weight_date)

#── Quick Summary by Source ──────────────────────────────────────────────────
cat("\n── Data Source Summary ──\n")
source_summary <- bmi_panel %>%
  group_by(source_type) %>%
  summarise(
    n_measurements = n(),
    n_persons = n_distinct(person_id),
    .groups = "drop"
  )
print(source_summary)

cat("\n── Detailed Source Summary ──\n")
detailed_source_summary <- bmi_panel %>%
  filter(!is.na(detailed_source)) %>%
  group_by(source_type, detailed_source) %>%
  summarise(
    n_measurements = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(n_measurements)) %>%
  head(10)
print(detailed_source_summary)

#── Quick Sanity Check ───────────────────────────────────────────────────────
bmi_panel %>%
  group_by(person_id) %>%
  summarize(
    n_obs     = n(),
    first_obs = min(weight_date),
    last_obs  = max(weight_date),
    has_ppi   = any(source_type == "PPI"),
    has_ehr   = any(source_type == "EHR"),
    .groups = "drop"
  ) %>%
  summarise(
    total_persons = n(),
    persons_with_ppi_only = sum(has_ppi & !has_ehr),
    persons_with_ehr_only = sum(!has_ppi & has_ehr),
    persons_with_both = sum(has_ppi & has_ehr),
    median_obs_per_person = median(n_obs)
  ) %>%
  print()

# define plausible ranges
min_weight <- 30
max_weight <- 250
min_height <- 130
max_height <- 210
min_bmi    <- 15
max_bmi    <- 70

# filter
bmi_panel_clean <- bmi_panel %>%
  filter(
    is.na(weight_kg) | (weight_kg >= min_weight & weight_kg <= max_weight),
    is.na(height_cm) | (height_cm >= min_height & height_cm <= max_height),
    bmi >= min_bmi & bmi <= max_bmi
  )

n_before <- nrow(bmi_panel)
n_after  <- nrow(bmi_panel_clean)

cat("\nRows before filtering:", n_before, "\n")
cat("Rows after filtering:",  n_after,  "\n")
cat("Percent kept:", round(100 * n_after / n_before, 1), "%\n")

cat("\n── Filtered Data by Source ──\n")
filtered_source_summary <- bmi_panel_clean %>%
  group_by(source_type) %>%
  summarise(
    n_measurements = n(),
    n_persons = n_distinct(person_id),
    mean_bmi = round(mean(bmi, na.rm = TRUE), 1),
    median_bmi = round(median(bmi, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(filtered_source_summary)

# This snippet assumes that you run setup first

# This code saves your dataframe into a csv file in a "data" folder in Google Bucket

# Replace df with THE NAME OF YOUR DATAFRAME
my_dataframe <- bmi_panel_clean

# Replace 'test.csv' with THE NAME of the file you're going to store in the bucket (don't delete the quotation marks)
destination_filename <- 'bmi_panel.csv'

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
