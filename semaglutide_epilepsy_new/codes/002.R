library(tidyverse)
library(bigrquery)

# This snippet assumes that you run setup first

# This code copies a file from your Google Bucket into a dataframe

# replace 'test.csv' with the name of the file in your google bucket (don't delete the quotation marks)
name_of_file_in_bucket <- 'mother_df.csv'

########################################################################
##
################# DON'T CHANGE FROM HERE ###############################
##
########################################################################

# Get the bucket name
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')

# Copy the file from current workspace to the bucket
system(paste0("gsutil cp ", my_bucket, "/data/", name_of_file_in_bucket, " ."), intern=T)

# Load the file into a dataframe
mother_df <- read_csv(name_of_file_in_bucket)

# This query represents dataset Basic Survey for domain "survey" and was generated for All of Us Controlled Tier Dataset v7
dataset_04958613_survey_sql <- paste("
    SELECT
        answer.person_id,
        answer.survey_datetime,
        answer.survey,
        answer.question_concept_id,
        answer.question,
        answer.answer_concept_id,
        answer.answer,
        answer.survey_version_concept_id,
        answer.survey_version_name  
    FROM
        `ds_survey` answer   
    WHERE
        (
            question_concept_id IN (SELECT
                DISTINCT concept_id                         
            FROM
                `cb_criteria` c                         
            JOIN
                (SELECT
                    CAST(cr.id as string) AS id                               
                FROM
                    `cb_criteria` cr                               
                WHERE
                    concept_id IN (1586134)                               
                    AND domain_id = 'SURVEY') a 
                    ON (c.path like CONCAT('%', a.id, '.%'))                         
            WHERE
                domain_id = 'SURVEY'                         
                AND type = 'PPI'                         
                AND subtype = 'QUESTION')
        )", sep="")

# Formulate a Cloud Storage destination path for the data exported from BigQuery.
# NOTE: By default data exported multiple times on the same day will overwrite older copies.
#       But data exported on a different days will write to a new location so that historical
#       copies can be kept as the dataset definition is changed.
survey_04958613_path <- file.path(
  Sys.getenv("WORKSPACE_BUCKET"),
  "bq_exports",
  Sys.getenv("OWNER_EMAIL"),
  strftime(lubridate::now(), "%Y%m%d"),  # Comment out this line if you want the export to always overwrite.
  "survey_04958613",
  "survey_04958613_*.csv")
message(str_glue('The data will be written to {survey_04958613_path}. Use this path when reading ',
                 'the data into your notebooks in the future.'))

# Perform the query and export the dataset to Cloud Storage as CSV files.
# NOTE: You only need to run `bq_table_save` once. After that, you can
#       just read data from the CSVs in Cloud Storage.
bq_table_save(
  bq_dataset_query(Sys.getenv("WORKSPACE_CDR"), dataset_04958613_survey_sql, billing = Sys.getenv("GOOGLE_PROJECT")),
  survey_04958613_path,
  destination_format = "CSV")

# Read the data directly from Cloud Storage into memory.
# NOTE: Alternatively you can `gsutil -m cp {survey_04958613_path}` to copy these files
#       to the Jupyter disk.
read_bq_export_from_workspace_bucket <- function(export_path) {
  col_types <- cols(survey = col_character(), question = col_character(), answer = col_character(), survey_version_name = col_character())
  bind_rows(
    map(system2('gsutil', args = c('ls', export_path), stdout = TRUE, stderr = TRUE),
        function(csv) {
          message(str_glue('Loading {csv}.'))
          chunk <- read_csv(pipe(str_glue('gsutil cat {csv}')), col_types = col_types, show_col_types = FALSE)
          if (is.null(col_types)) {
            col_types <- spec(chunk)
          }
          chunk
        }))
}
dataset_04958613_survey_df <- read_bq_export_from_workspace_bucket(survey_04958613_path)

# Rename and extract the date from survey_datetime
survey_date_summary <- dataset_04958613_survey_df %>%
  mutate(survey_date = as_date(survey_datetime)) %>% # Extract date
  group_by(person_id) %>%
  summarise(
    first_survey_date = min(survey_date, na.rm = TRUE),
    last_survey_date = max(survey_date, na.rm = TRUE)
  )

download_data <- function(query) {
    tb <- bq_project_query(Sys.getenv('GOOGLE_PROJECT'), query)
    bq_table_download(tb)
}

dataset = Sys.getenv('WORKSPACE_CDR')

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRMEA' AS EHRTYPE, MIN(measurement_datetime) AS MINehrDATE, 
    MAX(measurement_datetime) AS MAXehrDATE
FROM`{dataset}.measurement` AS m
LEFT JOIN `{dataset}.measurement_ext` AS mm on m.measurement_id = mm.measurement_id
WHERE 
LOWER(mm.src_id) LIKE 'ehr site%' 
GROUP BY person_id
")
                  
measdf <- suppressMessages(download_data(query))

measdf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRCON' AS EHRTYPE, MIN(condition_start_datetime) AS MINehrDATE, 
    MAX(condition_start_datetime) AS MAXehrDATE
FROM`{dataset}.condition_occurrence` AS m
LEFT JOIN `{dataset}.condition_occurrence_ext` AS mm on m.condition_occurrence_id = mm.condition_occurrence_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%' 
GROUP BY person_id
")
                  
condf <- suppressMessages(download_data(query))

condf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRDEV' AS EHRTYPE, MIN(device_exposure_start_date) AS MINehrDATE, 
    MAX(device_exposure_start_date) AS MAXehrDATE
FROM`{dataset}.device_exposure` AS m
LEFT JOIN `{dataset}.device_exposure_ext` AS mm on m.device_exposure_id = mm.device_exposure_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%'
GROUP BY person_id
")
                  
devdf <- suppressMessages(download_data(query))

devdf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRDRG' AS EHRTYPE, MIN(drug_exposure_start_date) AS MINehrDATE, 
    MAX(drug_exposure_start_date) AS MAXehrDATE
FROM`{dataset}.drug_exposure` AS m
LEFT JOIN `{dataset}.drug_exposure_ext` AS mm on m.drug_exposure_id = mm.drug_exposure_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%'
GROUP BY person_id
")
                  
drgdf <- suppressMessages(download_data(query))

drgdf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHROBS' AS EHRTYPE, MIN(observation_datetime) AS MINehrDATE, 
    MAX(observation_datetime) AS MAXehrDATE
FROM`{dataset}.observation` AS m
LEFT JOIN `{dataset}.observation_ext` AS mm on m.observation_id = mm.observation_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%'
GROUP BY person_id
")
                  
obsdf <- suppressMessages(download_data(query))

obsdf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRPRO' AS EHRTYPE, MIN(procedure_datetime) AS MINehrDATE, 
    MAX(procedure_datetime) AS MAXehrDATE
FROM`{dataset}.procedure_occurrence` AS m
LEFT JOIN `{dataset}.procedure_occurrence_ext` AS mm on m.procedure_occurrence_id = mm.procedure_occurrence_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%'
GROUP BY person_id
")
                  
prodf <- suppressMessages(download_data(query))

prodf$EHR <- 1

query <- str_glue("
SELECT
   DISTINCT person_id, 'EHRVIS' AS EHRTYPE, MIN(visit_start_datetime) AS MINehrDATE, 
    MAX(visit_start_datetime) AS MAXehrDATE
FROM`{dataset}.visit_occurrence` AS m
LEFT JOIN `{dataset}.visit_occurrence_ext` AS mm on m.visit_occurrence_id = mm.visit_occurrence_id
WHERE LOWER(mm.src_id) LIKE 'ehr site%'
GROUP BY person_id
")
                  
visdf <- suppressMessages(download_data(query))

visdf$EHR <- 1

ehrDTdf <- rbind.data.frame(measdf, condf, devdf, drgdf, obsdf, prodf, visdf)

ehrDTdf$EHRdtMIN <- substr(ehrDTdf$MINehrDATE, 1, 10)
ehrDTdf$EHRdtMIN <- as.Date(ehrDTdf$EHRdtMIN, format="%Y-%m-%d")

ehrDTdf$EHRdtMAX <- substr(ehrDTdf$MAXehrDATE, 1, 10)
ehrDTdf$EHRdtMAX <- as.Date(ehrDTdf$EHRdtMAX, format="%Y-%m-%d")

ehrDTdf1 <- ehrDTdf %>%
group_by(person_id) %>%
summarize(EHRTYPECOUNT=n(),
          EHR=1,
          EHRminDT=min(EHRdtMIN, na.rm=TRUE),
          EHRmaxDT=min(EHRdtMAX, na.rm=TRUE),
          EHRMEA=sum(EHRTYPE=='EHRMEA'),
          EHRCON=sum(EHRTYPE=='EHRCON'),
          EHRDEV=sum(EHRTYPE=='EHRDEV'),
          EHRDRG=sum(EHRTYPE=='EHRDRG'),
          EHROBS=sum(EHRTYPE=='EHROBS'),
          EHRPRO=sum(EHRTYPE=='EHRPRO'), 
          EHRVIS=sum(EHRTYPE=='EHRVIS')) 
colnames(ehrDTdf1)
ehrDTdf1$EHRMEA[is.na(ehrDTdf1$EHRMEA)] <- 0
ehrDTdf1$EHRCON[is.na(ehrDTdf1$EHRCON)] <- 0
ehrDTdf1$EHRDEV[is.na(ehrDTdf1$EHRDEV)] <- 0
ehrDTdf1$EHRDRG[is.na(ehrDTdf1$EHRDRG)] <- 0
ehrDTdf1$EHROBS[is.na(ehrDTdf1$EHROBS)] <- 0
ehrDTdf1$EHRPRO[is.na(ehrDTdf1$EHRPRO)] <- 0
ehrDTdf1$EHRVIS[is.na(ehrDTdf1$EHRVIS)] <- 0

# Merge df with ehrDTdf1, keeping all columns in df and adding EHRminDT and EHRmaxDT
merged_df <- mother_df %>%
  left_join(ehrDTdf1 %>% select(person_id, EHRminDT, EHRmaxDT), by = "person_id") %>%
  left_join(survey_date_summary %>% select(person_id, first_survey_date), by = "person_id")

# Save the filtered dataframe if needed
# write.csv(filtered_condition_summary, "filtered_condition_summary.csv", row.names = FALSE)

# This code saves your dataframe into a csv file in a "data" folder in Google Bucket

# Replace df with THE NAME OF YOUR DATAFRAME
my_dataframe <- merged_df

# Replace 'test.csv' with THE NAME of the file you're going to store in the bucket (don't delete the quotation marks)
destination_filename <- 'merged_df.csv'

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