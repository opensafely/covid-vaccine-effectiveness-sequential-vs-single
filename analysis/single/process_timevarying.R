######################################

# What this script does:
# imports data created by the `data_process.R` script
# converts hospitalisation and infection episodes into long format
# saves as a one-row-per-event dataset

# don't include probable or suspected covid in time since, as this is not included for seqtrials timesince

######################################

# Import libraries 
library('tidyverse')
library('here')
library('glue')

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# create output directories 
outdir <- here("output", "single", "process")
fs::dir_create(outdir)

# import timevarying extract 
data_extract <- arrow::read_feather(here("output", "single", "extract", "input_timevarying.feather")) %>%
  mutate(across(ends_with("_date"),  as.Date))

# import eligible data 
data_eligible <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds"))

# select and save outcomes
data_outcomes <- data_extract %>%
  select(patient_id, dereg_date, postest_date, covidadmitted_date, coviddeath_date, death_date)
write_rds(
  data_outcomes,
  file.path(outdir, "data_outcomes.rds"), 
  compress="gz"
  )

# join data_eligible and data extract so that we have the events
# for the timevarying variables that occur both before (data_eligible) 
# and after (data_timevarying) trial_date
data_timevarying <- data_extract %>%
  left_join(
    data_eligible %>% 
      select(
        patient_id,
        starts_with(c("positive", "admitted", "discharged", "covidemergency"))
        ), 
    by = "patient_id"
  )

## create one-row-per-event datasets ----
# for vaccination, positive test, hospitalisation/discharge, covid in primary care, death

data_admissions <- data_timevarying %>%
  select(
    patient_id,
    matches("^admitted\\_unplanned\\_\\d+\\_date"),
    matches("^discharged\\_unplanned\\_\\d+\\_date")
    ) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(".value", "index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_drop_na = TRUE
  ) %>%
  select(
    patient_id, index, 
    admitted_date=admitted_unplanned, 
    discharged_date = discharged_unplanned
    ) %>%
  arrange(patient_id, admitted_date)

data_admissions_infectious <- data_timevarying %>%
  select(
    patient_id, 
    matches("^admitted\\_unplanned\\_infectious\\_\\d+\\_date"), 
    matches("^discharged\\_unplanned\\_infectious\\_\\d+\\_date")
    ) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(".value", "index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_drop_na = TRUE
  ) %>%
  select(
    patient_id, index,
    admitted_date=admitted_unplanned_infectious, 
    discharged_date = discharged_unplanned_infectious
    ) %>%
  arrange(patient_id, admitted_date)

# remove infectious admissions from all admissions data
data_admissions_noninfectious <- anti_join(
  data_admissions,
  data_admissions_infectious,
  by = c("patient_id", "admitted_date", "discharged_date")
)


# data_pr_suspected_covid <- data_processed %>%
#   select(patient_id, matches("^primary_care_suspected_covid\\_\\d+\\_date")) %>%
#   pivot_longer(
#     cols = -patient_id,
#     names_to = c(NA, "suspected_index"),
#     names_pattern = "^(.*)_(\\d+)_date",
#     values_to = "date",
#     values_drop_na = TRUE
#   ) %>%
#   arrange(patient_id, date)
# 
# data_pr_probable_covid <- data_processed %>%
#   select(patient_id, matches("^primary_care_probable_covid\\_\\d+\\_date")) %>%
#   pivot_longer(
#     cols = -patient_id,
#     names_to = c(NA, "probable_index"),
#     names_pattern = "^(.*)_(\\d+)_date",
#     values_to = "date",
#     values_drop_na = TRUE
#   ) %>%
#   arrange(patient_id, date)
# 
data_postest <- data_timevarying %>%
  select(
    patient_id, 
    matches("^positive\\_test\\_\\d+\\_date")
    ) %>%
  pivot_longer(
    cols = -patient_id,
    names_to = c(NA, "postest_index"),
    names_pattern = "^(.*)_(\\d+)_date",
    values_to = "date",
    values_drop_na = TRUE
  ) %>%
  arrange(patient_id, date)

write_rds(
  data_admissions, 
  file.path(outdir, "data_long_admission_dates.rds"),
  compress="gz"
  )
write_rds(
  data_admissions_infectious, 
  file.path(outdir, "data_long_admission_infectious_dates.rds"), 
  compress="gz"
  )
write_rds(
  data_admissions_noninfectious, 
  file.path(outdir, "data_long_admission_noninfectious_dates.rds"),
  compress="gz"
  )
# write_rds(
#   data_pr_probable_covid, 
#   file.path(outdir, "data_long_pr_probable_covid_dates.rds"), 
#   compress="gz"
#   )
# write_rds(
#   data_pr_suspected_covid, 
#   file.path(outdir, "data_long_pr_suspected_covid_dates.rds"), 
#   compress="gz"
#   )
write_rds(
  data_postest, 
  file.path(outdir, "data_long_postest_dates.rds"), 
  compress="gz"
  )
