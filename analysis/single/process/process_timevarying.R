# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# imports data created by the `data_process.R` script
# converts hospitalisation episodes into long format
# saves as a one-row-per-event dataset
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# import libraries 
library('tidyverse')
library('here')
library('glue')
library('arrow')

# import local functions and parameters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# create output directories 
outdir <- here("output", "single", "process")
fs::dir_create(outdir)

# import eligible data 
data_eligible <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  studydef_path <- here("output", "single", "extract", "input_timevarying.feather")
  custom_path <- here("output", "single", "dummydata", "dummy_timevarying.feather")
  
  data_studydef_dummy <- read_feather(studydef_path) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
    # because of a bug in cohort extractor -- remove once pulled new version
    mutate(patient_id = as.integer(patient_id))
  
  data_custom_dummy <- read_feather(custom_path) 
  
  not_in_studydef <- names(data_custom_dummy)[!( names(data_custom_dummy) %in% names(data_studydef_dummy) )]
  not_in_custom  <- names(data_studydef_dummy)[!( names(data_studydef_dummy) %in% names(data_custom_dummy) )]
  
  if(length(not_in_custom)!=0) stop(
    paste(
      "These variables are in studydef but not in custom: ",
      paste(not_in_custom, collapse=", ")
    )
  )
  
  if(length(not_in_studydef)!=0) stop(
    paste(
      "These variables are in custom but not in studydef: ",
      paste(not_in_studydef, collapse=", ")
    )
  )
  
  # reorder columns
  data_studydef_dummy <- data_studydef_dummy[,names(data_custom_dummy)]
  
  unmatched_types <- cbind(
    map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")),
    map_chr(data_custom_dummy, ~paste(class(.), collapse=", "))
  )[ (map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")) != map_chr(data_custom_dummy, ~paste(class(.), collapse=", ")) ), ] %>%
    as.data.frame() %>% rownames_to_column()
  
  if(nrow(unmatched_types)>0) stop(
    #unmatched_types
    "inconsistent typing in studydef : dummy dataset\n",
    apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  )
  
  data_extract <- data_custom_dummy 
  
} else {
  
  data_extract <- read_feather(here("output", "single", "extract", "input_timevarying.feather")) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date)) %>%
    # only keep patients in data_eligible
    right_join(data_eligible %>% select(patient_id), by = "patient_id")
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# data processing ----

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

# create one-row-per-event datasets for hospital admissions and discharges
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

# save datasets
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
