# # # # # # # # # # # # # # # # # # # # #
# This script:
# creates metadata for aspects of the study design
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')

## create output directories ----
fs::dir_create(here("lib", "design"))



# number of matching rounds to perform

n_matching_rounds <- 2#4


# define key dates ----

study_dates <- lst(
  
  pfizer = lst( # pfizer dose 1
   start_date = "2020-12-08", #start of recruitment
   end_date = "2021-04-19", # end of recruitment
   # followupend_date = "2022-01-02", # end of follow-up
  ),
  
  az = lst( # az dose 1
    start_date = "2021-01-04", #start of recruitment
    end_date = "2021-04-19", # end of recruitment
    # followupend_date = "2022-07-10" # end of follow-up
  ),
  
  over80s = lst(
    start_date = "2020-12-08"
  ),
  
  in70s = lst(
    start_date = "2021-01-05"
  ),
  
  global = lst(
    
    index_date = "2020-12-08",
  
    studyend_date = "2021-04-19", # end of follow-up
  
    # dose 1 dates
    firstpfizer_date = "2020-12-08", # first pfizer vaccination in national roll-out
    firstaz_date = "2021-01-04", # first az vaccination in national roll-out
    firstmoderna_date = "2021-04-13", # first moderna vaccination in national roll-out
    firstpossiblevax_date = "2020-06-01", # used to catch "real" vaccination dates (eg not 1900-01-01)
  ),
)



extract_increment <- 14

study_dates$pfizer$control_extract_dates = as.Date(study_dates$pfizer$start_date) + (0:26)*extract_increment
study_dates$az$control_extract_dates = as.Date(study_dates$az$start_date) + (0:26)*extract_increment

jsonlite::write_json(study_dates, path = here("lib", "design", "study-dates.json"), auto_unbox=TRUE, pretty =TRUE)

# all as dates
study_dates <- map_depth(study_dates, 2, as.Date, .ragged=TRUE)


# define outcomes ----

events_lookup <- tribble(
  ~event, ~event_var, ~event_descr,

  # other
  "test", "covid_test_date", "SARS-CoV-2 test",

  # effectiveness
  "postest", "positive_test_date", "Positive SARS-CoV-2 test",
  "covidemergency", "covidemergency_date", "COVID-19 A&E attendance",
  "covidadmitted", "covidadmitted_date", "COVID-19 hospitalisation",
  "noncovidadmitted", "noncovidadmitted_date", "Non-COVID-19 hospitalisation",
  "covidadmittedproxy1", "covidadmittedproxy1_date", "COVID-19 hospitalisation (A&E proxy)",
  "covidadmittedproxy2", "covidadmittedproxy2_date", "COVID-19 hospitalisation (A&E proxy v2)",
  "covidcritcare", "covidcc_date", "COVID-19 critical care",
  "coviddeath", "coviddeath_date", "COVID-19 death",
  "noncoviddeath", "noncoviddeath_date", "Non-COVID-19 death",
  "death", "death_date", "Any death",

  # safety
  "admitted", "admitted_unplanned_1_date", "Unplanned hospitalisation",
  "emergency", "emergency_date", "A&E attendance",
)

model_outcomes <- c("postest", "covidadmitted", "death")

# define treatments ----

treatement_lookup <-
  tribble(
    ~treatment, ~treatment_descr,
    "pfizer", "BNT162b2",
    "az", "ChAdOx1-S",
  )

## lookups to convert coded variables to full, descriptive variables ----

recoder <-
  lst(
    subgroups = c(
      `Main` = "all",
      `Age band` = "ageband2"
    ),
    status = c(
      `Unmatched`= "unmatched",
      `Matched` = "matched"
    ),
    treated = c(
      `Unvaccinated` = "0",
      `Vaccinated` = "1"
    ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    ageband2 = c(
      `aged 80+` = "80+",
      `aged 70-79` = "70-79"
    ),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
  )


## follow-up time ----

# where to split follow-up time after recruitment
#postbaselinecuts <- c(14, 14 + ((1:6)*28))
postbaselinecuts <- c(3,7,14,21,28,35,70)


# maximum follow-up
maxfup <- max(postbaselinecuts)

# matching variables ----

# exact variables
exact_variables <- c(
  "jcvi_ageband",
  "cev_cv",
  "region",
  "sex",
  "timesince_covid_cat",
  "prior_covid_infection",
  NULL
)

# caliper variables
caliper_variables <- c(
  age = 3,
  NULL
)
matching_variables <- c(exact_variables, names(caliper_variables))

## important NOTE:

# these adjustment variables are different to those use for adjustment in the original MSM approach
# these are updated variables using the PRIMIS specification for identifying high risk individuals
# the older variable definitions used in the original MSM study have been superceded, although there is a great deal of overlap between the two
# some omissions, that would require re-extraction to include, are:
#   efi (electronic frailty index); flu vaccination; recent unplanned hospital admission (for any cause).
# the "shielded" indicator is accommodated by the matching,  as "cev_cv" captures all  shielded people.
# calendar-time does not need to be adjusted for as follow up is parallel on calendar time in each group. 

adjustment_variables <- c(
  "imd_Q5",
  "ethnicity_combined",
  "sev_obesity",
  "chronic_heart_disease",
  "chronic_kidney_disease",
  "diabetes",
  "chronic_liver_disease",
  "chronic_resp_disease", "asthma",
  "chronic_neuro_disease",
  "learndis",
  "sev_mental",
  "immunosuppressed", "asplenia",
  "multimorb",
  NULL
)
