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

n_matching_rounds <- 2


# define key dates ----

study_dates <- lst(
  
  pfizer = lst( # pfizer dose 3
   start_date = "2021-09-16", #start of recruitment thursday 16 september first pfizer booster jabs administered in england
   end_date = "2021-12-16", # end of recruitment (13 weeks later)
   # followupend_date = "2022-01-02", # end of follow-up
  ),
  
  moderna = lst( # moderna dose 3
    start_date = "2021-10-29", #start of recruitment friday 29 october first moderna booster jabs administered in england
    end_date = "2021-12-16", # end of recruitment (7 weeks later)
    # followupend_date = "2022-07-10" # end of follow-up
  ),
  
  studyend_date = "2021-12-31", # end of follow-up
  
  lastvax2_date = "2021-12-01", # don't recruit anyone with second vaccination after this date
  
  # dose 1 dates
  firstpfizer_date = "2020-12-08", # first pfizer vaccination in national roll-out
  firstaz_date = "2021-01-04", # first az vaccination in national roll-out
  firstmoderna_date = "2021-04-13", # first moderna vaccination in national roll-out
  firstpossiblevax_date = "2020-06-01", # used to catch "real" vaccination dates (eg not 1900-01-01)
)

study_dates$index_date = study_dates$pfizer$start_date

extract_increment <- 14

study_dates$pfizer$control_extract_dates = as.Date(study_dates$pfizer$start_date) + (0:26)*extract_increment
study_dates$moderna$control_extract_dates = as.Date(study_dates$moderna$start_date) + (0:26)*extract_increment

jsonlite::write_json(study_dates, path = here("lib", "design", "study-dates.json"), auto_unbox=TRUE, pretty =TRUE)

# all as dates
lens <- sapply(study_dates, length)
dates_general <- map(study_dates[lens==1], as.Date)
dates_cohort <- map(study_dates[lens==3], ~map(.x, as.Date))
study_dates <- splice(dates_general, dates_cohort)[names(study_dates)]

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
  "covidcc", "covidcc_date", "COVID-19 critical care",
  "coviddeath", "coviddeath_date", "COVID-19 death",
  "noncoviddeath", "noncoviddeath_date", "Non-COVID-19 death",
  "death", "death_date", "Any death",

  # safety
  "admitted", "admitted_unplanned_1_date", "Unplanned hospitalisation",
  "emergency", "emergency_date", "A&E attendance",
)

# define treatments ----

treatement_lookup <-
  tribble(
    ~treatment, ~treatment_descr,
    "pfizer", "BNT162b2",
    "az", "ChAdOx1-S",
    "moderna", "mRNA-1273",
    "pfizer-pfizer", "BNT162b2",
    "az-az", "ChAdOx1-S",
    "moderna-moderna", "mRNA-1273"
  )

## lookups to convert coded variables to full, descriptive variables ----

recoder <-
  lst(
    subgroups = c(
      `Main` = "all",
      `Prior SARS-CoV-2 infection` = "prior_covid_infection"
    ),
    status = c(
      `Unmatched`= "unmatched",
      `Matched` = "matched"
    ),
    treated = c(
      `Two doses` = "0",
      `Three doses` = "1"
    ),
    outcome = set_names(events_lookup$event, events_lookup$event_descr),
    all = c(` ` = "all"),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
  )


## follow-up time ----

# period width
postbaselinedays <- 28

# where to split follow-up time after recruitment
postbaselinecuts <- c(14, 14 + (1:6)*postbaselinedays)

# maximum follow-up
maxfup <- max(postbaselinecuts)

# matching variables ----

# exact variables
exact_variables <- c(
  
  "jcvi_ageband",
  "cev_cv",
  "vax12_type",
  #"vax2_week",
  "region",
  #"sex",
  #"cev_cv",
  
  #"multimorb",
  "prior_covid_infection",
  #"immunosuppressed",
  #"status_hospplanned"
  NULL
)

# caliper variables
caliper_variables <- c(
  age = 3,
  vax2_day = 7,
  NULL
)
matching_variables <- c(exact_variables, names(caliper_variables))

# cut-off for rolling 7 day average, that determines recruitment period
recruitment_period_cutoff <- 50
