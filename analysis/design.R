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

## round counts to nearest 6 for disclosure control
threshold <- 6

# number of matching rounds to perform

n_matching_rounds <- 2#4


# define key dates ----

study_dates <- lst(
  
  global = lst(
    
    # dose 1 dates
    firstpfizer_date = "2020-12-08", # first pfizer vaccination in national roll-out
    firstaz_date = "2021-01-04", # first az vaccination in national roll-out
    firstmoderna_date = "2021-04-13", # first moderna vaccination in national roll-out
    firstpossiblevax_date = "2020-06-01", # used to catch "real" vaccination dates (eg not 1900-01-01)
    
    index_date = firstpfizer_date,
    studyend_date = as.Date(firstmoderna_date) - 1, # end of follow-up
    
  ),
  
  pfizer = lst( # pfizer dose 1
   start_date = global$firstpfizer_date, #start of recruitment
   end_date = global$studyend_date, # end of recruitment and follow-up
  ),
  
  az = lst( # az dose 1
    start_date = global$firstaz_date, #start of recruitment
    end_date = global$studyend_date, # end of recruitment and follow-up
  ),
  
  over80s = lst(
    start_date = global$firstpfizer_date
  ),

  in70s = lst(
    start_date = "2021-01-05"
  ),
  
)



extract_increment <- 14

for (brand in c("pfizer", "az")) {
  study_dates[[brand]]$control_extract_dates = as.Date(study_dates[[brand]]$start_date) + (1:n_matching_rounds - 1)*extract_increment  
}

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
      `over80` = "80+",
      `in70s` = "70-79"
    ),
    prior_covid_infection = c(
      `No prior SARS-CoV-2 infection` = "FALSE",
      `Prior SARS-CoV-2 infection` = "TRUE"
    ),
  )

model_subgroups <- "all"

## follow-up time ----

# where to split follow-up time after recruitment
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

adjustment_variables_sequential <- c(
  "imd_Q5",
  "ethnicity_combined",
  "sev_obesity",
  "chronic_heart_disease",
  "chronic_kidney_disease",
  "diabetes",
  "chronic_liver_disease",
  "chronic_resp_disease", 
  "chronic_neuro_disease",
  "multimorb",
  "immunosuppressed", 
  "learndis",
  "sev_mental",
  "flu_vaccine",
  NULL
)

adjustment_variables_single <- c(
  "age",
  "sex",
  "region",
  "cev_cv",
  adjustment_variables_sequential
)

# define formulas for single trial
list_formula_single <- local({
  
  formula_exposure <- . ~ . + timesincevax_pw
  adjustment_variables_single_linear <- str_c(adjustment_variables_single[!(adjustment_variables_single%in%c("age", "region"))], collapse = " + ")
  formula_covars <- as.formula(glue(". ~ . + poly(age, degree=2, raw=TRUE) + {adjustment_variables_single_linear}"))
  # formula_secular <- . ~ . + ns(tstop, df=5)
  formula_secular_region <- . ~ . + ns(tstop, df=5)*region
  formula_timedependent <- . ~ . + timesince_hospinfectiousdischarge_pw + timesince_hospnoninfectiousdischarge_pw
  
  formula_all_rhsvars <- update(1 ~ 1, formula_exposure) %>%
    update(formula_covars) %>%
    # update(formula_secular) %>%
    update(formula_secular_region) %>%
    update(formula_timedependent)
  
  list_formula <- lst(
    formula_exposure,
    formula_covars,
    # formula_secular,
    formula_secular_region,
    formula_timedependent,
    formula_all_rhsvars,
  )
  
  return(list_formula)
  
})
