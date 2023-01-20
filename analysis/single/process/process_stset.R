# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# imports a cohort-specific processed dataset
# creates 3 datasets for that cohort:
# data_patients is one row per patient containing time-to-event data
# data_fixed is one row per patient containing baseline / time-invariant characteristics
# data_events is one row per patient event
# data_days is a (very large) one row per patient per day dataset
# creates additional survival variables for use in models (eg time to event from study start date)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# import libraries
library('tidyverse')
library('here')
library('glue')
library('survival')

# import local functions and parameters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "survival.R"))

# remove large objects when running on server
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  removeobs <- FALSE
} else {
  removeobs <- TRUE
}

# import processed data
data_eligible <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds"))
data_outcomes <- read_rds(here("output", "single", "process", "data_outcomes.rds"))

# create output directory
outdir <- here("output", "single", "stset")
fs::dir_create(outdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Generate data_fixed: time-invariant covariates ----

data_fixed <- data_eligible %>%
  select(
    patient_id,
    ageband2,
    region,
    all_of(unique(c(adjustment_variables_sequential, adjustment_variables_single)))
  )

# print dataset size
cat(" \n")
cat(glue("one-row-per-patient (time-independent) data size = ", nrow(data_fixed)), "\n")
cat(glue("memory usage = ", format(object.size(data_fixed), units="GB", standard="SI", digits=3L)), "\n")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Generate data_patients: one-row-per-patient data
data_patients <- data_eligible %>%
  left_join(data_outcomes, by = "patient_id") %>%
  transmute(
    patient_id,
    
    # since discrete dates are interpreted as the _end of_ the date, and we start follow up at the _start of_ the start date
    # the `minus 1` ensures that everybody has at least one discrete-time "day" that is event-free and vaccine-free
    # this is necessary for the propensity-type risk-of-vaccination model, essentially ensuring no events at time zero 
    start_date = study_dates$global$index_date - 1,
    end_date = study_dates$global$studyend_date,
    ageband2_start_date = case_when(
      ageband2 == "70-79" ~ study_dates$in70s$start_date - 1,
      ageband2 == "80+" ~ study_dates$over80s$start_date - 1,
    ),
    
    vax1_date,
    vax2_date,
    vax1_type,
    vax2_type,
    
    postest_date,
    covidadmitted_date,
    death_date,
    
    #composite of death, deregistration, end date, and no more than maxfup days after vaccination
    lastfup_date = pmin(vax1_date + maxfup, death_date, end_date, dereg_date, na.rm=TRUE),
    
    # Events are considered to occur the end of each day.
    # The study start date is the end of 7 december 2020 / start of 8 december 2020 = tstart=0.
    # The first possible vaccination date for individuals aged 80+ is 8 december 2020, 
    # ie between tstart=0 and tstop=1, so all patients are "unvaccinated" for at least 1 day of follow-up.
    # The first possible vaccination date for individuals aged 70-79 is 5 january 2021.
    # The first possible AZ vaccine date is 4 Jan 2021, ie between tstart=27 and tstop=28.
    # 
    # The first possible outcome date is 8 december 2020, ie between tstart=0 and tstop=1,
    # so all patients are event-free for at least 1 day of follow-up.
    # 
    # tte = "time to event"
    
    # time to end of study period
    tte_enddate = tte(start_date, end_date, end_date),
    
    # time to last follow up day
    tte_lastfup = tte(start_date, lastfup_date, lastfup_date),
    
    # time to deregistration
    tte_dereg = tte(start_date, dereg_date, dereg_date),
    
    # time to outcomes
    tte_postest = tte(start_date, postest_date, lastfup_date, na.censor=TRUE),
    tte_covidadmitted = tte(start_date, covidadmitted_date, lastfup_date, na.censor=TRUE),
    tte_death = tte(start_date, death_date, lastfup_date, na.censor=TRUE),

    # time to vaccination    
    tte_vaxany1 = tte(start_date, vax1_date, lastfup_date, na.censor=TRUE),
    tte_vaxany2 = tte(start_date, vax2_date, lastfup_date, na.censor=TRUE),
    
    tte_vaxpfizer1 = if_else(vax1_type == "pfizer", tte_vaxany1, NA_real_),
    tte_vaxpfizer2 = if_else(vax2_type == "pfizer", tte_vaxany2, NA_real_),
    
    tte_vaxaz1 = if_else(vax1_type == "az", tte_vaxany1, NA_real_),
    tte_vaxaz2 = if_else(vax2_type == "az", tte_vaxany2, NA_real_)
    
  ) %>%
  # convert tte variables to integer to save space, works since we know precision is to nearest day
  mutate(across(
    .cols = starts_with("tte_"),
    .fns = as.integer
  ))

stopifnot("vax1 time should not be same as vax2 time" = all(data_patients$tte_vaxany1 != data_patients$tte_vaxany2, na.rm=TRUE))

# print dataset size
cat(" \n")
cat(glue("one-row-per-patient (tte) data size = ", nrow(data_patients)), "\n")
cat(glue("memory usage = ", format(object.size(data_patients), units="MB", standard="SI", digits=3L)), "\n")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Generate data_events: one row per person per event dataset ----
# every time an event occurs or a covariate changes, a new row is generated

# import infectious hospitalisations data for time-updating "in-hospital" covariate
data_hospitalised_infectious <- 
  read_rds(here("output", "single", "process", "data_long_admission_infectious_dates.rds")) %>%
  pivot_longer(
    cols=c(admitted_date, discharged_date),
    names_to="status",
    values_to="date",
    values_drop_na = TRUE
  ) %>%
  inner_join(
    data_patients %>% select(patient_id, start_date, lastfup_date),
    .,
    by =c("patient_id")
  ) %>%
  mutate(
    tte = tte(start_date, date, lastfup_date, na.censor=TRUE) %>% as.integer(),
    hospinfectious_status = if_else(status=="admitted_date", 1L, 0L)
  )

# import non infectious hospitalisations data for time-updating "in-hospital" covariate
data_hospitalised_noninfectious <- 
  read_rds(here("output", "single", "process", "data_long_admission_noninfectious_dates.rds")) %>%
  pivot_longer(
    cols=c(admitted_date, discharged_date),
    names_to="status",
    values_to="date",
    values_drop_na = TRUE
  ) %>%
  inner_join(
    data_patients %>% select(patient_id, start_date, lastfup_date),
    .,
    by =c("patient_id")
  ) %>%
  mutate(
    tte = tte(start_date, date, lastfup_date, na.censor=TRUE) %>% as.integer(),
    hospnoninfectious_status = if_else(status=="admitted_date", 1L, 0L)
  )

# initial call based on events and vaccination status
data_events0 <- tmerge(
  
  data1 = data_patients %>% select(-starts_with("ind_"), -ends_with("_date")),
  data2 = data_patients,
  id = patient_id,
  
  vaxany_atrisk = tdc(ageband2_start_date-1-start_date),
  # brand specific = maximum of brand approval and age-based eligibility
  vaxpfizer_atrisk = tdc(pmax(study_dates$global$firstpfizer_date, ageband2_start_date)-1-start_date),
  vaxaz_atrisk = tdc(pmax(study_dates$global$firstaz_date, ageband2_start_date)-1-start_date),

  vaxany1_status = tdc(tte_vaxany1),
  vaxany2_status = tdc(tte_vaxany2),

  vaxpfizer1_status = tdc(tte_vaxpfizer1),
  vaxpfizer2_status = tdc(tte_vaxpfizer2),

  vaxaz1_status = tdc(tte_vaxaz1),
  vaxaz2_status = tdc(tte_vaxaz2),

  postest_status = tdc(tte_postest),
  covidadmitted_status = tdc(tte_covidadmitted),
  death_status = tdc(tte_death),
  dereg_status= tdc(tte_dereg),

  vaxany1 = event(tte_vaxany1),
  vaxany2 = event(tte_vaxany2),
  vaxpfizer1 = event(tte_vaxpfizer1),
  vaxpfizer2 = event(tte_vaxpfizer2),
  vaxaz1 = event(tte_vaxaz1),
  vaxaz2 = event(tte_vaxaz2),
  postest = event(tte_postest),
  covidadmitted = event(tte_covidadmitted),
  death = event(tte_death),
  dereg = event(tte_dereg),
  lastfup = event(tte_lastfup),
  
  tstart = 0L,
  tstop = tte_lastfup 
  
) 

# if not at risk of vax, not at risk of brand (use base R rather than dplyr to preserve tmerge class)
data_events0$vaxpfizer_atrisk <- if_else(data_events0$vaxany_atrisk == 0L, 0L, data_events0$vaxpfizer_atrisk)
data_events0$vaxaz_atrisk <- if_else(data_events0$vaxany_atrisk == 0L, 0L, data_events0$vaxaz_atrisk)

stopifnot("tstart should be  >= 0 in data_events0" = data_events0$tstart>=0)
stopifnot("tstop - tstart should be strictly > 0 in data_events0" = data_events0$tstop - data_events0$tstart > 0)

data_events <- data_events0 %>%
  tmerge(
    data1 = .,
    data2 = data_hospitalised_infectious,
    id = patient_id,
    hospinfectious_status = tdc(tte, hospinfectious_status),
    options = list(tdcstart = 0L)
  ) %>%
  tmerge(
    data1 = .,
    data2 = data_hospitalised_infectious %>% filter(status=="discharged_date"),
    id = patient_id,
    hospinfectiousdischarge = event(tte)
  ) %>%
  tmerge(
    data1 = .,
    data2 = data_hospitalised_noninfectious,
    id = patient_id,
    hospnoninfectious_status = tdc(tte, hospnoninfectious_status),
    options = list(tdcstart = 0L)
  ) %>%
  tmerge(
    data1 = .,
    data2 = data_hospitalised_noninfectious %>% filter(status=="discharged_date"),
    id = patient_id,
    hospnoninfectiousdischarge = event(tte)
  ) %>%
arrange(
  patient_id, tstart
  )

# do the following processing steps using base R rather than dplyr to preserve tmerge class
data_events$twidth <- data_events$tstop - data_events$tstart
data_events$vaxany_status <- data_events$vaxany1_status + data_events$vaxany2_status
data_events$vaxpfizer_status <- data_events$vaxpfizer1_status + data_events$vaxpfizer2_status
data_events$vaxaz_status <- data_events$vaxaz1_status + data_events$vaxaz2_status

cols_to_convert <-   c("vaxany1",
                       "vaxany2",
                       "vaxpfizer1",
                       "vaxpfizer2",
                       "vaxaz1",
                       "vaxaz2",
                       "postest",
                       "covidadmitted",
                       "death",
                       "hospinfectious_status",
                       "hospnoninfectious_status",
                       "hospinfectiousdischarge",
                       "hospnoninfectiousdischarge"
                       )

for (var in cols_to_convert) {
  data_events[[var]] <- as.integer(data_events[[var]])
}
 
# free up memory
if(removeobs){
  rm(data_events0, data_hospitalised_infectious, data_hospitalised_noninfectious)
}

stopifnot("tstart should be >= 0 in data_events" = data_events$tstart>=0)
stopifnot("tstop - tstart should be strictly > 0 in data_events" = data_events$tstop - data_events$tstart > 0)

# print dataset size
cat(" \n")
cat(glue("one-row-per-patient-per-event data size = ", nrow(data_events)), "\n")
cat(glue("memory usage = ", format(object.size(data_events), units="GB", standard="SI", digits=3L)), "\n")


# Save processed datasets ----
write_rds(data_fixed, file.path(outdir, "data_fixed.rds"), compress="gz")
write_rds(data_patients, file.path(outdir, "data_patients.rds"), compress="gz")
# don't compress data_events as this loses the tmerge class
write_rds(data_events, file.path(outdir, "data_events.rds"))
