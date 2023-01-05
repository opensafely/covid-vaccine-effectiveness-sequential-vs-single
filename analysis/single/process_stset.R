
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports a cohort-specific processed dataset
# creates 3 datasets for that cohort:
# data_tte is one row per patient containing time-to-event data
# data_fixed is one row per patient containing baseline / time-invariant characteristics
# data_pt is a (very large) one row per patient per day "person-time" dataset
# creates additional survival variables for use in models (eg time to event from study start date)
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "survival.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  
} else{
  removeobs <- TRUE
}

# Import processed data ----
data_eligible <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds"))
data_outcomes <- read_rds(here("output", "single", "process", "data_outcomes.rds"))

# create output directory ----
outdir <- here("output", "single", "stset")
fs::dir_create(outdir)

# Generate different data formats ----

## one-row-per-patient data ----

data_fixed <- data_eligible %>%
  select(
    patient_id,
    region,
    all_of(adjustment_variables)
  )

## print dataset size ----
cat(" \n")
cat(glue("one-row-per-patient (time-independent) data size = ", nrow(data_fixed)), "\n")
cat(glue("memory usage = ", format(object.size(data_fixed), units="GB", standard="SI", digits=3L)), "\n")

data_tte <- data_eligible %>%
  left_join(data_outcomes, by = "patient_id") %>%
  transmute(
    patient_id,
    
    # since discrete dates are interpreted as the _end of_ the date, and we start follow up at the _start of_ the start date
    # this ensures everybody has at least one discrete-time "day" event-free and vaccine-free
    start_date = study_dates$global$index_date - 1,
    end_date = study_dates$global$studyend_date,
    
    vax1_date,
    vax2_date,
    vax1_type,
    vax2_type,
    
    postest_date,
    covidadmitted_date,
    coviddeath_date,
    death_date,
    
    #composite of death, deregistration and end date
    lastfup_date = pmin(death_date, end_date, dereg_date, na.rm=TRUE),
    
    tte_enddate = tte(start_date, end_date, end_date),
    
    # events are considered to occur at midnight at the end of each day.
    # The study start date is the end of 7 december 2020 / start of 8 december 2020 = tstart=0
    # The first possible vaccination date is 8 december 2020, ie between tstart=0 and tstop=1, so all patients are "unvaccinated" for at least 1 day of follow-up
    # the first possible AZ vaccine date is 4 Jan 2021, ie between tstart=27 and tstop=28
    # the first possible outcome date is 8 december 2020, ie between tstart=0 and tstop=1, so all patients are event-free for at least 1 day of follow-up
    
    # time to last follow up day
    tte_lastfup = tte(start_date, lastfup_date, lastfup_date),
    
    # time to deregistration
    tte_dereg = tte(start_date, dereg_date, dereg_date),
    
    # time to outcomes
    tte_postest = tte(start_date, postest_date, lastfup_date, na.censor=TRUE),
    tte_covidadmitted = tte(start_date, covidadmitted_date, lastfup_date, na.censor=TRUE),
    tte_coviddeath = tte(start_date, coviddeath_date, lastfup_date, na.censor=TRUE),
    tte_death = tte(start_date, death_date, lastfup_date, na.censor=TRUE),

    # time to vaccination    
    tte_vaxany1 = tte(start_date, vax1_date, lastfup_date, na.censor=TRUE),
    tte_vaxany2 = tte(start_date, vax2_date, lastfup_date, na.censor=TRUE),
    
    tte_vaxpfizer1 = if_else(vax1_type == "pfizer", tte_vaxany1, NA_real_),
    tte_vaxpfizer2 = if_else(vax2_type == "pfizer", tte_vaxany2, NA_real_),
    
    tte_vaxaz1 = if_else(vax1_type == "az", tte_vaxany1, NA_real_),
    tte_vaxaz2 = if_else(vax2_type == "az", tte_vaxany2, NA_real_)
    
  ) %>%
  # convert tte variables to integer to save space. works since we know precision is to nearest day
  mutate(across(
    .cols = starts_with("tte_"),
    .fns = as.integer
  ))

stopifnot("vax1 time should not be same as vax2 time" = all(data_tte$tte_vaxany1 != data_tte$tte_vaxany2, na.rm=TRUE))

## print dataset size ----
cat(" \n")
cat(glue("one-row-per-patient (tte) data size = ", nrow(data_tte)), "\n")
cat(glue("memory usage = ", format(object.size(data_tte), units="MB", standard="SI", digits=3L)), "\n")

## create counting-process format dataset ----
# ie, one row per person per event
# every time an event occurs or a covariate changes, a new row is generated

# import infectious hospitalisations data for time-updating "in-hospital" covariate
data_hospitalised_infectious <- read_rds(here("output", "single", "process", "data_long_admission_infectious_dates.rds")) %>%
  pivot_longer(
    cols=c(admitted_date, discharged_date),
    names_to="status",
    values_to="date",
    values_drop_na = TRUE
  ) %>%
  inner_join(
    data_tte %>% select(patient_id, start_date, lastfup_date),
    .,
    by =c("patient_id")
  ) %>%
  mutate(
    tte = tte(start_date, date, lastfup_date, na.censor=TRUE) %>% as.integer(),
    hospinfectious_status = if_else(status=="admitted_date", 1L, 0L)
  )

# import non infectious hospitalisations data for time-updating "in-hospital" covariate
data_hospitalised_noninfectious <- read_rds(here("output", "single", "process", "data_long_admission_noninfectious_dates.rds")) %>%
  pivot_longer(
    cols=c(admitted_date, discharged_date),
    names_to="status",
    values_to="date",
    values_drop_na = TRUE
  ) %>%
  inner_join(
    data_tte %>% select(patient_id, start_date, lastfup_date),
    .,
    by =c("patient_id")
  ) %>%
  mutate(
    tte = tte(start_date, date, lastfup_date, na.censor=TRUE) %>% as.integer(),
    hospnoninfectious_status = if_else(status=="admitted_date", 1L, 0L)
  )


# data_suspected_covid <- read_rds(here("output", cohort, "data", "data_long_pr_suspected_covid_dates.rds")) %>%
#   inner_join(
#     data_tte %>% select(patient_id, start_date, lastfup_date),
#     .,
#     by =c("patient_id")
#   ) %>%
#   mutate(
#     tte = tte(start_date, date, lastfup_date, na.censor=TRUE)
#   )

# data_probable_covid <- read_rds(here("output", cohort, "data", "data_long_pr_probable_covid_dates.rds")) %>%
#   inner_join(
#     data_tte %>% select(patient_id, start_date, lastfup_date),
#     .,
#     by =c("patient_id")
#   ) %>%
#   mutate(
#     tte = tte(start_date, date, lastfup_date, na.censor=TRUE),
#   )

# data_postest <- read_rds(here("output", "single", "process", "data_long_postest_dates.rds")) %>%
#   inner_join(
#     data_tte %>% select(patient_id, start_date, lastfup_date),
#     .,
#     by =c("patient_id")
#   ) %>%
#   mutate(
#     tte = tte(start_date, date, lastfup_date, na.censor=TRUE),
#   )

# initial call based on events and vaccination status
data_tte_cp0 <- tmerge(
  
  data1 = data_tte %>% select(-starts_with("ind_"), -ends_with("_date")),
  data2 = data_tte,
  id = patient_id,
  
  vaxany_atrisk = tdc(start_date-start_date),
  vaxpfizer_atrisk = tdc(as.Date(study_dates$global$firstpfizer_date)-1-start_date),
  vaxaz_atrisk = tdc(as.Date(study_dates$global$firstaz_date)-1-start_date),

  vaxany1_status = tdc(tte_vaxany1),
  vaxany2_status = tdc(tte_vaxany2),

  vaxpfizer1_status = tdc(tte_vaxpfizer1),
  vaxpfizer2_status = tdc(tte_vaxpfizer2),

  vaxaz1_status = tdc(tte_vaxaz1),
  vaxaz2_status = tdc(tte_vaxaz2),

  # covidtest_status = tdc(tte_covidtest),
  postest_status = tdc(tte_postest),
  # emergency_status = tdc(tte_emergency),
  covidadmitted_status = tdc(tte_covidadmitted),
  # coviddeath_status = tdc(tte_coviddeath),
  # noncoviddeath_status = tdc(tte_noncoviddeath),
  death_status = tdc(tte_death),
  dereg_status= tdc(tte_dereg),
  lastfup_status = tdc(tte_lastfup),

  vaxany1 = event(tte_vaxany1),
  vaxany2 = event(tte_vaxany2),
  vaxpfizer1 = event(tte_vaxpfizer1),
  vaxpfizer2 = event(tte_vaxpfizer2),
  vaxaz1 = event(tte_vaxaz1),
  vaxaz2 = event(tte_vaxaz2),
  # covidtest = event(tte_covidtest),
  postest = event(tte_postest),
  # emergency = event(tte_emergency),
  covidadmitted = event(tte_covidadmitted),
  # coviddeath = event(tte_coviddeath),
  # noncoviddeath = event(tte_noncoviddeath),
  death = event(tte_death),
  dereg = event(tte_dereg),
  lastfup = event(tte_lastfup),
  
  tstart = 0L,
  tstop = tte_enddate # use enddate not lastfup because it's useful for status over time plots
  
) 

stopifnot("tstart should be  >= 0 in data_tte_cp0" = data_tte_cp0$tstart>=0)
stopifnot("tstop - tstart should be strictly > 0 in data_tte_cp0" = data_tte_cp0$tstop - data_tte_cp0$tstart > 0)

data_tte_cp <- data_tte_cp0 %>%
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
  # tmerge(
  #   data1 = .,
  #   data2 = data_suspected_covid,
  #   id = patient_id,
  #   suspectedcovid = event(tte)
  # ) %>%
  # tmerge(
  #   data1 = .,
  #   data2 = data_probable_covid,
  #   id = patient_id,
  #   probablecovid = event(tte)
  # ) %>%
  # tmerge(
  #   data1 = .,
  #   data2 = data_postest,
  #   id = patient_id,
  #   postesttdc = event(tte)
  # ) %>%
  arrange(
    patient_id, tstart
  ) %>%
  mutate(
    twidth = tstop - tstart,
    vaxany_status = vaxany1_status + vaxany2_status,
    vaxpfizer_status = vaxpfizer1_status + vaxpfizer2_status,
    vaxaz_status = vaxaz1_status + vaxaz2_status,
  ) %>%
  ungroup() %>%
  # for some reason tmerge converts event indicators to numeric. So convert back to save space
  mutate(across(
    .cols = c("vaxany1",
              "vaxany2",
              "vaxpfizer1",
              "vaxpfizer2",
              "vaxaz1",
              "vaxaz2",
              # "covidtest",
              "postest",
              # "emergency",
              "covidadmitted",
              # "coviddeath",
              # "noncoviddeath",
              "death",
              "lastfup",
              "hospinfectious_status",
              "hospnoninfectious_status",
              "hospinfectiousdischarge",
              "hospnoninfectiousdischarge"#,
              # "suspectedcovid",
              #"probablecovid",
              # "postesttdc"
    ),
    .fns = as.integer
  ))

# free up memory
if(removeobs){
  rm(
    data_tte_cp0, data_hospitalised_infectious, data_hospitalised_noninfectious, 
    # data_suspected_covid, data_probable_covid,
    data_postest
    )
}

stopifnot("tstart should be >= 0 in data_tte_cp" = data_tte_cp$tstart>=0)
stopifnot("tstop - tstart should be strictly > 0 in data_tte_cp" = data_tte_cp$tstop - data_tte_cp$tstart > 0)

### print dataset size ----
cat(" \n")
cat(glue("one-row-per-patient-per-event data size = ", nrow(data_tte_cp)), "\n")
cat(glue("memory usage = ", format(object.size(data_tte_cp), units="GB", standard="SI", digits=3L)), "\n")

## create person-time format dataset ----
# ie, one row per person per day (or per week or per month)
# this format has lots of redundancy but is necessary for MSMs
alltimes <- expand(
  data_tte, 
  patient_id, 
  times=as.integer(full_seq(c(1, tte_enddate),1))
  )

# do not use survSplit as this doesn't handle multiple events properly
# eg, a positive test will be expanded as if a tdc (eg c(0,0,1,1,1,..)) not an event (eg c(0,0,1,0,0,...))
# also, survSplit is slower!
data_tte_pt <- tmerge(
  data1 = data_tte_cp,
  data2 = alltimes,
  id = patient_id,
  alltimes = event(times, times)
) %>%
  arrange(patient_id, tstop) %>%
  group_by(patient_id) %>%
  mutate(
    hospinfectiousdischarge_time = if_else(hospinfectiousdischarge==1, tstop, NA_real_),
    hospnoninfectiousdischarge_time = if_else(hospnoninfectiousdischarge==1, tstop, NA_real_),
    # suspectedcovid_time = if_else(suspectedcovid==1, tstop, NA_real_),
    # probablecovid_time = if_else(probablecovid==1, tstop, NA_real_),
    postesttdc_time = if_else(postesttdc==1, tstop, NA_real_),
  ) %>%
  fill(
    hospinfectiousdischarge_time, hospnoninfectiousdischarge_time, 
    # suspectedcovid_time,
    # probablecovid_time,
    postesttdc_time
  ) %>%
  mutate(
    
    # define time since vaccination
    vaxany1_timesince = cumsum(vaxany1_status),
    vaxany2_timesince = cumsum(vaxany2_status),
    vaxpfizer1_timesince = cumsum(vaxpfizer1_status),
    vaxpfizer2_timesince = cumsum(vaxpfizer2_status),
    vaxaz1_timesince = cumsum(vaxaz1_status),
    vaxaz2_timesince = cumsum(vaxaz2_status),
    
    # postest_timesince = tstop - postesttdc_time,
    
    # define time since infectious hospitalisation
    timesince_hospinfectiousdischarge_pw = cut(
      tstop - hospinfectiousdischarge_time,
      breaks=c(0, 21, 28),
      labels=c( "1-21", "22-28"),
      right=TRUE
    ),
    timesince_hospinfectiousdischarge_pw = case_when(
      is.na(timesince_hospinfectiousdischarge_pw) & hospinfectious_status==0 ~ "Not in hospital",
      hospinfectious_status==1 ~ "In hospital",
      !is.na(timesince_hospinfectiousdischarge_pw) ~ as.character(timesince_hospinfectiousdischarge_pw),
      TRUE ~ NA_character_
    ) %>% factor(c("Not in hospital", "In hospital", "1-21", "22-28")),
    
    
    # define time since non infectious hospitalisation
    timesince_hospnoninfectiousdischarge_pw = cut(
      tstop - hospnoninfectiousdischarge_time,
      breaks=c(0, 21, 28),
      labels=c( "1-21", "22-28"),
      right=TRUE
    ),
    timesince_hospnoninfectiousdischarge_pw = case_when(
      is.na(timesince_hospnoninfectiousdischarge_pw) & hospnoninfectious_status==0 ~ "Not in hospital",
      hospnoninfectious_status==1 ~ "In hospital",
      !is.na(timesince_hospnoninfectiousdischarge_pw) ~ as.character(timesince_hospnoninfectiousdischarge_pw),
      TRUE ~ NA_character_
    ) %>% factor(c("Not in hospital", "In hospital", "1-21", "22-28")),
    
    # define time since covid primary care event
    # timesince_suspectedcovid_pw = cut(
    #   tstop - suspectedcovid_time,
    #   breaks=c(0, 21, 28, Inf),
    #   labels=c("1-21", "22-28", "29+"),
    #   right=TRUE
    # ) %>% fct_explicit_na(na_level="Not suspected") %>% factor(c("Not suspected", "1-21", "22-28", "29+")),
    
    # timesince_probablecovid_pw = cut(
    #   tstop - probablecovid_time,
    #   breaks=c(0, 21, 28, Inf),
    #   labels=c("1-21", "22-28", "29+"),
    #   right=TRUE
    # ) %>% fct_explicit_na(na_level="Not probable")  %>% factor(c("Not probable", "1-21", "22-28", "29+")),
    
    # define time since positive SGSS test
    timesince_postesttdc_pw = cut(
      tstop - postesttdc_time,
      breaks=c(0, 21, 28, Inf),
      labels=c("1-21", "22-28", "29+"),
      right=TRUE
    ) %>% fct_explicit_na(na_level="No positive test")  %>% factor(c("No positive test", "1-21", "22-28", "29+")),
    
  ) %>%
  ungroup() %>%
  select(
    -hospinfectiousdischarge_time,
    -hospnoninfectiousdischarge_time#,
    # -suspectedcovid_time,
    # -probablecovid_time,
    # -postesttdc_time,
  ) %>%
  # for some reason tmerge converts event indicators to numeric. So convert back to save space
  mutate(across(
    .cols = c("vaxany1",
              "vaxany2",
              "vaxpfizer1",
              "vaxpfizer2",
              "vaxaz1",
              "vaxaz2",
              # "covidtest",
              "postest",
              # "emergency",
              "covidadmitted",
              # "coviddeath",
              # "noncoviddeath",
              "death",
              "dereg",
              "lastfup",
              "hospinfectious_status",
              "hospnoninfectious_status",
              "hospinfectiousdischarge",
              "hospnoninfectiousdischarge"#,
              # "probablecovid",
              # "suspectedcovid",
              # "postesttdc",
              # "postest_timesince"
    ),
    .fns = as.integer
  ))

stopifnot("dummy 'alltimes' should be equal to tstop" = all(data_tte_pt$alltimes == data_tte_pt$tstop))


# remove unused columns
data_tte_pt <- data_tte_pt %>%
  mutate(
    vaxanyday1 = tte_vaxany1
  ) %>%
  select(
    -starts_with("tte_"),
    -ends_with("_date")
  )

### print dataset size ----
cat(" \n")
cat(glue("one-row-per-patient-per-time-unit data size = ", nrow(data_tte_pt)), "\n")
cat(glue("memory usage = ", format(object.size(data_tte_pt), units="GB", standard="SI", digits=3L)), "\n")

## Save processed tte data ----
write_rds(data_fixed, file.path(outdir, "data_fixed.rds"), compress="gz")
write_rds(data_tte, file.path(outdir, "data_tte.rds"), compress="gz")
write_rds(data_tte_pt, file.path(outdir, "data_pt.rds"), compress="gz")
