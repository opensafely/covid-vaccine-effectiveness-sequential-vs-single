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

# # import command-line arguments
# args <- commandArgs(trailingOnly=TRUE)
# if(length(args)==0){
#   removeobs <- FALSE
#   iteration <- 1L
# } else {
#   removeobs <- TRUE
#   iteration <- as.integer(args[[1]])
# }

# create output directory
outdir <- here("output", "single", "stset")
# fs::dir_create(outdir)


for (iteration in 1:process_data_days_n) {
  
  # import datasets from process_stset and restrict to iteration_ids  ----
  
  data_fixed <- read_rds(file.path(outdir, "data_fixed.rds")) 
  
  # run this script over `process_data_days_n` iterations (defined in analysis/design.R)
  # otherwise the script runs into memory issues
  iteration_size <- ceiling(nrow(data_fixed)/process_data_days_n)
  
  iteration_ids <- data_fixed %>%
    transmute(
      patient_id = as.integer(patient_id),
      rank_id = dense_rank(patient_id)
    ) %>%
    filter(
      (iteration - 1)*iteration_size < rank_id,
      rank_id <= iteration*iteration_size
    ) %>%
    distinct(patient_id) %>%
    pull(patient_id)
  
  data_fixed <- data_fixed %>%
    filter(as.integer(patient_id) %in% iteration_ids)
  
  data_patients <- read_rds(file.path(outdir, "data_patients.rds")) %>%
    filter(as.integer(patient_id) %in% iteration_ids)
  
  data_events <- read_rds(file.path(outdir, "data_events.rds")) %>%
    filter(as.integer(patient_id) %in% iteration_ids)
  
  rm(iteration_ids)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Generate data_days: one one row per person per day ----
  # this format has lots of redundancy but is necessary for MSMs
  alltimes <- expand(
    data_patients, 
    patient_id, 
    times=as.integer(full_seq(c(1, tte_enddate),1))
  )
  
  # do not use survSplit as this doesn't handle multiple events properly
  # eg, a positive test will be expanded as if a tdc (eg c(0,0,1,1,1,..)) not an event (eg c(0,0,1,0,0,...))
  # also, survSplit is slower!
  data_days <- tmerge(
    data1 = data_events,
    data2 = alltimes,
    id = patient_id,
    alltimes = event(times, times)
  ) %>%
    arrange(patient_id, tstop) %>%
    group_by(patient_id) %>%
    mutate(
      hospinfectiousdischarge_time = if_else(hospinfectiousdischarge==1, tstop, NA_real_),
      hospnoninfectiousdischarge_time = if_else(hospnoninfectiousdischarge==1, tstop, NA_real_)
    ) %>%
    fill(
      hospinfectiousdischarge_time, hospnoninfectiousdischarge_time 
    ) %>%
    mutate(
      
      # define time since vaccination
      vaxany1_timesince = cumsum(vaxany1_status),
      vaxany2_timesince = cumsum(vaxany2_status),
      vaxpfizer1_timesince = cumsum(vaxpfizer1_status),
      vaxpfizer2_timesince = cumsum(vaxpfizer2_status),
      vaxaz1_timesince = cumsum(vaxaz1_status),
      vaxaz2_timesince = cumsum(vaxaz2_status),
      
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
      ) %>% factor(c("Not in hospital", "In hospital", "1-21", "22-28"))
      
    ) %>%
    ungroup() %>%
    select(
      -hospinfectiousdischarge_time,
      -hospnoninfectiousdischarge_time
    ) %>%
    # tmerge converts event indicators to numeric - convert back to save space
    mutate(across(
      .cols = c("vaxany1",
                "vaxany2",
                "vaxpfizer1",
                "vaxpfizer2",
                "vaxaz1",
                "vaxaz2",
                "postest",
                "covidadmitted",
                "death",
                "dereg",
                "lastfup",
                "hospinfectious_status",
                "hospnoninfectious_status",
                "hospinfectiousdischarge",
                "hospnoninfectiousdischarge"
      ),
      .fns = as.integer
    ))
  
  stopifnot("dummy 'alltimes' should be equal to tstop" = all(data_days$alltimes == data_days$tstop))
  
  
  # remove unused columns
  data_days <- data_days %>%
    select(
      -starts_with("tte_"),
      -ends_with("_date")
    )
  
  # print dataset size
  cat(" \n")
  cat(glue("one-row-per-patient-per-time-unit data size = ", nrow(data_days)), "\n")
  cat(glue("memory usage = ", format(object.size(data_days), units="GB", standard="SI", digits=3L)), "\n")
  
  # Save processed dataset ----
  write_rds(data_days, file.path(outdir, glue("data_days_{iteration}.rds")), compress="gz")
  
  rm(data_fixed, data_events, data_patients, data_days, alltimes)
  
}
