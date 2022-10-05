
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports processed data
# chooses matching sets for each sequential trial
# outputs matching summary
#
# The script must be accompanied by two arguments:
# `agegroup` - over12s or under12s
# `matching_round` - the matching round (1,2,3,...)

# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('MatchIt')

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "pfizer"
  matching_round <- as.integer("1")
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
  matching_round <- as.integer(args[[2]])
}


## get cohort-specific parameters study dates and parameters ----

dates <- map(study_dates[[cohort]], as.Date)


matching_round_date <- dates$control_extract_dates[matching_round]



## create output directory ----
fs::dir_create(ghere("output", cohort, "matchround{matching_round}", "potential"))

# Import datasets ----

## import treated populations ----
data_alltreated <- read_rds(ghere("output", cohort, "treated", "data_treatedeligible.rds")) %>% mutate(treated=1L)

## import control populations ----
data_control <- read_rds(ghere("output", cohort, "matchround{matching_round}", "process", "data_controlpotential.rds")) %>% mutate(treated=0L)

# remove already-matched people from previous matching rounds
if(matching_round>1){
  
  data_matchstatusprevious <- read_rds(ghere("output", cohort, "matchround{matching_round-1L}", "actual", "data_matchstatus_allrounds.rds")) %>%
    select(patient_id, treated)
  
  # do not select treated people who have already been matched
  data_alltreated <- 
    data_alltreated %>%
    anti_join(
      data_matchstatusprevious, 
      by=c("patient_id", "treated")
    )
  
  # do not select untreated people who have already been matched
  data_control <- 
    data_control %>%
    anti_join(
      data_matchstatusprevious, 
      by=c("patient_id", "treated")
    )
}


## import matching variables ----

data_eligible <-
  bind_rows(data_alltreated, data_control) %>%
  mutate(
    treatment_date = if_else(treated==1L, vax3_date, as.Date(NA)),
  )


local({

  ## sequential trial matching routine is as follows:
  # each daily trial includes all n people who were vaccinated on that day (treated=1) and
  # a sample of n controls (treated=0) who:
  # - had not been vaccinated on or before that day (still at risk of treatment);
  # - still at risk of an outcome (not deregistered or dead); 
  # - had not already been selected as a control in a previous trial


  # set maximum number of daily trials
  # time index is relative to "start date"
  # trial index start at one, not zero. i.e., study start date is "day 1" (but the _time_ at the start of study start date is zero)
  start_trial_time <- 0
  end_trial_time <- as.integer(dates$end_date + 1 - dates$start_date)
  trials <- seq(start_trial_time+1, end_trial_time, 1) 
  
  # initialise list of candidate controls
  candidate_ids <- data_control$patient_id

  # initialise matching summary data
  data_treated <- NULL
  data_matched <- NULL

  #trial=1
  for(trial in trials){

    cat("matching trial ", trial, "\n")
    trial_time <- trial-1
    trial_date <- dates$start_date + trial_time

    
    # set of people vaccinated on trial day
    data_treated_i <-
      data_eligible %>%
      filter(
        # select treated
        treated == 1L,
        # select people vaccinated on trial day i
        treatment_date == trial_date
        ) %>% 
      transmute(
        patient_id,
        treated,
        trial_time=trial_time,
        trial_date=trial_date
      )
    
    # append total treated on trial day i to all previous treated people
    data_treated <- bind_rows(data_treated, data_treated_i)

    # set of people still eligible for control inclusion on trial day
    data_control_i <-
      data_eligible %>%
      filter(
        # select controls
        treated==0L,
        # remove anyone already vaccinated
        (vax3_date > trial_date) | is.na(vax3_date),
        # select only people not already selected as a control
        patient_id %in% candidate_ids
      ) %>%
      transmute(
        patient_id,
        treated,
        trial_time=trial_time,
        trial_date=trial_date
      )
    
    
    n_treated_all <- nrow(data_treated_i)
    
    if(n_treated_all<1) {
      message("Skipping trial ", trial, " - No treated people eligible for inclusion.")
      next
    }
  
    matching_candidates_i <- 
      bind_rows(data_treated_i, data_control_i) %>%
      left_join(
        data_eligible %>% 
          select(
            patient_id, 
            treated,
            all_of(c(
              exact_variables, 
              names(caliper_variables)
              )),
        ),
        by = c("patient_id", "treated")
      )
    

    safely_matchit <- purrr::safely(matchit)
    
    # run matching algorithm
    obj_matchit_i <-
      safely_matchit(
        formula = treated ~ 1,
        data = matching_candidates_i,
        method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
        replace = FALSE,
        estimand = "ATT",
        exact = exact_variables,
        caliper = caliper_variables, std.caliper=FALSE,
        m.order = "data", # data is sorted on (effectively random) patient ID
        #verbose = TRUE,
        ratio = 1L # 1:1 matching
      )[[1]]

    
    if(is.null(obj_matchit_i)) {
      message("Skipping trial ", trial, " - No exact matches found.")
      next
    }
    
    
    
    data_matchstatus_i <-
      if(is.null(obj_matchit_i)){
        tibble(
          patient_id = matching_candidates_i$patient_id,
          matched = FALSE,
          #thread_id = data_thread$thread_id,
          match_id = NA_integer_,
          treated = matching_candidates_i$treated,
          weight = 0,
          trial_time = trial_time,
          trial_date = trial_date,
        )
      } else {
          tibble(
            patient_id = matching_candidates_i$patient_id,
            matched = !is.na(obj_matchit_i$subclass),
            #thread_id = data_thread$thread_id,
            match_id = as.integer(as.character(obj_matchit_i$subclass)),
            treated = obj_matchit_i$treat,
            weight = obj_matchit_i$weights,
            trial_time = trial_time,
            trial_date = trial_date,
          ) 
      } %>%
      arrange(match_id, treated)
    
    
    
    # summary info for recruited people
    # - one row per person
    # - match_id is within matching_i
    data_matched_i <-
      data_matchstatus_i %>%
      filter(!is.na(match_id)) %>% # remove unmatched people. equivalent to weight != 0
      arrange(match_id, desc(treated)) %>%
      left_join(
        data_eligible %>% select(patient_id, treated, vax3_date),
        by = c("patient_id", "treated")
      ) %>%
      group_by(match_id) %>%
      mutate(
        controlistreated_date = vax3_date[treated==0], # this only works because of the group_by statement above! do not remove group_by statement!
      ) %>%
      ungroup()

    n_treated_matched <- sum(data_matched_i$treated)

    # append matched data to matches from previous trials
    data_matched <- bind_rows(data_matched, data_matched_i)
    
    # update list of candidate controls to those who have not already been recruited
    candidate_ids <- candidate_ids[!(candidate_ids %in% data_matched_i$patient_id)]

  }

  #remove trial_time and trial_date counters created by the loop
  trial_time <- NULL
  trial_date <- NULL

  data_matched <-
    data_matched %>%
    transmute(
      patient_id, 
      match_id, 
      matched=1L, 
      treated,
      control=1L-treated, 
      trial_time, 
      trial_date, 
      controlistreated_date
    )

  # matching status for all treated people and their controls (if matched).
  # includes: unmatched treated; matched treated; matched control
  data_matchstatus <<-
    data_treated %>%
    left_join(data_matched %>% filter(treated==1L, matched==1L), by=c("patient_id", "treated", "trial_time", "trial_date")) %>%
    mutate(
      matched = replace_na(matched, 0L), # 1 if matched, 0 if unmatched
      control = if_else(matched==1L, 0L, NA_integer_) # 1 if matched control, 0 if matched treated, NA if unmatched treated
    ) %>%
    bind_rows(
      data_matched %>% filter(control==1L) %>% mutate(treated=0L)
    )
  
  unmatched_control_ids <<- candidate_ids
})
 
# output matching status ----
write_rds(data_matchstatus, ghere("output", cohort, "matchround{matching_round}", "potential", "data_potential_matchstatus.rds"), compress="gz")

# number of treated/controls per trial
with(data_matchstatus %>% filter(matched==1), table(trial_time, treated))

# total matched pairs
with(data_matchstatus %>% filter(matched==1), table(treated))

# max trial date
print(paste0("max trial day is ", as.integer(max(data_matchstatus %>% filter(matched==1) %>% pull(trial_time), na.rm=TRUE))))


# output csv for subsequent study definition
data_matchstatus %>% 
  filter(control==1L, matched==1L) %>% 
  select(patient_id, trial_date, match_id) %>%
  mutate(
    trial_date=as.character(trial_date)
  ) %>%
  write_csv(ghere("output", cohort, "matchround{matching_round}", "potential", glue("potential_matchedcontrols.csv.gz")))


print(paste0("number of duplicate control IDs is ", data_matchstatus %>% filter(control==1L, matched==1L) %>% group_by(patient_id) %>% summarise(n=n()) %>% filter(n>1) %>% nrow() ))
# should be zero
# 
# 
# ## output dataset containing all matched pairs + matching factors
# data_matched <-
#   data_matchstatus %>%
#   filter(matched==1L) %>%
#   left_join(
#     data_eligible %>%
#     select(
#       patient_id,
#       treated,
#       all_of(
#         exact_variables#,
#         #names(caliper_variables)
#       ),
#     ),
#     by=c("patient_id", "treated")
#   ) %>%
#   arrange(trial_date, match_id, treated)
# 
# 
# write_rds(data_matched, fs::path(output_dir, glue("data_potential_matched{matching_round}.rds")), compress="gz")
