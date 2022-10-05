######################################

# This script:
# imports data extracted by the cohort extractor (or dummy data)
# fills in unknown ethnicity from GP records with ethnicity from SUS (secondary care)
# tidies missing values
# standardises some variables (eg convert to factor) and derives some new ones
# organises vaccination date data to "vax X type", "vax X date" (rather than "pfizer X date", "az X date", ...)
######################################

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('lubridate')
library('arrow')
library('here')
library('glue')

## import local functions and parameters ---

source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))

source(here("analysis", "process_functions.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  # use for interactive testing
  stage <- "treated"
  # stage <- "potential"
  # stage <- "actual"
  # stage <- "final"
  # cohort <- "pfizer"
  # matching_round <- as.integer("1")
} else {
  stage <- args[[1]]
  
  if (stage == "treated") {
    if (length(args) > 1) 
      stop("No additional args to be specified when `stage=\"treated\"")
  } else if (stage %in% c("potential", "actual")) {
    if (length(args) == 1) {
      stop("`cohort` and `matching_round` must be specified when `stage=\"potential\"` or \"actual\"")
    }
    
    cohort <- args[[2]] # NULL if treated
    matching_round <- as.integer(args[[3]]) # NULL if treated    
    
  } else if (stage == "final") {
    if (length(args) == 1) {
      stop("`cohort` must be specified when `stage=\"final\"`")
    }
    
    cohort <- args[[2]] # NULL if treated
    
  }
} 

## get cohort-specific parameters study dates and parameters ---- 
if (stage == "potential") {
  matching_round_date <- study_dates[[cohort]]$control_extract_dates[matching_round]
}

## create output directory ----
if (stage == "treated") {
  fs::dir_create(here("output", "pfizer", "treated"))
  fs::dir_create(here("output", "moderna", "treated"))
  fs::dir_create(here("output", "treated", "eligible"))
} else if (stage == "potential") {
  fs::dir_create(ghere("output", cohort, "matchround{matching_round}", "process"))
} else if (stage == "actual") {
  fs::dir_create(ghere("output", cohort, "matchround{matching_round}", "actual"))
} else if (stage == "final") {
  fs::dir_create(ghere("output", cohort, "match"))
}


# import data ----

if (stage == "actual") {
  ## trial info for potential matches in round X
  data_potential_matchstatus <- 
    read_rds(ghere("output", cohort, "matchround{matching_round}", "potential", "data_potential_matchstatus.rds")) %>% 
    filter(matched==1L)
}

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data
  
  if (stage == "treated") {
    studydef_path <- here("output", "treated", "extract", "input_treated.feather")
    custom_path <- here("lib", "dummydata", "dummy_treated.feather")
  } else if (stage %in% c("potential", "actual")) {
    studydef_path <- ghere("output", cohort, "matchround{matching_round}", "extract", "input_controlpotential.feather")
    custom_path <- here("lib", "dummydata", "dummy_control_potential1.feather")
  }  else if (stage == "final") {
    studydef_path <- ghere("output", cohort, "extract", "input_controlfinal.feather")
    custom_path <- ghere("output", cohort, "dummydata", "dummy_control_final.feather")
  }
  
  data_studydef_dummy <- read_feather(studydef_path) %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"), ~ as.Date(.))) %>%
    # because of a bug in cohort extractor -- remove once pulled new version
    mutate(patient_id = as.integer(patient_id))
  
  data_custom_dummy <- read_feather(custom_path) 
  
  if (stage != "final") {
    data_custom_dummy <- data_custom_dummy %>%
      mutate(
        msoa = sample(factor(c("1", "2")), size=n(), replace=TRUE) # override msoa so matching success more likely
      )
  }
  
  if (stage == "actual") {
    # reuse previous extraction for dummy run, dummy_control_potential1.feather
    data_custom_dummy <- data_custom_dummy %>%
      filter(patient_id %in% data_potential_matchstatus[(data_potential_matchstatus$treated==0L),]$patient_id) %>%
      # change a few variables to simulate new index dates
      mutate(
        region = if_else(runif(n())<0.05, sample(x=unique(region), size=n(), replace=TRUE), region),
      ) 
  }
  
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
  
  if (stage == "actual") {
    
    data_extract <- data_extract %>%
      # these variables are not included in the dummy data so join them on here
      # they're joined in the study def using `with_values_from_file`
      left_join(data_potential_matchstatus %>% filter(treated==0L), by=c("patient_id")) %>%
      # remove vax variables
      select(-starts_with("covid_vax"))
    
  }
  
} else {
  
  if (stage == "treated") {
    
    data_extract <- read_feather(ghere("output", "treated", "extract", "input_treated.feather")) %>%
      #because date types are not returned consistently by cohort extractor
      mutate(across(ends_with("_date"),  as.Date))
    
  } else if (stage == "potential") {
    
    data_extract <- read_feather(ghere("output", cohort, "matchround{matching_round}", "extract", "input_controlpotential.feather")) %>%
      # because date types are not returned consistently by cohort extractor
      mutate(across(ends_with("_date"), as.Date))
    
  } else if (stage == "actual") {
    
    data_extract <- read_feather(ghere("output", cohort, "matchround{matching_round}", "extract", glue("input_controlactual.feather"))) %>%
      #because date types are not returned consistently by cohort extractor
      mutate(across(ends_with("_date"),  as.Date)) %>% 
      mutate(treated=0L) %>%
      # these variables are not included in the dummy data so join them on here
      # they're joined in the study def using `with_values_from_file`
      left_join(data_potential_matchstatus %>% filter(treated==0L), by=c("patient_id", "treated", "trial_date", "match_id"))
    
  } else if (stage == "final") {
    
    data_extract <- read_feather(ghere("output", cohort, "extract", "input_controlfinal.feather")) %>%
      #because date types are not returned consistently by cohort extractor
      mutate(across(ends_with("_date"),  as.Date))
    
  }
  
}

# process the final dataset ----
if (stage == "final") {
  
  data_matchstatus <- read_rds(ghere("output", cohort, "matchround{n_matching_rounds}", "actual", "data_matchstatus_allrounds.rds"))
  
  # import data for treated group and select those who were successfully matched
  data_treatedeligible <- read_rds(ghere("output", cohort, "treated", "data_treatedeligible.rds"))
  
  data_treated <- 
    left_join(
      data_matchstatus %>% filter(treated==1L),
      data_treatedeligible,
      by="patient_id"
    ) 
  
  # import extracted data from controls
  
  
  # import final dataset of matched controls, including matching variables
  # alternative to this is re-extracting everything in the study definition
  data_control <- 
    data_matchstatus %>% filter(treated==0L) %>%
    left_join(
      map_dfr(
        seq_len(n_matching_rounds), 
        ~{read_rds(ghere("output", cohort, glue("matchround", .x), "actual", "data_successful_matchedcontrols.rds"))}
      ) %>% select(-match_id, -trial_date, -treated, -controlistreated_date), # remove to avoid clash with already-stored variables
      by=c("patient_id", "matching_round")
    ) %>%
    # merge with outcomes data
    left_join(
      data_extract,
      by=c("patient_id", "match_id", "trial_date")
    ) %>%
    mutate(
      treated=0L
    )
  
  # check final data agrees with matching status
  
  all(data_control$patient_id %in% (data_matchstatus %>% filter(treated==0L) %>% pull(patient_id)))
  all((data_matchstatus %>% filter(treated==0L) %>% pull(patient_id)) %in% data_control$patient_id)
  
  # merge treated and control groups
  data_matched <-
    bind_rows(
      data_treated,
      data_control %>% process_outcome() # process the post-baseline variables (done previously for data_treated)
    ) 
  
  write_rds(data_matched, here("output", cohort, "match", "data_matched.rds"), compress="gz")
  
  # matching status of all treated, eligible people ----
  
  data_treatedeligible_matchstatus <- 
    left_join(
      data_treatedeligible %>% select(patient_id, vax3_date),
      data_matchstatus %>% filter(treated==1L),
      by="patient_id"
    ) %>%
    mutate(
      matched = if_else(is.na(match_id), 0L, 1L),
      treated = if_else(is.na(match_id), 1L, treated),
    )
  
  print(
    glue(
      "all trial dates match vaccination dates for matched, treated people: ",
      data_treatedeligible_matchstatus %>% 
        filter(matched==1L) %>%
        mutate(
          agree = trial_date==vax3_date
        ) %>% pull(agree) %>% all()
    )
  )
  
  write_rds(data_treatedeligible_matchstatus, here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"), compress="gz")
  
} 

# script stops here when stage = "final"


# process data -----

## patient-level info ----

if (stage %in% c("treated", "potential", "actual")) {
  data_processed <- data_extract %>%
    process_jcvi() %>%
    process_demo() %>%
    process_pre() 
}

if (stage == "treated") {
  data_processed <- data_processed %>%
    process_outcome()
}

## process vaccination data ----

if (stage %in% c("treated", "potential")) {
  
  data_processed <- data_processed %>%
    process_vax(stage)
  
} else if (stage == "actual") {
  
  ### join to vax data 
  data_vax_wide <- 
    read_rds(ghere("output", cohort, "matchround{matching_round}", "process", "data_controlpotential.rds")) %>%
    select(patient_id, matches("^vax\\d"))
  
  data_processed <- data_processed %>%
    left_join(data_vax_wide, by = "patient_id")
  
}

####################################################################################

if (stage == "treated") {
  selection_stage <- rlang::quos(
    
    has_expectedvax3type = vax3_type %in% c("pfizer", "moderna"),
    
    has_vaxgap23 = vax3_date >= (vax2_date+17) | is.na(vax3_date), # at least 17 days between second and third vaccinations
    
    vax3_notbeforestartdate = case_when(
      (vax3_type=="pfizer") & (vax3_date < study_dates$pfizer$start_date) ~ FALSE,
      (vax3_type=="moderna") & (vax3_date < study_dates$moderna$start_date) ~ FALSE,
      TRUE ~ TRUE
    ),
    vax3_beforeenddate = case_when(
      (vax3_type=="pfizer") & (vax3_date <= study_dates$pfizer$end_date) & !is.na(vax3_date) ~ TRUE,
      (vax3_type=="moderna") & (vax3_date <= study_dates$moderna$end_date) & !is.na(vax3_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    index_date = vax3_date,
    
    c0 = is_adult & vax3_date <= study_dates$studyend_date,
    c1 = c0 & vax3_notbeforestartdate & vax3_beforeenddate & has_expectedvax3type & has_vaxgap23,
    
  )
  
} else if (stage %in% c("potential",  "actual")) {
  
  # define index_date
  if (stage == "potential") index_date <- "matching_round_date" else if (stage == "actual") index_date <- "trial_date"
  
  selection_stage <- rlang::quos(
    
    index_date = !! sym(index_date),
    
    vax3_notbeforeindexdate = case_when(
      is.na(vax3_date) | (vax3_date > index_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    c0 = is_adult,
    c1 = c0 & vax3_notbeforeindexdate,
    
  )
  
} 

# Define selection criteria ----
if (stage %in% c("treated", "potential", "actual")) {
  
data_criteria <- data_processed %>%
  transmute(
    
    patient_id,
    is_adult = age >= 18,
    has_age = !is.na(age),
    has_sex = !is.na(sex),
    has_imd = imd_Q5 != "Unknown",
    has_ethnicity = !is.na(ethnicity_combined),
    has_region = !is.na(region),
    #has_msoa = !is.na(msoa),
    isnot_hscworker = !hscworker,
    isnot_carehomeresident = !care_home_combined,
    isnot_endoflife = !endoflife,
    isnot_housebound = !housebound,
    
    vax1_afterfirstvaxdate = case_when(
      (vax1_type=="pfizer") & (vax1_date >= study_dates$firstpfizer_date) ~ TRUE,
      (vax1_type=="az") & (vax1_date >= study_dates$firstaz_date) ~ TRUE,
      (vax1_type=="moderna") & (vax1_date >= study_dates$firstmoderna_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    vax2_beforelastvaxdate = !is.na(vax2_date) & (vax2_date <= study_dates$lastvax2_date),
    
    has_knownvax1 = vax1_type %in% c("pfizer", "az"),
    has_knownvax2 = vax2_type %in% c("pfizer", "az"),
    
    vax12_homologous = vax1_type==vax2_type,
    has_vaxgap12 = vax2_date >= (vax1_date+17), # at least 17 days between first two vaccinations
    
    !!! selection_stage,
    
    no_recentcovid30 = is.na(anycovid_0_date) | ((index_date - anycovid_0_date) > 30),
    
    isnot_inhospital = is.na(admitted_unplanned_0_date) | (!is.na(discharged_unplanned_0_date) & discharged_unplanned_0_date < index_date),
    
    c2 = c1 & vax1_afterfirstvaxdate & vax2_beforelastvaxdate & has_vaxgap12 & has_knownvax1 & has_knownvax2 & vax12_homologous,
    c3 = c2 & isnot_hscworker,
    c4 = c3 & isnot_carehomeresident & isnot_endoflife & isnot_housebound,
    c5 = c4 & has_age & has_sex & has_imd & has_ethnicity & has_region,
    c6 = c5 & no_recentcovid30,
    c7 = c6 & isnot_inhospital,
    c8 = c7 & TRUE, # TODO define c8 (this will be TRUE when stage!=treated)
    
    include = c8,
    
  )

data_eligible <- data_criteria %>%
  filter(include) %>%
  select(patient_id) %>%
  left_join(data_processed, by="patient_id") %>%
  droplevels()

}

# save cohort-specific datasets ----
if (stage == "treated") {
  
  write_rds(data_eligible %>% filter(vax3_type == "pfizer"), 
            here("output", "pfizer", "treated", "data_treatedeligible.rds"),
            compress="gz")
  
  write_rds(data_eligible %>% filter(vax3_type == "moderna"), 
            here("output", "moderna", "treated", "data_treatedeligible.rds"), 
            compress="gz")
  
} else if (stage == "potential") {
  
  write_rds(data_eligible, 
            ghere("output", cohort, "matchround{matching_round}", "process", "data_controlpotential.rds"),
            compress = "gz")
  
}


# create flowchart (only when stage="treated") ----
if (stage == "treated") {
  
  data_flowchart <- data_criteria %>%
    summarise(
      across(matches("^c\\d"), .fns=sum)
    ) %>%
    pivot_longer(
      cols=everything(),
      names_to="criteria",
      values_to="n"
    ) %>%
    mutate(
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n),
      crit = str_extract(criteria, "^c\\d+"),
      criteria = fct_case_when(
        crit == "c0" ~ "Aged 18+ with 3rd dose on or before {study_dates$studyend_date}", 
        crit == "c1" ~ "  at least 17 days between 2nd and 3rd dose and 3rd dose of pfizer received between {study_dates$pfizer$start_date} and {study_dates$pfizer$end_date} or 3rd dose of moderna received between {study_dates$moderna$start_date} and {study_dates$moderna$end_date}",
        crit == "c2" ~ "  homologous primary vaccination course of pfizer or AZ and at least 17 days between doses",
        crit == "c3" ~ "  not a HSC worker",
        crit == "c4" ~ "  not a care/nursing home resident, end-of-life or housebound",
        crit == "c5" ~ "  no missing demographic information",
        crit == "c6" ~ "  no evidence of covid in 30 days before third dose",
        crit == "c7" ~ "  not in hospital (unplanned) during booster vaccination",
        crit == "c8" ~ "  did not received 3rd dose at unusual time given region, priority group, and 2nd dose date.", #TODO
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(across(criteria, factor, labels = sapply(levels(.$criteria), glue)))
  
  write_rds(data_flowchart, here("output", "treated", "eligible", "flowchart_treatedeligible.rds"))
  
}

# check matching (only when stage="actual") ----
if (stage == "actual") { 
  
  data_control <- data_eligible
  
  data_treated <- 
    left_join(
      data_potential_matchstatus %>% filter(treated==1L),
      read_rds(ghere("output", cohort, "treated", "data_treatedeligible.rds")) %>% 
        # only keep variables that are in data_control (this gets rid of outcomes and vax4 dates)
        select(any_of(names(data_control))),
      by="patient_id"
    )
  
  matching_candidates <- 
    bind_rows(data_treated, data_control) %>%
    arrange(treated, match_id, trial_date)
  
  #print missing values
  matching_candidates_missing <- map(matching_candidates, ~any(is.na(.x)))
  sort(names(matching_candidates_missing[unlist(matching_candidates_missing)]))
  
  # run matching algorithm ----
  obj_matchit <-
    MatchIt::matchit(
      formula = treated ~ 1,
      data = matching_candidates,
      method = "nearest", distance = "glm", # these two options don't really do anything because we only want exact + caliper matching
      replace = FALSE,
      estimand = "ATT",
      exact = c("match_id", "trial_date", exact_variables),
      caliper = caliper_variables, std.caliper=FALSE,
      m.order = "data", # data is sorted on (effectively random) patient ID
      #verbose = TRUE,
      ratio = 1L # irritatingly you can't set this for "exact" method, so have to filter later
    )
  
  
  data_matchstatus <-
    tibble(
      patient_id = matching_candidates$patient_id,
      matched = !is.na(obj_matchit$subclass)*1L,
      #thread_id = data_thread$thread_id,
      match_id = as.integer(as.character(obj_matchit$subclass)),
      treated = obj_matchit$treat,
      #weight = obj_matchit$weights,
      trial_time = matching_candidates$trial_time,
      trial_date = matching_candidates$trial_date,
      matching_round = matching_round,
      # controlistreated_date = matching_candidates$controlistreated_date
    ) %>%
    arrange(matched, match_id, treated) 
  
  
  ###
  matchstatus_vars <- c("patient_id", "match_id", "trial_date", "matching_round", "treated", "controlistreated_date")
  
  data_successful_matchstatus <- 
    data_matchstatus %>% 
    filter(matched) %>%
    left_join(
      # now joining all variables from the processed data as they are required for adjustments in the cox model
      matching_candidates %>% select(-all_of(c("trial_time", "trial_date", "match_id", "matched", "control", "controlistreated_date"))),
      by = c("patient_id", "treated")
    ) %>%
    group_by(match_id) %>%
    mutate(
      controlistreated_date = vax3_date[treated==0], # this only works because of the group_by statement above! do not remove group_by statement!
    ) %>%
    ungroup() %>%
    select(all_of(matchstatus_vars), everything())
  
  ## size of dataset
  print("data_successful_match treated/untreated numbers")
  table(treated = data_successful_matchstatus$treated, useNA="ifany")
  
  
  ## how many matches are lost?
  
  print(glue("{sum(data_successful_matchstatus$treated)} matched-pairs kept out of {sum(data_potential_matchstatus$treated)} 
           ({round(100*(sum(data_successful_matchstatus$treated) / sum(data_potential_matchstatus$treated)),2)}%)
           "))
  
  
  ## pick up all previous successful matches ----
  
  if(matching_round>1){
    
    data_matchstatusprevious <- 
      read_rds(ghere("output", cohort, "matchround{matching_round-1}", "actual", "data_matchstatus_allrounds.rds"))
    
    data_matchstatus_allrounds <- 
      data_successful_matchstatus %>% 
      select(all_of(matchstatus_vars)) %>%
      bind_rows(data_matchstatusprevious) 
    
  } else{
    data_matchstatus_allrounds <- 
      data_successful_matchstatus %>%
      select(all_of(matchstatus_vars))
  }
  
  write_rds(data_matchstatus_allrounds, ghere("output", cohort, "matchround{matching_round}", "actual", "data_matchstatus_allrounds.rds"), compress="gz")
  
  
  # output all control patient ids for finalmatched study definition
  data_matchstatus_allrounds %>%
    mutate(
      trial_date=as.character(trial_date)
    ) %>%
    filter(treated==0L) %>% #only interested in controls as all
    write_csv(ghere("output", cohort, "matchround{matching_round}", "actual", "cumulative_matchedcontrols.csv.gz"))
  
  ## size of dataset
  print("data_matchstatus_allrounds treated/untreated numbers")
  table(treated = data_matchstatus_allrounds$treated, useNA="ifany")
  
  
  
  ## duplicate IDs
  data_matchstatus_allrounds %>% group_by(treated, patient_id) %>%
    summarise(n=n()) %>% group_by(treated) %>% summarise(ndups = sum(n>1)) %>%
    print()
  
  
  write_rds(data_successful_matchstatus %>% filter(treated==0L), ghere("output", cohort, "matchround{matching_round}", "actual", "data_successful_matchedcontrols.rds"), compress="gz")
  
  ## size of dataset
  print("data_successful_match treated/untreated numbers")
  table(treated = data_successful_matchstatus$treated, useNA="ifany")
  
}

