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
  #stage <- "potential"
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
  fs::dir_create(here("output", "az", "treated"))
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
  
  ## set seed so results on dummy data are reproducible ---
  set.seed(10)
  
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
      data_treatedeligible %>% select(patient_id, vax1_date),
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
          agree = trial_date==vax1_date
        ) %>% pull(agree) %>% all()
    )
  )
  
  write_rds(data_treatedeligible_matchstatus, here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"), compress="gz")
  
} 

# script stops here when stage = "final"


# process data -----

## define index data ----
if (stage == "treated") {
  data_extract <- data_extract %>%
    mutate(index_date = covid_vax_disease_1_date) 
} else if(stage == "potential"){
  data_extract <- data_extract %>%
    mutate(index_date = matching_round_date) 
} else if(stage == "actual"){
  data_extract <- data_extract %>%
    mutate(index_date = trial_date) 
}

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
    left_join(data_vax_wide, by = "patient_id") %>%
    # the following line is needed for applying the eligibility criteria: covid_vax_disease_1_date_matches_vax1_date
    # it has already been checked that this is true in the process_potential stage, 
    # but `covid_vax_disease_1_date` is added to avoid having to add extra logic statements for the case when stage="actual"
    mutate(covid_vax_disease_1_date = vax1_date)
  
}



####################################################################################

if (stage == "treated") {
  selection_stage <- rlang::quos(
    
    has_expectedvax1type = vax1_type %in% c("pfizer", "az"),
    
    has_vaxgap12 = vax2_date >= (vax1_date+17) | is.na(vax2_date), # at least 17 days between first and second vaccinations. this is post-baseline conditioning but is essentially just removing a small number of people with unreliable vaccination data
    
    vax1_notbeforestartdate = case_when(
      (vax1_type=="pfizer") & (vax1_date < study_dates$pfizer$start_date) ~ FALSE,
      (vax1_type=="az") & (vax1_date < study_dates$az$start_date) ~ FALSE,
      TRUE ~ TRUE
    ),
    vax1_beforeenddate = case_when(
      (vax1_type=="pfizer") & (vax1_date <= study_dates$pfizer$end_date) & !is.na(vax1_date) ~ TRUE,
      (vax1_type=="az") & (vax1_date <= study_dates$az$end_date) & !is.na(vax1_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    vax1_notbeforeageeligible = case_when(
      jcvi_ageband %in% c("80+") & vax1_date < study_dates$over80s$start_date ~ FALSE,
      jcvi_ageband %in% c("70-74", "75-79") & vax1_date < study_dates$in70s$start_date ~ FALSE,
      TRUE ~ TRUE # ignore agebands under 70 years old as these are not being studied here
    ),
    
    c0 = vax1_notbeforestartdate & vax1_beforeenddate & vax1_notbeforeageeligible,
    c1 = c0 & has_expectedvax1type & has_vaxgap12  & covid_vax_disease_1_date_matches_vax1_date,
    
  )
  
} else if (stage %in% c("potential",  "actual")) {
  
  selection_stage <- rlang::quos(
    
    vax1_notbeforeindexdate = case_when(
      is.na(vax1_date) | (vax1_date > index_date) ~ TRUE,
      TRUE ~ FALSE
    ),
    
    vax1_notbeforeageeligible = case_when(
      ageband2 %in% c("80+") & vax1_date < study_dates$over80s$start_date ~ FALSE,
      ageband2 %in% c("70-79") & vax1_date < study_dates$in70s$start_date ~ FALSE,
      TRUE ~ TRUE # ignore agebands under 70 years old as these are not being studied here
    ),
    
    c0 = TRUE,
    c1 = c0 & vax1_notbeforeindexdate & vax1_notbeforeageeligible & covid_vax_disease_1_date_matches_vax1_date,
    
  )
  
} 

# Define selection criteria ----
if (stage %in% c("treated", "potential", "actual")) {
  
  data_criteria <- data_processed %>%
    left_join(
      data_extract %>% select(patient_id, matches("covid_vax_disease_\\d_date")),
      by = "patient_id"
      ) %>%
    transmute(
      
      patient_id,
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
      
      covid_vax_disease_1_date_matches_vax1_date = covid_vax_disease_1_date == vax1_date,
      
      !!! selection_stage,
      
      no_recentcovid30 = is.na(anycovid_0_date) | ((index_date - anycovid_0_date) > 30),
      
      isnot_inhospital = is.na(admitted_unplanned_0_date) | (!is.na(discharged_unplanned_0_date) & discharged_unplanned_0_date < index_date),
      
      c2 = c1 & isnot_hscworker,
      c3 = c2 & isnot_carehomeresident & isnot_endoflife & isnot_housebound,
      c4 = c3 & has_age & has_sex & has_imd & has_ethnicity & has_region,
      c5 = c4 & no_recentcovid30,
      c6 = c5 & isnot_inhospital,
      
      include = c6,
      
    )
  
  data_eligible <- data_criteria %>%
    filter(include) %>%
    select(patient_id) %>%
    left_join(data_processed, by="patient_id") %>%
    droplevels()
  
}

# save cohort-specific datasets ----
if (stage == "treated") {
  
  write_rds(data_eligible %>% filter(vax1_type == "pfizer"), 
            here("output", "pfizer", "treated", "data_treatedeligible.rds"),
            compress="gz")
  
  write_rds(data_eligible %>% filter(vax1_type == "az"), 
            here("output", "az", "treated", "data_treatedeligible.rds"), 
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
        crit == "c0" ~ "Aged 70+ with 1st dose between study dates", 
        crit == "c1" ~ "  no unreliable vaccination data",
        crit == "c2" ~ "  not a HSC worker",
        crit == "c3" ~ "  not a care/nursing home resident, end-of-life or housebound",
        crit == "c4" ~ "  no missing demographic information",
        crit == "c5" ~ "  no evidence of covid in 30 days before trial date",
        crit == "c6" ~ "  not in hospital (unplanned) on trial date",
        TRUE ~ NA_character_
      )
    ) %>%
    mutate(across(criteria, factor, labels = sapply(levels(.$criteria), glue)))
  
  write_rds(data_flowchart, here("output", "treated", "eligible", "flowchart_treatedeligible.rds"))
  
  data_flowchart %>%
    transmute(
      criteria, crit, 
      n = ceiling_any(n, to=7),
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n),
    ) %>%
    write_csv(here("output", "treated", "eligible", "flowchart_treatedeligible_rounded.csv")) 
  
  # distribution of vax1_date by jcvi_ageband
  vax1_date_plot <- data_eligible %>%
    select(patient_id, vax1_date, vax1_type, jcvi_ageband) %>%
    ggplot(aes(x = vax1_date, colour = vax1_type)) +
    geom_freqpoly(binwidth=1) +
    facet_grid(jcvi_ageband~., scales = "free_y")
  ggsave(
    filename = here("output", "treated", "eligible", "vax1_dates.png"),
    plot = vax1_date_plot,
    width=20, height=15, units="cm"
  )
  
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
  
  # rematch ----
  rematch <-
    # first join on exact variables + match_id + trial_date
    inner_join(
      x=data_treated %>% select(match_id, trial_date, all_of(c(names(caliper_variables), exact_variables))),
      y=data_control %>% select(match_id, trial_date, all_of(c(names(caliper_variables), exact_variables))),
      by = c("match_id", "trial_date", exact_variables)
    ) 
  
  
  if(length(caliper_variables) >0 ){
    # check caliper_variables are still within caliper
    rematch <- rematch %>%
      bind_cols(
        map_dfr(
          set_names(names(caliper_variables), names(caliper_variables)),
          ~ abs(rematch[[str_c(.x, ".x")]] - rematch[[str_c(.x, ".y")]]) <= caliper_variables[.x]
        )
      ) %>%
      # dplyr::if_all not in opensafely version of dplyr so use filter_at instead
      # filter(if_all(
      #   all_of(names(caliper_variables))
      # )) 
      filter_at(
        all_of(names(caliper_variables)),
        all_vars(.)
      )
    
    
  } 
  
  rematch <- rematch %>%
    select(match_id, trial_date) %>%
    mutate(matched=1)
  
  data_successful_match <-
    matching_candidates %>%
    inner_join(rematch, by=c("match_id", "trial_date", "matched")) %>%
    mutate(
      matching_round = matching_round
    ) %>%
    arrange(trial_date, match_id, treated)
  
  
  ###
  
  matchstatus_vars <- c("patient_id", "match_id", "trial_date", "matching_round", "treated", "controlistreated_date")
  
  data_successful_matchstatus <- 
    data_successful_match %>% 
    # keep all variables from the processed data as they are required for adjustments in the cox model
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

