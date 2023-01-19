# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# Does data processing, summarises data, and applies eligibility criteria
# Many steps are stage-specific, controlled by if statements
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# Import libraries
library('tidyverse')
library('lubridate')
library('arrow')
library('here')
library('glue')

# import local functions and parameters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "process", "process_functions.R"))

# import command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  # use for interactive testing
  stage <- "single"
  # stage <- "treated"
  # stage <- "potential"
  # stage <- "actual"
  # stage <- "final"
  if (stage %in% c("potential", "actual", "final")) {
    brand <- "pfizer"
    if (stage != "final") matching_round <- as.integer("1")
  }
} else {
  stage <- args[[1]]
  
  if (stage %in% c("treated", "single")) {
    if (length(args) > 1) 
      stop("No additional args to be specified when `stage`=\"treated\" or `stage`=\"single\"")
  } else if (stage %in% c("potential", "actual")) {
    if (length(args) == 1) {
      stop("`brand` and `matching_round` must be specified when `stage=\"potential\"` or \"actual\"")
    }
    
    # brand and matching_round undefined if stage=treated or stage=single
    brand <- args[[2]] 
    matching_round <- as.integer(args[[3]]) 
    
  } else if (stage == "final") {
    if (length(args) == 1) {
      stop("`brand` must be specified when `stage=\"final\"`")
    }
    
    brand <- args[[2]] # undefined if treated
    
  }
} 

# define brand and matching round when stage=single
if (stage == "single") {
  brand <- "pfizer"
  matching_round <- 1
}
# extract matching_round_date
if (stage %in% c("single", "potential")) {
  matching_round_date <- study_dates[[brand]]$control_extract_dates[matching_round]
}

# create output directories
if (stage == "single") {
  fs::dir_create(here("output", "single", "eligible"))
  fs::dir_create(here("output", "single", "process"))
} else if (stage == "treated") {
  fs::dir_create(here("output", "sequential", "pfizer", "treated"))
  fs::dir_create(here("output", "sequential", "az", "treated"))
  fs::dir_create(here("output", "sequential", "treated", "eligible"))
  fs::dir_create(here("output", "sequential", "treated", "process"))
} else if (stage == "potential") {
  fs::dir_create(ghere("output", "sequential", brand, "matchround{matching_round}", "process"))
  fs::dir_create(ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "potential"))
  fs::dir_create(ghere("output", "sequential", brand, "matchround{matching_round}", "potential"))
} else if (stage == "actual") {
  fs::dir_create(ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "actual"))
  fs::dir_create(ghere("output", "sequential", brand, "matchround{matching_round}", "actual"))
} else if (stage == "final") {
  fs::dir_create(ghere("output", "sequential", brand, "match"))
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# import data ----

if (stage == "actual") {
  # trial info for potential matches in matching_round
  data_potential_matchstatus <- 
    read_rds(ghere("output", "sequential", brand, "matchround{matching_round}", "potential", "data_potential_matchstatus.rds")) %>% 
    filter(matched==1L)
}

# use externally created dummy data if not running in the server
# check variables are as they should be
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){
  
  # set seed so results on dummy data are reproducible
  set.seed(10)
  
  # ideally in future this will check column existence and types from metadata,
  # rather than from a cohort-extractor-generated dummy data
  
  if (stage == "treated") {
    studydef_path <- here("output", "sequential", "treated", "extract", "input_treated.feather")
    custom_path <- here("lib", "dummydata", "dummy_treated.feather")
  } else if (stage %in% c("single", "potential")) {
    studydef_path <- ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "input_controlpotential.feather")
    custom_path <- here("lib", "dummydata", "dummy_control_potential1.feather")
  } else if (stage == "actual") {
    studydef_path <- ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "input_controlactual.feather")
    custom_path <- here("lib", "dummydata", "dummy_control_potential1.feather")
  } else if (stage == "final") {
    studydef_path <- ghere("output", "sequential", brand, "extract", "input_controlfinal.feather")
    custom_path <- ghere("output", "sequential", brand, "dummydata", "dummy_control_final.feather")
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
      # remove vaccine variables
      select(-starts_with("covid_vax_")) %>%
      # trial_date and match_id are not included in the dummy data so join them on here
      # they're joined in the study def using `with_values_from_file`
      left_join(
        data_potential_matchstatus %>% 
          filter(treated==0L) %>%
          select(patient_id, trial_date, match_id),
        by="patient_id"
      ) %>%
      # add some random noise to region so that matching fails for some individuals
      mutate(
        region = if_else(runif(n())<0.05, sample(x=unique(region), size=n(), replace=TRUE), region),
      ) 
  }

  # check all variables match in data_custom_dummy and data_studydef_dummy
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
  
  # check column classes are identical
  unmatched_types <- cbind(
    map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")),
    map_chr(data_custom_dummy, ~paste(class(.), collapse=", "))
  )[ (map_chr(data_studydef_dummy, ~paste(class(.), collapse=", ")) != map_chr(data_custom_dummy, ~paste(class(.), collapse=", ")) ), ] %>%
    as.data.frame() %>% 
    rownames_to_column()
  
  if(nrow(unmatched_types)>0) stop(
    #unmatched_types
    "inconsistent typing in studydef : dummy dataset\n",
    apply(unmatched_types, 1, function(row) paste(paste(row, collapse=" : "), "\n"))
  )
  
  data_extract <- data_custom_dummy 
  
} else {
  
  # if running in the server, simply read the extracted data
  if (stage == "treated") {
    data_extract <- read_feather(ghere("output", "sequential", "treated", "extract", "input_treated.feather")) 
  } else if (stage %in% c("single", "potential")) {
    data_extract <- read_feather(ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "input_controlpotential.feather")) 
  } else if (stage == "actual") {
    data_extract <- read_feather(ghere("output", "sequential", brand, "matchround{matching_round}", "extract", glue("input_controlactual.feather"))) 
  } else if (stage == "final") {
    data_extract <- read_feather(ghere("output", "sequential", brand, "extract", "input_controlfinal.feather")) 
  }
  
  data_extract <- data_extract %>%
    # because date types are not returned consistently by cohort extractor
    mutate(across(ends_with("_date"),  as.Date))
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# process the dataset ----
if (stage == "final") {
  
  # summarise extracted data
  my_skim(data_extract, path = ghere("output", "sequential", brand, "extract", "input_control{stage}_skim.txt"))
  
  data_matchstatus <- read_rds(ghere("output", "sequential", brand, "matchround{n_matching_rounds}", "actual", "data_matchstatus_allrounds.rds"))
  
  # import data for treated group and select those who were successfully matched
  data_treatedeligible <- read_rds(ghere("output", "sequential", brand, "treated", "data_treatedeligible.rds"))
  
  data_treated <- 
    left_join(
      data_matchstatus %>% filter(treated==1L),
      data_treatedeligible,
      by="patient_id"
    ) 
  
  # import final dataset of matched controls, including matching variables
  # alternative to this is re-extracting everything in the study definition
  data_control <- 
    data_matchstatus %>%
    filter(treated==0L) %>%
    left_join(
      map_dfr(
        seq_len(n_matching_rounds), 
        ~{read_rds(ghere("output", "sequential", brand, glue("matchround", .x), "actual", "data_successful_matchedcontrols.rds"))}
      ) %>% 
        # remove to avoid clash with already-stored variables
        select(-match_id, -trial_date, -treated, -controlistreated_date), 
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
      data_control 
    ) 
  write_rds(
    data_matched, 
    here("output", "sequential", brand, "match", "data_matched.rds"), 
    compress="gz"
    )
  
  # summarise matched data by treatment group
  data_matched %>% filter(treated==0) %>%
    my_skim(
      path = here("output", "sequential", brand, "match", "data_matched_control_skim.txt")
    )
  data_matched %>% filter(treated==1) %>%
    my_skim(
      path = here("output", "sequential", brand, "match", "data_matched_treated_skim.txt")
    )
  
  # matching status of all treated, eligible people
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
  # check trial dates and vaccination dates match
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
  # save dataset
  write_rds(
    data_treatedeligible_matchstatus, 
    here("output", "sequential", brand, "match", "data_treatedeligible_matchstatus.rds"),
    compress="gz"
    )
  
} 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# script stops here when stage = "final"
# ensure all code beyond this point wrapped in `if` statements conditional on `stage`
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# summarise and define index_date ----
if (stage == "single") {
  
  # no need to summarise as had already been summarised when stage="potential" and matching round=1
  # define index date
  data_extract <- data_extract %>%
    mutate(index_date = matching_round_date) 
  
} else if (stage == "treated") {
  
  # summarise extracted data
  my_skim(
    data_extract, 
    path = here("output", "sequential", "treated", "extract", "input_treated_skim.txt")
    )
  # define index date
  data_extract <- data_extract %>%
    mutate(index_date = covid_vax_disease_1_date) 
  
} else if(stage == "potential"){
  
  # summarise extracted data
  my_skim(
    data_extract, 
    path = ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "potential", "input_controlpotential_skim.txt")
    )
  # define index date
  data_extract <- data_extract %>%
    mutate(index_date = matching_round_date) 
  
} else if(stage == "actual") {
  
  # add certain matching variables when stage=actual
  data_extract <- data_extract %>%
    # add: treated 
    mutate(treated=0L) %>%
    # add: trial_time, matched, control, controlistreated_date to data_extract
    left_join(
      data_potential_matchstatus %>%
        filter(treated==0L),
      by=c("patient_id", "treated", "trial_date", "match_id")
    )
  
  # summarise extracted data
  my_skim(
    data_extract, 
    path = ghere("output", "sequential", brand, "matchround{matching_round}", "extract", "actual", "input_controlactual_skim.txt")
    )
  # define index date
  data_extract <- data_extract %>%
    mutate(index_date = trial_date) 
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# process jcvi, demo and pre variables ----
# the process_*() functions are defined in analysis/process/process_functions.R
if (stage %in% c("single", "treated", "potential", "actual")) {
  data_processed <- data_extract %>%
    process_jcvi() %>%
    process_demo() %>%
    process_pre() 
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# process vaccination data ----
if (stage %in% c("single", "treated", "potential")) {
  
  data_processed <- data_processed %>%
    process_vax(stage)
  
} else if (stage == "actual") {
  
  ### join to vax data 
  data_vax_wide <- 
    read_rds(ghere("output", "sequential", brand, "matchround{matching_round}", "process", "data_controlpotential.rds")) %>%
    select(patient_id, matches("^vax\\d"))
  
  data_processed <- data_processed %>%
    left_join(data_vax_wide, by = "patient_id") %>%
    # the following line is needed for applying the eligibility criteria: covid_vax_disease_1_date_matches_vax1_date
    # it has already been checked that this is true in the process_potential stage, 
    # but `covid_vax_disease_1_date` is added to avoid having to add extra logic statements for the case when stage="actual"
    mutate(covid_vax_disease_1_date = vax1_date)
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# summarise processed data ----
if (stage %in% c("single", "treated", "potential", "actual")) {
  if (stage %in% "single") {
    skim_path <- here("output", "single", "process", "data_processed_skim.txt")
  } else if (stage %in% "treated") {
    skim_path <- here("output", "sequential", "treated", "process", "data_processed_skim.txt")
  } else {
    skim_path <- ghere("output", "sequential", brand, "matchround{matching_round}", stage, "data_processed_skim.txt")
  }
  my_skim(data_processed, path = skim_path)
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# apply eligibility criteria ----
if (stage %in% c("single", "treated", "potential", "actual")) {
  
  # Define the eligibility criteria for the flowchart. 
  # These must match the definitions of c0:c8 in the data.
  criteria_descr <- character()
  if (stage == "single") {
    criteria_descr["Aged 70+"] <- "c0"
  } else if (stage == "treated") {
    criteria_descr["Aged 70+ with 1st dose between study dates and eligible for single trial"] <- "c0"
  } else if (stage %in% c("potential", "actual")) {
    # c0 not straightforward to define when stage%in%c(potential,actual)
    # but we don't create flowcharts at those stages, so ok not to define c0 here
  }
  criteria_descr["  no unreliable vaccination data"] <- "c1"
  criteria_descr["  at least 1 year continuous registration"] <- "c2"
  criteria_descr["  not a HSC worker"] <- "c3"
  criteria_descr["  not a care/nursing home resident or housebound"] <- "c4"
  criteria_descr["  not receiving end-of-life care"] <- "c5"
  criteria_descr["  no missing demographic information"] <- "c6"
  if (stage=="single") {
    criteria_descr["  no evidence of covid before eligible for vaccination"] <- "c7"
    include <- "c7"
  } else {
    criteria_descr["  no evidence of covid before trial date"] <- "c7"
    criteria_descr["  not in hospital (unplanned) on trial date"] <- "c8"
    include <- "c8"
  }
  criteria_descr <- sort(criteria_descr)
  
  # specify stage-dependent criteria
  if (stage == "treated") {
    
    selection_stage <- rlang::quos(
      
      # Define the brands of interest for dose 1.
      has_expectedvax1type = vax1_type %in% c("pfizer", "az"),
      
      # At least 17 days between first and second vaccinations. 
      # This is post-baseline conditioning but is essentially just removing a small number of people with unreliable vaccination data.
      has_vaxgap12 = vax2_date >= (vax1_date+17) | is.na(vax2_date),
      
      # Check the vaccination was not recorded before the brand start date.
      vax1_notbeforestartdate = case_when(
        (vax1_type=="pfizer") & (vax1_date < study_dates$pfizer$start_date) ~ FALSE,
        (vax1_type=="az") & (vax1_date < study_dates$az$start_date) ~ FALSE,
        TRUE ~ TRUE
      ),
      # Check the vaccination was not recorded after the study end date.
      vax1_beforeenddate = case_when(
        (vax1_type=="pfizer") & (vax1_date <= study_dates$pfizer$end_date) & !is.na(vax1_date) ~ TRUE,
        (vax1_type=="az") & (vax1_date <= study_dates$az$end_date) & !is.na(vax1_date) ~ TRUE,
        TRUE ~ FALSE
      ),
      
      # Define criteria c0 and c1 when stage=treated
      c0 = vax1_notbeforestartdate & vax1_beforeenddate & vax1_notbeforeageeligible,
      c1 = c0 & has_expectedvax1type & has_vaxgap12  & covid_vax_disease_1_date_matches_vax1_date,
      
    )
    
  } 
  if (stage %in% c("single", "potential",  "actual")) {
    
    # function so that when stage=single, include people who were vaccinated *on or after* index_date,
    # but when stage%in%c(potential,actual), include people who were vaccinated *after* index_date
    vax1_notbeforeindexdate_fun <- function(stage, vax1_date, index_date) {
      
      if (stage == "single") {
        # on or after index_date
        out <- case_when(
          is.na(vax1_date) | (vax1_date >= index_date) ~ TRUE,
          TRUE ~ FALSE
        )
      } 
      if (stage %in% c("potential",  "actual")) {
        # after index date
        out <- case_when(
          is.na(vax1_date) | (vax1_date > index_date) ~ TRUE,
          TRUE ~ FALSE
        )
      }
      
      return(out)
      
    }
    
    selection_stage <- rlang::quos(
      
      # Define vax1_notbeforeindexdate using vax1_notbeforeindexdate_fun function defined above.
      vax1_notbeforeindexdate = vax1_notbeforeindexdate_fun(stage, vax1_date, index_date),
      # Define criteria c0 and c1 when stage=treated
      c0 = TRUE,
      c1 = c0 & vax1_notbeforeindexdate & vax1_notbeforeageeligible & covid_vax_disease_1_date_matches_vax1_date,
      
    )
    
  } 
  
  # define the eligibility criteria in the data
  data_criteria <- data_processed %>%
    left_join(
      data_extract %>% select(patient_id, matches("covid_vax_disease_\\d_date")),
      by = "patient_id"
      ) %>%
    transmute(
      
      patient_id,
      vax1_type,
      has_age = !is.na(age),
      has_sex = !is.na(sex),
      has_imd = imd_Q5 != "Unknown",
      has_ethnicity = !is.na(ethnicity_combined),
      has_region = !is.na(region),
      isnot_hscworker = !hscworker,
      isnot_carehomeresident = !care_home_combined,
      isnot_endoflife = !endoflife,
      isnot_housebound = !housebound,
      
      covid_vax_disease_1_date_matches_vax1_date = case_when(
        # FALSE if only one is missing
        is.na(covid_vax_disease_1_date) & !is.na(vax1_date) ~ FALSE,
        !is.na(covid_vax_disease_1_date) & is.na(vax1_date) ~ FALSE,
        # TRUE if both missing
        is.na(covid_vax_disease_1_date) & is.na(vax1_date) ~ TRUE,
        # TRUE if equal
        covid_vax_disease_1_date == vax1_date ~ TRUE,
        # FALSE otherwise
        TRUE ~ FALSE
        ),
      
      # Check the individual was not vaccinated before eligible.
      vax1_notbeforeageeligible = case_when(
        is.na(vax1_date) ~ TRUE,
        ageband2 %in% c("80+") & vax1_date < study_dates$over80s$start_date ~ FALSE,
        ageband2 %in% c("70-79") & vax1_date < study_dates$in70s$start_date ~ FALSE,
        # ignore agebands under 70 years old as these are not being studied here
        TRUE ~ TRUE 
      ),
      
      # stage-specific variables defined above
      !!! selection_stage,
      
      isnot_inhospital = is.na(admitted_unplanned_0_date) | (!is.na(discharged_unplanned_0_date) & discharged_unplanned_0_date < index_date),
      
      # these must match the descriptions in criteria_descr
      c2 = c1 & has_follow_up_previous_year,
      c3 = c2 & isnot_hscworker,
      c4 = c3 & isnot_carehomeresident & isnot_housebound,
      c5 = c4 & isnot_endoflife,
      c6 = c5 & has_age & has_sex & has_imd & has_ethnicity & has_region,
      c7 = c6 & !prior_covid_infection,
      c8 = c7 & isnot_inhospital,
      
      # stage-dependent, `include` defined above
      include = !! sym(include),
      
    )
  
  data_eligible <- data_criteria %>%
    filter(include) %>%
    select(patient_id) %>%
    left_join(data_processed, by="patient_id") %>%
    droplevels()
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# save cohort-specific datasets ----
if (stage == "single") {
  
  my_skim(
    data_eligible,
    path = ghere("output", "single", "eligible", "data_eligible_skim.txt")
    )
  write_rds(
    data_eligible, 
    here("output", "single", "eligible", "data_singleeligible.rds"),
    compress = "gz"
    )
  # also save as csv.gz for reading into study_definition
  write_csv(
    data_eligible %>% distinct(patient_id),  
    here("output", "single", "eligible", "data_singleeligible.csv.gz")
    )
  
} else if (stage == "treated") {
  
  data_eligible %>% 
    filter(vax1_type == "pfizer") %>%
    my_skim(path = here("output", "sequential", "treated", "eligible", "data_eligible_pfizer_skim.txt"))
  
  write_rds(
    data_eligible %>% filter(vax1_type == "pfizer"), 
    here("output", "sequential", "pfizer", "treated", "data_treatedeligible.rds"),
    compress="gz"
    )
  
  data_eligible %>% 
    filter(vax1_type == "az") %>%
    my_skim(path = here("output", "sequential", "treated", "eligible", "data_eligible_az_skim.txt"))
  
  write_rds(
    data_eligible %>% filter(vax1_type == "az"), 
    here("output", "sequential", "az", "treated", "data_treatedeligible.rds"), 
    compress="gz"
    )
  
} else if (stage == "potential") {
  
  my_skim(
    data_eligible, 
    path = ghere("output", "sequential", brand, "matchround{matching_round}", "process", "data_controlpotential_skim.txt")
    )
  
  write_rds(
    data_eligible, 
    ghere("output", "sequential", brand, "matchround{matching_round}", "process", "data_controlpotential.rds"),
    compress = "gz"
    )
  
  if (matching_round == 1) {
      write_csv(
        data_eligible %>% distinct(patient_id),
        ghere("output", "sequential", brand, "matchround{matching_round}", "process", "data_controlpotential.csv.gz")
        )
  }
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# create flowchart (only when stage="treated" or "single") ----
if (stage %in% c("single", "treated")) {
  
  create_flowchart <- function(brand="any") {
    
    if (brand == "any") {
      data_flowchart <- data_criteria 
    } else {
      data_flowchart <- data_criteria %>%
        filter(vax1_type == brand)
    }
    
    data_flowchart <- data_flowchart %>%
      summarise(
        across(unname(criteria_descr), .fns=sum)
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
        criteria = fct_recoderelevel(crit, criteria_descr)
      ) %>%
      mutate(across(starts_with("pct_"), ~round(.x,4)))
    
    if (stage == "single") {
      flowchart_path <- here("output", "single", "eligible")
    } else if (stage == "treated") {
      flowchart_path <- here("output", "sequential", "treated", "eligible")
    }
    
    write_rds(
      data_flowchart, 
      file.path(flowchart_path, glue("flowchart_{stage}eligible_{brand}_unrounded.rds"))
      )
    
    # apply rounding
    data_flowchart %>%
      transmute(
        criteria, crit, 
        n = ceiling_any(n, to=threshold),
        n_exclude = lag(n) - n,
        pct_exclude = n_exclude/lag(n),
        pct_all = n / first(n),
        pct_step = n / lag(n),
      ) %>%
      mutate(across(starts_with("pct_"), ~round(.x, 4))) %>%
      write_csv(file.path(flowchart_path, glue("flowchart_{stage}eligible_{brand}_rounded.csv"))) 
    
  }
  
  # create flowchart when stage=single, don't split brands
  if (stage == "single") {
    create_flowchart("any")
  }
  # create flowchart when stage=treated, split brands
  if (stage == "treated") {
    create_flowchart("pfizer")
    create_flowchart("az")
  }
  
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# check matches made when stage=potiential still hold when stage=actual ----
if (stage == "actual") { 
  
  data_control <- data_eligible
  
  data_treated <- 
    left_join(
      data_potential_matchstatus %>% filter(treated==1L),
      read_rds(ghere("output", "sequential", brand, "treated", "data_treatedeligible.rds")) %>% 
        # only keep variables that are in data_control (this gets rid of outcomes and vax4 dates)
        select(any_of(names(data_control))),
      by="patient_id"
    )
  
  matching_candidates <- 
    bind_rows(data_treated, data_control) %>%
    arrange(treated, match_id, trial_date)
  
  # print missing values
  matching_candidates_missing <- map(matching_candidates, ~any(is.na(.x)))
  sort(names(matching_candidates_missing[unlist(matching_candidates_missing)]))
  
  # rematch
  # first join on exact variables + match_id + trial_date
  rematch <-
    inner_join(
      x=data_treated %>% select(match_id, trial_date, all_of(c(names(caliper_variables), exact_variables))),
      y=data_control %>% select(match_id, trial_date, all_of(c(names(caliper_variables), exact_variables))),
      by = c("match_id", "trial_date", exact_variables)
    ) 
  
  # check caliper_variables are still within caliper
  if (length(caliper_variables) > 0) {
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
        vars(names(caliper_variables)),
        all_vars(.)
      )
  } 
  
  rematch <- rematch %>%
    select(match_id, trial_date) %>%
    mutate(matched=1)
  # only keep individuals who were successfully rematched
  data_successful_match <-
    matching_candidates %>%
    inner_join(rematch, by=c("match_id", "trial_date", "matched")) %>%
    mutate(matching_round = matching_round) %>%
    arrange(trial_date, match_id, treated)
  
  # move match status cols to front of dataset
  matchstatus_vars <- c("patient_id", "match_id", "trial_date", "matching_round", "treated", "controlistreated_date")
  data_successful_matchstatus <- 
    data_successful_match %>% 
    # keep all variables from the processed data as they are required for adjustments in the cox model
    select(all_of(matchstatus_vars), everything())
  
  ## size of dataset
  print("data_successful_match treated/untreated numbers")
  table(treated = data_successful_matchstatus$treated, useNA="ifany")
  
  ## how many matches are lost?
  print(glue(
    "{sum(data_successful_matchstatus$treated)} matched-pairs kept out of {sum(data_potential_matchstatus$treated)} 
    ({round(100*(sum(data_successful_matchstatus$treated) / sum(data_potential_matchstatus$treated)),2)}%)
    "
  ))
  
  # pick up all previous successful matches
  if (matching_round > 1) {
    
    data_matchstatusprevious <- 
      read_rds(ghere("output", "sequential", brand, "matchround{matching_round-1}", "actual", "data_matchstatus_allrounds.rds"))
    
    data_matchstatus_allrounds <- 
      data_successful_matchstatus %>% 
      select(all_of(matchstatus_vars)) %>%
      bind_rows(data_matchstatusprevious) 
    
  } else {
    
    data_matchstatus_allrounds <- 
      data_successful_matchstatus %>%
      select(all_of(matchstatus_vars))
    
  }
  
  write_rds(
    data_matchstatus_allrounds, 
    ghere("output", "sequential", brand, "matchround{matching_round}", "actual", "data_matchstatus_allrounds.rds"),
    compress="gz"
    )
  
  
  # output all control patient ids for finalmatched study definition
  data_matchstatus_allrounds %>%
    # to preserve date format when read into study definition
    mutate(trial_date=as.character(trial_date)) %>%
    # only interested in controls
    filter(treated==0L) %>%
    write_csv(
      ghere("output", "sequential", brand, "matchround{matching_round}", "actual", "cumulative_matchedcontrols.csv.gz")
      )
  
  # size of dataset
  print("data_matchstatus_allrounds treated/untreated numbers")
  table(treated = data_matchstatus_allrounds$treated, useNA="ifany")
  
  # check for duplicate IDs
  print("check for duplicate ids")
  data_matchstatus_allrounds %>% 
    group_by(treated, patient_id) %>%
    summarise(n=n()) %>% 
    group_by(treated) %>% 
    summarise(ndups = sum(n>1)) %>%
    print()
  
  my_skim(
    data_eligible, 
    path = ghere("output", "sequential", brand, "matchround{matching_round}", "actual", "data_successful_matchedcontrols_skim.txt")
    )
  write_rds(
    data_successful_matchstatus %>% filter(treated==0L),
    ghere("output", "sequential", brand, "matchround{matching_round}", "actual", "data_successful_matchedcontrols.rds"), 
    compress="gz"
    )
  
  # size of dataset
  print("data_successful_match treated/untreated numbers")
  table(treated = data_successful_matchstatus$treated, useNA="ifany")
  
}

