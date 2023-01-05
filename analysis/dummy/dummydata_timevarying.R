# create dummy data for timevarying and outcome variables in single trial approach ----

library('tidyverse')
library('arrow')
library('here')
library('glue')

fs::dir_create(here("output", "single", "dummydata"))

# no indents to make it easier to compare diff
if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("", "expectations")){ 
  
  source(here("analysis", "functions", "utility.R"))
  
  source(here("analysis", "design.R"))
  
  # import all datasets of matched controls, including matching variables
  data_singleeligible <- 
    read_rds(here("output", "single", "eligible", "data_singleeligible.rds")) %>%
    select(
      patient_id, index_date,
      admitted_unplanned_0_date, discharged_unplanned_0_date
      ) %>%
    mutate(across(ends_with("_date"), ~as.integer(.x - study_dates$global$index_date))) %>%
    rename_with(.fn = ~ str_replace(.x, "_date", "_day"), .cols = ends_with("_date"))
    
    
  
  missing <- function(x, rate){
    missing_index <- seq_len(length(x))[rbinom(length(x), 1, rate)==1]
    x[missing_index] <- NA
    x
  }
  
  set.seed(10)
  
  dummydata <- data_singleeligible %>%
    mutate(
      dereg_day = missing(as.integer(runif(n=n(), index_day, index_day+120)), 0.99),
      postest_day = missing(as.integer(runif(n=n(), index_day, index_day+100)), 0.7),
      covidadmitted_day = missing(as.integer(runif(n=n(), index_day, index_day+100)), 0.7),
      death_day = missing(as.integer(runif(n=n(), index_day, index_day+100)), 0.9),
      coviddeath_day = missing(death_day, 0.7),
      admitted_unplanned_infectious_0_day = missing(as.integer(runif(n=n(), index_day, index_day+100)), 0.9),
      discharged_unplanned_infectious_0_day = missing(as.integer(runif(n=n(), admitted_unplanned_infectious_0_day+1, admitted_unplanned_infectious_0_day+30)), 0.5)
    )
  
  # add the hospitalisation vairables
  for (x in c("", "_infectious")) {
    for (i in 1:6) { # check in variables_timevarying.py to see what n is set to
      dummydata <- dummydata %>%
        mutate(
          !! sym(glue("admitted_unplanned{x}_{i}_day")) := pmax(
            missing(as.integer(runif(n=n(), !! sym(glue("discharged_unplanned{x}_{i-1}_day"))+1, !! sym(glue("discharged_unplanned{x}_{i-1}_day"))+30)), 0.5), 0 
            ),
          !! sym(glue("discharged_unplanned{x}_{i}_day")) :=  pmax(
            missing(as.integer(runif(n=n(), !! sym(glue("admitted_unplanned{x}_{i}_day"))+1, !! sym(glue("admitted_unplanned{x}_{i}_day"))+100)), 0.5), 0
          )
        )
    }
  }
  
  dummydata %>%
    select(-index_day, -admitted_unplanned_0_day, -discharged_unplanned_0_day) %>%
    #convert integer days to dates since index date and rename vars
    mutate(across(ends_with("_day"), ~ as.Date(as.character(study_dates$global$index_date + .x)))) %>%
    rename_with(~str_replace(., "_day", "_date"), ends_with("_day")) %>%
    write_feather(sink = here("output", "single", "dummydata", "dummy_timevarying.feather"))
  
  
} else {
  
  # save empty output to save space if running on real data
  tibble() %>%
    write_feather(sink = here("output", "single", "dummydata", "dummy_timevarying.feather"))
  
}
