library(tidyverse)
library(glue)

# load and define functions
source(here::here("lib", "functions", "utility.R"))

calculate_flow_stats <- function(.data) {
  .data %>%
    mutate(
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n)
    )
}

# derive data for flowcharts with brands combined and specified
for (brand in c("any", "pfizer", "az")) {
  
  # read data for eligible treated
  treated_eligible <- readr::read_rds(here::here("output", "treated", "eligible", glue("flowchart_treatedeligible_{brand}_unrounded.rds"))) %>%
    mutate(section = "eligibility")
  
  # define function for extracting values for the flow chart
  match_flow_fun <- function(
    brand, 
    n_eligible = treated_eligible[treated_eligible$crit == "c6",][["n"]]
  ) {
    
    # read matching data
    data_matched <- readr::read_rds(here::here("output", brand, "match", "data_matched.rds"))
    
    # unique ids included
    included_ids <- data_matched %>%
      distinct(patient_id)
    
    # unique treated ids
    treated_matched <- data_matched %>%
      filter(treated==1) %>%
      distinct(patient_id)
    
    # unique control ids
    control_matched <- data_matched %>%
      filter(treated==0) %>%
      distinct(patient_id)
    
    # control-before-treated ids
    control_before_treated <- treated_matched %>%
      inner_join(control_matched, by = "patient_id")
    
    # fill in data for flow chart
    out <- tibble(
      criteria = "Eligible", crit=NA_character_, n = n_eligible
    ) %>%
      add_row(
        criteria = glue("  vaccinated with {brand} and matched"),
        crit = glue("c7-{brand}"), 
        n = nrow(treated_matched)
      ) %>%
      add_row(
        criteria = glue("  matched as a control prior to vaccination with {brand}"),
        crit = glue("c8-{brand}"),
        n = nrow(control_before_treated)
      ) %>%
      mutate(across(starts_with("crit"), as.character)) %>%
      calculate_flow_stats() %>%
      mutate(section = brand)
    
    return(list(included_ids = included_ids, flow = out))
    
  }
  
  # apply the function for each brand
  if (brand == "any") {
    
    match_flow_pfizer <- match_flow_fun("pfizer")
    match_flow_az <- match_flow_fun("az")
    
    included_ids <- bind_rows(
      match_flow_pfizer$included_ids,
      match_flow_az$included_ids
    )
    
    flow <- bind_rows(
      match_flow_pfizer$flow,
      match_flow_az$flow
    )
    
  } else {
    
    match_flow_brand <- match_flow_fun(brand)
    included_ids <- match_flow_brand$included_ids
    flow <- match_flow_brand$flow
    
  }
  
  # calculate total unique ids inluded across both brands
  total_included <- included_ids %>% 
    distinct() %>% 
    nrow()
  
  # bind to treated_eligible
  flow_unrounded <- treated_eligible %>%
    bind_rows(flow) %>%
    add_row(
      criteria = "Total unique ids included", n = total_included, section = "total"
    ) 
  
  flow_unrounded %>%
    print(n=Inf)
  
  flow_rounded <- flow_unrounded %>%
    mutate(across(n, ceiling_any, to = 7)) %>%
    group_by(section) %>%
    calculate_flow_stats() %>%
    ungroup()
  
  flow_rounded %>%
    print(n=Inf)
  
}




