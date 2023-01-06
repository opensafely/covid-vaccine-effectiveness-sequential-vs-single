library(tidyverse)
library(glue)
library(here)

# load and define functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

calculate_flow_stats <- function(.data) {
  .data %>%
    mutate(
      n_exclude = lag(n) - n,
      pct_exclude = n_exclude/lag(n),
      pct_all = n / first(n),
      pct_step = n / lag(n)
    )
}

# create output directories
outdir <- here("output", "report", "flowchart")
fs::dir_create(outdir)

data_singleeligible <- readr::read_rds(here("output", "single", "eligible", "data_singleeligible.rds")) %>%
  select(patient_id, vax1_date)

flow_categories <- tribble(
  ~crit, ~criteria,
  # those who are vaccinated on day 1 of recruitment
  "A", "vaccinated and unmatched",
  "B", "vaccinated and matched",
  # those who are vaccinated during the recruitment period but not on day 1
  "C", "unvaccinated and unmatched then vaccinated and unmatched",
  "D", "unvaccinated and unmatched then vaccinated and matched",
  "E", "unvaccinated and matched then vaccinated and unmatched",
  "F", "unvaccinated and matched then vaccinated and matched",
  # those who remain unvaccinated at the end of the recruitment period
  "G", "unvaccinated and unmatched",
  "H", "unvaccinated and matched"
) 

# TODO use fuzzy join to match on crits????
flow_boxes <- tribble(
  ~crits, ~box,
  "ABCDEF", "Vaccinated with {brand} during recruitment period",
  "AC", "Vaccinated with {brand} during recruitment period, unmatched as treated, unmatched as control",
  "BDF", "Vaccinated with {brand} during recruitment period, matched as treated",
  "E", "Vaccinated with {brand} during recruitment period, matched as treated, matched as control",
  "F", "Vaccinated with {brand} during recruitment period, matched as treated, matched as control",
  "EFH", "Unvaccinated, matched as control in {brand} trial",
  "GH", "Unvaccinated up to recruitment end",
  "G", "Unvaccinated up to recruitment end, unmatched as control",
  "H", "Unvaccinated up to recruitment end, matched as control in {brand} trial"
)

# derive data for flowcharts with brands combined and specified
flowchart_matching_function <- function(brand) {
  
  # read data for eligible treated
  treated_eligible <- readr::read_rds(here("output", "sequential", "treated", "eligible", glue("flowchart_treatedeligible_{brand}_unrounded.rds"))) 
  assign(
    glue("flowchart_treatedeligible_{brand}_unrounded"),
    treated_eligible
  )
  
  # reshape so one row per patient, and logical columns to indicate if matched as treated, control or both
  data_matched <- readr::read_rds(here("output", "sequential", brand, "match", "data_matched.rds")) %>%
    select(patient_id, treated) %>%
    mutate(matched = 1) %>%
    pivot_wider(
      names_from = treated,
      values_from = matched
    ) %>%
    rename("treated" = "1", "control" = "0") 
  
  # categorise individuals
  data_match_flow  <- data_singleeligible %>%
    left_join(data_matched, by = "patient_id") %>%
    mutate(across(c(treated, control), ~ if_else(is.na(.x), FALSE, TRUE))) %>%
    mutate(
      crit = case_when(
        vax1_date == study_dates[[brand]]$start_date & !treated & !control ~ "A",
        vax1_date == study_dates[[brand]]$start_date & treated & !control ~ "B",
        vax1_date <= study_dates$global$studyend_date & !control & !treated ~ "C",
        vax1_date <= study_dates$global$studyend_date & !control & treated ~ "D",
        vax1_date <= study_dates$global$studyend_date & control & !treated ~ "E",
        vax1_date <= study_dates$global$studyend_date & control & treated ~ "F",
        (is.na(vax1_date) | vax1_date > study_dates$global$studyend_date) & !control & !treated ~ "G",
        (is.na(vax1_date) | vax1_date > study_dates$global$studyend_date) & control & !treated ~ "H"
      )
    )
  
  # count number in each category
  flowchart_matching <- data_match_flow %>% group_by(crit) %>% count() %>%
    right_join(flow_categories, by = "crit") %>%
    arrange(crit) %>% 
    mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
    # mutate(across(criteria, ~str_c(brand, .x, sep = ": ")))
    mutate(stage = brand)
  
  # store patient_ids for those in category G, as there will be some overlap between pfizer and az
  GH_ids <- data_match_flow %>% filter(crit %in% c("G", "H")) %>% select(patient_id, crit)
  
  return(
    list(
      flowchart_matching = flowchart_matching,
      GH_ids = GH_ids
    )
  )
  
}

flowchart_matching_pfizer <- flowchart_matching_function("pfizer")
flowchart_matching_az <- flowchart_matching_function("az")

# distinct patients meeting criteria G
n_G <- bind_rows(flowchart_matching_pfizer$GH_ids, flowchart_matching_az$GH_ids) %>% 
  filter(crit == "G") %>%
  distinct(patient_id) %>% 
  nrow()

n_H <- bind_rows(flowchart_matching_pfizer$GH_ids, flowchart_matching_az$GH_ids) %>% 
  filter(crit == "H") %>%
  distinct(patient_id) %>% 
  nrow()

flowchart_matching_pfizer <- flowchart_matching_pfizer$flowchart_matching
flowchart_matching_az <- flowchart_matching_az$flowchart_matching

flowchart_singleeligible_unrounded <- readr::read_rds(here("output", "single", "eligible", "flowchart_singleeligible_any_unrounded.rds"))
  
flowchart_full <- bind_rows(
  flowchart_singleeligible_unrounded %>% mutate(stage = "single"),
  flowchart_matching_pfizer %>% filter(crit != "G"),
  flowchart_matching_az %>% filter(crit != "G")
) %>%
  add_row(
    crit = "G",
    criteria = "unvaccinated and unmatched",
    stage = "pfizer and az",
    n = n_G
  ) %>%
  add_row(
    crit = "H",
    criteria = "unvaccinated and matched",
    stage = "pfizer and az",
    n = n_H
  )



# flowchart_treatedeligible_pfizer_unrounded
# 
# flowchart_treatedeligible_az_unrounded
# 
# flowchart_matching_pfizer_unrounded
# 
# flowchart_matching_az_unrounded
# 
# flowchart_matching_any_unrounded

# flowchart:
c_i_0 <- flowchart_singleeligible_unrounded %>% filter(crit=="c0") %>% select(criteria, n)
c_e <- flowchart_singleeligible_unrounded %>% filter(crit%in%c("c1", "c2", "c3", "c4", "c5")) %>% select(criteria, n, n_exclude)
c_i_5 <- flowchart_singleeligible_unrounded %>% filter(crit=="c5") %>% transmute(criteria = "Included in single trial", n)
c_i_5_pf <- flowchart_treatedeligible_pfizer_unrounded %>% filter(crit=="c5") %>% transmute(criteria = "Vaccinated with pfizer during recruitment period", n)
c_i_5_az <- flowchart_treatedeligible_az_unrounded %>% filter(crit=="c5") %>% transmute(criteria = "Vaccinated with az during recruitment period", n)

c_i_7_pf <- flowchart_matching_pfizer_unrounded %>% filter(crit=="c7-pfizer") %>% transmute(criteria = "Treated and matched (pfizer)", n)
c_i_7_az <- flowchart_matching_az_unrounded %>% filter(crit=="c7-az") %>% transmute(criteria = "Treated and matched (az)", n)

c_i_8_pf <- flowchart_matching_pfizer_unrounded %>% filter(crit=="c8-pfizer") %>% transmute(criteria = "Control then treated (pfizer)", n)
c_i_8_az <- flowchart_matching_az_unrounded %>% filter(crit=="c8-az") %>% transmute(criteria = "Control then treated (az)", n)


tribble(
  ~id, ~name, ~description,
  # those who are vaccinated on day 1 of recruitment
  "A", "vax_unmatched", "vaccinated and unmatched",
  "B", "vax_matched", "vaccinated and matched",
  # those who are vaccinated during the recruitment period but not on day 1
  "C", "unvax_unmatched_then_vax_unmatched", "unvaccinated and unmatched then vaccinated and unmatched",
  "D", "unvax_unmatched_then_vax_matched", "unvaccinated and unmatched then vaccinated and matched",
  "E", "unvax_matched_then_vax_unmatched", "unvaccinated and matched then vaccinated and unmatched",
  "F", "unvax_matched_then_vax_matched", "unvaccinated and matched then vaccinated and matched",
  # those who remain unvaccinated at the end of the recruitment period
  "G", "unvax_unmatched", "unvaccinated and unmatched",
  "H", "unvax_matched", "unvaccinated and matched"
)


# TODO work out why total_unique_ids so high???
total_unique_ids <- flowchart_matching_any_unrounded %>% filter(section == "total") %>% pull(n)

A_un <- pull(c_i_5, n) - pull(c_i_5_pf, n) - pull(c_i_5_az, n)
B_pf <- pull(c_i_5_pf, n) - pull(c_i_7_pf, n)
B_az <- pull(c_i_5_az, n) - pull(c_i_7_az, n)
# B_un <- A_un - (total_unique_ids - pull(c_i_7_pf, n) - pull(c_i_7_az, n))
C_pf <- pull(c_i_7_pf, n) - pull(c_i_8_pf, n)
C_az <- pull(c_i_7_az, n) - pull(c_i_8_az, n)

bind_rows(
  c_i_0,
  c_e,
  c_i_5,
  c_i_5_pf,
  c_i_5_az
) %>%
  add_row(
    criteria = "Unvaccinated on final date of recruitment period",
    n = A_un
  ) %>%
  bind_rows(
    c_i_7_pf %>% mutate(n_exclude = B_pf),
    c_i_7_az %>% mutate(n_exclude = B_az),
    c_i_8_pf,
    c_i_8_az
  ) %>%
  add_row(
    criteria = "Unvaccinated and matched (pfizer)",
    n = C_pf
  ) %>%
  add_row(
    criteria = "Unvaccinated and matched (az)",
    n = C_az
  )# %>%
  # add_row(
  #   criteria = "Unvaccinated and unmatched",
  #   n_exclude = B_un
  # )
  




