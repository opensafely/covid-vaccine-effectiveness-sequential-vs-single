# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# calculates the counts for the flowchart in the manuscript
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# import libraries
library(tidyverse)
library(glue)
library(here)

# load functions and paramters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "fuzzy_join.R"))

# create output directories
outdir <- here("output", "report", "flowchart")
fs::dir_create(outdir)

# import flowchart for single trial cohort
data_singleeligible <- readr::read_rds(here("output", "single", "eligible", "data_singleeligible.rds")) %>%
  select(patient_id, vax1_date, vax1_type)

# define all flow categories for sequential trial cohorts
flow_categories <- tribble(
  ~crit, ~criteria,
  # those who are vaccinated with the brand on day 1 of recruitment
  "A", "vaccinated and unmatched",
  "B", "vaccinated and matched",
  # those who are vaccinated with the brand during the recruitment period but not on day 1
  "C", "unvaccinated and unmatched then vaccinated and unmatched",
  "D", "unvaccinated and unmatched then vaccinated and matched",
  "E", "unvaccinated and matched then vaccinated and unmatched",
  "F", "unvaccinated and matched then vaccinated and matched",
  # those who remain unvaccinated at the end of the recruitment period
  "G", "unvaccinated and unmatched",
  "H", "unvaccinated and matched",
  # those who are vaccinated with the other brand
  "I", "vaccinated with the other brand and matched as control"
) 

# define boxes for sequential trial flow
flow_boxes <- tribble(
  ~box_crit, ~box_descr,
  "ABCDEF", "Vaccinated with {brand} during recruitment period",
  "AC", "Vaccinated with {brand} during recruitment period, unmatched as treated, unmatched as control",
  "BDF", "Vaccinated with {brand} during recruitment period, matched as treated",
  "E", "Vaccinated with {brand} during recruitment period, matched as treated, matched as control",
  "F", "Vaccinated with {brand} during recruitment period, matched as treated, matched as control",
  "EFHI", "Unvaccinated, matched as control in {brand} trial",
  "GH", "Unvaccinated up to recruitment end",
  "G", "Unvaccinated up to recruitment end, unmatched as control",
  "H", "Unvaccinated up to recruitment end, matched as control in {brand} trial",
  "I", "Vaccinated with the other brand during recruitment period, matched as control"
)

# derive data for flowcharts with brands combined and specified
flowchart_matching_function <- function(brand) {
  
  cat(" \n")
  cat(glue("{brand}: \n"))
  cat(" \n")
  
  # reshape so one row per patient, and logical columns to indicate if matched as treated, control or both
  data_matched <- readr::read_rds(here("output", "sequential", brand, "match", "data_matched.rds")) %>%
    select(patient_id, treated, vax1_date, vax1_type) %>%
    mutate(matched = 1) %>%
    pivot_wider(
      names_from = treated,
      values_from = matched
    ) %>%
    rename("treated" = "1", "control" = "0") %>%
    mutate(sequential = TRUE)
  
  cat("Check there are the same number of treated and control:\n")
  data_matched %>%
    summarise(
      treated = sum(treated, na.rm = TRUE),
      control = sum(control, na.rm = TRUE)
    ) %>%
    print()
  
  # categorise individuals
  data_match_flow  <- data_singleeligible %>%
    select(patient_id, vax1_date, vax1_type) %>%
    mutate(single = TRUE) %>%
    full_join(
      data_matched,
      by = "patient_id",
      suffix = c("", "_sequential")
      ) %>%
    mutate(across(c(treated, control, single, sequential), ~ replace_na(as.logical(.x), replace=FALSE)))  %>%
    mutate(
      crit = case_when(
        # those who are vaccinated on day 1 of recruitment
        vax1_type == brand & vax1_date == study_dates[[brand]]$start_date & !control & !treated ~ "A",
        vax1_type == brand & vax1_date == study_dates[[brand]]$start_date & !control & treated ~ "B",
        # those who are vaccinated during the recruitment period but not on day 1
        vax1_type == brand & vax1_date <= study_dates$global$studyend_date & !control & !treated ~ "C",
        vax1_type == brand & vax1_date <= study_dates$global$studyend_date & !control & treated ~ "D",
        vax1_type == brand & vax1_date <= study_dates$global$studyend_date & control & !treated ~ "E",
        vax1_type == brand & vax1_date <= study_dates$global$studyend_date & control & treated ~ "F",
        # those who remain unvaccinated at the end of the recruitment period
        (is.na(vax1_date) | vax1_date > study_dates$global$studyend_date) & !control & !treated ~ "G",
        (is.na(vax1_date) | vax1_date > study_dates$global$studyend_date) & control & !treated ~ "H",
        # those who are vaccinated with the other brand
        vax1_date <= study_dates$global$studyend_date & control & !treated ~ "I",
        TRUE ~ NA_character_
      )
    )
  
  cat("Check `vax1_date` and `vax_type` match in single and sequential data:\n")
  data_match_flow %>%
    filter(vax1_date != vax1_date_sequential) %>%
    select(vax1_date, vax1_date_sequential, vax1_type, vax1_type_sequential) %>% 
    print()
  data_match_flow %>%
    filter(vax1_type != vax1_type_sequential) %>%
    select(vax1_date, vax1_date_sequential, vax1_type, vax1_type_sequential) %>%
    print()
  
  # check NAs (capture output as too wide to print in log)
  capture.output(
    data_match_flow %>%
      filter(is.na(crit)) %>%
      mutate(
        vax1_date_cat = case_when(
          is.na(vax1_date) ~ "no vax",
          vax1_type == brand & vax1_date < study_dates[[brand]]$start_date ~ as.character(glue("vax with {brand} before start")),
          vax1_type == brand & vax1_date <= study_dates$global$studyend_date ~ as.character(glue("vax with {brand} between start and end (inclusive)")),
          vax1_type != brand & vax1_date < study_dates[[brand]]$start_date ~ "vax with other brand before start",
          vax1_type != brand & vax1_date <= study_dates$global$studyend_date ~ "vax with other brand between start and end (inclusive)",
          TRUE ~ "vax after end"
        )) %>% 
      group_by(vax1_date_cat, treated, control) %>%
      summarise(
        min_date = min(vax1_date, na.rm = TRUE),
        max_date = max(vax1_date, na.rm = TRUE),
        na_date = sum(is.na(vax1_date)),
        single_only = sum(single & !sequential),
        sequential_only = sum(sequential & !single),
        both = sum(single & sequential),
        total = n(),
        .groups = "keep"
      ) %>%
      ungroup() %>%
      knitr::kable(format = "pipe"),
    file = file.path(outdir, glue("check_NAs_{brand}.txt"))
  )
  
  
  
  # count number in each category
  flowchart_matching <- data_match_flow %>% group_by(crit) %>% count() %>%
    right_join(flow_categories, by = "crit") %>%
    arrange(crit) %>% 
    mutate(across(n, ~if_else(is.na(.x), 0L, .x))) %>%
    mutate(brand = brand)
  
  # store patient_ids for those in category G, as there will be some overlap between pfizer and az
  GH_ids <- data_match_flow %>% filter(crit %in% c("G", "H")) %>% select(patient_id, crit)
  
  return(
    list(
      # data_match_flow = data_match_flow,
      flowchart_matching = flowchart_matching,
      GH_ids = GH_ids
    )
  )
  
}

flowchart_matching_pfizer <- flowchart_matching_function("pfizer")
flowchart_matching_az <- flowchart_matching_function("az")

# distinct patients meeting criteria H for *either* pfizer or az
n_H <- bind_rows(
  flowchart_matching_pfizer$GH_ids, 
  flowchart_matching_az$GH_ids
  ) %>% 
  filter(crit == "H") %>%
  distinct(patient_id) %>% 
  nrow()

# distinct patients meeting criteria G for *both* pfizer and az
n_G <- inner_join(
  flowchart_matching_pfizer$GH_ids, 
  flowchart_matching_az$GH_ids,
  by = c("patient_id", "crit")
  ) %>% 
  filter(crit == "G") %>%
  distinct(patient_id) %>% 
  nrow()


flowchart_matching_pfizer <- flowchart_matching_pfizer$flowchart_matching
flowchart_matching_az <- flowchart_matching_az$flowchart_matching

# single trial flowchart
flowchart_singleeligible_rounded <- readr::read_csv(here("output", "single", "eligible", "flowchart_singleeligible_any_rounded.csv"))
  
# brand-specific
flow_boxes_brand <- flow_boxes %>%
  # get rid of boxes with G as these have duplicates across the brands
  filter(!str_detect(box_crit, "G")) %>%
  # join to the counts for each criteria
  fuzzy_join(
    bind_rows(flowchart_matching_pfizer, flowchart_matching_az), 
    by = c("box_crit" = "crit"), 
    match_fun = str_detect,
    mode = "left"
  ) %>%
  # sum across all criteria in each box
  group_by(brand, box_crit, box_descr) %>%
  summarise(n = sum(n), .groups = "keep")  %>%
  ungroup()

# unvaccinated counts (not brand-specific)
flow_boxes_unvax <- flow_boxes %>%
  # get rid of boxes with G as these have duplicates across the brands
  filter(str_detect(box_crit, "G")) %>%
  fuzzy_join(
    flow_categories %>%
      filter(crit %in% c("G", "H")) %>%
      mutate(
        n = case_when(
          crit == "G" ~ n_G,
          crit == "H" ~ n_H
        )
      ), 
    by = c("box_crit" = "crit"), 
    match_fun = str_detect,
    mode = "left"
  ) %>%
  # sum across all criteria in each box
  group_by(box_crit, box_descr) %>%
  summarise(n = sum(n), .groups = "keep")  %>% 
  ungroup() %>%
  mutate(brand = "unvax")

# bind together and process final flowchart for report
flowchart_final <- bind_rows(
  flowchart_singleeligible_rounded %>% mutate(brand = "single"),
  bind_rows(
    flow_boxes_brand, 
    flow_boxes_unvax 
  ) %>%
    mutate(across(n, ~roundmid_any(n, to = threshold))) %>%
    rename(criteria = box_descr, crit = box_crit) 
) %>%
  rowwise() %>%
  mutate(across(criteria,~glue(.x)))
  
# save flowchart_final  
write_csv(
  flowchart_final,
  file.path(outdir, "flowchart_final_rounded.csv")
)

cat("Check that c7 = ABCDEF_pfizer + ABCDEF_az + GH\n")
cat("c7:\n")
c7 <- flowchart_singleeligible_rounded %>% filter(crit=="c7") %>% pull(n) 
print(c7)
cat("ABCDEF_pfizer:\n")
ABCDEF_pfizer <- flow_boxes_brand %>% filter(box_crit == "ABCDEF", brand=="pfizer") %>% pull(n)
print(ABCDEF_pfizer)
cat("ABCDEF_az:\n")
ABCDEF_az <- flow_boxes_brand %>% filter(box_crit == "ABCDEF", brand=="az") %>% pull(n)
print(ABCDEF_az)
cat("GH:\n")
GH <- flow_boxes_unvax %>% filter(box_crit == "GH") %>% pull(n)
print(GH)

cat("ceiling_any:\n")
c7 == ceiling_any(ABCDEF_pfizer + ABCDEF_az + GH, to=threshold)
cat("roundmid_any:\n")
c7 == roundmid_any(ABCDEF_pfizer + ABCDEF_az + GH, to=threshold)
