# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe matching results
# imports matching data
# reports on matching coverage, matching flowcharts, creates a "table 1", etc
# # # # # # # # # # # # # # # # # # # # #


# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('arrow')

## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "pfizer"
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
}

## get cohort-specific parameters study dates and parameters ----
dates <- study_dates[[cohort]]

## create output directory ----

output_dir <- here("output", "sequential", cohort, "match", "report")
fs::dir_create(output_dir)

# matching coverage on each day of recruitment period ----

data_treatedeligible_matchstatus <- read_rds(
  here("output", "sequential", cohort, "match", "data_treatedeligible_matchstatus.rds")
  )

# matching coverage for boosted people
data_coverage <-
  data_treatedeligible_matchstatus %>%
  group_by(vax1_date) %>%
  summarise(
    n_eligible = n(),
    n_matched = sum(matched, na.rm=TRUE),
  ) %>%
  mutate(
    n_unmatched = n_eligible - n_matched,
  ) %>%
  pivot_longer(
    cols = c(n_unmatched, n_matched),
    names_to = "status",
    names_prefix = "n_",
    values_to = "n"
  ) %>%
  arrange(vax1_date, status) %>%
  group_by(vax1_date, status) %>%
  summarise(
    n = sum(n),
  ) %>%
  group_by(status) %>%
  complete(
    vax1_date = full_seq(c(dates$start_date, dates$end_date), 1), # go X days before to
    fill = list(n=0)
  ) %>%
  mutate(
    cumuln = cumsum(n)
  ) %>%
  ungroup() %>%
  mutate(
    status = factor(status, levels=c("unmatched", "matched")),
    status_descr = fct_recoderelevel(status, recoder$status)
  ) %>%
  arrange(status_descr, vax1_date)

data_coverage_rounded <-
  data_coverage %>%
  group_by(status) %>%
  mutate(
    cumuln = roundmid_any(cumuln, to = threshold),
    n = diff(c(0,cumuln)),
  )

write_csv(data_coverage_rounded, fs::path(output_dir, "coverage.csv"))

# table 1 style baseline characteristics ----

library('gt')
library('gtsummary')

data_matched <- read_rds(ghere("output", "sequential", cohort, "match", "data_matched.rds")) 

var_labels <- list(
  N  ~ "Total N",
  treated ~ "Status",
  jcvi_ageband ~ "JCVI ageband",
  sex ~ "Sex",
  ethnicity_combined ~ "Ethnicity",
  imd_Q5 ~ "Deprivation",
  region ~ "Region",
  
  cev_cv ~ "Clinically vulnerable",
  
  sev_obesity ~ "Body Mass Index > 40 kg/m^2",
  chronic_heart_disease ~ "Chronic heart disease",
  chronic_kidney_disease ~ "Chronic kidney disease",
  diabetes ~ "Diabetes",
  chronic_liver_disease ~ "Chronic liver disease",
  chronic_resp_disease ~ "Chronic respiratory disease",
  chronic_neuro_disease ~ "Chronic neurological disease",
  
  multimorb ~ "Morbidity count",
  immunosuppressed ~ "Immunosuppressed",
  learndis ~ "Learning disabilities",
  sev_mental ~ "Serious mental illness"
  
  #prior_tests_cat ~ "Number of SARS-CoV-2 tests",
  #prior_covid_infection ~ "Prior documented SARS-CoV-2 infection"
) %>%
set_names(., map_chr(., all.vars))

map_chr(var_labels[-c(1,2)], ~last(as.character(.)))


# use gtsummary to obtain stnadardised table 1 data
tab_summary_baseline <-
  data_matched %>%
  mutate(
    N = 1L,
    #treated_descr = fct_recoderelevel(as.character(treated), recoder$treated),
    age = factor(age, levels=sort(unique(age)))
  ) %>%
  select(
    treated,
    all_of(names(var_labels)),
  ) %>%
  tbl_summary(
    by = treated,
    label = unname(var_labels[names(.)]),
    statistic = list(N = "{N}")
  ) 

raw_stats <- tab_summary_baseline$meta_data %>%
  select(var_label, df_stats) %>%
  unnest(df_stats)


raw_stats_redacted <- raw_stats %>%
  mutate(
    n=roundmid_any(n, threshold),
    N=roundmid_any(N, threshold),
    p=n/N,
    N_miss = roundmid_any(N_miss, threshold),
    N_obs = roundmid_any(N_obs, threshold),
    p_miss = N_miss/N_obs,
    N_nonmiss = roundmid_any(N_nonmiss, threshold),
    p_nonmiss = N_nonmiss/N_obs,
    var_label = factor(var_label, levels=map_chr(var_labels[-c(1,2)], ~last(as.character(.)))),
    variable_levels = replace_na(as.character(variable_levels), "")
  ) 

write_csv(raw_stats_redacted, fs::path(output_dir, "table1.csv"))
