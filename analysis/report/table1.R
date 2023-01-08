# table 1 style baseline characteristics ----

# # # # # # # # # # # # # # # # # # # # #
# Purpose: describe baseline characteristics
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
# library('lubridate')
library('here')
library('glue')
# library('arrow')
library('gt')
library('gtsummary')

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
  # cohort <- "az"
  # cohort <- "single"
} else {
  #FIXME replace with actual eventual action variables
  removeobjects <- TRUE
  cohort <- args[[1]]
}

# create output directory ----
outdir <- here("output", "report", "table1")
fs::dir_create(outdir)

# read in data ----
if (cohort == "single") {
  data_in <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds")) %>%
    mutate(treated = "single")
} else {
  data_in <- read_rds(ghere("output", "sequential", cohort, "match", "data_matched.rds")) 
}

# set variable names
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
  sev_mental ~ "Serious mental illness",
  flu_vaccine ~ "Influenza vaccination in previous 5 years"
  
) %>%
  set_names(., map_chr(., all.vars))

map_chr(var_labels[-c(1,2)], ~last(as.character(.)))


# use gtsummary to obtain stnadardised table 1 data
tab_summary_baseline <-
  data_in %>%
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

write_csv(raw_stats_redacted, fs::path(outdir, glue("table1_{cohort}.csv")))
