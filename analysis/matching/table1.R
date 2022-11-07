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
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "redaction.R"))

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


## create output directories ----

output_dir <- here("output", cohort, "table1")
fs::dir_create(output_dir)

## Import data and derive some variables ----

data_matched0 <- read_rds(ghere("output", cohort, "match", "data_matched.rds")) 

data_treatedeligible_matchstatus <- read_rds(here("output", cohort, "match", "data_treatedeligible_matchstatus.rds"))


## apply additional exclusion criteria ----
## see https://github.com/opensafely/covid-vaccine-effectiveness-seqtrial/pull/10/files#diff-e05cc930e240857c06e4ea714fcabf609e0782ab827e11524be6f9af186dc80c


# define additional criteria in the treated, elgible population
data_treatedeligible_exclusion <- 
  read_rds(ghere("output", cohort, "treated", "data_treatedeligible.rds")) %>%
  transmute(
    patient_id, 
    treatedeligible_nopriorcovid = (
      (is.na(positive_test_0_date) | positive_test_0_date > study_dates[[cohort]][["start_date"]]) &
        (is.na(primary_care_covid_case_0_date) | primary_care_covid_case_0_date > study_dates[[cohort]][["start_date"]]) &
        (is.na(admitted_covid_0_date) | admitted_covid_0_date > study_dates[[cohort]][["start_date"]])
    ),
  )

# define additional criteria in the treated, matched population, ie, also removing matches where the control doesn't meet the criteria
data_matched_preexclusion <- 
  data_matched0 %>%
  group_by(match_id, trial_date, matching_round) %>% 
  mutate(uniquematch_id = cur_group_id()) %>% 
  ungroup() %>%
  group_by(uniquematch_id) %>%
  mutate(
    patient_id, 
    treated,
    nopriorcovid = (
      (is.na(positive_test_0_date) | positive_test_0_date > study_dates[[cohort]][["start_date"]]) &
        (is.na(primary_care_covid_case_0_date) | primary_care_covid_case_0_date > study_dates[[cohort]][["start_date"]]) &
        (is.na(admitted_covid_0_date) | admitted_covid_0_date > study_dates[[cohort]][["start_date"]])
    ),
    nopriorcovid_pair = all(nopriorcovid),
  ) %>%
  ungroup() 

data_treatedmatched_exclusion <- 
  data_matched_preexclusion %>%
  select(patient_id, treated, nopriorcovid, nopriorcovid_pair) %>%
  filter(treated==1L) 

# combine criteria
data_exclusion <- data_treatedeligible_exclusion %>%
  left_join(data_treatedmatched_exclusion, by="patient_id") %>%
  transmute(
    patient_id,
    include_eligible = treatedeligible_nopriorcovid,
    include_matched = (treatedeligible_nopriorcovid & nopriorcovid_pair)
  )

# redefine matching success in the matchstatus dataset
data_treatedeligible_matchstatus <- 
  data_treatedeligible_matchstatus %>%
  left_join(data_exclusion, by="patient_id") %>%
  # remove all eligible treated people originally eligible but no longer eligible due to new criteria
  filter(include_eligible) %>%
  # remove all matched, treated people who are no longer matched because of new inclusion criteria
  mutate(matched = (matched & include_matched)*1L)


data_matched <- data_matched_preexclusion %>%
  filter(nopriorcovid_pair) %>%
  select(-uniquematch_id, -nopriorcovid_pair)


# matching coverage on each day of recruitment period ----


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



## round to nearest 6 for disclosure control
threshold <- 6

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

var_labels <- list(
  N  ~ "Total N",
  treated ~ "Status",
  jcvi_ageband ~ "JCVI ageband",
  sex ~ "Sex",
  #ethnicity_combined ~ "Ethnicity",
  imd_Q5 ~ "Deprivation",
  region ~ "Region",
  
  cev_cv ~ "Clinically vulnerable",
  
  #prior_tests_cat ~ "Number of SARS-CoV-2 tests",
  prior_covid_infection ~ "Prior documented SARS-CoV-2 infection"
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
    n=roundmid_any(n, 6),
    N=roundmid_any(N, 6),
    p=n/N,
    N_miss = roundmid_any(N_miss, 6),
    N_obs = roundmid_any(N_obs, 6),
    p_miss = N_miss/N_obs,
    N_nonmiss = roundmid_any(N_nonmiss, 6),
    p_nonmiss = N_nonmiss/N_obs,
    var_label = factor(var_label, levels=map_chr(var_labels[-c(1,2)], ~last(as.character(.)))),
    variable_levels = replace_na(as.character(variable_levels), "")
  ) 

write_csv(raw_stats_redacted, fs::path(output_dir, "table1.csv"))



# flowchart ----

# data_flowchart_match <-
#   read_rds(here("output", "data", "data_inclusioncriteria.rds")) %>%
#   left_join(
#     data_matchstatus %>% select(patient_id, matched),
#     by="patient_id"
#   ) %>%
#   mutate(
#     c7 = c6 & matched,
#   ) %>%
#   select(-patient_id, -matched) %>%
#   group_by(vax3_type) %>%
#   summarise(
#     across(.fns=sum)
#   ) %>%
#   pivot_longer(
#     cols=-vax3_type,
#     names_to="criteria",
#     values_to="n"
#   ) %>%
#   group_by(vax3_type) %>%
#   mutate(
#     n_exclude = lag(n) - n,
#     pct_exclude = n_exclude/lag(n),
#     pct_all = n / first(n),
#     pct_step = n / lag(n),
#     crit = str_extract(criteria, "^c\\d+"),
#     criteria = fct_case_when(
#       crit == "c0" ~ "Aged 18+ and received booster dose of BNT162b2 or mRNA-1273 between 29 October 2021 and 31 January 2022", # paste0("Aged 18+\n with 2 doses on or before ", format(study_dates$lastvax2_date, "%d %b %Y")),
#       crit == "c1" ~ "  with no missing demographic information",
#       crit == "c2" ~ "  with homologous primary vaccination course of BNT162b2 or ChAdOx1",
#       crit == "c3" ~ "  and not a health and social care worker",
#       crit == "c4" ~ "  and not a care/nursing home resident, end-of-life or housebound",
#       crit == "c5" ~ "  and no COVID-19-related events within 90 days",
#       crit == "c6" ~ "  and not admitted in hospital at time of booster",
#       crit == "c7" ~ "  and successfully matched",
#       TRUE ~ NA_character_
#     )
#   )

## follow-up / outcomes table

source(here("lib", "functions", "survival.R"))


fuptable <- function(data, ...) {
  
  
  data %>%
    mutate(
          # follow-up time is up to and including censor date
      censor_date = pmin(
        dereg_date,
        #vax4_date-1, # -1 because we assume vax occurs at the start of the day
        death_date,
        study_dates[["global"]]$studyend_date,
        trial_date -1 + maxfup,
        na.rm=TRUE
      ),
      matchcensor_date = pmin(censor_date, controlistreated_date -1, na.rm=TRUE), # new censor date based on whether control gets treated or not
      
      tte_postest = tte(trial_date - 1, postest_date, matchcensor_date, na.censor=FALSE),
      tte_covidadmitted = tte(trial_date - 1, covidadmitted_date, matchcensor_date, na.censor=FALSE),
      tte_death = tte(trial_date - 1, death_date, matchcensor_date, na.censor=FALSE),
      
      ind_postest = censor_indicator(postest_date, matchcensor_date),
      ind_covidadmitted = censor_indicator(covidadmitted_date, matchcensor_date),
      ind_death = censor_indicator(death_date, matchcensor_date),
    ) %>%
    group_by( ... ) %>%
    summarise(
      N = roundmid_any(n(), 6),
      
      events_postest = roundmid_any(sum(ind_postest),6),
      events_covidadmitted = roundmid_any(sum(ind_covidadmitted),6),
      events_death = roundmid_any(sum(ind_death),6),
      
      fup_sum_postest = roundmid_any(sum(tte_postest),6),
      fup_sum_covidadmitted = roundmid_any(sum(tte_covidadmitted),6),
      fup_sum_death = roundmid_any(sum(tte_death),6),
      
      fup_mean_postest = fup_sum_postest/N,
      fup_mean_covidadmitted = fup_sum_covidadmitted/N,
      fup_mean_death = fup_sum_death/N,
      
      fup_median_postest = median(tte_postest),
      fup_median_covidadmitted = median(tte_covidadmitted),
      fup_median_death = median(tte_death),
      
      fup_q1_postest = quantile(tte_postest, 0.25),
      fup_q1_covidadmitted = quantile(tte_covidadmitted, 0.25),
      fup_q1_death = quantile(tte_death, 0.25),
      
      fup_q3_postest = quantile(tte_postest, 0.75),
      fup_q3_covidadmitted = quantile(tte_covidadmitted, 0.75),
      fup_q3_death = quantile(tte_death, 0.75),
    )
}


fup_table_all <- data_matched %>% fuptable()
fup_table_treated <- data_matched %>% fuptable(treated)
fup_table_ageband2 <- data_matched %>% fuptable(ageband2)
fup_table_treated_ageband2 <- data_matched %>% fuptable(treated, ageband2)


fup_table <- 
  bind_rows(
    fup_table_all,
    fup_table_treated,
    fup_table_ageband2,
    fup_table_treated_ageband2
  ) %>%
  mutate(
    treated = coalesce(as.character(treated), "all"),
    ageband2 = coalesce(ageband2, "all"),
  ) %>%
  select(
    treated, ageband2, everything()
  )

write_csv(fup_table, fs::path(output_dir, "fup.csv"))
