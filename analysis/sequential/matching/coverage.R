# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: describe matching results
# imports matching data
# reports on matching coverage
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

# Import libraries 
library('tidyverse')
library('lubridate')
library('here')
library('glue')
library('arrow')

# import local functions and parameters 
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))

# import command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  # use for interactive testing
  brand <- "pfizer"
} else {
  brand <- args[[1]]
}

# get brand-specific parameters study dates
dates <- study_dates[[brand]]

# create output directory
output_dir <- here("output", "report", "coverage")
fs::dir_create(output_dir)

# matching coverage on each day of recruitment period ----

data_treatedeligible_matchstatus <- read_rds(
  here("output", "sequential", brand, "match", "data_treatedeligible_matchstatus.rds")
  )

# matching coverage for treated people
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
    vax1_date = full_seq(c(dates$start_date, dates$end_date), 1), 
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

# don't compress as this file is being released
write_csv(data_coverage_rounded, fs::path(output_dir, glue("coverage_{brand}.csv")))
