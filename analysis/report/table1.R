# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: describe baseline characteristics
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# Import libraries
library('tidyverse')
library('here')
library('glue')
library('gt')
library('gtsummary')

# import local functions and parameters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))

# import command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  # use for interactive testing
  approach <- "single"
  brand <- "any"
  # approach <- "sequential"
  # brand <- "pfizer"
} else {
  approach <- args[[1]]
  brand <- args[[2]]
}

# create output directory
outdir <- here("output", "report", "table1")
fs::dir_create(outdir)

# read data
if (approach == "single") {
  data_in <- read_rds(here("output", "single", "eligible", "data_singleeligible.rds")) %>%
    mutate(treated = "single")
} else {
  data_in <- read_rds(ghere("output", "sequential", brand, "match", "data_matched.rds")) 
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# create table 1 ----
# set variable names as factor levels
var_levels <- map_chr(
  var_lookup[-which(names(var_lookup) %in% c("N", "treated"))],
  ~last(as.character(.))
  )

# use gtsummary to obtain standardised table 1 data
tab_summary_baseline <-
  data_in %>%
  mutate(
    N = 1L,
    age_factor = factor(age, levels=sort(unique(age)))
  ) %>%
  select(
    treated, age_factor,
    any_of(names(var_lookup)),
  ) %>%
  tbl_summary(
    by = treated,
    statistic = list(N = "{N}")
  )

raw_stats <- tab_summary_baseline$meta_data %>%
  select(var_label, df_stats) %>%
  unnest(df_stats)

raw_stats_rounded <- raw_stats %>%
  mutate(
    n=roundmid_any(n, threshold),
    N=roundmid_any(N, threshold),
    p=n/N,
    N_miss = roundmid_any(N_miss, threshold),
    N_obs = roundmid_any(N_obs, threshold),
    p_miss = N_miss/N_obs,
    N_nonmiss = roundmid_any(N_nonmiss, threshold),
    p_nonmiss = N_nonmiss/N_obs,
    median, p25, p75,
    var_label = factor(var_label, levels=var_levels),
    variable_levels = replace_na(as.character(variable_levels), "")
  ) 

write_csv(raw_stats_rounded, fs::path(outdir, glue("table1_{approach}_{brand}_rounded.csv")))
