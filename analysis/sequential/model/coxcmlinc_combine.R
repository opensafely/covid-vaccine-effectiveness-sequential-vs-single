
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Combine km estimates from different outcomes
#  - The script must be accompanied by three arguments:
#    `cohort` - the cohort used
#    `subgroup` - the subgroup variable
#    `outcome` - the outcome variable
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')

## Import custom user functions from lib
source(here("lib", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "pfizer"
} else {
  cohort <- args[[1]]
}



output_dir <- ghere("output", cohort, "models", "coxcmlinc", "combined")
fs::dir_create(output_dir)


metaparams <-
  expand_grid(
    outcome = factor(c("postest", "emergency", "covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath", "death")),
    #outcome = factor(c("postest", "covidadmitted")),
    subgroup = factor(recoder$subgroups),
  ) %>%
  mutate(
    #subgroup_level = map(as.character(subgroup), ~unname(recoder[[.x]])),
    outcome_descr = fct_recoderelevel(outcome,  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroups),
    #subgroup_level_descr = map(as.character(subgroup), ~names(recoder[[.x]])),
  )


# combine cox estimates ----
coxcmlinc_cuts <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "coxcmlinc", subgroup, outcome, "coxcmlinc_cuts.rds"))
      dat %>%
        add_column(
          subgroup_level = as.character(.[[subgroup]]),
          subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
          .before=1
        ) %>%
        select(-all_of(subgroup))
    }
    )
  ) %>%
  unnest(data)

write_csv(coxcmlinc_cuts, fs::path(output_dir, "coxcmlinc_cuts.csv"))


coxcmlinc_overall <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "coxcmlinc", subgroup, outcome, "coxcmlinc_overall.rds"))
      dat %>%
        add_column(
          subgroup_level = as.character(.[[subgroup]]),
          subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
          .before=1
        ) %>%
        select(-all_of(subgroup))
    }
    )
  ) %>%
  unnest(data)

write_csv(coxcmlinc_overall, fs::path(output_dir, "coxcmlinc_overall.csv"))

