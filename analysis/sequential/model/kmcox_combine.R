
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
source(here("analysis", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))

## define input and output directories
outdir <- here("output", "sequential", "combine")
fs::dir_create(outdir)

# define metaparams
metaparams <-
  expand_grid(
    cohort = model_brands,
    outcome = model_outcomes,
    subgroup = model_subgroups
  )  

# define function for combining files
combine_files <- function(filename) {
  
  metaparams %>%
    mutate(
      data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome) {
        subgroup <- as.character(subgroup)
        dat <- read_rds(here("output", "sequential", cohort, "model", subgroup, outcome, glue("{filename}.rds")))
        dat %>%
          ungroup() %>%
          add_column(
            subgroup_level = as.character(.[[subgroup]]),
            .before=1
          ) %>%
          select(-all_of(subgroup))
      })
    ) %>%
    unnest(data) %>%
    write_csv(file.path(outdir, glue("{filename}.csv")))
  
}

combine_files("km_estimates_rounded")
combine_files("contrasts_km_cuts_rounded")
combine_files("contrasts_cox_cuts")

## move km plots to single folder ----
move_plots <- function(filename) {
  
  metaparams %>%
    rowwise() %>%
    mutate(
      plotdir = here("output", "sequential", cohort, "model", subgroup, outcome, glue("{filename}.png")),
      plotnewdir = file.path(outdir, glue("{filename}_{cohort}_{subgroup}_{outcome}.png")),
    ) %>%
    {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}
  
  
}

move_plots("km_plot_rounded")
move_plots("km_plot_unrounded")
