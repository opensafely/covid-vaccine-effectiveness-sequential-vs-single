# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Purpose: Combine km estimates from different outcomes
#  - The script must be accompanied by three arguments:
#    `brand` - the brand used
#    `subgroup` - the subgroup variable
#    `outcome` - the outcome variable
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

# import libraries
library('tidyverse')
library('here')
library('glue')
library('survival')

# import custom user functions and paramters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# create output directory
outdir <- here("output", "sequential", "combine")
fs::dir_create(outdir)

# define metaparams
metaparams <-
  expand_grid(
    brand = model_brands,
    outcome = model_outcomes,
    subgroup = model_subgroups
  )  

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# combine files ----

combine_files <- function(filename) {
  
  metaparams %>%
    mutate(
      data = pmap(list(brand, subgroup, outcome), function(brand, subgroup, outcome) {
        subgroup <- as.character(subgroup)
        dat <- read_rds(here("output", "sequential", brand, "model", subgroup, outcome, glue("{filename}.rds")))
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
combine_files("km_estimates_unrounded")
combine_files("contrasts_km_cuts_rounded")
combine_files("contrasts_cox_cuts")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## move km plots to single folder ----

move_plots <- function(filename) {
  
  metaparams %>%
    rowwise() %>%
    mutate(
      plotdir = here("output", "sequential", brand, "model", subgroup, outcome, glue("{filename}.png")),
      plotnewdir = file.path(outdir, glue("{filename}_{brand}_{subgroup}_{outcome}.png")),
    ) %>%
    {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}
  
  
}

move_plots("km_plot_rounded")
move_plots("km_plot_unrounded")
