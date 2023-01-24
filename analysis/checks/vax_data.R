library(tidyverse)
library(here)
library(glue)

# create output directory
outdir <- here("output", "checks")
fs::dir_create(outdir)

process <- function(.data) {
  .data %>%
    select(patient_id, starts_with("covid_vax_disease")) %>%
    mutate(across(ends_with("_date"), as.Date))
    
}

extract_treated <- arrow::read_feather(here("output", "sequential", "treated", "extract", "input_treated.feather")) %>%
  process() 

extract_controlpotential <- arrow::read_feather(here("output", "sequential", "pfizer", "matchround1", "extract", "input_controlpotential.feather")) %>%
  process() 

capture.output(
  inner_join(
    extract_treated,
    extract_controlpotential,
    by = "patient_id",
    suffix = c("_treated", "_controlpotential")
  ) %>%
    summarise(
      dose_1_mismatch = sum(covid_vax_disease_1_date_treated != covid_vax_disease_1_date_controlpotential, na.rm = TRUE),
      dose_2_mismatch = sum(covid_vax_disease_2_date_treated != covid_vax_disease_2_date_controlpotential, na.rm = TRUE)
    ) %>%
    print(),
  file = file.path(outdir, "vax_data.txt")
)

