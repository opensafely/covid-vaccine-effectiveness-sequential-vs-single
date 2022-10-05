# create dummy data for treated and potential control population ----

library('tidyverse')
library('arrow')
library('here')
library('glue')

source(here("lib", "functions", "utility.R"))

# remotes::install_github("https://github.com/wjchulme/dd4d")
library('dd4d')
source('/Users/eh1415/Documents/dd4d/R/bn_simulate.R')

source(here("analysis", "design.R"))

source(here("analysis", "dummy", "sim_lst.R"))

population_size <- 20000

# import globally defined repo variables from
study_dates <- jsonlite::read_json(
  path=here("lib", "design", "study-dates.json")
)

# create dummy data for variables defined before baseline ----
pfizerstart_date <- as.Date(study_dates$pfizer$start_date)
modernastart_date <- as.Date(study_dates$moderna$start_date)

firstpfizer_date <- as.Date(study_dates$firstpfizer_date)
firstaz_date <- as.Date(study_dates$firstaz_date)
firstmoderna_date <- as.Date(study_dates$firstmoderna_date)

index_date <- as.Date(study_dates$index_date)

index_day <- 0L
pfizerstart_day <- as.integer(pfizerstart_date - index_date)
modernastart_day <- as.integer(modernastart_date - index_date)

firstpfizer_day <- as.integer(firstpfizer_date - index_date)
firstaz_day <- as.integer(firstaz_date - index_date)
firstmoderna_day <- as.integer(firstmoderna_date - index_date)


known_variables <- c(
  "index_date", "pfizerstart_date", "modernastart_date", "firstpfizer_date", "firstaz_date", "firstmoderna_date",
  "index_day", "pfizerstart_day", "modernastart_day", "firstpfizer_day", "firstaz_day", "firstmoderna_day"
)

sim_list <- splice(
  sim_list_vax,
  sim_list_jcvi,
  sim_list_demographic,
  sim_list_pre,
  sim_list_outcome
)

bn <- bn_create(sim_list, known_variables = known_variables)

# bn_plot(bn)
# bn_plot(bn, connected_only = TRUE)

set.seed(10)

dummydata <- bn_simulate(bn, pop_size = population_size, keep_all = FALSE, .id = "patient_id")

# create covid_vax_disease variables, as the dependencies are difficult to specify using bn_node
dummydata_vax <- dummydata %>% 
  select(patient_id, starts_with("covid_vax")) %>% 
  pivot_longer(
    cols = -patient_id,
    values_drop_na = TRUE
  ) %>%
  group_by(patient_id) %>%
  mutate(sequence = rank(value, ties = "random")) %>%
  ungroup() %>%
  filter(sequence<=4) %>%
  select(-name) %>%
  pivot_wider(
    names_from = sequence,
    names_glue = "covid_vax_disease_{sequence}_day",
    values_from = value
  )


dummydata_processed <- dummydata  %>%
  left_join(dummydata_vax, by = "patient_id") %>%
  # convert logical to integer as study defs output 0/1 not TRUE/FALSE
  # mutate(across(where(is.logical), ~ as.integer(.))) %>%
  # re-index outcomes on cavid_vax_disease_3_day
  mutate(across(all_of(names(sim_list_outcome)), ~ covid_vax_disease_3_day + .)) %>%
  # convert integer days to dates since index date and rename vars
  mutate(across(ends_with("_day"), ~ as.Date(as.character(index_date + .)))) %>%
  rename_with(~ str_replace(., "_day", "_date"), ends_with("_day"))


# save dummy data files ----

fs::dir_create(here("lib", "dummydata"))

# dummy_treated
dummydata_processed %>% 
  filter(!is.na(covid_vax_disease_3_date)) %>% 
  write_feather(sink = here("lib", "dummydata", "dummy_treated.feather"))

# dummy_control_potential1 (reused for actual)
dummydata_processed %>% 
  select(-all_of(str_replace(names(sim_list_outcome), "_day", "_date"))) %>%
  select(-matches("covid_vax_\\w+_4_date")) %>%
  write_feather(sink = here("lib", "dummydata", "dummy_control_potential1.feather"))

