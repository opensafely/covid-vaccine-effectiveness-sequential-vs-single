
# # # # # # # # # # # # # # # # # # # # #
# Purpose: To gather level 4 files ("moderately sensitive") place in a single directory for easy review and release
# # # # # # # # # # # # # # # # # # # # #

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')


## Import custom user functions from lib
source(here("lib", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))


## post-matching ----

output_dir <- here("output", "release")
fs::dir_create(output_dir)

for(cohort in c("over12", "under12")){
#for(cohort in c("over12")){

  input_dir <- ghere("output", cohort, "models", "km", "combined")
  
  ## table1 ----

  fs::file_copy(here("output", cohort, "table1", "coverage.csv"), fs::path(output_dir, glue("{cohort}_coverage.csv")), overwrite = TRUE)
  fs::file_copy(here("output", cohort, "table1", "table1.csv"), fs::path(output_dir, glue("{cohort}_table1.csv")), overwrite = TRUE)
  # fs::file_copy(here("output", cohort, "table1", "flowchart.csv"), fs::path(output_dir, glue("{cohort}_flowchart.csv")), overwrite = TRUE)

  ## KM ----
  fs::file_copy(fs::path(input_dir, "km_estimates_rounded.csv"), fs::path(output_dir, glue("{cohort}_km_estimates_rounded.csv")), overwrite = TRUE)
  fs::file_copy(fs::path(input_dir, "contrasts_daily_rounded.csv"), fs::path(output_dir, glue("{cohort}_contrasts_daily_rounded.csv")), overwrite = TRUE)
  fs::file_copy(fs::path(input_dir, "contrasts_cuts_rounded.csv"), fs::path(output_dir, glue("{cohort}_contrasts_cuts_rounded.csv")), overwrite = TRUE)
  fs::file_copy(fs::path(input_dir, "contrasts_overall_rounded.csv"), fs::path(output_dir, glue("{cohort}_contrasts_overall_rounded.csv")), overwrite = TRUE)
}


fs::dir_create(here("output", "meta-release"))

## create text for output review issue ----
fs::dir_ls(output_dir, type="file", recurse=TRUE) %>%
  map_chr(~str_remove(., fixed(here()))) %>%
  map_chr(~paste0("- [ ] ", str_remove(.,fixed("/")))) %>%
  paste(collapse="\n") %>%
  writeLines(here("output", "meta-release", "files-for-release.txt"))


## create command for releasing using osrelease ----
fs::dir_ls(output_dir, type="file", recurse=TRUE) %>%
  map_chr(~str_remove(., fixed(here()))) %>%
  #map_chr(~paste0("'",. ,"'")) %>%
  paste(., collapse=" ") %>%
  paste("osrelease", .) %>%
  writeLines(here("output", "meta-release", "osrelease-command.txt"))

