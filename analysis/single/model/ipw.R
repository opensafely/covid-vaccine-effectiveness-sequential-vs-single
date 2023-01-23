# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# fits marginal structural models (msm) for vaccine effectiveness, with different adjustment sets
# saves model summaries (tables and figures)
# "tte" = "time-to-event"
#
# The script should be run via an action in the project.yaml
# The script must be accompanied by 5 arguments:
# 1. the name of the brand
# 2. the subgroup variable. Use "all" if no subgroups
# 3. the name of the outcome
# 4. the sample size for the vaccination models (a completely random sample of participants)
# 5. the sample size for those who did not experience the outcome for the main MSM models 
#    (all those who did experience an outcome are included)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

# import libraries
library('tidyverse')
library('here')
library('glue')
library('survival')
library('splines')
library('parglm')
library('gtsummary')
library('gt')

# import custom user functions from lib
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))
source(here("analysis", "functions", "survival.R"))
source(here("analysis", "single", "process", "process_data_days_function.R"))
source(here("analysis", "single", "model", "get_ipw_weights_function.R"))

# import command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  removeobs <- FALSE
  brand <- "pfizer"
  subgroup <- "all"
  outcome <- "postest"
} else {
  removeobs <- TRUE
  brand <- args[[1]]
  subgroup <- args[[2]]
  outcome <- args[[3]]
}

# define parglm optimisation parameters
parglmparams <- parglm.control(
  method = "LINPACK",
  nthreads = 8,
  maxit = 40 # default = 25
)


# import formulae 
## these are created in data_define_cohorts.R script
list2env(list_formula_single, globalenv())
formula_1 <- outcome ~ 1
# remove subgroup variable from covariate set
formula_remove_subgroup <- as.formula(paste0(". ~ . - ", subgroup))


# create output directory
outdir <- here("output", "single", brand, subgroup, outcome, "ipw")
fs::dir_create(outdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loop over all subgroup_levels ----

subgroup_levels <- recoder[[subgroup]]

for(subgroup_level in subgroup_levels){
  
  cat("  \n")
  cat(subgroup_level, "  \n")
  cat("  \n")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # data processing ----
  
  # import processed data
  data_fixed <- read_rds(here("output", "single", "stset", "data_fixed.rds")) %>%
    mutate(all = factor("all",levels=c("all"))) %>%
    filter(
      # select patients in current subgroup_level
      .[[subgroup]] == subgroup_level
    )
  
  # generate data_samples
  data_samples <- read_rds(here("output", "single", "stset", "data_patients.rds")) %>%
    # right join because data_fixed filtered on subgroup
    right_join(data_fixed, by="patient_id") %>%
    mutate(
      all = factor("all",levels=c("all")),
      tte_outcome = .[[glue("tte_{outcome}")]]
    ) %>%
    transmute(
      patient_id,
      sample_outcome = sample_nonoutcomes_n(!is.na(tte_outcome), patient_id, msm_sample_nonoutcomes_n),
      sample_weights = sample_weights(!is.na(tte_outcome), sample_outcome),
    )
  
  ## read and process data_days (one row per person day)
  # see analysis/single/process/process_data_days_function.R for the function process_data_days_function
  # (do this _within_ loop so that it can be deleted just before models are run, to reduce RAM use)
  cat("Start `process_data_days_function` for stage=vaccine\n")
  data_days_sub <- bind_rows(
    lapply(
      1:process_data_days_n,
      function(i) 
        process_data_days_function(
          file = "model",
          stage = "vaccine",
          iteration = i
        )
    )
  )
  cat("End `process_data_days_function`\n")
  
  if(removeobs) rm(data_samples, data_fixed)
  
  #print dataset size
  cat(glue("data_days_sub data size = ", nrow(data_days_sub)), "\n  ")
  cat(glue("data_days_sub patient size = ", n_distinct(data_days_sub$patient_id)), "\n  ")
  cat(glue("memory usage = ", format(object.size(data_days_sub), units="GB", standard="SI", digits=3L)), "\n  ")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # IPW model ----
  
  
  # IPW model for pfizer / az vaccination
  # these models are shared across brands (one is treatment model, one is censoring model)
  # these could be separated out and run only once, but it complicates the remaining workflow so leaving as is
  weights_vaxpfizer1 <- get_ipw_weights(
    
    data_days_sub, "vaxpfizer1", "vaxpfizer1_status", "vaxpfizer1_atrisk",
    
    # select no more than n non-outcome samples
    # sample_type = "random_n", 
    sample_amount=ipw_sample_random_n, 
    
    ipw_formula = update(vaxpfizer1 ~ 1, formula_covars) %>% 
      update(formula_secular_region) %>%
      update(formula_timedependent) %>%
      update(formula_remove_subgroup),
    
    ipw_formula_fxd = update(vaxpfizer1 ~ 1, formula_covars) %>%
      update(formula_secular_region) %>% 
      update(formula_remove_subgroup),
    
    subgroup_level = subgroup_level
    
  )
  
  weights_vaxaz1 <- get_ipw_weights(
    
    data_days_sub, "vaxaz1", "vaxaz1_status", "vaxaz1_atrisk",
    
    # sample_type="random_n", 
    sample_amount=ipw_sample_random_n,
    
    ipw_formula = update(vaxaz1 ~ 1, formula_covars) %>% 
      update(formula_secular_region) %>% 
      update(formula_timedependent) %>% 
      update(formula_remove_subgroup),
    
    ipw_formula_fxd = update(vaxaz1 ~ 1, formula_covars) %>% 
      update(formula_secular_region) %>% 
      update(formula_remove_subgroup),
    
    subgroup_level = subgroup_level
    
  )
  
  
  data_weights <- data_days_sub %>%
    filter(
      # select all patients who experienced the outcome, and a proportion (determined in data_sample action) of those who don't
      sample_outcome==1L 
    ) %>%
    left_join(weights_vaxpfizer1, by=c("patient_id", "tstart", "tstop")) %>%
    left_join(weights_vaxaz1, by=c("patient_id", "tstart", "tstop")) %>%
    replace_na(list( 
      # weight is 1 if patient is not yet at risk or has already been vaccinated / censored
      ipweight_stbl_vaxpfizer1 = 1,
      ipweight_stbl_vaxaz1 = 1
    )) %>%
    arrange(patient_id, tstop) %>%
    group_by(patient_id) %>%
    mutate(
      cmlipweight_stbl_vaxpfizer1 = cumprod(ipweight_stbl_vaxpfizer1),
      cmlipweight_stbl_vaxaz1 = cumprod(ipweight_stbl_vaxaz1)
    ) %>%
    ungroup() %>%
    mutate(
      ## COMBINE WEIGHTS
      # take product of all weights
      ipweight_stbl = ipweight_stbl_vaxpfizer1 * ipweight_stbl_vaxaz1,
      ipweight_stbl_sample = ipweight_stbl * sample_weights,
      
      cmlipweight_stbl = cmlipweight_stbl_vaxpfizer1 * cmlipweight_stbl_vaxaz1,
      cmlipweight_stbl_sample = cmlipweight_stbl * sample_weights,
    )
  
  if(removeobs) rm(weights_vaxpfizer1, weights_vaxaz1, data_days_sub)
  
  # report weights ----
  summarise_weights <-
    data_weights %>%
    select(contains("ipweight")) %>%
    map(redacted_summary_num) %>%
    enframe()
  
  capture.output(
    walk2(summarise_weights$value, summarise_weights$name, print_num),
    file = file.path(outdir, glue("weights_table_{subgroup_level}.txt")),
    append=FALSE
  )
  
  # plot and save distribution of weights
  ipweight_histogram <- data_weights %>%
    filter(vax_atrisk==1) %>%
    ggplot() +
    geom_histogram(aes(x=ipweight_stbl)) +
    scale_x_log10() +
    theme_bw()
  
  ggsave(
    file.path(outdir, glue("weights_prob_histogram_{subgroup_level}.svg")),
    plot = ipweight_histogram
  )
  
  cmlipweight_histogram <- data_weights %>%
    ggplot() +
    geom_histogram(aes(x=cmlipweight_stbl)) +
    scale_x_log10() +
    theme_bw()
  
  ggsave(
    file.path(outdir, glue("weights_cmlprob_histogram_{subgroup_level}.svg")), 
    plot = cmlipweight_histogram
  )
  
  if(removeobs) rm(ipweight_histogram, cmlipweight_histogram)
  
  # output weight distribution file
  data_weights <- data_weights %>%
    select(
      "patient_id",
      "tstart", "tstop",
      any_of(all.vars(formula_all_rhsvars)),
      "sample_weights",
      "ipweight_stbl",
      "ipweight_stbl_sample",
      "cmlipweight_stbl",
      "cmlipweight_stbl_sample",
      "outcome",
    )
  
  cat("  \n")
  cat(glue("data_weights data size = ", nrow(data_weights)), "  \n")
  cat(glue("memory usage = ", format(object.size(data_weights), units="GB", standard="SI", digits=3L)), "  \n")
  
  write_rds(
    data_weights, 
    file.path(outdir, glue("data_weights_{subgroup_level}.rds")), 
    compress="gz"
  )
  
  if(removeobs) rm(data_weights)
  
}

