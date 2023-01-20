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

# import command-line arguments
args <- commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  removeobs <- FALSE
  brand <- "pfizer"
  subgroup <- "all"
  outcome <- "postest"
  # vax models use less follow up time because median time to vaccination (=outcome) is ~ 30 days
  ipw_sample_random_n <- 150000 
  # outcome models use more follow up time because longer to outcome, and much fewer outcomes than vaccinations
  msm_sample_nonoutcomes_n <- 5000
  
} else {
  removeobs <- TRUE
  brand <- args[[1]]
  subgroup <- args[[2]]
  outcome <- args[[3]]
  ipw_sample_random_n <- as.integer(args[[4]])
  msm_sample_nonoutcomes_n <- as.integer(args[[5]])
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
outdir <- here("output", "single", brand, subgroup, outcome, "msm")
fs::dir_create(outdir)

# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # 
# function to calculate weights for treatment model ----
## if exposure is any vaccine, then create model for vaccination by any brand + model for death for censoring weights
## if exposure is pfizer vaccine, then create model for vaccination by pfizer + model for az and model for death for censoring weights
## if exposure is az vaccine, then create model for vaccination by az + model for pfizer and model death for censoring weights

get_ipw_weights <- function(
  data,
  event,
  event_status,
  event_atrisk,
  
  sample_type,
  sample_amount,
  
  ipw_formula,
  ipw_formula_fxd,
  
  subgroup_level
){
  
  stopifnot(sample_type %in% c("random_prop", "random_n", "nonoutcomes_n"))
  
  name <- str_remove(event_atrisk, "_atrisk")
  
  data_atrisk <- data %>%
    mutate(
      event = data[[event]],
      event_status = data[[event_status]],
      event_atrisk = data[[event_atrisk]],
    ) %>%
    filter(event_atrisk)
  
  if(sample_type=="nonoutcomes_prop"){
    data_sample <- data_atrisk %>%
      group_by(patient_id) %>%
      summarise(
        had_event = any(event>0)
      ) %>%
      ungroup() %>%
      transmute(
        patient_id,
        sample_event = sample_nonoutcomes_prop(had_event, patient_id, sample_amount),
        sample_weights_event = sample_weights(had_event, sample_event),
      )
  }
  
  if(sample_type=="nonoutcomes_n"){
    data_sample <- data_atrisk %>%
      group_by(patient_id) %>%
      summarise(
        had_event = any(event>0)
      ) %>%
      ungroup() %>%
      transmute(
        patient_id,
        sample_event = sample_nonoutcomes_n(had_event, patient_id, sample_amount),
        sample_weights_event = sample_weights(had_event, sample_event),
      )
  }
  
  if(sample_type=="random_prop"){
    data_sample <- data_atrisk %>%
      group_by(patient_id) %>%
      summarise() %>%
      ungroup() %>%
      transmute(
        patient_id,
        sample_event = sample_random_prop(patient_id, sample_amount),
        sample_weights_event = sample_event*1L,
      )
  }
  
  
  if(sample_type=="random_n"){
    data_sample <- data_atrisk %>%
      group_by(patient_id) %>%
      summarise() %>%
      ungroup() %>%
      transmute(
        patient_id,
        sample_event = sample_random_n(patient_id, sample_amount),
        sample_weights_event = sample_event*1L,
      )
  }
  
  data_atrisk_sample <- data_atrisk %>%
    left_join(data_sample, by="patient_id") %>%
    filter(sample_event)
  
  rm("data_sample")
  
  # with time-updating covariates
  cat("  \n")
  cat(glue("{event}  \n"))
  event_model <- parglm(
    formula = ipw_formula,
    data = data_atrisk_sample,
    family = binomial,
    weights = sample_weights_event,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  
  cat(glue("{event} data size = ", length(event_model$y)), "\n")
  cat(glue("memory usage = ", format(object.size(event_model), units="GB", standard="SI", digits=3L)), "\n")
  cat("warnings: ", "\n")
  print(warnings())
  
  # without time-updating covariates
  cat("  \n")
  cat(glue("{event}_fxd  \n"))
  event_model_fxd <- parglm(
    formula = ipw_formula_fxd,
    data = data_atrisk_sample,
    family = binomial,
    weights = sample_weights_event,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  
  cat(glue("{event}_fxd data size = ", length(event_model_fxd$y)), "\n")
  cat(glue("memory usage = ", format(object.size(event_model_fxd), units="GB", standard="SI", digits=3L)), "\n")
  cat("warnings: ", "\n")
  print(warnings())
  
  # save outputs
  write_rds(
    event_model, 
    file.path(outdir, glue("model_{name}_{subgroup_level}.rds")), 
    compress="gz"
    )
  write_rds(
    ipw_formula, 
    file.path(outdir, glue("model_formula_{name}_{subgroup_level}.rds")), 
    compress="gz"
    )
  
  rm("data_atrisk_sample")
  
  # get predictions from model
  data_atrisk <- data_atrisk %>%
    transmute(
      patient_id,
      tstart, tstop,
      event,
      event_status,
      # get predicted probabilities from ipw models
      pred_event=predict(event_model, type="response", newdata=data_atrisk),
      pred_event_fxd=predict(event_model_fxd, type="response", newdata=data_atrisk)
    ) %>%
    arrange(patient_id, tstop) %>%
    group_by(patient_id) %>%
    mutate(
      
      # get probability of occurrence of realised event status (time varying model)
      probevent_realised = case_when(
        event!=1L ~ 1 - pred_event,
        event==1L ~ pred_event,
        TRUE ~ NA_real_
      ),
      
      # get probability of occurrence of realised event status (non-time varying model)
      probevent_realised_fxd = case_when(
        event!=1L ~ 1 - pred_event_fxd,
        event==1L ~ pred_event_fxd,
        TRUE ~ NA_real_
      ),
      
      # stabilised inverse probability weights
      ipweight_stbl = probevent_realised_fxd/probevent_realised
      
    ) %>%
    ungroup()
  
  stopifnot("probs should all be non-null" = all(!is.na(data_atrisk$probevent_realised)))
  stopifnot("probs (fxd) should all be non-null" = all(!is.na(data_atrisk$probevent_realised_fxd)))
  
  weights <- data_atrisk %>%
    select(
      patient_id,
      tstart, tstop,
      ipweight_stbl
    )
  
  weights[[glue("ipweight_stbl_{name}")]] <- weights$ipweight_stbl
  weights$ipweight_stbl <- NULL
  
  return(weights)
  
}

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
  data_fixed <- read_rds(here("output", "single", "stset", "data_fixed.rds"))
  
  # generate data_samples
  data_samples <- read_rds(here("output", "single", "stset", "data_patients.rds")) %>%
    left_join(data_fixed, by="patient_id") %>%
    mutate(
      all = factor("all",levels=c("all")),
      tte_outcome = .[[glue("tte_{outcome}")]]
    ) %>%
    filter(
      # select patients in current subgroup_level
      .[[subgroup]] == subgroup_level 
    ) %>%
    transmute(
      patient_id,
      sample_outcome = sample_nonoutcomes_n(!is.na(tte_outcome), patient_id, msm_sample_nonoutcomes_n),
      sample_weights = sample_weights(!is.na(tte_outcome), sample_outcome),
    )
  
  ## read and process data_days (one row per person day)
  # see analysis/single/process/process_data_days_function.R for the function process_data_days_function
  # (do this _within_ loop so that it can be deleted just before models are run, to reduce RAM use)
  data_days_sub <- process_data_days_function("msm")
  
  if(removeobs) rm(data_samples, data_fixed)
  
  #print dataset size
  cat(glue("data_days_sub data size = ", nrow(data_days_sub)), "\n  ")
  cat(glue("data_days_sub patient size = ", n_distinct(data_days_sub$patient_id)), "\n  ")
  cat(glue("memory usage = ", format(object.size(data_days_sub), units="GB", standard="SI", digits=3L)), "\n  ")
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # IPW model ----
  
  if(brand=="any"){
    
    # IPW model for any vaccination
    weights_vaxany1 <- get_ipw_weights(
      data_days_sub, "vaxany1", "vaxany1_status", "vaxany1_atrisk",
      sample_type = "random_n", sample_amount=ipw_sample_random_n,
      ipw_formula = update(vaxany1 ~ 1, formula_covars) %>% update(formula_secular_region) %>% update(formula_timedependent) %>% update(formula_remove_subgroup),
      ipw_formula_fxd = update(vaxany1 ~ 1, formula_covars) %>% update(formula_secular_region) %>% update(formula_remove_subgroup),
      subgroup_level = subgroup_level
    )
    
  }
  
  if(brand!="any"){
    
    # IPW model for pfizer / az vaccination
    # these models are shared across brands (one is treatment model, one is censoring model)
    # these could be separated out and run only once, but it complicates the remaining workflow so leaving as is
    weights_vaxpfizer1 <- get_ipw_weights(
      data_days_sub, "vaxpfizer1", "vaxpfizer1_status", "vaxpfizer1_atrisk",
      # select no more than n non-outcome samples
      sample_type = "random_n", 
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
      sample_type="random_n", sample_amount=ipw_sample_random_n,
      ipw_formula = update(vaxaz1 ~ 1, formula_covars) %>% 
        update(formula_secular_region) %>% 
        update(formula_timedependent) %>% 
        update(formula_remove_subgroup),
      ipw_formula_fxd = update(vaxaz1 ~ 1, formula_covars) %>% 
        update(formula_secular_region) %>% 
        update(formula_remove_subgroup),
      subgroup_level = subgroup_level
    )
    
  }
  
  
  if(brand=="any"){
    
    data_weights <- data_days_sub %>%
      filter(
        # select all patients who experienced the outcome, and a proportion (determined in data_sample action) of those who don't
        sample_outcome==1L 
      ) %>%
      left_join(weights_vaxany1, by=c("patient_id", "tstart", "tstop")) %>%
      replace_na(list(
        # weight is 1 if patient is not yet at risk or has already been vaccinated / censored
        ipweight_stbl_vaxany1 = 1
      )) %>%
      arrange(patient_id, tstop) %>%
      group_by(patient_id) %>%
      mutate(
        cmlipweight_stbl_vaxany1 = cumprod(ipweight_stbl_vaxany1),
      ) %>%
      ungroup() %>%
      mutate(
        ipweight_stbl = ipweight_stbl_vaxany1,
        ipweight_stbl_sample = ipweight_stbl * sample_weights,
        cmlipweight_stbl = cmlipweight_stbl_vaxany1,
        cmlipweight_stbl_sample = cmlipweight_stbl * sample_weights,
      )
    
    if(removeobs) rm(weights_vaxany1)
    
  }
  
  if(brand != "any"){
    
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
    
    if(removeobs) rm(weights_vaxpfizer1, weights_vaxaz1)
    
  }
  
  
  if(removeobs) rm(data_days_sub)
  
  # report weights
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
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # MSM model ----
  
  # model 4 
  # baseline, comorbs, secular trend adjusted vaccination effect model + IP-weighted + do not use time-dependent covariates ----
  cat("  \n")
  cat("msmmod4 \n")
  msmmod4_par <- parglm(
    formula = formula_1 %>% 
      update(formula_exposure) %>% 
      update(formula_covars) %>% 
      update(formula_secular_region) %>% 
      update(formula_remove_subgroup),
    data = data_weights,
    weights = cmlipweight_stbl_sample,
    family = binomial,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  msmmod4_par$data <- NULL
  print(jtools::summ(msmmod4_par, digits =3))
  cat("warnings: ", "\n")
  print(warnings())
  
  cat(glue("msmmod4_par data size = ", length(msmmod4_par$y)), "\n")
  cat(glue("memory usage = ", format(object.size(msmmod4_par), units="GB", standard="SI", digits=3L)), "\n")
  write_rds(
    msmmod4_par, 
    file.path(outdir, glue("model4_{subgroup_level}.rds")),
    compress="gz"
    )
  if(removeobs) rm(msmmod4_par)
  
  ## print warnings
  cat("warnings: ", "\n")
  print(warnings())
  cat("  \n")
  print(gc(reset=TRUE))
  
  data_weights %>%
    summarise(
      obs = n(),
      patients = n_distinct(patient_id),
      outcomes = sum(outcome),
      incidence_prop = outcomes/patients,
      incidence_rate = outcomes/obs
    ) %>%
    write_csv(
      file.path(outdir, glue("summary_substantive_{subgroup_level}.csv"))
      )
  
  if(removeobs) rm(data_weights)
  
}

