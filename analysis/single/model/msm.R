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
outdir <- here("output", "single", brand, subgroup, outcome, "msm")
fs::dir_create(outdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Loop over all subgroup_levels ----

subgroup_levels <- recoder[[subgroup]]

for(subgroup_level in subgroup_levels){
  
  cat("  \n")
  cat(subgroup_level, "  \n")
  cat("  \n")
  
  # import IPWs 
  data_weights <- read_rds(
    here("output", "single", brand, subgroup, outcome, "ipw", glue("data_weights_{subgroup_level}.rds")) 
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
  
}

