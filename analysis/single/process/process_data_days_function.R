# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# this function applies some data processing to avoid code replication
# is to be used in:
# analysis/model/preflight.R with stage = preflight
# analysis/model/msm.R with stage = msm
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

process_data_days_function <- function(
  file, # preflight or model
  stage, # vaccine or outcome
  iteration # 1:process_data_days_n
  ) {
  
  cat(glue("Iteration {iteration}:"), "\n")
  
  cat("Import data_days:\n")
  data_days0 <- read_rds(here("output", "single", "stset", glue("data_days_{iteration}.rds")))
  
  if (file == "preflight") {
    
    stage_vars <- c("vax", "death", "dereg")
    
  } else if (file == "model") {
    
    cat("Join data_days to data_samples:\n")
    data_days0 <- data_days0 %>%
      left_join(data_samples, by="patient_id")
    
    stage_vars <- c("sample_weights", "sample_outcome")
    
  } 
  
  cat("Process step 1:\n")
  data_days1 <- data_days0 %>%
    filter(
      # stop follow-up maxfup days after vax 1
      vaxany1_timesince <= maxfup,
      # follow up ends at (day after) occurrence of outcome, ie where status not >0
      .[[glue("{outcome}_status")]] == 0, 
      # if brand-specific, follow up ends at (day after) occurrence of competing vaccination, 
      # ie where vax{competingbrand}_status not >0
      vaxany1_status == .[[glue("vax{brand}1_status")]], 
      # # censor at second dose
      # vaxany2_status == 0, 
      # select follow-up time where vax brand is being administered
      .[[glue("vax{brand}_atrisk")]] == 1, 
    ) %>%
    # join fixed covariates
    # inner_join as lhs data filtered on vax status, and rhs data (data_fixed)
    # contains only patients in the relevant subgroup
    inner_join(data_fixed, by="patient_id") %>%
    mutate( 
      # this step converts logical to integer so that model coefficients print nicely in gtsummary methods
      across(where(is.logical), ~.x*1L)
    ) 
  
  cat("Process step 2:\n")
  postest_when_unvax <- data_days1 %>%
    # identify those who had a positive test while unvaccinated
    filter(vaxany_status == 0 & postest_status == 1) %>%
    distinct(patient_id) %>%
    mutate(postest_when_unvax = TRUE)
    
  cat("Process step 3:\n")
  data_days <- data_days1 %>%
    # join those who had a positive test while unvaccinated
    left_join(postest_when_unvax, by = "patient_id") %>%
    replace_na(list(postest_when_unvax = FALSE)) %>%
    mutate(
      # vax*1_atrisk FALSE after a positive test:
      # (important to keep these as logical as used for filtering when calculating in `get_ipw_weights`)
      vaxany1_atrisk = (vaxany1_status==0 & vaxany_atrisk==1 & postest_status==0),
      vaxpfizer1_atrisk = (vaxany1_status==0 & vaxpfizer_atrisk==1 & postest_status==0),
      vaxaz1_atrisk = (vaxany1_status==0 & vaxaz_atrisk==1 & postest_status==0),
      death_atrisk = (death_status==0),
    ) %>%
    # update vax1 and vax variables to be always zero for patients who have a positive test when unvaccinated
    # this enables estimation of the "modified estimand":
    # compare all vaccinated people with no previous record of SARS-CoV-02 infection with everybody else (including people vaccinated after infection)
    # this is implemented by ensuring that the vaccination status never changes after first documented infection
    mutate(
      across(
        c(
          vaxany1, vaxpfizer1, vaxaz1, 
          matches(c("^vax[[:alpha:]]+1_status$", "^vax[[:alpha:]]+_status", "^vax[[:alpha:]]+1_timesince$"))
          ),
        ~ if_else(postest_when_unvax, 0L, .x)
      )
    ) %>%
    mutate(
      timesincevax_pw = timesince_cut(vaxany1_timesince, c(0,postbaselinecuts), "pre-vax"),
      outcome = .[[outcome]],
      vax = .[[glue("vax{brand}1")]],
      vax_atrisk = .[[glue("vax{brand}1_atrisk")]]
    ) %>%
    select(
      "patient_id",
      "all",
      "tstart", 
      "tstop",
      "outcome",
      "timesincevax_pw",
      any_of(all.vars(formula_all_rhsvars)),
      all_of(stage_vars),
      "vaxany1_atrisk",
      "vaxpfizer1_atrisk",
      "vaxaz1_atrisk",
      "death_atrisk",
      "vax_atrisk",
      "vaxany1",
      "vaxpfizer1",
      "vaxaz1",
      "vaxany1_status",
      "vaxpfizer1_status",
      "vaxaz1_status",
    )
  
  if (file == "preflight") {
    
    if (stage == "vaccine") {
      
      # vaccination models
      data_days_vax <- data_days %>%
        # select follow-up time where vax brand is being administered
        filter(vax_atrisk) 
      
      return(data_days_vax)
      
    } else if (stage == "outcome") {
      
      data_days_outcome <- data_days %>%
        group_by(patient_id) %>%
        mutate(
          had_outcome = any(outcome>0),
        ) %>%
        ungroup()
      
      return(data_days_outcome)
      
    }
    
  } else if (file == "model") {
    
    return(data_days)
    
  }
  
  
}



