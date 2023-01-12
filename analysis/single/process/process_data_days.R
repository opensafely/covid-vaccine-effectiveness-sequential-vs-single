# this function is to be used in:
# analysis/model/preflight.R with stage = preflight
# analysis/model/msm.R with stage = msm

process_data_days <- function(stage) {
  
  data_days0 <- read_rds(here("output", "single", "stset", "data_days.rds")) # one row per patient per day
  
  if (stage == "msm") {
    
    data_days0 <- data_days0 %>%
      left_join(data_samples, by="patient_id")
    
    stage_vars <- c("sample_weights", "sample_outcome")
    
  } else if (stage == "preflight") {
    
    stage_vars <- c("vax", "death", "dereg")
    
  }
  
  data_days1 <- data_days0 %>%
    mutate(all = factor("all",levels=c("all"))) %>%
    filter(
      .[[glue("{outcome}_status")]] == 0, # follow up ends at (day after) occurrence of outcome, ie where status not >0
      vaxany1_status == .[[glue("vax{brand}1_status")]], # if brand-specific, follow up ends at (day after) occurrence of competing vaccination, ie where vax{competingbrand}_status not >0
      vaxany2_status == 0, # censor at second dose
      .[[glue("vax{brand}_atrisk")]] == 1, # select follow-up time where vax brand is being administered
    ) %>%
    left_join(data_fixed, by="patient_id") %>%
    filter(
      .[[subgroup]] == subgroup_level # select patients in current subgroup_level
    ) %>%
    mutate( 
      # this step converts logical to integer so that model coefficients print nicely in gtsummary methods
      across(where(is.logical), ~.x*1L)
    ) 
  
  postest_when_unvax <- data_days1 %>%
    # identify those who had a positive test while unvaccinated
    filter(vaxany_status == 0 & postest_status == 1) %>%
    distinct(patient_id) %>%
    mutate(postest_when_unvax = TRUE)
    
  data_days1 %>%
    left_join(postest_when_unvax, by = "patient_id") %>%
    replace_na(list(postest_when_unvax = FALSE)) %>%
    # vax*1_atrisk FALSE after a positive test:
    # (important to keep these as logical as used for filtering when calculating in `get_ipw_weights`)
    mutate(
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
    # # update vaxanyday1 to be always missing for patients who have a positive test when unvaccinated
    # mutate(
    #   across(
    #     vaxanyday1, 
    #     ~if_else(postest_when_unvax, NA_integer_, .x)
    #     )
    #   ) %>%
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
  
}



