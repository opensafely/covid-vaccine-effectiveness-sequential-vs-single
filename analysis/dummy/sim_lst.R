# sim list vax ----
sim_list_vax <- lst(
  
  first_vax_type = bn_node(~rcat(n=..n, c("pfizer","az","moderna",""), c(0.49,0.4,0.1,0.01)), keep=FALSE),
  covid_vax_pfizer_1_day = bn_node(
    ~as.integer(runif(n=..n, firstpfizer_day, firstpfizer_day+60)),
    missing_rate = ~1-(first_vax_type=="pfizer")
  ),
  covid_vax_pfizer_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_pfizer_1_day+30, covid_vax_pfizer_1_day+60)),
    needs = c("covid_vax_pfizer_1_day"),
    missing_rate = ~0.01
  ),
  covid_vax_pfizer_3_day = bn_node(
    ~as.integer(runif(n=..n, max(covid_vax_pfizer_2_day+15,pfizerstart_day), max(covid_vax_pfizer_2_day, pfizerstart_day)+100)),
    needs = c("covid_vax_pfizer_2_day"),
    missing_rate = ~0.5
  ),
  covid_vax_pfizer_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_pfizer_3_day+120, covid_vax_pfizer_3_day+200)),
    needs = c("covid_vax_pfizer_3_day"),
    missing_rate = ~0.99
  ),
  covid_vax_az_1_day = bn_node(
    ~as.integer(runif(n=..n, firstaz_day, firstaz_day+60)),
    missing_rate = ~1-(first_vax_type=="az")
  ),
  covid_vax_az_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_az_1_day+30, covid_vax_az_1_day+60)),
    needs = c("covid_vax_az_1_day"),
    missing_rate = ~0.01
  ),
  covid_vax_az_3_day = bn_node(
    ~as.integer(runif(n=..n, max(covid_vax_az_2_day+15,pfizerstart_day), max(covid_vax_az_2_day,pfizerstart_day)+100)),
    needs = c("covid_vax_az_2_day"),
    missing_rate = ~0.99
  ),
  covid_vax_az_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_az_3_day+120, covid_vax_az_3_day+200)),
    needs = c("covid_vax_az_3_day"),
    missing_rate = ~0.99
  ),
  covid_vax_moderna_1_day = bn_node(
    ~as.integer(runif(n=..n, firstmoderna_day, firstmoderna_day+60)),
    missing_rate = ~1-(first_vax_type=="moderna")
  ),
  covid_vax_moderna_2_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_moderna_1_day+30, covid_vax_moderna_1_day+60)),
    needs = c("covid_vax_moderna_1_day"),
    missing_rate = ~0.01
  ),
  covid_vax_moderna_3_day = bn_node(
    ~as.integer(runif(n=..n, max(covid_vax_moderna_2_day+15, modernastart_day), max(covid_vax_moderna_2_day,modernastart_day)+100)),
    needs = c("covid_vax_moderna_2_day"),
    missing_rate = ~0.5
  ),
  covid_vax_moderna_4_day = bn_node(
    ~as.integer(runif(n=..n, covid_vax_moderna_3_day+120, covid_vax_moderna_3_day+200)),
    needs = c("covid_vax_moderna_3_day"),
    missing_rate = ~0.99
  ),
  
)

# sim list jcvi ----
sim_list_jcvi <- lst(
  age_aug2021 = bn_node(~age),
  
  bmi = bn_node(
    ~rfactor(n=..n, levels = c("Not obese", "Obese I (30-34.9)", "Obese II (35-39.9)", "Obese III (40+)"), p = c(0.5, 0.2, 0.2, 0.1)),
  ),
  
  care_home_type = bn_node(
    ~rfactor(n=..n, levels=c("Carehome", "Nursinghome", "Mixed", ""), p = c(0.01, 0.01, 0.01, 0.97))
  ),
  
  care_home_tpp = bn_node(
    ~care_home_type!=""
  ),
  
  care_home_code = bn_node(
    ~rbernoulli(n=..n, p = 0.01)
  ),
  
  asthma = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_neuro_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_resp_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  sev_obesity = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  diabetes = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  sev_mental = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_heart_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_kidney_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  chronic_liver_disease = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  immunosuppressed = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  asplenia = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  learndis = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  
  cev_ever = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  cev = bn_node( ~rbernoulli(n=..n, p = 0.02)),
  
  endoflife = bn_node( ~rbernoulli(n=..n, p = 0.001)),
  housebound = bn_node( ~rbernoulli(n=..n, p = 0.001)),
  
)

sim_list_demographic <- lst(
  
  has_follow_up_previous_6weeks = bn_node(
    ~rbernoulli(n=..n, p=0.999)
  ),
  
  hscworker = bn_node(
    ~rbernoulli(n=..n, p=0.01)
  ),
  
  age = bn_node(
    ~as.integer(rnorm(n=..n, mean=60, sd=15))
  ), 
  
  sex = bn_node(
    ~rfactor(n=..n, levels = c("F", "M"), p = c(0.51, 0.49)),
    missing_rate = ~0.001 # this is shorthand for ~(rbernoulli(n=1, p = 0.2))
  ),
  
  ethnicity = bn_node(
    ~rfactor(n=..n, levels = c(1,2,3,4,5), p = c(0.8, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0.25
  ),
  
  ethnicity_6_sus = bn_node(
    ~rfactor(n=..n, levels = c(0,1,2,3,4,5), p = c(0.1, 0.7, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0
  ),
  
  practice_id = bn_node(
    ~as.integer(runif(n=..n, 1, 200))
  ),
  
  msoa = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 100)), levels=1:100),
    missing_rate = ~ 0.005
  ),
  
  stp = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 36)), levels=1:36)
  ),
  
  region = bn_node(
    variable_formula = ~rfactor(n=..n, levels=c(
      "North East",
      "North West",
      "Yorkshire and The Humber",
      "East Midlands",
      "West Midlands",
      "East",
      "London",
      "South East",
      "South West"
    ), p = c(0.2, 0.2, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))
  ),
  
  imd = bn_node(
    ~factor(plyr::round_any(runif(n=..n, 1, 32000), 100), levels=seq(0,32000,100)),
    missing_rate = ~0.02,
    keep = FALSE
  ),
  
  imd_integer = bn_node(
    ~as.integer(as.character(imd)),
    keep=FALSE
  ),
  
  imd_Q5 = bn_node(
    ~factor(
      case_when(
        (imd_integer >= 0) & (imd_integer < 32844*1/5) ~ "1 (most deprived)",
        (imd_integer >= 32844*1/5) & (imd_integer < 32844*2/5) ~ "2",
        (imd_integer >= 32844*2/5) & (imd_integer < 32844*3/5) ~ "3",
        (imd_integer >= 32844*3/5) & (imd_integer < 32844*4/5) ~ "4",
        (imd_integer >= 32844*4/5) & (imd_integer <= 32844*5/5) ~ "5 (least deprived)",
        TRUE ~ "Unknown"
      ),
      levels= c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")
    ),
    missing_rate = ~0
  ),
  
  rural_urban = bn_node(
    ~rfactor(n=..n, levels = 1:9, p = rep(1/9, 9)),
    missing_rate = ~ 0.1
  ),
  
  
)

# sim list demographic ----
sim_list_demographic <- lst(
  
  has_follow_up_previous_6weeks = bn_node(
    ~rbernoulli(n=..n, p=0.999)
  ),
  
  hscworker = bn_node(
    ~rbernoulli(n=..n, p=0.01)
  ),
  
  age = bn_node(
    ~as.integer(rnorm(n=..n, mean=60, sd=15))
  ), 
  
  sex = bn_node(
    ~rfactor(n=..n, levels = c("F", "M"), p = c(0.51, 0.49)),
    missing_rate = ~0.001 # this is shorthand for ~(rbernoulli(n=1, p = 0.2))
  ),
  
  ethnicity = bn_node(
    ~rfactor(n=..n, levels = c(1,2,3,4,5), p = c(0.8, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0.25
  ),
  
  ethnicity_6_sus = bn_node(
    ~rfactor(n=..n, levels = c(0,1,2,3,4,5), p = c(0.1, 0.7, 0.05, 0.05, 0.05, 0.05)),
    missing_rate = ~ 0
  ),
  
  practice_id = bn_node(
    ~as.integer(runif(n=..n, 1, 200))
  ),
  
  msoa = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 100)), levels=1:100),
    missing_rate = ~ 0.005
  ),
  
  stp = bn_node(
    ~factor(as.integer(runif(n=..n, 1, 36)), levels=1:36)
  ),
  
  region = bn_node(
    variable_formula = ~rfactor(n=..n, levels=c(
      "North East",
      "North West",
      "Yorkshire and The Humber",
      "East Midlands",
      "West Midlands",
      "East",
      "London",
      "South East",
      "South West"
    ), p = c(0.2, 0.2, 0.3, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05))
  ),
  
  imd = bn_node(
    ~factor(plyr::round_any(runif(n=..n, 1, 32000), 100), levels=seq(0,32000,100)),
    missing_rate = ~0.02,
    keep = FALSE
  ),
  
  imd_integer = bn_node(
    ~as.integer(as.character(imd)),
    keep=FALSE
  ),
  
  imd_Q5 = bn_node(
    ~factor(
      case_when(
        (imd_integer >= 0) & (imd_integer < 32844*1/5) ~ "1 (most deprived)",
        (imd_integer >= 32844*1/5) & (imd_integer < 32844*2/5) ~ "2",
        (imd_integer >= 32844*2/5) & (imd_integer < 32844*3/5) ~ "3",
        (imd_integer >= 32844*3/5) & (imd_integer < 32844*4/5) ~ "4",
        (imd_integer >= 32844*4/5) & (imd_integer <= 32844*5/5) ~ "5 (least deprived)",
        TRUE ~ "Unknown"
      ),
      levels= c("1 (most deprived)", "2", "3", "4", "5 (least deprived)", "Unknown")
    ),
    missing_rate = ~0
  ),
  
  rural_urban = bn_node(
    ~rfactor(n=..n, levels = 1:9, p = rep(1/9, 9)),
    missing_rate = ~ 0.1
  ),
  
  
)

# sim list pre ----
sim_list_pre = lst(
  
  covid_test_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.7
  ),
  
  primary_care_covid_case_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.99
  ),
  
  
  prior_covid_test_frequency = bn_node(
    ~as.integer(rpois(n=..n, lambda=3)),
    missing_rate = ~0
  ),
  
  positive_test_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.9
  ),
  
  admitted_unplanned_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.9
  ),
  
  discharged_unplanned_0_day = bn_node(
    ~as.integer(runif(n=..n, admitted_unplanned_0_day+1, admitted_unplanned_0_day+20)),
    needs="admitted_unplanned_0_day"
  ),
  
  admitted_planned_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.9
  ),
  
  discharged_planned_0_day = bn_node(
    ~as.integer(runif(n=..n, admitted_planned_0_day+1, admitted_planned_0_day+20)),
    needs="admitted_planned_0_day"
  ),
  
  covidemergency_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.99
  ),
  
  admitted_covid_0_day = bn_node(
    ~as.integer(runif(n=..n, index_day-100, index_day-1)),
    missing_rate = ~0.99
  ),
  
)

# sim list post ----
sim_list_outcome = lst(
  
  # ## post-baseline events (outcomes)
  dereg_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 120)),
    missing_rate = ~0.99
  ),
  primary_care_covid_case_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.7
  ),
  covid_test_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.7
  ),
  postest_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.7
  ),
  emergency_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 200)),
    missing_rate = ~0.8
  ),
  emergencyhosp_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 200)),
    missing_rate = ~0.85
  ),
  covidemergency_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 200)),
    missing_rate = ~0.8
  ),
  covidemergencyhosp_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 200)),
    missing_rate = ~0.85
  ),
  
  # respemergency_day = bn_node(
  #   ~as.integer(runif(n=..n, index_day, index_day+100)),
  #   missing_rate = ~0.95
  # ),
  #
  # respemergencyhosp_day = bn_node(
  #   ~as.integer(runif(n=..n, index_day, index_day+100)),
  #   missing_rate = ~0.95
  # ),
  
  covidadmitted_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.7
  ),
  
  # placeholder for single criticalcare variable ---
  covidcritcare_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.8
  ),
  admitted_unplanned_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.7
  ),
  # admitted_planned_day = bn_node(
  #   ~as.integer(runif(n=..n, index_day, index_day+100)),
  #   missing_rate = ~0.7
  # ),
  
  coviddeath_day = bn_node(
    ~death_day,
    missing_rate = ~0.7,
    needs = "death_day"
  ),
  death_day = bn_node(
    ~ as.integer(runif(n = ..n, index_day, index_day + 100)),
    missing_rate = ~0.90
  ),
  
)

