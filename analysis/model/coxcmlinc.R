
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Get cumulative incidence (kaplan meier) estimates for specified outcome, and derive risk differences
#  - import matched data
#  - adds outcome variable and restricts follow-up
#  - gets KM estimates, with covid and non covid death as competing risks
#  - The script must be accompanied by three arguments:
#    `cohort` - pfizer or moderna
#    `subgroup` - prior_covid_infection
#    `outcome` - the dependent variable

# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----


## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')


## import local functions and parameters ---

source(here("analysis", "design.R"))
source(here("lib", "functions", "utility.R"))
source(here("lib", "functions", "survival.R"))


# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  removeobjects <- FALSE
  cohort <- "pfizer"
  subgroup <- "all"
  outcome <- "postest"
  
} else {
  removeobjects <- TRUE
  cohort <- args[[1]]
  subgroup <- args[[2]]
  outcome <- args[[3]]
}

# derive symbolic arguments for programming with

cohort_sym <- sym(cohort)
subgroup_sym <- sym(subgroup)

# create output directories ----

output_dir <- ghere("output", cohort, "models", "coxcmlinc", subgroup, outcome)
fs::dir_create(output_dir)


data_matched <- read_rds(ghere("output", cohort, "match", "data_matched.rds"))

## import baseline data, restrict to matched individuals and derive time-to-event variables
data_matched <- 
  data_matched %>%
  group_by(match_id, trial_date, matching_round) %>% 
  mutate(uniquematch_id = cur_group_id()) %>% 
  ungroup() %>%
  group_by(uniquematch_id) %>%
  mutate(
    # nopriorcovid = (
    #   (is.na(positive_test_0_date) | positive_test_0_date > study_dates[[cohort]][["start_date"]]) &
    #   (is.na(primary_care_covid_case_0_date) | primary_care_covid_case_0_date > study_dates[[cohort]][["start_date"]]) &
    #   (is.na(admitted_covid_0_date) | admitted_covid_0_date > study_dates[[cohort]][["start_date"]])
    # ),
    # nopriorcovid_pair = all(nopriorcovid),
    nopriorcovid_pair = !any(prior_covid_infection),
  ) %>%
  ungroup() %>%
  filter(nopriorcovid_pair) %>%
  select(-uniquematch_id) %>%
  mutate(all="all") %>%
  group_by(patient_id, match_id, matching_round, treated) %>% 
  mutate(new_id = cur_group_id()) %>% 
  ungroup() %>%
  select(
    # select only variables needed for models to save space
    patient_id, treated, trial_date, match_id, new_id,
    controlistreated_date,
    vax1_date,
    death_date, dereg_date, coviddeath_date, noncoviddeath_date, vax2_date,
    all_of(c(glue("{outcome}_date"), subgroup, adjustment_variables))
  ) %>%
  
  mutate(

    #trial_date,
    outcome_date = .data[[glue("{outcome}_date")]],
    
    # follow-up time is up to and including censor date
    censor_date = pmin(
      dereg_date,
      vax2_date-1, # -1 because we assume vax occurs at the start of the day
      death_date,
      study_dates[["global"]]$studyend_date,
      na.rm=TRUE
    ),
    
    matchcensor_date = pmin(censor_date, controlistreated_date -1, na.rm=TRUE), # new censor date based on whether control gets treated or not

    tte_outcome = tte(trial_date - 1, outcome_date, matchcensor_date, na.censor=FALSE), # -1 because we assume vax occurs at the start of the day, and so outcomes occurring on the same day as treatment are assumed "1 day" long
    ind_outcome = censor_indicator(outcome_date, matchcensor_date),
    
  )

# outcome frequency
outcomes_per_treated <- table(outcome=data_matched$ind_outcome, treated=data_matched$treated)

table(
  data_matched$treated,
  cut(data_matched$tte_outcome, c(-Inf, 0, 1, Inf), right=FALSE, labels= c("<0", "0", ">0")), 
  useNA="ifany"
)
# should be c(0, 0, nrow(data_matched)) in each row


## competing risks cumulative risk differences ----

## cox models ----

coxcmlinc <- function(data, cuts=NULL){
 
  # data <- data_matched
  # cuts <- c(0,postbaselinecuts)
 
  # if(is.null(cuts)){cuts <- unique(c(0,data$time))}
  if(is.null(cuts)){stop("Specify cuts.")}
  
  data <- data %>% 
    # create variable for cuts[1] for tstart in tmerge
    mutate(time0 = cuts[1])
  
    fup_split <-
      data %>%
      select(new_id, treated) %>%
      uncount(weights = length(cuts)-1, .id="period_id") %>%
      mutate(
        fup_time = cuts[period_id],
        fup_period = paste0(cuts[period_id], "-", cuts[period_id+1]-1)
      ) %>%
      droplevels() %>%
      select(
        new_id, period_id, fup_time, fup_period
      )
    
    data_split <-
      tmerge(
        data1 = data,
        data2 = data,
        id = new_id,
        tstart = time0,
        tstop = tte_outcome,
        ind_outcome = event(if_else(ind_outcome, tte_outcome, NA_real_))
      ) %>%
      # add post-treatment periods
      tmerge(
        data1 = .,
        data2 = fup_split,
        id = new_id,
        period_id = tdc(fup_time, period_id)
      ) %>%
      mutate(
        period_start = cuts[period_id],
        period_end = cuts[period_id+1],
      )
    
    data_split_hotcode <- 
      data_split %>%
      mutate(
        treated_period = if_else(treated==1, period_id, 0L)
      ) %>%
      fastDummies::dummy_cols(select_columns = c("treated_period")) 
    
    treatment_term_hotcode <- paste("treated_period_", seq_len(length(cuts)-1), sep="",  collapse=" + ")
  
  
  if(length(cuts)>2){
    strata_term <-"strata(period_id) + "
  } else if(length(cuts==2)){
    strata_term <- ""
  } else
    stop("cuts must be >1")

  
  # adjustment_formula <- as.formula(
  #   eval(
  #     paste(
  #       "Surv(tstart, tstop, ind_outcome) ~", treatment_term, "+",
  #       paste0(adjustment_variables, collapse=" + ")
  #     )
  #   )
  # )
  
  adjustment_formula_hotcode <- as.formula(
    eval(
      paste(
        "Surv(tstart, tstop, ind_outcome) ~", strata_term , treatment_term_hotcode, "+",
        paste0(adjustment_variables, collapse=" + ")
      )
    )
  )
  
  
  
  
  # test equality of models under alternative specification of model formula
  # test <- coxph(
  #   adjustment_formula,
  #   data = data_split_hotcode, y=FALSE, robust=TRUE, id=new_id, na.action="na.fail"
  # )
  # test_hotcode <- coxph(
  #   adjustment_formula_hotcode,
  #   data = data_split_hotcode, y=FALSE, robust=TRUE, id=new_id, na.action="na.fail"
  # )
  # survexp(~1, data=data_split_hotcode, ratetable = test_hotcode, method= "ederer", times = times)
  
  times <- seq_len(last(postbaselinecuts))
  
  data_cmlinc <-
    data_split_hotcode %>%
    group_by(!!subgroup_sym) %>%
    nest() %>%
    mutate(
      
      # cox model
      cox_obj = map(data, ~{
        coxph(
          adjustment_formula_hotcode, 
          data = .x, x=FALSE, model=TRUE, y=FALSE, robust=TRUE, id=new_id, na.action="na.fail" # keep model (model=TRUE) in object to avoid scoping issues
        )
      }),
      # tidy model object
      cox_obj_tidy = map(cox_obj, ~broom::tidy(.x)),
      
      # marginal survival prob if no-one treated (set all treated:period_id values to zero)
      surv0 = map2(data, cox_obj, ~{
        cf_data <- mutate(.x,
            across(starts_with("treated_period_"), ~0L)
          )
        survexp(~1, data=cf_data, ratetable = .y, method="ederer", times=times)
      }),
      
      # marginal survival prob if everyone treated (set all treated:period_id values to 1 if in appropriate period)
      surv1 = map2(data, cox_obj, ~{
        cf_data <- .x %>%
          select(-starts_with("treated_period")) %>%
          mutate(treated_period = period_id) %>%
          fastDummies::dummy_cols(select_columns = c("treated_period")) 
        survexp(~1, data=cf_data, ratetable = .y, method="ederer", times=times)
      }),
      cumulrisk = map2(surv0, surv1, ~{
        tibble(
          times,
          surv1 = .y$surv,
          surv0 = .x$surv,
          #diff = surv1 - surv0
        ) %>%
        pivot_longer(
          cols=starts_with("surv"),
          names_to="treated",
          values_to="surv",
          names_prefix= "surv"
        ) %>%
        mutate(
          treated_descr = fct_recode(treated, `Vaccinated`="1", `Unvaccinated`="0"),
          treated = as.integer(treated)
        )
      })
    ) %>%
    select(!!subgroup_sym, cumulrisk) %>%
    unnest(cumulrisk)
}



# no rounding necessary as HRs are a safe statistic
coxcmlinc_cuts <- coxcmlinc(data_matched, c(0,postbaselinecuts))
coxcmlinc_overall <- coxcmlinc(data_matched, c(0,max(postbaselinecuts)))

write_rds(coxcmlinc_cuts, fs::path(output_dir, "coxcmlinc_cuts.rds"))
write_rds(coxcmlinc_overall, fs::path(output_dir, "coxcmlinc_overall.rds"))


plot_cuts <- 
  ggplot(coxcmlinc_cuts)+
  geom_line(aes(x=times, y=1-surv, group=treated)) 
ggsave(filename = fs::path(output_dir, "coxcmlinc_cuts.png"), plot = plot_cuts)

plot_overall <- 
  ggplot(coxcmlinc_overall)+
  geom_line(aes(x=times, y=1-surv, group=treated))

ggsave(filename=fs::path(output_dir, "coxcmlinc_overall.png"), plot = plot_overall)

