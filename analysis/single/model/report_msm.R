
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports fitted MSMs from `model_msm.R`
# calculates robust CIs taking into account patient-level clustering
# outputs plots for the primary vaccine-outcome relationship
# outputs plots showing model-estimated spatio-temporal trends
#
# The script should only be run via an action in the project.yaml only
# The script must be accompanied by five arguments: cohort, outcome, brand, recentpostest_period, and stratum
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('lubridate')
library('survival')
library('splines')
library('parglm')
library('gtsummary')
library("sandwich")
library("lmtest")

## Import custom user functions from lib
source(here("lib", "utility_functions.R"))
source(here("lib", "redaction_functions.R"))
source(here("lib", "survival_functions.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)


if(length(args)==0){
  # use for interactive testing
  cohort <- "over80s"
  strata_var <- "all"
  recentpostest_period <- as.numeric("Inf")
  brand <- "any"
  outcome <- "postest"
  removeobs <- FALSE
} else {
  cohort <- args[[1]]
  strata_var <- args[[2]]
  recentpostest_period <- as.numeric(args[[3]])
  brand <- args[[4]]
  outcome <- args[[5]]
  removeobs <- TRUE
}



# import global vars ----
gbl_vars <- jsonlite::fromJSON(
  txt="./analysis/global-variables.json"
)

# Import metadata for outcome ----
## these are created in data_define_cohorts.R script

metadata_outcomes <- read_rds(here("output", "metadata", "metadata_outcomes.rds"))
stopifnot("outcome does not exist" = (outcome %in% metadata_outcomes[["outcome"]]))
metadata_outcomes <- metadata_outcomes[metadata_outcomes[["outcome"]]==outcome, ]

list2env(metadata_outcomes, globalenv())

### import outcomes, exposures, and covariate formulae ----
## these are created in data_define_cohorts.R script

list_formula <- read_rds(here("output", "metadata", "list_formula.rds"))
list2env(list_formula, globalenv())

formula_1 <- outcome ~ 1
formula_remove_strata_var <- as.formula(paste0(". ~ . - ",strata_var))

##  Create big loop over all categories

strata <- read_rds(here("output", "metadata", "list_strata.rds"))[[strata_var]]
strata_descr <- read_rds(here("output", "metadata", "list_strata_descr.rds"))[[strata_var]]
strata_names <- paste0("strata_",strata)
summary_list <- vector("list", length(strata))
names(summary_list) <- strata_names

for(stratum in strata){
  stratum_name <- strata_names[which(strata==stratum)]
  stratum_descr <- strata_descr[which(strata==stratum)]
  # Import processed data ----
  
  data_weights <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("data_weights_{stratum}.rds")))
  
  # import models ----
  
  #msmmod0 <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("model0_{stratum}.rds")))
  msmmod1 <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("model1_{stratum}.rds")))
  msmmod2 <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("model2_{stratum}.rds")))
  #msmmod3 <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("model3_{stratum}.rds")))
  msmmod4 <- read_rds(here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("model4_{stratum}.rds")))
  
  ## report models ----
  
  #robust0 <- tidy_plr(msmmod0, cluster=data_weights$patient_id)
  robust1 <- tidy_plr(msmmod1, cluster=data_weights$patient_id)
  robust2 <- tidy_plr(msmmod2, cluster=data_weights$patient_id)
  #robust3 <- tidy_plr(msmmod3, cluster=data_weights$patient_id)
  robust4 <- tidy_plr(msmmod4, cluster=data_weights$patient_id)
  
  robust_summary <- bind_rows(
    list(
      #"0"=mutate(robust0,, model_descr="unadjusted"),
      "1"=mutate(robust1, model_descr="Region-stratified Cox model, with no further adjustment"),
      "2"=mutate(robust2, model_descr="Region-stratified Cox model, with adjustment for baseline confounders"),
      #"3"=mutate(robust3, model_descr="Region-stratified Cox model, with adjustment for baseline and time-varying confounders"),
      "4"=mutate(robust4, model_descr="Region-stratified marginal structural Cox model, with adjustment for baseline and time-varying confounders")
    ),
    .id = "model"
  ) %>%
    mutate(
      stratum=stratum_descr,
    )
  
  summary_list[[stratum_name]] <- robust_summary
  
}


summary_df <- summary_list %>% bind_rows %>%
  mutate(
    model_descr = fct_reorder(model_descr, as.numeric(model)),
    model_descr_wrap = fct_inorder(str_wrap(model_descr, 30)),
  ) %>%
  select(
    stratum, model, model_descr, model_descr_wrap, term, estimate, conf.low, conf.high, std.error, statistic, p.value, or, or.ll, or.ul
  )

write_csv(summary_df, path = here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("estimates.csv")))

# create plot
msmmod_effect_data <- summary_df %>%
  filter(str_detect(term, "timesincevax")) %>%
  mutate(
    term=str_replace(term, pattern="timesincevax\\_pw", ""),
    term=fct_inorder(term),
    term_left = as.numeric(str_extract(term, "\\d+"))-1,
    term_right = as.numeric(str_extract(term, "\\d+$")),
    term_right = if_else(is.na(term_right), max(term_left)+7, term_right),
    term_midpoint = term_left + (term_right-term_left)/2,
  ) %>%
  group_by(model, term) %>%
  mutate(
    #Vaccine effectiveness
    ve = 1-or,
    ve.ll = 1-or.ul,
    ve.ul = 1-or.ll,
    # conduct z test for difference in effects between strata
    diff = estimate - first(estimate),
    diff.std.error = if_else(row_number()!=1, sqrt((std.error^2) + (first(std.error))^2), NA_real_),
    z = diff/diff.std.error,
    diff.ll = diff + qnorm(0.025)*diff.std.error,
    diff.ul = diff + qnorm(0.975)*diff.std.error,
    diff.p.value = 2 * pmin(pnorm(z), pnorm(-z))
  ) %>%
  ungroup()

write_csv(msmmod_effect_data, path = here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("estimates_timesincevax.csv")))

msmmod_effect <-
  ggplot(data = msmmod_effect_data, aes(colour=as.factor(stratum))) +
  geom_point(aes(y=or, x=term_midpoint), position = position_dodge(width = 0.6))+
  geom_linerange(aes(ymin=or.ll, ymax=or.ul, x=term_midpoint), position = position_dodge(width = 0.6))+
  geom_hline(aes(yintercept=1), colour='grey')+
  facet_grid(rows=vars(model_descr_wrap), switch="y")+
  scale_y_log10(
    #breaks=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    #sec.axis = sec_axis(~(1-.), name="Effectiveness", breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99), labels = scales::label_percent(1))
  )+
  scale_x_continuous(breaks=unique(msmmod_effect_data$term_left))+
  scale_colour_brewer(type="qual", palette="Set2")+#, guide=guide_legend(reverse = TRUE))+
  #coord_cartesian(ylim=c(max(c(0.005, min(msmmod_effect_data$or.ll))), max(c(1, msmmod_effect_data$or.ul)))) +
  labs(
    y="Hazard ratio, versus no vaccination",
    x="Time since first dose",
    colour=NULL#,
    #title=glue("{outcome_descr} by time since first {brand} vaccine"),
    #subtitle=cohort_descr
  ) +
  theme_bw()+
  theme(
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "right"
  )

## save plot
ggsave(filename=here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("VE_plot.svg")), msmmod_effect, width=20, height=18, units="cm")
ggsave(filename=here("output", cohort, strata_var, recentpostest_period, brand, outcome, glue("VE_plot.png")), msmmod_effect, width=20, height=18, units="cm")