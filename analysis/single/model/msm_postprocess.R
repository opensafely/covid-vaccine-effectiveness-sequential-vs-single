
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports fitted IPW models from `msm.R`
# calculates robust CIs taking into account patient-level clustering
# 
# imports fitted MSMs from `msm.R`
# calculates robust CIs taking into account patient-level clustering
# outputs plots for the primary vaccine-outcome relationship
# outputs plots showing model-estimated spatio-temporal trends
#
# The script should only be run via an action in the project.yaml only
# The script must be accompanied by five arguments: brand, outcome, subgroup
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
# library('gtsummary')
# library('gt')
library("sandwich")
library("lmtest")

## Import custom user functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))
source(here("analysis", "functions", "survival.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
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

# create output directory ----
indir <- here("output", "single", brand, subgroup, outcome, "msm") 
outdir <- here("output", "single", brand, subgroup, outcome, "postprocess") 
fs::dir_create(outdir)

# extract formulas
list2env(list_formula_single, globalenv())
formula_1 <- outcome ~ 1
formula_remove_subgroup <- as.formula(paste0(". ~ . - ", subgroup))

var_lookup$age <- `age, degree = 2` ~ "Age"
covar_labels <- unname(var_lookup)

## define functions ----
broom_model_summary <- function(model, cluster, subgroup_level) {
  
  tbl_reg <- broom.helpers::tidy_plus_plus(
    model = model,
    tidy_fun = partial(tidy_plr, cluster = cluster),
    include = -contains("ns(tstop")#,
    # variable_labels = covar_labels
  ) %>%
    add_column(
      subgroup=subgroup_level,
      .before=1
    )
  
}

##  Create big loop over all categories

subgroup_levels <- recoder[[subgroup]]

summary_list <- vector("list", length(subgroup_levels))
names(summary_list) <- subgroup_levels

# subgroup_level <- "all"

for(subgroup_level in subgroup_levels) {
  
  cat("  \n")
  cat(subgroup_level, "  \n")
  cat("  \n")
  
  # IPW
  # pfizer
  # import models
  model_vaxbrand1 <- read_rds(file.path(indir, glue("model_vax{brand}1_{subgroup_level}.rds")))
  ipw_formula <- read_rds(file.path(indir, glue("model_formula_vax{brand}1_{subgroup_level}.rds")))
  assign(as.character(model_vaxbrand1$call$data), model_vaxbrand1$data)
  # calculate robust CIs taking into account patient-level clustering
  broom_vaxbrand1 <- broom_model_summary(model_vaxbrand1, model_vaxbrand1$data$patient_id, subgroup_level)
  # warnings()
  # save outputs
  write_rds(broom_vaxbrand1, file.path(outdir, glue("broom_vax{brand}1_{subgroup_level}.rds")))
  write_csv(broom_vaxbrand1, file.path(outdir, glue("broom_vax{brand}1_{subgroup_level}.csv")))
  
  if(removeobs) rm(list= c(as.character(model_vaxbrand1$call$data), "ipw_formula", "model_vaxbrand1", "broom_vaxbrand1"))
  
  
  # MSM
  data_weights <- read_rds(file.path(indir, glue("data_weights_{subgroup_level}.rds")))
  msmmod4 <- read_rds(file.path(indir, glue("model4_{subgroup_level}.rds")))
  robust4 <- tidy_plr(msmmod4, cluster=data_weights$patient_id)
  robust_summary <- bind_rows(
    list("4" = mutate(robust4, model_descr="Region-stratified marginal structural Cox model, with adjustment for baseline and time-varying confounders")),
    .id = "model"
  ) %>%
    mutate(
      subgroup_level=subgroup_level,
    )
  
  summary_list[[subgroup_level]] <- robust_summary
  
}

# tidy and save output
summary_df <- summary_list %>% bind_rows %>%
  mutate(
    model_descr = fct_reorder(model_descr, as.numeric(model)),
    model_descr_wrap = fct_inorder(str_wrap(model_descr, 30)),
  ) %>%
  select(
    subgroup_level, model, model_descr, model_descr_wrap, term, estimate, conf.low, conf.high, std.error, statistic, p.value, or, or.ll, or.ul
  )

write_csv(summary_df, file.path(outdir, "estimates.csv"))

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

write_csv(msmmod_effect_data, file.path(outdir, "estimates_timesincevax.csv"))

msmmod_effect <-
  ggplot(data = msmmod_effect_data, aes(colour=as.factor(subgroup_level))) +
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
# ggsave(
#   filename=file.path(outdir, glue("VE_plot.svg")),
#   msmmod_effect, 
#   width=20, height=18, units="cm"
# )
ggsave(
  filename=here(outdir, glue("VE_plot.png")), 
  msmmod_effect, 
  width=20, height=18, units="cm"
)


