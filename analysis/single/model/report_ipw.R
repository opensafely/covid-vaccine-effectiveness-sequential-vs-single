
# # # # # # # # # # # # # # # # # # # # #
# This script:
# imports fitted IPW models from `model_msm.R`
# calculates robust CIs taking into account patient-level clustering
# outputs forest plots for these models
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
library('gtsummary')
library('gt')
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
  outcome <- "postest"
  subgroup <- "all"
} else {
  removeobs <- TRUE
  brand <- args[[1]]
  outcome <- args[[2]]
  subgroup <- args[[3]]
}

# create output directory ----
indir <- here("output", "single", brand, subgroup, outcome, "msm") 
outdir <- file.path(indir, "report_ipw")
fs::dir_create(outdir)

# Import metadata for outcome ----
## these are created in data_define_cohorts.R script
metadata_outcomes <- events_lookup %>% 
  filter(event %in% model_outcomes) %>%
  pull(event_descr)

# list2env(metadata_outcomes, globalenv())

# reweight censored deaths or not?
# reweight_death defined in design.R

# extract formulas
list2env(list_formula_single, globalenv())
formula_remove_subgroup <- as.formula(paste0(". ~ . - ", subgroup))

var_lookup$age <- `age, degree = 2` ~ "Age"
# characteristics <- read_rds(here("output", "metadata", "baseline_characteristics.rds"))
# characteristics$age <- `age, degree = 2` ~ "Age"
# characteristics[[subgroup]] <- NULL
covar_labels <- unname(var_lookup)

## table and plot functions ----

gt_model_summary <- function(model, cluster) {
  
  # covar_labels = append(
  #   characteristics,
  #   list(
  #     timesince_hospinfectiousdischarge_pw ~ "Time since discharge from infectious hosp admission",
  #     timesince_hospnoninfectiousdischarge_pw ~ "Time since discharge from non-infectious hosp admission",
  #     #timesince_probablecovid_pw ~ "Time since probable COVID",
  #     timesince_suspectedcovid_pw ~ "Time since suspected COVID",
  #     timesince_postesttdc_pw ~ "Time since positive SARS-CoV-2 test"
  #   ) %>% set_names(., map_chr(., all.vars))
  # )
  
  tbl_reg <- tbl_regression(
    x = model,
    pvalue_fun = ~style_pvalue(.x, digits=3),
    tidy_fun = partial(tidy_plr, cluster = cluster),
    include = -contains("ns(tstop"),
    label = covar_labels
  )
  
  tbl_reg$model_obj <- NULL
  tbl_reg$inputs <- NULL
  tbl_reg$table_header$fmt_fun <- NULL
  
}


broom_model_summary <- function(model, cluster, subgroup_level) {
  
  # covar_labels = append(
  #   characteristics,
  #   list(
  #     timesince_hospinfectiousdischarge_pw ~ "Time since discharge from infectious hosp admission",
  #     timesince_hospnoninfectiousdischarge_pw ~ "Time since discharge from non-infectious hosp admission",
  #     #timesince_probablecovid_pw ~ "Time since probable COVID",
  #     timesince_suspectedcovid_pw ~ "Time since suspected COVID",
  #     timesince_postesttdc_pw ~ "Time since positive SARS-CoV-2 test"
  #   ) %>% set_names(., map_chr(., all.vars))
  # )
  
  # covar_labels <- unname(covar_labels)
  
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

#broom_model_summary(model_vaxany1, model_vaxany1$data$patient_id)


# gt_from_broom <- function(broom_obj, title){
#   
#   gt_data <- broom_obj %>%
#     filter(
#       !str_detect(variable,fixed("ns(tstop")),
#       !str_detect(variable,fixed("region")),
#       !is.na(term)
#     ) %>%
#     group_by(var_label) %>%
#     mutate(
#       characteristic = if_else(row_number()==1, var_label, ""),
#       df = n(),
#       label = if_else(df!=1, label, ""),
#     ) %>%
#     ungroup() %>%
#     transmute(
#       characteristic,
#       label,
#       or,
#       or.ll,
#       or.ul,
#       p.value
#     ) %>%
#     gt(
#       #groupname_col="var_label"
#     ) %>%
#     cols_label(
#       characteristic = "Characteristic",
#       or = "HR",
#       or.ll = "95% CI",
#       p.value = "P value",
#     ) %>%
#     fmt_number(
#       columns = starts_with(c("or")),
#       decimals = 2
#     ) %>%
#     fmt(
#       columns = all_of("p.value"),
#       fns = function(x){gtsummary::style_pvalue(x, digits=3)}
#     ) %>%
#     cols_merge_range("or.ll", "or.ul", sep = "--", autohide = TRUE) %>%
#     fmt_missing(
#       everything(),
#       missing_text="--"
#     ) %>%
#     cols_align(
#       align = "right",
#       columns = everything()
#     ) %>%
#     cols_align(
#       align = "left",
#       columns = "characteristic"
#     )
#   
#   
#   gt_data
#   
# }
# 
# forest_from_broom <- function(broom_obj, title){
#   
#   #jtools::plot_summs(ipwvaxany1)
#   #modelsummary::modelplot(ipwvaxany1, coef_omit = 'Interc|tstop', conf.type="wald", exponentiate=TRUE)
#   #sjPlot::plot_model(ipwvaxany1)
#   #all these methods use broom::tidy to get the coefficients. but tidy.glm only uses profile CIs, not Wald. (yTHO??)
#   #profile CIs will take forever on large datasets.
#   #so need to write custom function for plotting wald CIs. grr
#   
#   
#   plot_data <- broom_obj %>%
#     filter(
#       !str_detect(variable,fixed("ns(tstop")),
#       !str_detect(variable,fixed("region")),
#       !is.na(term)
#     ) %>%
#     mutate(
#       var_label = if_else(var_class=="integer", "", var_label),
#       label = if_else(reference_row %in% TRUE, paste0(label, " (ref)"),label),
#       or = if_else(reference_row %in% TRUE, 1, or),
#       variable = fct_inorder(variable),
#       variable_card = as.numeric(variable)%%2,
#     ) %>%
#     group_by(variable) %>%
#     mutate(
#       variable_card = if_else(row_number()!=1, 0, variable_card),
#       level = fct_rev(fct_inorder(paste(variable, label, sep="__"))),
#       level_label = label
#     ) %>%
#     ungroup() %>%
#     droplevels()
#   
#   var_lookup <- plot_data$var_label
#   names(var_lookup) <- plot_data$variable
#   
#   level_lookup <- plot_data$level
#   names(level_lookup) <- plot_data$level_label
#   
#   ggplot(plot_data) +
#     geom_point(aes(x=or, y=level)) +
#     geom_linerange(aes(xmin=or.ll, xmax=or.ul, y=level)) +
#     geom_vline(aes(xintercept=1), colour='black', alpha=0.8)+
#     facet_grid(rows=vars(variable), scales="free_y", switch="y", space="free_y", labeller = labeller(variable = var_lookup))+
#     scale_x_log10(breaks=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5))+
#     scale_y_discrete(breaks=level_lookup, labels=names(level_lookup))+
#     geom_rect(aes(alpha = variable_card), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, fill='grey', colour="transparent") +
#     scale_alpha_continuous(range=c(0,0.3), guide=FALSE)+
#     labs(
#       y="",
#       x="Hazard ratio",
#       colour=NULL,
#       title=title
#       #subtitle=cohort_descr
#     ) +
#     theme_minimal() +
#     theme(
#       strip.placement = "outside",
#       strip.background = element_rect(fill="transparent", colour="transparent"),
#       strip.text.y.left = element_text(angle = 0, hjust=1),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor.y = element_blank(),
#       panel.spacing = unit(0, "lines")
#     )
# }


##  Create big loop over all categories

subgroup_levels <- recoder[[subgroup]]
subgroup_level <- "all"

for(subgroup_level in subgroup_levels){
  
  # import models ----
  if(brand=="any"){
    
    model_vaxany1 <- read_rds(file.path(indir, glue("model_vaxany1_{subgroup_level}.rds")))
    ipw_formula <- read_rds(file.path(indir, glue("model_formula_vaxany1_{subgroup_level}.rds")))
    assign(as.character(model_vaxany1$call$data), model_vaxany1$data) # alternative to `data_atrisk <- model_vaxany1$data` that ensures the right model name is used
    
    ## output model coefficients
    broom_vaxany1 <- broom_model_summary(model_vaxany1, model_vaxany1$data$patient_id, subgroup_level)
    write_rds(broom_vaxany1, file.path(outdir, glue("broom_vaxany1_{subgroup_level}.rds")))
    write_csv(broom_vaxany1, file.path(outdir, glue("broom_vaxany1_{subgroup_level}.csv")))
    
    if(removeobs)
      rm(list= c(
        as.character(model_vaxany1$call$data),
        "ipw_formula", "model_vaxany1",
        "plot_vaxany1", "ggsecular_vaxany1"
      ))
    
  }
  
  if(brand!="any"){
    
    # pfizer
    model_vaxpfizer1 <- read_rds(file.path(indir, glue("model_vaxpfizer1_{subgroup_level}.rds")))
    ipw_formula <- read_rds(file.path(indir, glue("model_formula_vaxpfizer1_{subgroup_level}.rds")))
    assign(as.character(model_vaxpfizer1$call$data), model_vaxpfizer1$data)
    
    broom_vaxpfizer1 <- broom_model_summary(model_vaxpfizer1, model_vaxpfizer1$data$patient_id, subgroup_level)
    write_rds(broom_vaxpfizer1, file.path(outdir, glue("broom_vaxpfizer1_{subgroup_level}.rds")))
    write_csv(broom_vaxpfizer1, file.path(outdir, glue("broom_vaxpfizer1_{subgroup_level}.csv")))
    
    if(removeobs)
      rm(list= c(
        as.character(model_vaxpfizer1$call$data),
        "ipw_formula", "model_vaxpfizer1",
        "plot_vaxpfizer1", "ggsecular_vaxpfizer1"
      ))
    
    # AZ
    model_vaxaz1 <- read_rds(here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("model_vaxaz1_{subgroup_level}.rds")))
    ipw_formula <- read_rds(here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("model_formula_vaxaz1_{subgroup_level}.rds")))
    assign(as.character(model_vaxaz1$call$data), model_vaxaz1$data)
    
    broom_vaxaz1 <- broom_model_summary(model_vaxaz1, model_vaxaz1$data$patient_id, subgroup_level)
    write_rds(broom_vaxaz1, here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("broom_vaxaz1_{subgroup_level}.rds")))
    write_csv(broom_vaxaz1, here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("broom_vaxaz1_{subgroup_level}.csv")))
    
    if(removeobs)
      rm(list= c(
        as.character(model_vaxaz1$call$data),
        "ipw_formula", "model_vaxaz1",
        "plot_vaxaz1", "ggsecular_vaxaz1"
      ))
    
  }
  
  
  # if(outcome!="death" & reweight_death){
  #   model_death <- read_rds(here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("model_death_{subgroup_level}.rds")))
  #   ipw_formula <- read_rds(here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("model_formula_death_{subgroup_level}.rds")))
  #   assign(as.character(model_death$call$data), model_death$data)
  #   
  #   
  #   ## output model coefficients
  #   broom_death <- broom_model_summary(model_death, model_death$data$patient_id, subgroup_level)
  #   write_rds(broom_death, here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("broom_death_{subgroup_level}.rds")))
  #   write_csv(broom_death, here("output", cohort, subgroup, recentpostest_period, brand, outcome, glue("broom_death_{subgroup_level}.csv")))
  #   
  #   if(removeobs)
  #     rm(list= c(
  #       as.character(model_death$call$data),
  #       "ipw_formula", "model_death",
  #       "tab_death", "plot_death", "ggsecular_death"
  #     ))
  # }
  
}

warnings()