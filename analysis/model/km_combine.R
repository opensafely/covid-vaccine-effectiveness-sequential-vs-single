
# # # # # # # # # # # # # # # # # # # # #
# Purpose: Combine km estimates from different outcomes
#  - The script must be accompanied by three arguments:
#    `cohort` - the cohort used
#    `subgroup` - the subgroup variable
#    `outcome` - the outcome variable
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('survival')

## Import custom user functions from lib
source(here("lib", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "over12"
} else {
  cohort <- args[[1]]
}



output_dir <- ghere("output", cohort, "models", "km", "combined")
fs::dir_create(output_dir)


metaparams <-
  expand_grid(
    outcome = factor(c("postest", "emergency", "covidemergency", "covidadmitted", "covidcritcare", "coviddeath", "noncoviddeath")),
    #outcome = factor(c("postest", "covidadmitted")),
    subgroup = factor(recoder$subgroups),
  ) %>%
  mutate(
    #subgroup_level = map(as.character(subgroup), ~unname(recoder[[.x]])),
    outcome_descr = fct_recoderelevel(outcome,  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroups),
    #subgroup_level_descr = map(as.character(subgroup), ~names(recoder[[.x]])),
  )

km_estimates <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome) {
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "km", subgroup, outcome, "km_estimates_rounded.rds"))
      dat %>%
      add_column(
        subgroup_level = as.character(.[[subgroup]]),
        subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
        .before=1
      ) %>%
      select(-all_of(subgroup))
    })
  ) %>%
  unnest(data)

write_csv(km_estimates, fs::path(output_dir, "km_estimates_rounded.csv"))


contrasts_daily <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "km", subgroup, outcome, "contrasts_daily_rounded.rds"))
      dat %>%
        add_column(
          subgroup_level = as.character(.[[subgroup]]),
          subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
          .before=1
        ) %>%
        select(-all_of(subgroup))
      }
    )
  ) %>%
  unnest(data)

write_csv(contrasts_daily, fs::path(output_dir, "contrasts_daily_rounded.csv"))


contrasts_cuts <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "km", subgroup, outcome, "contrasts_cuts_rounded.rds"))
      dat %>%
        add_column(
          subgroup_level = as.character(.[[subgroup]]),
          subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
          .before=1
        ) %>%
        select(-all_of(subgroup))
    }
    )
  ) %>%
  unnest(data)

write_csv(contrasts_cuts, fs::path(output_dir, "contrasts_cuts_rounded.csv"))


contrasts_overall <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(here("output", cohort, "models", "km", subgroup, outcome, "contrasts_overall_rounded.rds"))
      dat %>%
        add_column(
          subgroup_level = as.character(.[[subgroup]]),
          subgroup_level_descr = fct_recoderelevel(.[[subgroup]], recoder[[subgroup]]),
          .before=1
        ) %>%
        select(-all_of(subgroup))
      }
    )
  ) %>%
  unnest(data)

write_csv(contrasts_overall, fs::path(output_dir, "contrasts_overall_rounded.csv"))


## move km plots to single folder ----
fs::dir_create(here("output", cohort, "models", "km", "combined"))

metaparams %>%
  mutate(
    plotdir = here("output", cohort, "models", "km", subgroup, outcome, "km_plot_rounded.png"),
    plotnewdir = glue("output", cohort, "models", "km", "combined", "km_plot_rounded_{subgroup}_{outcome}.png", .sep="/"),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}

metaparams %>%
  mutate(
    plotdir = here("output", cohort, "models", "km", subgroup, outcome, "km_plot_unrounded.png"),
    plotnewdir = glue("output", cohort, "models", "km", "combined", "km_plot_unrounded_{subgroup}_{outcome}.png", .sep="/"),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}

## plot overall estimates for inspection ----

plot_estimates <- function(estimate, estimate.ll, estimate.ul, name){

  plot_temp <-
    contrasts_overall %>%
    group_by(outcome_descr) %>%
    mutate(
      outcome_descr = fct_relabel(outcome_descr, str_wrap, width=10),
      subgroup_level_descr = fct_rev(subgroup_level_descr),

    ) %>%
    ggplot(aes(y=subgroup_level_descr)) +
    geom_vline(aes(xintercept=0), linetype="dotted", colour="darkgrey")+
    geom_point(aes(x={{estimate}}), position=position_dodge(width=-0.3))+
    geom_linerange(aes(xmin={{estimate.ll}}, xmax={{estimate.ul}}), position=position_dodge(width=-0.3))+
    facet_grid(rows=vars(subgroup_descr), cols=vars(outcome_descr), scales="free", space="free_y", switch="y")+
    scale_x_continuous(expand = expansion(mult=c(0,0.01)))+
    labs(y=NULL)+
    theme_minimal()+
    theme(
      legend.position="bottom",
      axis.text.x.top=element_text(hjust=0),

      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement="outside",
      #strip.text.y.left = element_text(angle=0),
      strip.text.y.left = element_blank(),

      panel.border = element_blank(),
      panel.spacing = unit(0.3, "lines"),
    )


  ggsave(
    filename=glue("output", cohort, "models", "km", "combined", "overall_plot_rounded_{name}.png", .sep="/"),
    plot_temp,
    width=20, height=15, units="cm"
  )

  plot_temp
}

plot_estimates(rd, rd.ll, rd.ul, "rd")
plot_estimates(rr, rr.ll, rr.ul, "rr")

