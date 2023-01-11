
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
source(here("analysis", "functions", "utility.R"))

## Import design elements
source(here("analysis", "design.R"))

## import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  # use for interactive testing
  cohort <- "pfizer"
} else {
  cohort <- args[[1]]
}

indir <- here("output", "sequential", cohort, "models", "km")
outdir <- file.path(indir, "combined")
fs::dir_create(outdir)

metaparams <-
  expand_grid(
    outcome = factor(model_outcomes),
    subgroup = factor(recoder$subgroups),
  ) %>%
  mutate(
    outcome_descr = fct_recoderelevel(outcome,  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroups),
  )


# combine km estimates ----
km_estimates <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome) {
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "km_estimates_rounded.rds"))
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

write_csv(km_estimates, fs::path(outdir, "km_estimates_rounded.csv"))


contrasts_km_daily <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "contrasts_km_daily_rounded.rds"))
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

write_csv(contrasts_km_daily, fs::path(outdir, "contrasts_km_daily_rounded.csv"))


contrasts_km_cuts <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "contrasts_km_cuts_rounded.rds"))
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

write_csv(contrasts_km_cuts, fs::path(outdir, "contrasts_km_cuts_rounded.csv"))


contrasts_km_overall <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "contrasts_km_overall_rounded.rds"))
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

write_csv(contrasts_km_overall, fs::path(outdir, "contrasts_km_overall_rounded.csv"))

# combine cox estimates ----
contrasts_cox_cuts <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "contrasts_cox_cuts.rds"))
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

write_csv(contrasts_cox_cuts, fs::path(outdir, "contrasts_cox_cuts.csv"))


contrasts_cox_overall <- metaparams %>%
  mutate(
    data = pmap(list(cohort, subgroup, outcome), function(cohort, subgroup, outcome){
      subgroup <- as.character(subgroup)
      dat <- read_rds(file.path(indir, subgroup, outcome, "contrasts_cox_overall.rds"))
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

write_csv(contrasts_cox_overall, fs::path(outdir, "contrasts_cox_overall.csv"))

## move km plots to single folder ----
fs::dir_create(here("output", cohort, "models", "km", "combined"))

metaparams %>%
  mutate(
    plotdir = file.path(indir, subgroup, outcome, "km_plot_rounded.png"),
    plotnewdir = file.path(outdir, "km_plot_rounded_{subgroup}_{outcome}.png"),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}

metaparams %>%
  mutate(
    plotdir = file.path(indir, subgroup, outcome, "km_plot_unrounded.png"),
    plotnewdir = file.path(outdir, "km_plot_unrounded_{subgroup}_{outcome}.png"),
  ) %>%
  {walk2(.$plotdir, .$plotnewdir, ~fs::file_copy(.x, .y, overwrite = TRUE))}

## plot overall km estimates for inspection ----

plot_estimates <- function(estimate, estimate.ll, estimate.ul, name){

  plot_temp <-
    contrasts_km_overall %>%
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
    filename=file.path(outdir, "overall_plot_rounded_{name}.png"),
    plot_temp,
    width=20, height=15, units="cm"
  )

  plot_temp
}

plot_estimates(rd, rd.ll, rd.ul, "rd")
plot_estimates(rr, rr.ll, rr.ul, "rr")

