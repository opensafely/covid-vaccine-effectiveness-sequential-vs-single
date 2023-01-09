
# # # # # # # # # # # # # # # # # # # # #
# This script:
# combines outputs from the `report_msm.R` scripts across different outcomes and brands
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('glue')
library('here')
library('lubridate')
library('survival')
library('splines')
library('parglm')
library('gtsummary')
library("sandwich")
library("lmtest")
library('gt')

## Import custom user functions from lib
source(here("lib", "utility_functions.R"))
source(here("lib", "redaction_functions.R"))
source(here("lib", "survival_functions.R"))

# import command-line arguments ----

args <- commandArgs(trailingOnly=TRUE)



if(length(args)==0){
  # use for interactive testing
  removeobs <- FALSE
  cohort <- "over80s"
  strata_var <- "all"
  recent_postestperiod <- as.numeric("0")
  
} else{
  removeobs <- TRUE
  cohort <- args[[1]]
  strata_var <- args[[2]]
  recent_postestperiod <- as.numeric(args[[3]])
}



# import global vars ----
gbl_vars <- jsonlite::fromJSON(
  txt="./analysis/global-variables.json"
)

# Import metadata for outcomes ----
## these are created in data_define_cohorts.R script

metadata_outcomes <- read_rds(here("output", "metadata", "metadata_outcomes.rds"))


fs::dir_create(here("output", cohort, strata_var, recent_postestperiod, "combined"))

##  Create big loop over all categories

strata <- read_rds(here("output", "metadata", "list_strata.rds"))[[strata_var]]
strata_descr <- read_rds(here("output", "metadata", "list_strata_descr.rds"))[[strata_var]]
summary_list <- vector("list", length(strata_descr))
names(summary_list) <- strata_descr

# import models ----

# select outcomes
if(recent_postestperiod==0){
  outcomes <- c(
    "postest",
    "covidadmitted",
    "death"
  )
  
  deathtypes <- c("alldeath")
}
if(recent_postestperiod>0){
  outcomes <- c(
    "postest",
    "covidadmitted",
    "coviddeath",
    "noncoviddeath",
    "death"
  )
  
  deathtypes <- c("alldeath", "specificdeath")
}

estimates <-
  metadata_outcomes %>%
  filter(outcome %in% outcomes) %>%
  mutate(
    outcome = fct_inorder(outcome),
    outcome_descr = fct_inorder(map_chr(outcome_descr, ~paste(stringi::stri_wrap(., width=14, simplify=TRUE, whitespace_only=TRUE), collapse="\n")))
  ) %>%
  crossing(
    tibble(
      brand = fct_inorder(c("any", "pfizer", "az")),
      brand_descr = fct_inorder(c("Any vaccine", "BNT162b2", "ChAdOx1"))
    )
  ) %>%
  mutate(
    brand = fct_inorder(brand),
    brand_descr = fct_inorder(brand_descr),
    estimates = map2(brand, outcome, ~read_csv(here("output", cohort, strata_var, recent_postestperiod, .x, .y, glue("estimates_timesincevax.csv"))))
  ) %>%
  unnest(estimates) %>%
  mutate(
    model_descr = fct_inorder(model_descr),
  )



estimates_formatted <- estimates %>%
  transmute(
    outcome_descr,
    brand_descr,
    stratum,
    model,
    model_descr,
    term=str_replace(term, pattern="timesincevax\\_pw", ""),
    HR =scales::label_number(accuracy = .01, trim=TRUE)(or),
    HR_CI = paste0("(", scales::label_number(accuracy = .01, trim=TRUE)(or.ll), "-", scales::label_number(accuracy = .01, trim=TRUE)(or.ul), ")"),
    VE = scales::label_number(accuracy = .1, trim=FALSE, scale=100)(ve),
    VE_CI = paste0("(", scales::label_number(accuracy = .1, trim=TRUE, scale=100)(ve.ll), "-", scales::label_number(accuracy = .1, trim=TRUE, scale=100)(ve.ul), ")"),
    
    HR_ECI = paste0(HR, " ", HR_CI),
    VE_ECI = paste0(VE, " ", VE_CI),
  )

if(strata_var!="all"){
  estimates_formatted <- estimates_formatted %>%
    mutate(
      diff = scales::label_number(accuracy = .01, trim=TRUE)(diff),
      diff_CI = paste0("(", scales::label_number(accuracy = .01, trim=TRUE)(diff.ll), "-", scales::label_number(accuracy = .01, trim=TRUE)(diff.ul), ")"),
      diff_ECI = paste0(diff, " ", diff_CI),
    )
}

estimates_formatted_wide <- estimates_formatted %>%
  select(outcome_descr, brand_descr, stratum, model, term, HR_ECI, VE_ECI) %>%
  pivot_wider(
    id_cols=c(outcome_descr, brand_descr, term, stratum),
    names_from = model,
    values_from = c(HR_ECI, VE_ECI),
    names_glue = "{model}_{.value}"
  )

write_csv(estimates, path = here("output", cohort, strata_var, recent_postestperiod, "combined", glue("estimates_timesincevax.csv")))
write_csv(estimates_formatted, path = here("output", cohort, strata_var, recent_postestperiod, "combined", glue("estimates_formatted_timesincevax.csv")))
write_csv(estimates_formatted_wide, path = here("output", cohort, strata_var, recent_postestperiod, "combined", glue("estimates_formatted_wide_timesincevax.csv")))



formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}



# create forest plot
msmmod_effect_data <- estimates %>%
  mutate(
    term=str_replace(term, pattern="timesincevax\\_pw", ""),
    term=fct_inorder(term),
    term_left = as.numeric(str_extract(term, "\\d+"))-1,
    term_right = as.numeric(str_extract(term, "\\d+$")),
    term_right = if_else(is.na(term_right), 63, term_right),
    term_midpoint = term_left + (term_right-term_left)/2,
    #stratum = if_else(stratum=="all", "", stratum)
  )

if(strata_var=="all"){
  for(deathtype in deathtypes){
    
    if(deathtype=="alldeath") outcomes <- c("postest","covidadmitted","death")
    if(deathtype=="specificdeath") outcomes <- c("postest","covidadmitted","coviddeath","noncoviddeath")
    
    msmmod_effect_data_plot <- msmmod_effect_data %>%
      filter(
        or !=Inf, or !=0,
        outcome %in% outcomes
      )
    
    msmmod_effect <-
      msmmod_effect_data_plot %>%
      ggplot(aes(colour=model_descr)) +
      geom_hline(aes(yintercept=1), colour='grey')+
      geom_point(aes(y=or, x=term_midpoint), position = position_dodge(width = 1.5), size=0.8)+
      geom_linerange(aes(ymin=or.ll, ymax=or.ul, x=term_midpoint), position = position_dodge(width = 1.5))+
      facet_grid(rows=vars(outcome_descr), cols=vars(brand_descr), switch="y")+
      scale_y_log10(
        breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
        limits = c(0.05, max(c(1, msmmod_effect_data_plot$or.ul))),
        oob = scales::oob_keep,
        sec.axis = sec_axis(
          ~(1-.),
          name="Effectiveness",
          breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99),
          labels = function(x){formatpercent100(x, 1)}
        )
      )+
      scale_x_continuous(breaks=unique(msmmod_effect_data_plot$term_left), expand=expansion(mult=c(0, NA)), limits=c(0,NA))+
      scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=1))+
      coord_cartesian() +
      labs(
        y="Hazard ratio, versus no vaccination",
        x="Days since first dose",
        colour=NULL#,
        #title=glue("Outcomes by time since first {brand} vaccine"),
        #subtitle=cohort_descr
      ) +
      theme_bw(base_size=12)+
      theme(
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #strip.text.y.left = element_text(angle = 0),
        
        panel.spacing = unit(1, "lines"),
        
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0, face= "italic"),
        
        legend.position = "bottom"
      )
    
    ## save plot
    ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_{deathtype}.svg")), msmmod_effect, width=20, height=20, units="cm")
    ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_{deathtype}.png")), msmmod_effect, width=20, height=20, units="cm")
    
    
    msmmod_effect_free <-
      msmmod_effect_data_plot %>%
      ggplot(aes(colour=model_descr)) +
      geom_hline(aes(yintercept=0), colour='grey')+
      geom_point(aes(y=log(or), x=term_midpoint), position = position_dodge(width = 1.5), size=0.8)+
      geom_linerange(aes(ymin=log(or.ll), ymax=log(or.ul), x=term_midpoint), position = position_dodge(width = 1.5))+
      facet_grid(rows=vars(outcome_descr), cols=vars(brand_descr), switch="y", scales="free_y")+
      scale_y_continuous(
        labels = function(x){scales::label_number(0.001)(exp(x))},
        breaks = function(x){log(scales::breaks_log(n=6, base=10)(exp(x)))},
        sec.axis = sec_axis(
          ~(1-exp(.)),
          name="Effectiveness",
          breaks = function(x){1-(scales::breaks_log(n=6, base=10)(1-x))},
          labels = function(x){formatpercent100(x, 1)}
        )
      )+
      scale_x_continuous(breaks=unique(msmmod_effect_data_plot$term_left), expand=expansion(mult=c(0, NA)), limits=c(0,NA))+
      scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=1))+
      labs(
        y="Hazard ratio, versus no vaccination",
        x="Days since first dose",
        colour=NULL#,
        #title=glue("Outcomes by time since first {brand} vaccine"),
        #subtitle=cohort_descr
      ) +
      theme_bw(base_size=12)+
      theme(
        panel.border = element_blank(),
        axis.line.y = element_line(colour = "black"),
        
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.background = element_blank(),
        strip.placement = "outside",
        #strip.text.y.left = element_text(angle = 0),
        
        panel.spacing = unit(1, "lines"),
        
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0, face= "italic"),
        
        legend.position = "bottom"
      )
    msmmod_effect_free
    
    ## save plot
    ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_free_{deathtype}.svg")), msmmod_effect_free, width=20, height=20, units="cm")
    ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_free_{deathtype}.png")), msmmod_effect_free, width=20, height=20, units="cm")
  }
  
}

if(strata_var!="all"){
  
  msmmod_effect_data_plot <- msmmod_effect_data %>%
    filter(
      model==max(model),
      or !=Inf, or !=0
    )
  
  msmmod_effect <-
    ggplot(data = msmmod_effect_data_plot, aes(colour=stratum)) +
    geom_hline(aes(yintercept=1), colour='grey')+
    geom_point(aes(y=or, x=term_midpoint), position = position_dodge(width = 1.5), size=0.8)+
    geom_linerange(aes(ymin=or.ll, ymax=or.ul, x=term_midpoint), position = position_dodge(width = 1.5))+
    facet_grid(rows=vars(outcome_descr), cols=vars(brand_descr), switch="y")+
    scale_y_log10(
      breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
      limits = c(0.009, max(c(1, msmmod_effect_data_plot$or.ul))),
      oob = scales::oob_keep,
      sec.axis = sec_axis(
        ~(1-.),
        name="Effectiveness",
        breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99),
        labels = function(x){formatpercent100(x, 1)}
      )
    )+
    scale_x_continuous(breaks=unique(msmmod_effect_data_plot$term_left))+
    scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=1))+
    coord_cartesian() +
    labs(
      y="Hazard ratio, versus no vaccination",
      x="Days since first dose",
      colour=NULL
    ) +
    theme_bw(base_size=12)+
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      #strip.text.y.left = element_text(angle = 0),
      
      panel.spacing = unit(1, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom"
    )
  
  ## save plot
  ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot.svg")), msmmod_effect, width=20, height=20, units="cm")
  ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot.png")), msmmod_effect, width=20, height=20, units="cm")
  
  
  msmmod_effect_free <-
    ggplot(data = msmmod_effect_data_plot, aes(colour=stratum)) +
    geom_hline(aes(yintercept=0), colour='grey')+
    geom_point(aes(y=log(or), x=term_midpoint), position = position_dodge(width = 1.5), size=0.8)+
    geom_linerange(aes(ymin=log(or.ll), ymax=log(or.ul), x=term_midpoint), position = position_dodge(width = 1.5))+
    facet_grid(rows=vars(outcome_descr), cols=vars(brand_descr), switch="y", scales="free_y")+
    scale_y_continuous(
      labels = function(x){scales::label_number(0.001)(exp(x))},
      breaks = function(x){log(scales::breaks_log(n=6, base=10)(exp(x)))},
      sec.axis = sec_axis(
        ~(1-exp(.)),
        name="Effectiveness",
        breaks = function(x){1-(scales::breaks_log(n=6, base=10)(1-x))},
        labels = function(x){formatpercent100(x, 1)}
      )
    )+
    scale_x_continuous(breaks=unique(msmmod_effect_data_plot$term_left))+
    scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=1))+
    labs(
      y="Hazard ratio, versus no vaccination",
      x="Days since first dose"
    ) +
    theme_bw(base_size=12)+
    theme(
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black"),
      
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      #strip.text.y.left = element_text(angle = 0),
      
      panel.spacing = unit(1, "lines"),
      
      plot.title = element_text(hjust = 0),
      plot.title.position = "plot",
      plot.caption.position = "plot",
      plot.caption = element_text(hjust = 0, face= "italic"),
      
      legend.position = "bottom"
    )
  msmmod_effect_free
  
  ## save plot
  ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_free.svg")), msmmod_effect_free, width=20, height=20, units="cm")
  ggsave(filename=here("output", cohort, strata_var, recent_postestperiod, "combined", glue("VE_plot_free.png")), msmmod_effect_free, width=20, height=20, units="cm")
  
  
  
}