
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
library("sandwich")
library("lmtest")

## Import custom user functions
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))
source(here("analysis", "functions", "survival.R"))

# create output directory 
indir <- here("output", "single") 
outdir <- here("output", "single", "combine") 
fs::dir_create(outdir)


model_key <- expand_grid(
  brand=model_brands,
  subgroup=model_subgroups,
  outcome=model_outcomes
)

add_descr <- function(.data) {
  
  brand_descr <- brand_lookup %>% filter(brand %in% model_brands) %>% pull(brand_descr)
  subgroup_descr <- names(recoder$subgroups[recoder$subgroups %in% model_subgroups])
  outcome_descr <- events_lookup %>% filter(event %in% model_outcomes) %>% pull(event_descr)
  
  .data %>%
    mutate( 
      brand_descr = factor(brand, levels = model_brands, labels = brand_descr),
      subgroup_descr = factor(subgroup, levels = model_subgroups, labels = subgroup_descr),
      outcome_descr = factor(outcome, levels = model_outcomes, labels = outcome_descr),
    ) %>%
    mutate(across(brand, factor, levels = model_brands)) %>%
    mutate(across(subgroup, factor, levels = model_subgroups)) %>%
    mutate(across(outcome, factor, levels = model_outcomes)) 
  
}


## IPW ----
forest_from_broomstack <- function(broomstack, title){
  
  #jtools::plot_summs(ipwvaxany1)
  #modelsummary::modelplot(ipwvaxany1, coef_omit = 'Interc|tstop', conf.type="wald", exponentiate=TRUE)
  #sjPlot::plot_model(ipwvaxany1)
  #all these methods use broom::tidy to get the coefficients. but tidy.glm only uses profile CIs, not Wald. (yTHO??)
  #profile CIs will take forever on large datasets.
  #so need to write custom function for plotting wald CIs. grr
  
  
  plot_data <- broomstack %>%
    filter(
      !str_detect(variable,fixed("ns(tstop")),
      !str_detect(variable,fixed("region")),
      !is.na(term)
    ) %>%
    mutate(
      var_label = if_else(var_class=="integer", "", var_label),
      label = if_else(reference_row %in% TRUE, paste0(label, " (ref)"),label),
      or = if_else(reference_row %in% TRUE, 1, or),
      variable = fct_inorder(variable),
      variable_card = as.numeric(variable)%%2,
    ) %>%
    group_by(variable) %>%
    mutate(
      variable_card = if_else(row_number()!=1, 0, variable_card),
      level = fct_rev(fct_inorder(paste(variable, label, sep="__"))),
      level_label = label
    ) %>%
    ungroup() %>%
    droplevels()
  
  var_lookup <- str_wrap(plot_data$var_label, 20)
  names(var_lookup) <- plot_data$variable
  
  level_lookup <- plot_data$level
  names(level_lookup) <- plot_data$level_label
  
  ggplot(plot_data) +
    geom_point(aes(x=or, y=level, colour=subgroup_level), position = position_dodge(width = 0.3)) +
    geom_linerange(aes(xmin=or.ll, xmax=or.ul, y=level, colour=subgroup_level), position = position_dodge(width = 0.3)) +
    geom_vline(aes(xintercept=1), colour='black', alpha=0.8)+
    facet_grid(
      rows=vars(variable),
      cols=vars(brand_descr),
      scales="free_y", switch="y", space="free_y", labeller = labeller(variable = var_lookup)
    )+
    scale_x_log10(
      breaks=c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
      limits=c(0.01, 3)
    )+
    scale_y_discrete(breaks=level_lookup, labels=names(level_lookup))+
    geom_rect(aes(alpha = variable_card), xmin = -Inf,xmax = Inf, ymin = -Inf, ymax = Inf, fill='grey', colour="transparent") +
    scale_alpha_continuous(range=c(0,0.3), guide=FALSE)+
    labs(
      y="",
      x="Hazard ratio",
      colour=NULL,
      title=title
      #subtitle=cohort_descr
    ) +
    theme_minimal() +
    theme(
      strip.placement = "outside",
      strip.background = element_rect(fill="transparent", colour="transparent"),
      strip.text.y.left = element_text(angle = 0, hjust=1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.spacing = unit(0, "lines"),
      legend.position = "bottom"
    )
}

# broomstack
broomstack <- model_key %>%
  filter(outcome=="postest", brand=="pfizer") %>%
  mutate(
    subgroup_level =  pmap(
      list(.$subgroup), 
      function(subgroup) unname(recoder[[subgroup]])
      )
  ) %>%
  unnest(subgroup_level) %>%
  mutate(
    broom = pmap(
      ., 
      ~read_rds(file.path(indir, brand, subgroup, outcome, "postprocess", glue("broom_vax{brand}1_{subgroup_level}.rds"))) %>%
        select(-subgroup_level)
    )
  ) %>%
  unnest(broom) %>%
  add_descr()
  

broomstack_formatted <- broomstack %>%
  transmute(
    subgroup_level,
    brand_descr,
    var_label,
    label,
    HR = scales::label_number(accuracy = .01, trim=TRUE)(or),
    HR = if_else(is.na(HR), "1", HR),
    CI = paste0("(", scales::label_number(accuracy = .01, trim=TRUE)(or.ll), "-", scales::label_number(accuracy = .01, trim=TRUE)(or.ul), ")"),
    CI = if_else(is.na(HR), "", CI),
    HR_ECI = paste0(HR, " ", CI)
  )

broomstack_formatted_wide <- broomstack_formatted %>%
  select(
    subgroup_level, brand_descr, var_label, label, HR_ECI
  ) %>%
  pivot_wider(
    id_cols = c(subgroup_level, var_label, label),
    names_from = brand_descr,
    values_from = HR_ECI,
    names_glue = "{brand_descr}_{.value}"
  )

write_csv(broomstack_formatted, file.path(outdir, "tab_vax1.csv"))
write_csv(broomstack_formatted_wide, file.path(outdir, "tab_vax1_wide.csv"))

# plot broomstack
plot_vax <- forest_from_broomstack(broomstack, "Vaccination model")
ggsave(
  file.path(outdir, "plot_vax1.svg"),
  plot_vax,
  units="cm", width=30, height=25
)

## MSM ----

estimates <- model_key %>%
  filter(outcome=="postest", brand=="pfizer") %>%
  mutate(
    estimates = pmap(
      .,
      ~read_csv(file.path(indir, brand, subgroup, outcome, "postprocess", glue("estimates_timesincevax.csv"))))
  ) %>%
  unnest(estimates) %>%
  mutate(
    outcome = factor(outcome, levels = model_outcomes),
  ) %>%
  add_descr() 
  
  # metadata_outcomes %>%
  # mutate(
  #   outcome = fct_inorder(outcome),
  #   outcome_descr = fct_inorder(map_chr(outcome_descr, ~paste(stringi::stri_wrap(., width=14, simplify=TRUE, whitespace_only=TRUE), collapse="\n")))
  # ) %>%
  # crossing(
  #   tibble(
  #     brand = fct_inorder(c("any", "pfizer", "az")),
  #     brand_descr = fct_inorder(c("Any vaccine", "BNT162b2", "ChAdOx1"))
  #   )
  # ) %>%
  # mutate(
  #   brand = fct_inorder(brand),
  #   brand_descr = fct_inorder(brand_descr),
  #   estimates = map2(brand, outcome, ~read_csv(here("output", cohort, strata_var, recent_postestperiod, .x, .y, glue("estimates_timesincevax.csv"))))
  # ) %>%
  # unnest(estimates) %>%
  # mutate(
  #   model_descr = fct_inorder(model_descr),
  # )

estimates_formatted <- estimates %>%
  transmute(
    outcome_descr,
    brand_descr,
    subgroup_descr,
    subgroup_level,
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

# if(subgroup!="all"){
#   estimates_formatted <- estimates_formatted %>%
#     mutate(
#       diff = scales::label_number(accuracy = .01, trim=TRUE)(diff),
#       diff_CI = paste0("(", scales::label_number(accuracy = .01, trim=TRUE)(diff.ll), "-", scales::label_number(accuracy = .01, trim=TRUE)(diff.ul), ")"),
#       diff_ECI = paste0(diff, " ", diff_CI),
#     )
# }

estimates_formatted_wide <- estimates_formatted %>%
  select(outcome_descr, brand_descr, subgroup_descr, model, term, HR_ECI, VE_ECI) %>%
  pivot_wider(
    id_cols=c(outcome_descr, brand_descr, term, subgroup_descr),
    names_from = model,
    values_from = c(HR_ECI, VE_ECI),
    names_glue = "{model}_{.value}"
  )

write_csv(estimates, file.path(outdir, glue("estimates_timesincevax.csv")))
write_csv(estimates_formatted, file.path(outdir, glue("estimates_formatted_timesincevax.csv")))
write_csv(estimates_formatted_wide, file.path(outdir, glue("estimates_formatted_wide_timesincevax.csv")))

# create forest plot for each subgroup_level
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

msmmod_effect_plot <- function(subgroup_level, scales_option) {
  
  msmmod_effect_data_plot <- msmmod_effect_data %>%
    filter(
      or !=Inf, or !=0, 
      subgroup_level %in% subgroup_level
      ) %>%
    mutate(across(model_descr, ~str_wrap(.x, width = 120)))
  
  plot_function <- function(scales_option) {
    
    if (scales_option == "fixed") {
      
      y_intercept <- 1
      
      scale_y_custom <- function() {
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
          )
      }
      
    }
    
    if (scales_option == "free_y") {
      
      msmmod_effect_data_plot <- msmmod_effect_data_plot %>%
        mutate(across(c(or, or.ll, or.ul), log))
      
      y_intercept <- 0
      
      scale_y_custom <- function() {
        scale_y_continuous(
          labels = function(x){scales::label_number(0.001)(exp(x))},
          breaks = function(x){log(scales::breaks_log(n=6, base=10)(exp(x)))},
          sec.axis = sec_axis(
            ~(1-exp(.)),
            name="Effectiveness",
            breaks = function(x){1-(scales::breaks_log(n=6, base=10)(1-x))},
            labels = function(x){formatpercent100(x, 1)}
          )
        )
      }
      
    }
    
    msmmod_effect <-
      msmmod_effect_data_plot %>%
      ggplot(aes(colour = model_descr)) +
      geom_hline(aes(yintercept = y_intercept), colour='grey') +
      geom_point(
        aes(y = or, x = term_midpoint), 
        position = position_dodge(width = 1.5), 
        size = 0.8
        ) +
      geom_linerange(
        aes(ymin = or.ll, ymax = or.ul, x = term_midpoint),
        position = position_dodge(width = 1.5)
        ) +
      facet_grid(
        rows = vars(outcome_descr), 
        cols = vars(brand_descr), 
        switch = "y",
        scales = scales_option
        ) +
      scale_y_custom() +
      scale_x_continuous(
        breaks=c(0, postbaselinecuts), 
        expand=expansion(mult=c(0, NA)),
        limits=c(0,max(postbaselinecuts))
      ) +
      scale_colour_brewer(
        type="qual", 
        palette="Set2", 
        guide=guide_legend(ncol=1)
        ) +
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
        
        panel.spacing = unit(1, "lines"),
        
        plot.title = element_text(hjust = 0),
        plot.title.position = "plot",
        plot.caption.position = "plot",
        plot.caption = element_text(hjust = 0, face= "italic"),
        
        legend.position = "bottom"
      )
    
    ## save plot
    ggsave(
      filename=file.path(outdir, glue("VE_plot_{subgroup_level}_{scales_option}.svg")), 
      msmmod_effect, 
      width=20, height=20, units="cm"
    )
    
    return(msmmod_effect)
    
  }
  
  plot_function(scales_option)
  
}

for (subgroup_level in unique(msmmod_effect_data$subgroup_level)) {
  msmmod_effect_plot(subgroup_level, "fixed")
  msmmod_effect_plot(subgroup_level, "free_y")
}


  
