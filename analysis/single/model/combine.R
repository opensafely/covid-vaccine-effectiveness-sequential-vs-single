# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# 
# imports fitted IPW models from `msm.R`
# calculates robust CIs taking into account patient-level clustering
# 
# imports fitted MSMs from `msm.R`
# calculates robust CIs taking into account patient-level clustering
# outputs plots for the primary vaccine-outcome relationship
# outputs plots showing model-estimated spatio-temporal trends
#
# The script should only be run via an action in the project.yaml only
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Preliminaries ----

# import libraries
library('tidyverse')
library('here')
library('glue')
library('lubridate')
library('survival')
library('splines')
library('parglm')
library("sandwich")
library("lmtest")

# import custom user functions and metadata
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))
source(here("analysis", "functions", "redaction.R"))
source(here("analysis", "functions", "survival.R"))

# create output directory 
indir <- here("output", "single") 
outdir <- here("output", "single", "combine") 
fs::dir_create(outdir)

# define metaparams
metaparams <- expand_grid(
  brand=model_brands,
  subgroup=model_subgroups,
  outcome=model_outcomes
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# IPW ----

# define forest_from_broomstack function
forest_from_broomstack <- function(broomstack, title){
  
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
  
  var_lookup <- str_trunc(plot_data$var_label, width = 20, side = "right")
  names(var_lookup) <- plot_data$variable
  
  level_lookup <- plot_data$level
  names(level_lookup) <- plot_data$level_label
  
  ggplot(plot_data) +
    geom_point(
      aes(x=or, y=level, colour=outcome),
      position = position_dodge(width = 0.3)
      ) +
    geom_linerange(
      aes(xmin=or.ll, xmax=or.ul, y=level, colour=outcome),
      position = position_dodge(width = 0.3)
      ) +
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

# generate broomstack ----
broomstack <- metaparams %>%
  mutate(
    subgroup_level =  pmap(
      list(.$subgroup), 
      function(subgroup) unname(recoder[[subgroup]])
      )
  ) %>%
  unnest(subgroup_level) %>%
  pmap(
    function(brand, subgroup, outcome, subgroup_level) 
      read_rds(file.path(indir, brand, subgroup, outcome, "postprocess", glue("broom_vax{brand}1_{subgroup_level}.rds"))) %>%
      add_column(brand=brand, subgroup=subgroup, outcome=outcome, .before=1)
  ) %>%
  bind_rows() 
  
# plot broomstack
plot_vax <- forest_from_broomstack(
  broomstack %>% add_descr(), 
  "Vaccination model"
  )
ggsave(
  file.path(outdir, "plot_vax1.svg"),
  plot_vax,
  units="cm", width=30, height=25
)

# selected required columns and save
broomstack %>%
  select(
    brand, subgroup, subgroup_level, outcome,
    var_label, label, starts_with("or")
    ) %>%
  write_csv(file.path(outdir, "tab_vax1.csv"))
  
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
## MSM ----

# combine all estimates
estimates <- metaparams %>%
  pmap(
    function(brand, subgroup, outcome)
      read_csv(file.path(indir, brand, subgroup, outcome, "postprocess", "estimates_timesincevax.csv")) %>%
      add_column(subgroup=subgroup, brand=brand, outcome=outcome, .before=1)
      ) %>%
  bind_rows() %>%
  add_descr() 

# only save the required columns
estimates %>%
  select(brand, subgroup, subgroup_level, model, term, outcome, starts_with("or")) %>%
  write_csv(file.path(outdir, "estimates_timesincevax.csv"))
  
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
