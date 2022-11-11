
library('tidyverse')
library('here')
library('glue')
library('lubridate')
library('gt')
library('patchwork')
library('scales')


## Import design elements
source(here("analysis", "design.R"))

source(here("lib", "functions", "utility.R"))


# where are the outputs (ie the inputs for this manuscript!) saved?
output_dir_os <- here("release")



# where should we put the objects created within this rmd script?
output_dir_rmd <- here("write-up", "figures")
fs::dir_create(output_dir_rmd)


## sequential trial outputs ----
# from:  https://jobs.opensafely.org/datalab/covid-19-vaccine-effectiveness/covid-vaccine-effectiveness-seqtrial-main/


brands <- tibble(brand=c("pfizer", "az")) %>% group_by(brand)

# wrap outcome descriptions over two lines for nicer plots
new_names <- sapply(
  names(recoder$outcome),
  function(x)
    if_else(nchar(x) > 10, str_replace(x, " ", "\\\n"), x)
)

names(recoder$outcome) <- new_names

## KM estimates and contrasts are commented out, as not used for now

# km_estimates <- brands %>%
#   mutate(
#     dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "models", "km", "combined", "km_estimates_rounded.csv")))
#   ) %>% 
#   unnest(dat) %>%
#   mutate(
#     treated_descr = fct_recoderelevel(as.character(treated),  recoder$treated),
#     outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
#     subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
#     subgroup_level_descr = replace_na(subgroup_level_descr, "")
#   ) %>%
#   group_by(treated_descr, outcome_descr, subgroup_descr, subgroup_level_descr) %>%
#   mutate(
#     lagrisk=lag(risk,1,0),
#     lagsurv=lag(surv,1,1),
#   )
# 
# 
# contrasts_km_daily <- brands %>%
#   mutate(
#     dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "models", "km", "combined", "contrasts_km_daily_rounded.csv")))
#   ) %>% 
#   unnest(dat) %>%
#   mutate(
#     outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
#     subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
#     subgroup_level_descr = replace_na(subgroup_level_descr, "")
#   )
# 
# 
# contrasts_km_cuts <- brands %>%
#   mutate(
#     dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "models", "km", "combined", "contrasts_km_cuts_rounded.csv")))
#   ) %>% 
#   unnest(dat) %>%
#   mutate(
#     outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
#     subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
#     subgroup_level_descr = replace_na(subgroup_level_descr, "")
#   )
# 
# contrasts_km_overall <- 
#   brands %>%
#   mutate(
#     dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "models", "km", "combined", "contrasts_km_overall_rounded.csv")))
#   ) %>% 
#   unnest(dat) %>%
#   mutate(
#     outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
#     subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
#     subgroup_level_descr = replace_na(subgroup_level_descr, ""),
#   )



st_cox_cuts <- brands %>%
  mutate(
    dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "contrasts_cox_cuts.csv")))
  ) %>% 
  unnest(dat) %>%
  mutate(
    outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
    subgroup_level_descr = replace_na(subgroup_level_descr, "")
  )

st_cox_overall <- brands %>%
  mutate(
    dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "contrasts_cox_overall.csv")))
  ) %>% 
  unnest(dat) %>%
  mutate(
    outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
    subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
    subgroup_level_descr = replace_na(subgroup_level_descr, ""),
  )


## meta-analyse across age groups


st_cox_overall_MA <-
  st_cox_overall %>%
  filter(subgroup=="ageband2") %>%
  group_by(brand, outcome_descr, period_start, period_end) %>%
  mutate(loghr = log(coxhazr)) %>%
  summarise(
    # NOTE: coxhr.se is the standard error of the log hazard ratio, not the hazard ratio
    loghr = weighted.mean(loghr, coxhr.se^-2),
    coxhr = exp(loghr),
    loghr.se = sqrt(1/sum(coxhr.se^-2)),
    statistic = loghr/loghr.se,
    p.value = pchisq(statistic^2, df=1, lower.tail=FALSE),
    conf.low = exp(loghr + qnorm(0.025)*loghr.se),
    conf.high = exp(loghr + qnorm(0.975)*loghr.se),
  )


st_cox_cuts_MA <-
  st_cox_cuts %>%
  filter(subgroup=="ageband2") %>%
  group_by(brand, outcome_descr,  outcome,period_start, period_end) %>%
  mutate(
    loghr = log(coxhazr)
  ) %>%
  summarise(
    # NOTE: coxhr.se is the standard error of the log hazard ratio, not the hazard ratio
    loghr = weighted.mean(loghr, coxhr.se^-2),
    hr = exp(loghr),
    loghr.se = sqrt(1/sum(coxhr.se^-2)),
    statistic = loghr/loghr.se,
    p.value = pchisq(statistic^2, df=1, lower.tail=FALSE),
    conf.low = exp(loghr + qnorm(0.025)*loghr.se),
    conf.high = exp(loghr + qnorm(0.975)*loghr.se),
    approach="st",
  )
  


## MSM outputs ----

# from: https://jobs.opensafely.org/datalab/covid-19-vaccine-effectiveness/covid-vaccine-effectiveness-research_positivity/



#msm_cox_cuts <- read_csv(fs::path(output_dir_os, "msm", "msmvstcox_estimates_timesincevax.csv")) %>%
#filter(approach=="msm", model==4) %>%
# mutate(
#   subgroup="ageband2",
#   subgroup_level = case_when(
#     cohort == "over80s" ~ "80+",
#     cohort == "in70s" ~ "70-79",
#   ),
#   subgroup_level_descr = case_when(
#     cohort == "over80s" ~ "aged 80+",
#     cohort == "in70s" ~ "aged 70-79",
#   ),
#   period_start = as.integer(str_extract(term, "^\\d+"))-1,
#   period_end = as.integer(str_extract(term, "\\d+$")),
# ) %>%
#   mutate(
#     outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
#     subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
#     subgroup_level_descr = replace_na(subgroup_level_descr, "")
#   )

# msm_cox_cuts_MA <- 
#   msm_cox_cuts %>%
#   filter(subgroup=="ageband2") %>%
#   group_by(brand, outcome_descr, outcome, period_start, period_end) %>%
#   mutate(
#     loghr = log(or),
#     period_end = if_else(is.na(period_end), 63L, period_end),
#   ) %>%
#   summarise(
#     # NOTE: std.error is the standard error of the log hazard ratio, not the hazard ratio
#     loghr = weighted.mean(loghr, std.error^-2),
#     hr = exp(loghr),
#     loghr.se = sqrt(1/sum(std.error^-2)),
#     statistic = loghr/loghr.se,
#     p.value = pchisq(statistic^2, df=1, lower.tail=FALSE),
#     conf.low = exp(loghr + qnorm(0.025)*loghr.se),
#     conf.high = exp(loghr + qnorm(0.975)*loghr.se),
#     approach="msm",
#   ) %>%
#   ungroup()


msm_cox_cuts_MA <- read_csv(fs::path(output_dir_os, "msm", "meta_estimates.csv")) %>%
  filter(recent_postestperiod == Inf, model==4, brand != "any") %>%  
  mutate(
    # subgroup="ageband2",
    # subgroup_level = case_when(
    #   cohort == "over80s" ~ "80+",
    #   cohort == "in70s" ~ "70-79",
    # ),
    # subgroup_level_descr = case_when(
    #   cohort == "over80s" ~ "aged 80+",
    #   cohort == "in70s" ~ "aged 70-79",
    # ),
    period_start = as.integer(str_extract(term, "^\\d+"))-1,
    period_end = as.integer(str_extract(term, "\\d+$")),
  ) %>%
  mutate(
    outcome_descr = fct_recoderelevel(as.character(outcome),  recoder$outcome),
    #subgroup_descr = fct_recoderelevel(subgroup,  recoder$subgroup),
    #subgroup_level_descr = replace_na(subgroup_level_descr, "")
    period_end = if_else(is.na(period_end), 63L, period_end),
    hr=or,
    conf.low=or.ll,
    conf.high=or.ul,
    approach="msm"
  )




## combine and plot ----

cox_cuts_MA <-
  bind_rows(
    msm_cox_cuts_MA,
    st_cox_cuts_MA
  ) %>%
  filter(
    outcome %in% c("postest", "covidadmitted", "death")
  ) %>%
  mutate(
    period_end = if_else(period_end>35, 63, period_end),
    midpoint = (period_end+period_start)/2,
  ) %>%
  mutate(across(approach, factor, levels = c("msm", "st"), labels = c("Single trial", "Sequential trial"))) %>%
  ungroup() %>%
  mutate(across(brand, factor, levels = c("az", "pfizer"), labels = c("ChAdOx1", "BNT162b2"))) %>%
  arrange(
    approach, brand, outcome_descr, period_start
  )

## table of hazard ratios ----
doc <- officer::read_docx() 

ftab <- cox_cuts_MA %>%
  transmute(
    approach,
    Brand = brand, 
    Outcome = str_replace(as.character(outcome_descr), "\\n", " "),
    `Days since first dose` = if_else(period_end>35, glue("{period_start+1}+*"), glue("{period_start+1}-{period_end}")),
    value = glue("{format(round(hr,3),nsmall=3)} ({format(round(conf.low,3),nsmall=3)}, {format(round(conf.high,3),nsmall=3)})")
    ) %>%
  pivot_wider(
    names_from = approach,
    values_from = value
  ) %>%
  flextable::flextable() %>%
  flextable::merge_v(j=1:2, part = "body") %>%
  flextable::fontsize(size = 8, part = "all") %>%
  flextable::theme_booktabs() 

doc <- doc %>%
  flextable::body_add_flextable(value = ftab, split = FALSE) %>%
  print(target = file.path(output_dir_os, glue::glue("hr_table.docx")))
# add footnote manually:
# *final period is 36-xx days for marginal structural model and 35-70 days for sequential trial??

## effect plot ----
formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

effect_plot <-
  cox_cuts_MA %>%
  ggplot(aes(colour=approach)) +
  geom_hline(aes(yintercept=1), colour='black')+
  geom_vline(aes(xintercept=0), colour='black')+
  geom_point(aes(y=hr, x=midpoint), position = position_dodge(width = 2), size=0.8)+
  geom_linerange(aes(ymin=conf.low, ymax=conf.high, x=midpoint), position = position_dodge(width = 2))+
  facet_grid(rows=vars(outcome_descr), cols=vars(brand), switch="y")+
  scale_y_log10(
    breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5),
    limits = c(0.05, max(c(1, cox_cuts_MA$conf.high))),
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name="Effectiveness",
      breaks = c(-4, -1, 0, 0.5, 0.80, 0.9, 0.95, 0.98, 0.99),
      labels = function(x){formatpercent100(x, 1)}
    )
  )+
  scale_x_continuous(
    breaks = c(0,postbaselinecuts),
    expand = expansion(mult=c(0), add=c(0,7)), limits=c(0,NA)
  )+
  scale_colour_brewer(type="qual", palette="Set2", guide=guide_legend(ncol=2))+
  coord_cartesian() +
  labs(
    y="Hazard ratio, versus no vaccination",
    x="Days since first dose",
    colour=NULL
  ) +
  theme_bw(base_size=12)+
  theme(
    panel.border = element_blank(),
    
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
    
    axis.title.y.right = element_text(margin = margin(t=0,r=0,b=0,l=10)),
    axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
    
    legend.position = "bottom"
  )

effect_plot
ggsave(
  file.path(output_dir_os, "ve-plot.png"),
  effect_plot,
  width=20, height=15, units="cm"
)

## coverage plots ----

data_coverage <- brands %>%
  mutate(
    dat = map(brand, ~read_csv(fs::path(output_dir_os, "st", brand, "coverage.csv")))
  ) %>% 
  unnest(dat) %>%
  ungroup() %>%
  mutate(across(brand, factor)) %>%
  mutate(across(status, factor, levels = c("unmatched", "matched")))

xmin <- min(data_coverage$vax1_date )
xmax <- max(data_coverage$vax1_date )+1

# daily numbers
plot_coverage_n <-
  data_coverage %>%
  mutate(
    brand = case_when(
      brand == "pfizer" ~ "BNT162b2",
      brand == "az" ~ "ChAdOx1",
      TRUE ~ NA_character_
    ),
  ) %>%
  ggplot()+
  geom_col(
    aes(
      x=vax1_date+0.5,
      y=n,
      group=paste0(brand,status),
      fill=brand,
      alpha=fct_rev(status),
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    #alpha=0.8,
    width=1
  )+
  #geom_rect(xmin=xmin, xmax= xmax+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  geom_hline(yintercept = 0, colour="black")+
  facet_grid(rows = "brand")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax1_date, "1 month")),
    limits = c(xmin-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.x)),
    expand = expansion(c(0, NA))
  )+
  scale_fill_brewer(type="qual", palette="Set2")+
  scale_colour_brewer(type="qual", palette="Set2")+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Count",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_bw()+
  theme(
    
    panel.border = element_blank(),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    
    panel.spacing = unit(1, "lines"),
    
    axis.title.y = element_text(margin = margin(t=0,r=10,b=0,l=0)),
    axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
    
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    strip.text.y.right = element_text(angle = 0),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )

ggsave(
  filename=file.path(output_dir_os, "coverage_n.png"),
  plot_coverage_n,
  width=15, height=20, units="cm"
)  

# cumulative numbers
plot_coverage_cumuln <-
  data_coverage %>%
  mutate(
    brand = case_when(
      brand == "pfizer" ~ "BNT162b2",
      brand == "az" ~ "ChAdOx1",
      TRUE ~ NA_character_
    ),
  ) %>%
  ggplot()+
  geom_col(
    aes(
      x=vax1_date+0.5,
      y=cumuln,
      group=paste0(brand,status),
      fill=brand,
      alpha=fct_rev(status),
      colour=NULL
    ),
    position=position_stack(reverse=TRUE),
    width=1
  )+
  geom_rect(xmin=xmin, xmax= xmax+1, ymin=-6, ymax=6, fill="grey", colour="transparent")+
  facet_grid(rows = "brand")+
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(data_coverage$vax1_date, "1 month")),
    limits = c(xmin-1, NA),
    labels = scales::label_date("%d/%m"),
    expand = expansion(add=1),
  )+
  scale_y_continuous(
    labels = ~scales::label_number(accuracy = 1, big.mark=",")(abs(.)),
    expand = expansion(c(0, NA))
  )+
  #scale_fill_brewer(type="qual", palette="Set2")+
  #scale_colour_brewer(type="qual", palette="Set2")+
  scale_fill_manual(values = brewer_pal(type="qual", palette="Set2")(3)[c(3,1)])+
  scale_colour_manual(values = brewer_pal(type="qual", palette="Set2")(3)[c(3,1)])+
  scale_alpha_discrete(range= c(0.8,0.4))+
  labs(
    x="Date",
    y="Cumulative count",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_bw()+
  theme(
    
    panel.border = element_blank(),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    
    panel.spacing = unit(1, "lines"),
    
    axis.title.y = element_text(margin = margin(t=0,r=10,b=0,l=0)),
    axis.title.x = element_text(margin = margin(t=10,r=0,b=0,l=0)),
    
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    #strip.text.y.right = element_text(angle = 0),
    strip.text.y.right = element_blank(),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )

ggsave(
  filename=file.path(output_dir_os, "coverage_cumuln.png"),
  plot_coverage_cumuln,
  width=15, height=20, units="cm"
)  
