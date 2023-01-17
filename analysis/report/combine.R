# table of summary stats and coefficients from vaccination model

library(tidyverse)
library(here)
library(glue)

outdir <- here("output", "release") 
fs::dir_create(outdir)

# import custom user functions and metadata
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# cohort characteristics (single and sequential)

table1 <- bind_rows(
  read_csv(here("output", "report", "table1", "table1_sequential_pfizer_rounded.csv")) %>%
    mutate(across(by, factor, levels = c(0,1), labels = c("unvaccinated", "pfizer"))) %>%
    add_column(brand="pfizer", .before=1),
  read_csv(here("output", "report", "table1", "table1_sequential_az_rounded.csv")) %>% 
    mutate(across(by, factor, levels = c(0,1), labels = c("unvaccinated", "az"))) %>%
    add_column(brand = "az", .before = 1),
  read_csv(here("output", "report", "table1", "table1_single_any_rounded.csv"))
) %>%
  filter(variable != "age_factor") %>%
  rowwise() %>%
  mutate(
    across(
      c(p, p_miss, p_nonmiss),
      ~if_else(.x < 10,
               format(round(100*.x, 1), nsmall=1),
               format(round(100*.x, 0), nsmall=0)
      )
    )
  ) %>%
  mutate(across(stat_display, glue)) %>%
  select(brand, by, variable, variable_levels, stat_display) %>%
  pivot_wider(
    names_from = c("brand", "by"),
    values_from = stat_display
  ) %>%
  rename(
    pfizer = pfizer_pfizer,
    az = az_az,
    single = NA_single
  )

# file to release (table for manuscript)
write_csv(
  table1,
  file.path(outdir, "table1_rounded.csv")
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# vax1 model coefficients (single only)

# file to release
tab_vax1 <- read_csv(here("output", "single", "combine", "tab_vax1.csv"))

# post release processing  (table for manuscript)
tab_vax1_wide <- tab_vax1 %>%
  filter(
    subgroup=="all", 
    # use death as the outcome, as it includes all uncensored follow-up time
    outcome=="death" 
  ) %>%
  add_descr() %>%
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
  ) %>%
  select(
    -HR, -CI
    # subgroup_level, brand_descr, var_label, label, HR_ECI
  ) %>%
  pivot_wider(
    id_cols = c(subgroup_level, var_label, label),
    names_from = brand_descr,
    values_from = HR_ECI,
    names_glue = "{brand_descr}"
  )

# write_csv(here(outdir, "tab_vax1_wide.csv"))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# outcome model coefficients (single and sequential)

# import sequential
estimates_sequential <- read_csv(here("output", "sequential", "combine", "contrasts_cox_cuts.csv")) %>%
  transmute(
    approach = "Sequential trial",
    brand, subgroup, subgroup_level, outcome, 
    estimate = coxhazr,
    lci = coxhr.ll,
    uci = coxhr.ul,
    period_start,
    period_end,
    fup_period
  )

# import single
estimates_single <- read_csv(here("output", "single", "combine", "estimates_timesincevax.csv")) %>%
  transmute(
    approach = "Single trial",
    brand, subgroup, subgroup_level, outcome, 
    estimate = or,
    lci = or.ll,
    uci = or.ul,
    period_start = as.numeric(str_extract_all(term, "^\\d+")) - 1,
    period_end = as.numeric(str_extract_all(term, "\\d+$")),
    fup_period = str_c(period_start, period_end, sep = "-")
  ) 

# create table (supplementary table xxx)
bind_rows(
  estimates_sequential,
  estimates_single
) %>%
  add_descr() %>%
  mutate(across(fup_period, factor, levels = paste0(c(0,postbaselinecuts[-length(postbaselinecuts)]), "-", postbaselinecuts))) %>%
  mutate(
    estimate = scales::label_number(accuracy = .01, trim=TRUE)(estimate),
    ci = paste0("(", scales::label_number(accuracy = .01, trim=TRUE)(lci), "-", scales::label_number(accuracy = .01, trim=TRUE)(uci), ")")
  ) %>%
  transmute(
    brand_descr, outcome_descr, approach, fup_period,
    value = paste0(estimate, " ", ci),
  ) %>%
  pivot_wider(
    names_from = approach,
    values_from = value
  ) %>% 
  arrange(outcome_descr, brand_descr, fup_period) %>%
  print(n=Inf)

# figure xxx

plot_data <- bind_rows(
  estimates_sequential,
  estimates_single
) %>%
  add_descr() %>%
  rowwise() %>%
  mutate(mid_point = (period_start + period_end)/2) 


position_dodge_val <- 1

primary_vax_y1 <- lst(
  breaks = c(0.05, 0.1, 0.2, 0.5, 1), 
  limits = c(min(breaks), max(breaks))
)
primary_vax_y2 <- list(
  breaks = c(0, 0.5, 0.8, 0.9, 0.95)
)

shape_palette <- c(17,16)
names(shape_palette) <- brand_lookup$brand_descr[1:2]

outcome_descr_long <- levels(plot_data$outcome_descr)
outcome_descr_wrap <- str_replace(outcome_descr_long, "\\s", "\\\n")

plot_data %>%
  mutate(
    across(
      outcome_descr, 
      factor, 
      levels = outcome_descr_long,
      labels = outcome_descr_wrap
    )
  ) %>%
  ggplot(aes(
    x = mid_point, 
    colour = approach, 
    shape = approach
  )) +
  geom_hline(aes(yintercept=1), colour='grey') +
  geom_linerange(
    aes(ymin = lci, ymax = uci),
    position = position_dodge(width = position_dodge_val)
  ) +
  geom_point(
    aes(y = estimate),
    alpha = 0.8,
    position = position_dodge(width = position_dodge_val)
  ) +
  facet_grid(
    rows = vars(outcome_descr), cols = vars(brand_descr), 
    switch = "y", scales = "free", space = "free_x"
  ) +
  scale_x_continuous(
    breaks = c(0, postbaselinecuts),
    labels = c(0, postbaselinecuts),
    limits = c(0, max(postbaselinecuts))
  ) +
  scale_y_log10(
    name = "Hazard ratio, versus no vaccination",
    breaks = primary_vax_y1[["breaks"]],
    limits = primary_vax_y1[["limits"]],
    oob = scales::oob_keep,
    sec.axis = sec_axis(
      ~(1-.),
      name= "Vaccine effectiveness",
      breaks = primary_vax_y2[["breaks"]],
      labels = function(x){formatpercent100(x, 1)}
    )
  ) +
  labs(
    x = "Days since first dose"
  ) +
  scale_color_brewer(palette = "Dark2") +
  scale_shape(guide = "none") +
  guides(
    color = guide_legend(
      title = NULL,
      override.aes = list(shape = shape_palette[1:2])
    )) +
  theme_bw() +
  theme(
    
    panel.border = element_blank(),
    axis.line.y = element_line(colour = "black"),
    
    axis.text = element_text(size=10),
    
    axis.title.x = element_text(size=10, margin = margin(t = 10)),
    axis.title.y = element_text(size=10, margin = margin(r = 10)),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    strip.text = element_text(size=10),
    
    panel.spacing = unit(0.8, "lines"),
    
    plot.title = element_text(hjust = 0),
    plot.title.position = "plot",
    plot.caption.position = "plot",
    plot.caption = element_text(hjust = 0, face= "italic"),
    
    legend.position = "bottom",
    legend.text = element_text(size=10)
    
  ) 

# ggsave()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# sequential matching coverage (sequential only)

coverage <- bind_rows(
  read_csv(here("output", "report", "coverage", "coverage_pfizer.csv")) %>%
    mutate(brand="pfizer"),
  read_csv(here("output", "report", "coverage", "coverage_az.csv")) %>%
    mutate(brand = "az")
) 

# to release
write_csv(
  coverage,
  file.path(outdir, "coverage_rounded.csv")
)

# figure (supplementary figure xxx)

colour_palette <- c(
  "BNT162b2, matched" = "#e7298a", # dark pink / dark grey
  "BNT162b2, unmatched" = "#e78ac3", # medium pink / medium grey
  "ChAdOx1, matched" = "#7570b3", # dark purple / dark grey
  "ChAdOx1, unmatched" = "#8da0cb" # medium purple / medium grey
)

# this is necessary because there is an older version of a package in opensafely 
# and I think it requires breaks to be unique
# if(Sys.getenv("OPENSAFELY_BACKEND") %in% c("")) {
#   y_labels <- ~scales::label_number(accuracy = 1, big.mark=",")(.x)
# } else {
#   y_labels <- waiver()
# }

plot_coverage_cumuln <-
  coverage %>%
  mutate(brand_descr = factor(brand, levels = brand_lookup$brand, labels = brand_lookup$brand_descr)) %>%
  mutate(
    colour_var = factor(
      paste0(brand_descr, ", ", status),
      levels = names(colour_palette)
    )
  ) %>%
  droplevels() %>%
  ggplot() +
  geom_col(
    aes(
      x = vax1_date,
      y = cumuln,
      group = colour_var,
      fill = colour_var,
      colour = NULL
    ),
    position = position_stack(reverse=TRUE),
    width = 1
  ) +
  facet_wrap(
    facets = vars(brand_descr),
    nrow = 2
  ) +
  # facet_grid(
  #   rows = vars(brand_descr)
  # ) +
  scale_x_date(
    breaks = unique(lubridate::ceiling_date(coverage$vax1_date, "1 month")),
    # limits = c(xmin-1, NA),
    labels = scales::label_date("%b %Y"),
    expand = expansion(add=7),
  ) +
  scale_y_continuous(
    labels = ~scales::label_number(accuracy = 1, big.mark=",")(.x),
    expand = expansion(c(0, NA))
  ) +
  scale_fill_manual(values = colour_palette) +
  labs(
    x = "Date",
    y = "Cumulative count per day",
    colour = NULL,
    fill = NULL,
    alpha = NULL
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size=10, margin = margin(t = 10)),
    axis.title.y = element_text(size=10, margin = margin(r = 10)),
    axis.line.x.bottom = element_line(),
    axis.text.x.top=element_text(hjust=0),
    # strip.text.y.right = element_text(angle = 90),
    strip.text = element_text(hjust=0),
    axis.ticks.x=element_line(),
    legend.position = "bottom"
  )

# ggsave()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# sequential KM cumulative incidence (sequential only)

km_estimates_rounded <- read_csv(here("output", "sequential", "combine", "km_estimates_rounded.csv"))

# supplementary figure xxx

colour_palette <- c(
  "Unvaccinated" = "#616161", # light grey
  "BNT162b2" = "#e7298a", # dark pink / dark grey
  "ChAdOx1" = "#7570b3" # dark purple / dark grey
)

linetype_palette <- c(
  "Unvaccinated" = "dashed",
  "BNT162b2" = "solid",
  "ChAdOx1" = "solid" 
)

km_estimates_rounded %>%
  add_descr() %>%
  mutate(
    treated_descr = factor(
      treated,
      levels = c(0,1),
      labels = c("Unvaccinated", "Vaccinated")
    ),
    colour_var = factor(
      if_else(treated_descr %in% "Unvaccinated", as.character(treated_descr), as.character(brand_descr)),
      levels = names(colour_palette)
    )
  ) %>%
  mutate(
    across(
      outcome_descr, 
      factor, 
      levels = outcome_descr_long,
      labels = outcome_descr_wrap
    )
  ) %>%
  group_by(brand_descr, outcome_descr, treated_descr) %>%
  group_modify(
    ~add_row(
      .x,
      time=0,
      lagtime=0,
      leadtime=1,
      #interval=1,
      surv=1,
      surv.ll=1,
      surv.ul=1,
      risk=0,
      risk.ll=0,
      risk.ul=0,
      .before=0
    )
  ) %>%
  ungroup() %>%
  ggplot(aes(
    group = colour_var, 
    colour = colour_var, 
    fill = colour_var, 
    linetype = colour_var
  )) +
  geom_step(
    aes(x=time, y=risk), 
    direction="vh"
  ) +
  geom_rect(
    aes(xmin=lagtime, xmax=time, ymin=risk.ll, ymax=risk.ul), 
    alpha=0.1, colour="transparent"
  ) +
  facet_grid(
    rows = vars(outcome_descr),
    cols = vars(brand_descr),
    switch = "y",
    scales = "free_y"
  ) +
  scale_color_manual(name = NULL, values = colour_palette) +
  scale_fill_manual(name = NULL, values = colour_palette) +
  scale_linetype_manual(name = NULL, values = linetype_palette) +
  scale_x_continuous(breaks = c(postbaselinecuts)) +
  scale_y_continuous(expand = expansion(mult=c(0,0.01))) +
  # coord_cartesian(xlim=c(0, NA)) +
  labs(
    x = "Days since first dose",
    y = NULL,#"Cumulative incidence",
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size=10, margin = margin(t = 10)),
    # axis.title.y = element_text(size=10, margin = margin(r = 10)),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    strip.background = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90),
    # strip.text = element_text(size=8),
    axis.line.x = element_line(colour = "black"),
    legend.box = "vertical",
    legend.position="bottom"
  )

# ggsave()

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# event counts (single and sequential; text only)

# TODO

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# counts of vaccination following a positive test (single only; text only)

# TODO

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Of xxx person-years of follow-up xxx (xxx%) were after vaccination.

