# # # # # # # # # # # # # # # # # # # # #
# This script:
# plots the cumultive incidence of doses1 and 2 over time
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
# library('lubridate')
# library('survival')

## Import custom user functions from lib
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

## create output directories ----
fs::dir_create(here("output", "report", "brand12counts"))

## Import data_days and carry out some further processing ----
data_days <- read_rds(here("output", "single", "stset", "data_days.rds")) %>%
  transmute(
    
    patient_id,
    
    # this was tstop but I think should be tstart?
    date = as.Date(study_dates$global$index_date) + tstart,
    
    lastfup_status,
    dereg_status,
    death_status,

    vaxbrand12_status = fct_case_when(
      vaxpfizer_status==0 & vaxaz_status==0  ~ "Not vaccinated",
      vaxpfizer_status==1  ~ "BNT162b2\ndose 1",
      vaxpfizer_status==2  ~ "BNT162b2\ndose 2",
      vaxaz_status==1  ~ "ChAdOx1\ndose 1",
      vaxaz_status==2  ~ "ChAdOx1\ndose 2",
      TRUE ~ NA_character_
    ),
    
  )

plot_data <- data_days %>%
  filter(dereg_status==0) %>%
  group_by(date, vaxbrand12_status, lastfup_status, .drop=FALSE) %>%
  summarise(
    n = n(),
    .groups = "keep"
  ) %>%
  group_by(date) %>%
  mutate(
    n_per_10000 = (n/sum(n))*10000
  ) %>%
  ungroup() %>%
  arrange(date, vaxbrand12_status, lastfup_status) %>%
  mutate(
    lastfup_status = if_else(lastfup_status %in% 1, "Died / deregistered", "At-risk"),
    group = factor(
      paste0(lastfup_status,":",vaxbrand12_status),
      levels= map_chr(
        cross2(c("At-risk", "Died / deregistered"), levels(vaxbrand12_status)),
        paste, sep = ":", collapse = ":"
      )
    )
  )

#################################################################################

colour_palette <- c(
  "Not vaccinated" = "#f7f7f7", # light grey
  "BNT162b2\ndose 1" = "#e78ac3", # medium pink / medium grey
  "BNT162b2\ndose 2" = "#e7298a", # dark pink / dark grey
  "ChAdOx1\ndose 1" = "#8da0cb", # medium purple / medium grey
  "ChAdOx1\ndose 2" = "#7570b3" # dark purple / dark grey
)


annotation_data <- tribble(
  ~vaxbrand12_status, ~x_pos, ~y_pos,
  "Not vaccinated", "2020-12-28", 7500,
  "BNT162b2\ndose 1", "2021-03-15", 5000,
  "BNT162b2\ndose 2", "2021-04-05", 2500,
  "ChAdOx1\ndose 1", "2021-03-15", 1250,
  "ChAdOx1\ndose 2", "2021-04-05", 500,
) %>%
  mutate(across(x_pos, as.Date)) 

plot <- plot_data %>%
  group_by(date, vaxbrand12_status) %>%
  summarise(across(n_per_10000, sum)) %>%
  ungroup() %>%
  ggplot() +
  geom_area(
    aes(
      x=date, y=n_per_10000,
      group=vaxbrand12_status,
      fill=vaxbrand12_status,
      # alpha=lastfup_status,
      ),
  colour="black"
  # outline.type = "full"
  )+
  geom_text(
    data=annotation_data,
    aes(x=x_pos, y=y_pos,label=vaxbrand12_status)
    )+
  scale_x_date(
    date_breaks = "1 week", 
    labels = scales::date_format("%Y-%m-%d"),
    expand = c(0,0)
    )+
  scale_y_continuous(expand=c(0,0)) +
  scale_fill_manual(values=colour_palette)+
  # scale_alpha_manual(values=c(0.9,0.1), breaks=c(0.1))+
  labs(
    x="Date",
    y="Status per 10,000 people",
    colour=NULL,
    fill=NULL,
    alpha=NULL
  ) +
  theme_bw(base_size = 12)+
  theme(
    legend.position = "left",
    
    panel.border=element_rect(colour='black', fill=NA),
    
    strip.background = element_blank(),
    # strip.text.y.right = element_text(angle = 0),
    
    axis.line.x = element_line(colour = "black"),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    axis.ticks.x = element_line(colour = 'black')
  ) +
  theme(legend.position = "none") 



ggsave(
  plot = plot,
  filename = "brandcounts12.png",
  path=here("output", "report"),
  units="cm",
  width = 25,
  height = 25
)