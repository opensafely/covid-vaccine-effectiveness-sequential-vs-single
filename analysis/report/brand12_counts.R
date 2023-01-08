# # # # # # # # # # # # # # # # # # # # #
# This script:
# plots the cumultive incidence of doses1 and 2 over time
# # # # # # # # # # # # # # # # # # # # #

# Preliminaries ----

## Import libraries ----
library('tidyverse')
library('here')
library('glue')
library('lubridate')
library('survival')

## Import custom user functions from lib
source("analysis", "design.R")
# source(here("lib", "utility_functions.R"))
# source(here("lib", "redaction_functions.R"))


## custom functions ----

# function to extract total plot height minus panel height
plotHeight <- function(plot, unit){
  grob <- ggplot2::ggplotGrob(plot)
  grid::convertHeight(gtable::gtable_height(grob), unitTo=unit, valueOnly=TRUE)
}

# function to extract total plot width minus panel height
plotWidth <- function(plot, unit){
  grob <- ggplot2::ggplotGrob(plot)
  grid::convertWidth(gtable::gtable_width(grob), unitTo=unit, valueOnly=TRUE)
}

# function to extract total number of bars plot (strictly this is the number of rows in the build of the plot data)
plotNbars <- function(plot){
  length(unique(ggplot2::ggplot_build(plot)$data[[1]]$x))
}

# function to extract total number of bars plot (strictly this is the number of rows in the build of the plot data)
plotNfacetrows <- function(plot){
  length(levels(ggplot2::ggplot_build(plot)$data[[1]]$PANEL))
}

plotNyscales <- function(plot){
  length(ggplot2::ggplot_build(plot)$layout$panel_scales_y[[1]]$range$range)
}

# function to extract total number of panels
plotNpanelrows <- function(plot){
  length(unique(ggplot2::ggplot_build(plot)$layout$layout$ROW))
}

## create output directories ----
# fs::dir_create(here("output", cohort, "descriptive", "plots"))
# fs::dir_create(here("output", cohort, "descriptive", "tables"))

## define theme ----

plot_theme <-
  theme_minimal()+
  theme(
    legend.position = "left",
    panel.border=element_rect(colour='black', fill=NA),
    strip.text.y.right = element_text(angle = 0),
    axis.line.x = element_line(colour = "black"),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_line(colour = 'black')
  )

## Import processed data ----

data_fixed <- read_rds(here("output", cohort, "data", "data_fixed.rds"))
data_tte <- read_rds(here("output", cohort, "data", "data_tte.rds"))
data_pt <- read_rds(here("output", cohort, "data", "data_pt.rds"))

# create plots ----
data_by_day <-
  data_pt %>%
  mutate(
    
    day = tstop,
    date = as.Date(start_date) + day,
    week = lubridate::floor_date(date, unit="week", week_start=2), #week commencing tuesday (since index date is a tuesday)
    #date = week,
    
    vaxany_status_onedose = vaxany_status!=0,
    vaxany_status = fct_case_when(
      vaxany_status==0 & death_status==0 & dereg_status==0 ~ "Not vaccinated",
      vaxany_status==1 ~ "One dose",
      vaxany_status==2 ~ "Two doses",
      death_status==1 | dereg_status==1 ~ "Died/deregistered",
      TRUE ~ NA_character_
    ),
    
    vaxbrand1_status = fct_case_when(
      vaxpfizer_status==0 & vaxaz_status==0  ~ "Not vaccinated",
      vaxpfizer_status>0  ~ "BNT162b2",
      vaxaz_status>0  ~ "ChAdOx1",
      TRUE ~ NA_character_
    ),
    
    vaxbrand12_status = fct_case_when(
      vaxpfizer_status==0 & vaxaz_status==0  ~ "Not vaccinated",
      vaxpfizer_status==1  ~ "BNT162b2\ndose 1",
      vaxpfizer_status==2  ~ "BNT162b2\ndose 2",
      vaxaz_status==1  ~ "ChAdOx1\ndose 1",
      vaxaz_status==2  ~ "ChAdOx1\ndose 2",
      TRUE ~ NA_character_
    ),
    
  )




plot_brand1_counts <- function(var, var_descr){
  
  data1 <- data_by_day %>%
    mutate(
      variable = data_by_day[[var]]
    ) %>%
    filter(dereg_status==0) %>%
    group_by(date, variable, vaxbrand1_status, lastfup_status, .drop=FALSE) %>%
    summarise(
      n = n(),
    ) %>%
    group_by(date, variable) %>%
    mutate(
      n_per_10000 = (n/sum(n))*10000
    ) %>%
    ungroup() %>%
    arrange(date, variable, vaxbrand1_status, lastfup_status) %>%
    mutate(
      lastfup_status = if_else(lastfup_status %in% 1, "Died / deregistered", "At-risk"),
      group= factor(
        paste0(lastfup_status,":",vaxbrand1_status),
        levels= map_chr(
          cross2(c("At-risk", "Died / deregistered"), levels(vaxbrand1_status)),
          paste, sep = ":", collapse = ":"
        )
      )
    )
  
  
  #colorspace::lighten("#1b9e77", 0.25)
  
  plot <- data1 %>%
    ggplot() +
    geom_area(aes(x=date, y=n_per_10000,
                  group=group,
                  fill=vaxbrand1_status,
                  alpha=lastfup_status
    )
    )+
    facet_grid(rows=vars(variable))+
    scale_x_date(date_breaks = "1 week", labels = scales::date_format("%Y-%m-%d"))+
    scale_fill_manual(values=c("#d95f02", "#7570b3", "#1b9e77"))+
    scale_alpha_manual(values=c(0.8,0.3), breaks=c(0.3))+
    labs(
      x="Date",
      y="Status per 10,000 people",
      colour=NULL,
      fill=NULL,
      alpha=NULL
    ) +
    plot_theme+
    theme(legend.position = "bottom")
  
  plot
}



plot_brand12_counts <- function(var, var_descr){
  
  data1 <- data_by_day %>%
    mutate(
      variable = data_by_day[[var]]
    ) %>%
    filter(dereg_status==0) %>%
    group_by(date, variable, vaxbrand12_status, lastfup_status, .drop=FALSE) %>%
    summarise(
      n = n(),
    ) %>%
    group_by(date, variable) %>%
    mutate(
      n_per_10000 = (n/sum(n))*10000
    ) %>%
    ungroup() %>%
    arrange(date, variable, vaxbrand12_status, lastfup_status) %>%
    mutate(
      lastfup_status = if_else(lastfup_status %in% 1, "Died / deregistered", "At-risk"),
      group= factor(
        paste0(lastfup_status,":",vaxbrand12_status),
        levels= map_chr(
          cross2(c("At-risk", "Died / deregistered"), levels(vaxbrand12_status)),
          paste, sep = ":", collapse = ":"
        )
      )
    )
  
  #colorspace::lighten("#1b9e77", 0.25)
  
  plot <- data1 %>%
    ggplot() +
    geom_area(aes(x=date, y=n_per_10000,
                  group=group,
                  fill=vaxbrand12_status,
                  alpha=lastfup_status
    )
    )+
    facet_grid(rows=vars(variable))+
    scale_x_date(date_breaks = "1 week", labels = scales::date_format("%Y-%m-%d"))+
    scale_fill_manual(values=c("#d95f02", "#7570b3", "#9590D3", "#1b9e77",  "#4BBA93"))+
    scale_alpha_manual(values=c(0.9,0.1), breaks=c(0.1))+
    labs(
      x="Date",
      y="Status per 10,000 people",
      colour=NULL,
      fill=NULL,
      alpha=NULL
    ) +
    plot_theme+
    theme(legend.position = "bottom")
  
  plot
}


#################################################################################




plot_data <- bind_rows(
  plot_data_over80s,
  plot_data_in70s
)


plot_theme <-
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
  )

colour_palette <- c(
  "Not vaccinated" = "#f7f7f7", # light grey
  "BNT162b2\ndose 1" = "#e78ac3", # medium pink / medium grey
  "BNT162b2\ndose 2" = "#e7298a", # dark pink / dark grey
  "ChAdOx1\ndose 1" = "#8da0cb", # medium purple / medium grey
  "ChAdOx1\ndose 2" = "#7570b3" # dark purple / dark grey
)


annotation_data <- tribble(
  ~cohort_descr, ~vaxbrand12_status, ~x_pos, ~y_pos,
  "80+ years", "Not vaccinated", "2020-12-28", 7500,
  "80+ years", "BNT162b2\ndose 1", "2021-03-15", 5000,
  "80+ years", "BNT162b2\ndose 2", "2021-04-05", 2500,
  "80+ years", "ChAdOx1\ndose 1", "2021-03-15", 1250,
  "80+ years", "ChAdOx1\ndose 2", "2021-04-05", 500,
  "70-79 years", "Not vaccinated", "2021-01-25", 7500,
  "70-79 years", "BNT162b2\ndose 1","2021-03-15", 5000,
  "70-79 years", "BNT162b2\ndose 2","2021-04-05",2500,
  "70-79 years", "ChAdOx1\ndose 1","2021-03-15", 1250,
  "70-79 years", "ChAdOx1\ndose 2","2021-04-05",500
) %>%
  mutate(across(x_pos, as.Date)) 

plot <- plot_data %>%
  group_by(date, vaxbrand12_status, cohort_descr) %>%
  summarise(across(n_per_10000, sum)) %>%
  ungroup() %>%
  mutate(across(cohort_descr, ~str_c(.x, " years"))) %>%
  ggplot() +
  geom_area(aes(x=date, y=n_per_10000,
                group=vaxbrand12_status,
                fill=vaxbrand12_status,
                # alpha=lastfup_status,
  ),
  colour="black"
  # outline.type = "full"
  )+
  geom_text(data=annotation_data, aes(x=x_pos,y=y_pos,label=vaxbrand12_status))+
  facet_wrap(vars(fct_rev(cohort_descr)), ncol = 1)+
  # facet_grid(rows=vars(fct_rev(cohort_descr)))+
  scale_x_date(date_breaks = "1 week", labels = scales::date_format("%Y-%m-%d"), expand = c(0,0))+
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
  plot_theme+
  theme(legend.position = "none")

fs::dir_create(here("output", "combined"))

ggsave(
  plot = plot,
  filename = paste0("brandcounts12.svg"),
  path=here("output", "combined"),
  units="cm",
  width = 25,
  height = 25
)


ggsave(
  plot = plot,
  filename = paste0("brandcounts12.png"),
  path=here("output", "combined"),
  units="cm",
  width = 25,
  height = 25
)