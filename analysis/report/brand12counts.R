# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# This script:
# plots the cumulative incidence of doses1 and 2 over time
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# Preliminaries ----

# import libraries
library('tidyverse')
library('here')
library('glue')

# import custom user functions and parameters
source(here("analysis", "design.R"))
source(here("analysis", "functions", "utility.R"))

# create output directory
outdir <- here("output", "report", "brand12counts")
fs::dir_create(outdir)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# import data_days ----
cat("Start import:\n")
data_days <- lapply(
  1:process_data_days_n,
  function(iteration) read_rds(here("output", "single", "stset", glue("data_days_{iteration}.rds"))) 
) %>%
  bind_rows() 

# process data_days ----
cat("Start processing:\n")
data_days <- data_days %>%
  transmute(
    
    patient_id,
    
    date = as.Date(study_dates$global$index_date) + tstart,
    
    dereg_status,
    death_status,

    vaxbrand12_status = fct_case_when(
      vaxpfizer_status==0 & vaxaz_status==0  ~ "Not vaccinated",
      vaxpfizer_status==1  ~ "BNT162b2\ndose 1",
      vaxpfizer_status==2  ~ "BNT162b2\ndose 2",
      vaxaz_status==1  ~ "ChAdOx1\ndose 1",
      vaxaz_status==2  ~ "ChAdOx1\ndose 2",
      TRUE ~ NA_character_
    )
    
  )

cat("Generate plot_data:\n")
plot_data <- data_days %>%
  group_by(date, vaxbrand12_status) %>%
  summarise(
    n = roundmid_any(n()),
    .groups = "keep"
  ) %>%
  group_by(date) %>%
  mutate(
    n_per_10000 = (n/sum(n))*10000
  ) %>%
  ungroup() %>%
  arrange(date, vaxbrand12_status) %>%
  select(-n)

# save data for plotting
write_csv(
  plot_data,
  file.path(outdir, "brand12counts_rounded.csv")
)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# create plot ----
cat("Create plot:\n")
colour_palette <- c(
  "Not vaccinated" = "#f7f7f7", # light grey
  "BNT162b2\ndose 1" = "#e78ac3", # medium pink / medium grey
  "BNT162b2\ndose 2" = "#e7298a", # dark pink / dark grey
  "ChAdOx1\ndose 1" = "#8da0cb", # medium purple / medium grey
  "ChAdOx1\ndose 2" = "#7570b3" # dark purple / dark grey
)

# update these values when final results released
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
  ggplot() +
  geom_area(
    aes(
      x=date, y=n_per_10000,
      group=vaxbrand12_status,
      fill=vaxbrand12_status,
      ),
  colour="black",
  position = "stack"
  ) +
  geom_text(
    data=annotation_data,
    aes(x=x_pos, y=y_pos,label=vaxbrand12_status)
    ) +
  scale_x_date(
    date_breaks = "1 week", 
    labels = scales::date_format("%Y-%m-%d"),
    expand = c(0,0)
    ) +
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
    
    panel.border=element_rect(colour='black', fill=NA),
    
    strip.background = element_blank(),
    
    axis.line.x = element_line(colour = "black"),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust=1),
    
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    
    axis.ticks.x = element_line(colour = 'black'),
    
    legend.position = "none"
    
  ) 

ggsave(
  plot = plot,
  filename = "brandcounts12_rounded.png",
  path=outdir,
  units="cm",
  width = 25,
  height = 25
)