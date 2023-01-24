process_vax <- function(.data, stage) {
  
  data_vax <- local({
    
    data_vax_pfizer <- .data %>%
      select(patient_id, matches("covid\\_vax\\_pfizer\\_\\d+\\_date")) %>%
      pivot_longer(
        cols = -patient_id,
        names_to = c(NA, "vax_pfizer_index"),
        names_pattern = "^(.*)_(\\d+)_date",
        values_to = "date",
        values_drop_na = TRUE
      ) %>%
      arrange(patient_id, date)
    
    data_vax_az <- .data %>%
      select(patient_id, matches("covid\\_vax\\_az\\_\\d+\\_date")) %>%
      pivot_longer(
        cols = -patient_id,
        names_to = c(NA, "vax_az_index"),
        names_pattern = "^(.*)_(\\d+)_date",
        values_to = "date",
        values_drop_na = TRUE
      ) %>%
      arrange(patient_id, date)
    
    data_vax <-
      data_vax_pfizer %>%
      full_join(data_vax_az, by=c("patient_id", "date")) %>%
      mutate(
        type = fct_case_when(
          (!is.na(vax_az_index)) & is.na(vax_pfizer_index) ~ "az",
          is.na(vax_az_index) & (!is.na(vax_pfizer_index)) ~ "pfizer",
          (!is.na(vax_az_index)) + (!is.na(vax_pfizer_index)) > 1 ~ "duplicate",
          TRUE ~ NA_character_
        )
      ) %>%
      arrange(patient_id, date) %>%
      group_by(patient_id) %>%
      mutate(
        vax_index=row_number()
      ) %>%
      ungroup()
    
    data_vax
    
  })
  
  data_vax_wide = data_vax %>%
    pivot_wider(
      id_cols= patient_id,
      names_from = c("vax_index"),
      values_from = c("date", "type"),
      names_glue = "covid_vax_{vax_index}_{.value}"
    )
  
  # only add variables corresponding to 2nd dose if stage = single or treated
  if (stage %in% c("single", "treated")) {
    vax2_vars <- rlang::quos(
      
      vax2_type = covid_vax_2_type,
      
      vax2_type_descr = fct_case_when(
        vax2_type == "pfizer" ~ "BNT162b2",
        vax2_type == "az" ~ "ChAdOx1",
        TRUE ~ NA_character_
      ),
      
      vax2_date = covid_vax_2_date,
      
      # day 0 is the day before "start_date"
      vax2_day = as.integer(floor((vax2_date - study_dates$global$index_date)) + 1)
      
    )
  } else if (stage == "potential") {
    vax2_vars <- rlang::quos(
      NULL
    )
  }
  
  .data %>%
    left_join(data_vax_wide, by ="patient_id") %>%
    mutate(
      vax1_type = covid_vax_1_type,
      
      vax1_type_descr = fct_case_when(
        vax1_type == "pfizer" ~ "BNT162b2",
        vax1_type == "az" ~ "ChAdOx1",
        TRUE ~ NA_character_
      ),
      vax1_date = covid_vax_1_date,
      
      # day 0 is the day before "start_date"
      vax1_day = as.integer(floor((vax1_date - study_dates$global$index_date)) + 1),
      
      !!! vax2_vars
      
    ) %>%
    select(
      -starts_with("covid_vax_"),
    ) 
  
}