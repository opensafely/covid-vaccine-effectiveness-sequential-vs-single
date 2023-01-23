# # # # # # # # # # # # # # # # # # # # ## # # # # # # # # # # # # # # # # # # # 
# function to calculate weights for treatment model ----
## if exposure is pfizer vaccine, then create model for vaccination by pfizer + model for az for censoring weights
## if exposure is az vaccine, then create model for vaccination by az + model for pfizer for censoring weights

get_ipw_weights <- function(
  # data,
  event,
  event_status,
  event_atrisk,
  
  # sample_type,
  sample_amount,
  
  ipw_formula,
  ipw_formula_fxd,
  
  subgroup_level
){
  
  # stopifnot(sample_type %in% c("random_prop", "random_n", "nonoutcomes_n"))
  
  ## read and process data_days (one row per person day)
  # see analysis/single/process/process_data_days_function.R for the function process_data_days_function
  # (do this _within_ loop so that it can be deleted just before models are run, to reduce RAM use)
  cat("Start `process_data_days_function` for stage=vaccine\n")
  data_atrisk <- bind_rows(
    lapply(
      1:process_data_days_n,
      function(i) 
        process_data_days_function(
          file = "model",
          stage = "vaccine",
          iteration = i
        )
    )
  ) 
  cat("End `process_data_days_function`\n")
  
  data_atrisk <- data_atrisk %>%
    mutate(
      event = data_atrisk[[event]],
      event_status = data_atrisk[[event_status]],
      event_atrisk = data_atrisk[[event_atrisk]],
    ) %>%
    filter(event_atrisk)
  
  # if(removeobs) rm(data_samples, data_fixed)
  
  name <- str_remove(event_atrisk, "_atrisk")
  
  # if(sample_type=="random_n"){
    data_sample <- data_atrisk %>%
      distinct(patient_id) %>%
      transmute(
        patient_id,
        sample_event = sample_random_n(patient_id, sample_amount),
        sample_weights_event = sample_event*1L,
      ) %>%
      filter(sample_event)
  # }
  
  data_atrisk_sample <- data_atrisk %>%
    right_join(data_sample, by="patient_id") 
  
  rm(data_sample)
  
  # with time-updating covariates
  cat("  \n")
  cat(glue("{event}  \n"))
  event_model <- parglm(
    formula = ipw_formula,
    data = data_atrisk_sample,
    family = binomial,
    weights = sample_weights_event,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  
  cat(glue("{event} data size = ", length(event_model$y)), "\n")
  cat(glue("memory usage = ", format(object.size(event_model), units="GB", standard="SI", digits=3L)), "\n")
  cat("warnings: ", "\n")
  print(warnings())
  
  # without time-updating covariates
  cat("  \n")
  cat(glue("{event}_fxd  \n"))
  event_model_fxd <- parglm(
    formula = ipw_formula_fxd,
    data = data_atrisk_sample,
    family = binomial,
    weights = sample_weights_event,
    control = parglmparams,
    na.action = "na.fail",
    model = FALSE
  )
  
  cat(glue("{event}_fxd data size = ", length(event_model_fxd$y)), "\n")
  cat(glue("memory usage = ", format(object.size(event_model_fxd), units="GB", standard="SI", digits=3L)), "\n")
  cat("warnings: ", "\n")
  print(warnings())
  
  # save outputs
  write_rds(
    event_model, 
    file.path(outdir, glue("model_{name}_{subgroup_level}.rds")), 
    compress="gz"
  )
  write_rds(
    ipw_formula, 
    file.path(outdir, glue("model_formula_{name}_{subgroup_level}.rds")), 
    compress="gz"
  )
  
  rm("data_atrisk_sample")
  
  # split the dataset into `process_data_days_n` (defined in analysis/design.R)
  id_split <- data_atrisk %>%
    distinct(patient_id) %>%
    mutate(id_group = sample(1:process_data_days_n, size = n(), replace=TRUE))
  
  data_atrisk <- data_atrisk %>%
    left_join(id_split, by = "patient_id") %>%
    group_split(id_group) %>%
    as.list()
  
  # loop so that only calculating model predictions in 1/10th of dataset at a time
  for (i in 1:process_data_days_n) {
    
    # get predictions from model
    data_atrisk[[i]] <- data_atrisk[[i]] %>%
      transmute(
        patient_id,
        tstart, tstop,
        event,
        event_status,
        # get predicted probabilities from ipw models
        pred_event=predict(event_model, type="response", newdata=data_atrisk[[i]]),
        pred_event_fxd=predict(event_model_fxd, type="response", newdata=data_atrisk[[i]])
      ) %>%
      arrange(patient_id, tstop) %>%
      group_by(patient_id) %>%
      mutate(
        
        # get probability of occurrence of realised event status (time varying model)
        probevent_realised = case_when(
          event!=1L ~ 1 - pred_event,
          event==1L ~ pred_event,
          TRUE ~ NA_real_
        ),
        
        # get probability of occurrence of realised event status (non-time varying model)
        probevent_realised_fxd = case_when(
          event!=1L ~ 1 - pred_event_fxd,
          event==1L ~ pred_event_fxd,
          TRUE ~ NA_real_
        ),
        
        # stabilised inverse probability weights
        ipweight_stbl = probevent_realised_fxd/probevent_realised
        
      ) %>%
      ungroup()
    
  }
  
  data_atrisk <- bind_rows(data_atrisk)
  
  stopifnot("probs should all be non-null" = all(!is.na(data_atrisk$probevent_realised)))
  stopifnot("probs (fxd) should all be non-null" = all(!is.na(data_atrisk$probevent_realised_fxd)))
  
  weights <- data_atrisk %>%
    select(
      patient_id,
      tstart, tstop,
      ipweight_stbl
    )
  
  weights[[glue("ipweight_stbl_{name}")]] <- weights$ipweight_stbl
  weights$ipweight_stbl <- NULL
  
  return(weights)
  
}