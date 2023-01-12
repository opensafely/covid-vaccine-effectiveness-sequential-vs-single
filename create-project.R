library("tidyverse")
library("yaml")
library("here")
library("glue")
#library("rlang")
 
source(here("analysis", "design.R"))

# create action functions ----

## create comment function ----
comment <- function(...) {
  list_comments <- list(...)
  comments <- map(list_comments, ~paste0("## ", ., " ##"))
  comments
}


## create function to convert comment "actions" in a yaml string into proper comments
convert_comment_actions <-function(yaml.txt) {
  yaml.txt %>%
    str_replace_all("\\\n(\\s*)\\'\\'\\:(\\s*)\\'", "\n\\1")  %>%
    #str_replace_all("\\\n(\\s*)\\'", "\n\\1") %>%
    str_replace_all("([^\\'])\\\n(\\s*)\\#\\#", "\\1\n\n\\2\\#\\#") %>%
    str_replace_all("\\#\\#\\'\\\n", "\n")
}


## generic action function ----
action <- function(
  name,
  run,
  arguments=NULL,
  needs=NULL,
  highly_sensitive=NULL,
  moderately_sensitive=NULL,
  ... # other arguments / options for special action types
) {

  outputs <- list(
    highly_sensitive = highly_sensitive,
    moderately_sensitive = moderately_sensitive
  )
  outputs[sapply(outputs, is.null)] <- NULL

  action <- list(
    run = paste(c(run, arguments), collapse=" "),
    needs = needs,
    outputs = outputs,
    ... = ...
  )
  action[sapply(action, is.null)] <- NULL

  action_list <- list(name = action)
  names(action_list) <- name

  action_list
}

namelesslst <- function(...){
  unname(lst(...))
}

## actions for a single matching round ----




action_1matchround <- function(cohort, matching_round){
  
  control_extract_date <- study_dates[[cohort]][[glue("control_extract_dates")]][matching_round]
  
  if (matching_round == 1) {
    process_controlpotential_hsoutputs <- lst(
      rds = glue("output/sequential/{cohort}/matchround{matching_round}/process/*.rds"),
      csv = glue("output/sequential/{cohort}/matchround{matching_round}/process/*.csv.gz")
    )
  } else {
    process_controlpotential_hsoutputs <- lst(
      rds = glue("output/sequential/{cohort}/matchround{matching_round}/process/*.rds")
    )
  }
  
  splice(
    action(
      name = glue("extract_controlpotential_{cohort}_{matching_round}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlpotential", 
        " --output-file output/sequential/{cohort}/matchround{matching_round}/extract/input_controlpotential.feather", 
        " --param cohort={cohort}",
        " --param matching_round={matching_round}",
        " --param index_date={control_extract_date}"
      ),
      needs = c(
        "design",
        if(matching_round>1) {
          c(glue("process_controlpotential_{cohort}_{matching_round-1}"),glue("process_controlactual_{cohort}_{matching_round-1}"))
        } else {
            NULL
          }
      ) %>% as.list,
      highly_sensitive = lst(
        cohort = glue("output/sequential/{cohort}/matchround{matching_round}/extract/input_controlpotential.feather")
      )
    ),
    
    action(
      name = glue("process_controlpotential_{cohort}_{matching_round}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("potential", cohort, matching_round),
      needs = namelesslst(
        glue("extract_controlpotential_{cohort}_{matching_round}"),
      ),
      highly_sensitive = process_controlpotential_hsoutputs,
      moderately_sensitive = lst(
        input_controlpotential_skim = glue("output/sequential/{cohort}/matchround{matching_round}/extract/potential/*.txt"),
        data_processed_skim = glue("output/sequential/{cohort}/matchround{matching_round}/potential/*.txt"),
        data_controlpotential_skim = glue("output/sequential/{cohort}/matchround{matching_round}/process/*.txt")
      )
    ),
    
    action(
      name = glue("match_potential_{cohort}_{matching_round}"),
      run = glue("r:latest analysis/sequential/matching/match_potential.R"),
      arguments = c(cohort, matching_round),
      needs = c(
        glue("process_treated"), 
        glue("process_controlpotential_{cohort}_{matching_round}"),
        if(matching_round>1) {glue("process_controlactual_{cohort}_{matching_round-1}")} else {NULL}
      ) %>% as.list,
      highly_sensitive = lst(
        rds = glue("output/sequential/{cohort}/matchround{matching_round}/potential/*.rds"),
        csv = glue("output/sequential/{cohort}/matchround{matching_round}/potential/*.csv.gz"),
      )
    ),
    
    action(
      name = glue("extract_controlactual_{cohort}_{matching_round}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlactual", 
        " --output-file output/sequential/{cohort}/matchround{matching_round}/extract/input_controlactual.feather", 
        " --param cohort={cohort}",
        " --param matching_round={matching_round}",
      ),
      needs = namelesslst(
        "design",
        glue("match_potential_{cohort}_{matching_round}"), 
      ),
      highly_sensitive = lst(
        cohort = glue("output/sequential/{cohort}/matchround{matching_round}/extract/input_controlactual.feather")
      )
    ),
    
    
    action(
      name = glue("process_controlactual_{cohort}_{matching_round}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("actual", cohort, matching_round),
      needs = c(
        glue("process_treated"),
        glue("match_potential_{cohort}_{matching_round}"), 
        glue("extract_controlpotential_{cohort}_{matching_round}"),  # this is only necessary for the dummy data
        glue("process_controlpotential_{cohort}_{matching_round}"), # this is necessary for the vaccine data
        glue("extract_controlactual_{cohort}_{matching_round}"),
        if(matching_round>1){glue("process_controlactual_{cohort}_{matching_round-1}")} else {NULL}
      ) %>% as.list,
      highly_sensitive = lst(
        rds = glue("output/sequential/{cohort}/matchround{matching_round}/actual/*.rds"),
        csv = glue("output/sequential/{cohort}/matchround{matching_round}/actual/*.csv.gz"),
      ),
      moderately_sensitive = lst(
        input_controlactual_skim = glue("output/sequential/{cohort}/matchround{matching_round}/extract/actual/*.txt"),
        data_actual_skim = glue("output/sequential/{cohort}/matchround{matching_round}/actual/*.txt"),
      )
    )

  )
}

# test function
#action_1matchround("pfizer", 2)

# create all necessary actions for n matching rounds
action_extract_and_match <- function(cohort, n_matching_rounds){
  
  allrounds <- map(seq_len(n_matching_rounds), ~action_1matchround(cohort, .x)) %>% flatten
  
  splice(
    
    allrounds,
    
    
    action(
      name = glue("extract_controlfinal_{cohort}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlfinal", 
        " --output-file output/sequential/{cohort}/extract/input_controlfinal.feather",
        " --param cohort={cohort}",
        " --param n_matching_rounds={n_matching_rounds}",
      ),
      needs = namelesslst(
        "design",
        glue("process_controlactual_{cohort}_{n_matching_rounds}")
      ),
      highly_sensitive = lst(
        extract = glue("output/sequential/{cohort}/extract/input_controlfinal.feather")
      )
    ),
    
    action(
      name = glue("dummydata_controlfinal_{cohort}"),
      run = glue("r:latest analysis/dummy/dummydata_controlfinal.R"),
      arguments = c(cohort),
      needs =map(
        seq_len(n_matching_rounds),
        ~glue("process_controlactual_{cohort}_",.x)
      ),
      highly_sensitive = lst(
        dummydata_controlfinal = glue("output/sequential/{cohort}/dummydata/dummy_control_final.feather")
      ),
    ),
    
    action(
      name = glue("process_controlfinal_{cohort}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("final", cohort),
      needs = c(
        map(
          seq_len(n_matching_rounds),
          ~glue("process_controlactual_{cohort}_",.x)
        ),
        glue("extract_controlfinal_{cohort}"),
        glue("process_treated"),
        glue("dummydata_controlfinal_{cohort}")
      ),
      highly_sensitive = lst(
        extract = glue("output/sequential/{cohort}/match/*.rds")
      ),
      moderately_sensitive = lst(
        input_controlfinal_skim = glue("output/sequential/{cohort}/extract/*.txt"),
        data_matched_skim = glue("output/sequential/{cohort}/match/*.txt")
      )
    )
  )
  
}

# test action
# action_extract_and_match("pfizer", 2)


action_kmcox <- function(cohort, subgroup, outcome){
  action(
    name = glue("kmcox_{cohort}_{subgroup}_{outcome}"),
    run = glue("r:latest analysis/sequential/model/kmcox.R"),
    arguments = c(cohort, subgroup, outcome),
    needs = namelesslst(
      glue("process_controlfinal_{cohort}"),
    ),
    moderately_sensitive= lst(
      #csv= glue("output/sequential/{cohort}/models/km/{subgroup}/{outcome}/*.csv"),
      rds= glue("output/sequential/{cohort}/model/{subgroup}/{outcome}/*.rds"),
      png= glue("output/sequential/{cohort}/model/{subgroup}/{outcome}/*.png"),
    )
  )
}

## model action function ----
# action_kmcox_combine <- function(
#     cohort
# ){
# 
#   action(
#     name = glue("combine_kmcox_{cohort}"),
#     run = glue("r:latest analysis/sequential/model/kmcox_combine.R"),
#     arguments = c(cohort),
#     needs = splice(
#       as.list(
#         glue_data(
#           .x=expand_grid(
#             subgroup=model_subgroups,
#             outcome=model_outcomes,
#           ),
#           "kmcox_{cohort}_{subgroup}_{outcome}"
#         )
#       )
#     ),
#     moderately_sensitive = lst(
#       rds = glue("output/sequential/{cohort}/model/combine/*.csv"),
#       png = glue("output/sequential/{cohort}/model/combine/*.png"),
#     )
#   )
# }

cohort_seqtrial <- function(cohort) {
  
  splice(
    comment("# # # # # # # # # # # # # # # # # # #",
            glue("{cohort} cohort"),
            "# # # # # # # # # # # # # # # # # # #"),
    
    comment("# # # # # # # # # # # # # # # # # # #",
            "Extract, process and match control data",
            "# # # # # # # # # # # # # # # # # # #"),
    
    action_extract_and_match(cohort, n_matching_rounds),
    
    action(
      name = glue("coverage_{cohort}"),
      run = glue("r:latest analysis/sequential/matching/coverage.R"),
      arguments = c(cohort),
      needs = namelesslst(
        "process_treated",
        glue("process_controlfinal_{cohort}"),
      ),
      moderately_sensitive= lst(
        coverage = glue("output/report/coverage/coverage_{cohort}.csv")
      )
    ),
    
    action(
      name = glue("table1_{cohort}"),
      run = glue("r:latest analysis/report/table1.R"),
      arguments = c(cohort),
      needs = namelesslst(
        "process_treated",
        glue("process_controlfinal_{cohort}"),
      ),
      moderately_sensitive= lst(
        coverage = glue("output/report/table1/table1_{cohort}.csv")
      )
    ),
    
    comment("# # # # # # # # # # # # # # # # # # #",
            "Model",
            "# # # # # # # # # # # # # # # # # # #"),
    
    
    expand_grid(
      subgroup=model_subgroups,
      outcome=model_outcomes,
    ) %>%
      pmap(
        function(brand, subgroup, outcome) action_kmcox(cohort, subgroup, outcome)
      ) %>%
      unlist(recursive = FALSE)#,
    
    # action_kmcox_combine(cohort)
    
  )
  
}


model_single <- function(brand, subgroup, outcome, ipw_sample_random_n, msm_sample_nonoutcomes_n) {
  
  splice(
    
    comment("# # # # # # # # # # # # # # # # # # #", 
            glue("Model: {brand}; {subgroup}; {outcome};"), 
            "# # # # # # # # # # # # # # # # # # #"),
    
    action(
      name = glue("msm_preflight_{brand}_{subgroup}_{outcome}_{ipw_sample_random_n}_{msm_sample_nonoutcomes_n}"),
      run = "r:latest analysis/single/model/msm_preflight.R",
      arguments = c(brand, subgroup, outcome, ipw_sample_random_n, msm_sample_nonoutcomes_n),
      needs = namelesslst(
        "process_stset"
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/*.csv"),
        html = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/*.html")
      )
    ),
    
    action(
      name = glue("msm_{brand}_{subgroup}_{outcome}_{ipw_sample_random_n}_{msm_sample_nonoutcomes_n}"),
      run = "r:latest analysis/single/model/msm.R",
      arguments = c(brand, subgroup, outcome, ipw_sample_random_n, msm_sample_nonoutcomes_n),
      needs = namelesslst(
        "process_stset",
        glue("msm_preflight_{brand}_{subgroup}_{outcome}_{ipw_sample_random_n}_{msm_sample_nonoutcomes_n}")
      ),
      highly_sensitive = lst(
        rds = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.rds")
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.csv"),
        png = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.svg"),
        txt = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.txt")
      )
    ),
    
    action(
      name = glue("msm_postprocess_{brand}_{subgroup}_{outcome}"),
      run = "r:latest analysis/single/model/msm_postprocess.R",
      arguments = c(brand, subgroup, outcome),
      needs = namelesslst(
        glue("msm_{brand}_{subgroup}_{outcome}_{ipw_sample_random_n}_{msm_sample_nonoutcomes_n}")
      ),
      highly_sensitive = lst(
        rds = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.rds")
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.csv"),
        png = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.svg")
      )
    )

  )
  
}


# specify project ----

## defaults ----
defaults_list <- lst(
  version = "3.0",
  expectations= lst(population_size=100000L)
)

## actions ----
actions_list <- splice(

  comment("# # # # # # # # # # # # # # # # # # #",
          "DO NOT EDIT project.yaml DIRECTLY",
          "This file is created by create-project.R",
          "Edit and run create-project.R to update the project.yaml",
          "# # # # # # # # # # # # # # # # # # #",
           " "
          ),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "Preliminaries", 
          "# # # # # # # # # # # # # # # # # # #"),
  action(
    name = "design",
    run = glue("r:latest analysis/design.R"),
    moderately_sensitive = lst(
      lib = glue("lib/design/study-dates.json")
    ),
  ),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "SEQUENTIAL TRIAL APPROACH", 
          "# # # # # # # # # # # # # # # # # # #"),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "Extract and process treated data", 
          "# # # # # # # # # # # # # # # # # # #"),
  # all treated people
  action(
    name = "extract_treated",
    run = glue(
      "cohortextractor:latest generate_cohort", 
      " --study-definition study_definition_treated", 
      " --output-file output/sequential/treated/extract/input_treated.feather",
    ),
    needs = namelesslst(
      "design"
    ),
    highly_sensitive = lst(
      extract = "output/sequential/treated/extract/input_treated.feather"
    ),
  ),
  
  # all treated people
  action(
    name = "process_treated",
    run = "r:latest analysis/process/process_data.R",
    arguments = "treated",
    needs = namelesslst(
      "extract_treated"
    ),
    highly_sensitive = lst(
      eligible = "output/sequential/treated/eligible/*.rds",
      pfizer = "output/sequential/pfizer/treated/*.rds",
      az = "output/sequential/az/treated/*.rds"
    ),
    moderately_sensitive = lst(
      eligiblecsv = "output/sequential/treated/eligible/*.csv",
      input_treated_skim = "output/sequential/treated/extract/*.txt",
      data_processed_skim = "output/sequential/treated/process/*.txt",
      data_eligible_skim = "output/sequential/treated/eligible/*.txt"
    )
  ),

  cohort_seqtrial("pfizer"),
  cohort_seqtrial("az"),
  
  action(
    name = glue("combine_kmcox"),
    run = glue("r:latest analysis/sequential/model/kmcox_combine.R"),
    needs = splice(
      as.list(
        glue_data(
          .x=expand_grid(
            cohort=model_brands,
            subgroup=model_subgroups,
            outcome=model_outcomes,
          ),
          "kmcox_{cohort}_{subgroup}_{outcome}"
        )
      )
    ),
    moderately_sensitive = lst(
      rds = "output/sequential/combine/*.csv",
      png = "output/sequential/combine/*.png"
    )
  ),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "SINGLE TRIAL APPROACH", 
          "# # # # # # # # # # # # # # # # # # #"),
  comment("Extract and process data", 
          "# # # # # # # # # # # # # # # # # # #"),
  
  action(
    name = "process_single",
    run = "r:latest analysis/process/process_data.R",
    arguments = "single",
    needs = namelesslst(
      "extract_treated",
      "extract_controlpotential_pfizer_1"
    ),
    highly_sensitive = lst(
      eligible = "output/single/eligible/*.rds"
    ),
    moderately_sensitive = lst(
      eligiblecsv = "output/single/eligible/*.csv",
      eligiblecsvgz = "output/single/eligible/*.csv.gz",
      data_processed_skim = "output/single/process/*.txt",
      data_eligible_skim = "output/single/eligible/*.txt"
    )
  ),
  
  action(
    name = "table1_single",
    run = glue("r:latest analysis/report/table1.R"),
    arguments = c("single"),
    needs = namelesslst(
      "process_single"
    ),
    moderately_sensitive= lst(
      coverage = "output/report/table1/table1_single.csv"
    )
  ),
  
  # extract outcome and timevarying variables for the single trial approach
  action(
    name = "extract_timevarying",
    run = glue(
      "cohortextractor:latest generate_cohort", 
      " --study-definition study_definition_timevarying", 
      " --output-file output/single/extract/input_timevarying.feather",
    ),
    needs = namelesslst(
      "design", "process_single"
    ),
    highly_sensitive = lst(
      extract = "output/single/extract/input_timevarying.feather"
    )
  ),
  
  action(
    name = "dummydata_timevarying",
    run = "r:latest analysis/dummy/dummydata_timevarying.R",
    needs = namelesslst(
      "process_single",
      "extract_timevarying"
    ),
    highly_sensitive = lst(
      dummydata = "output/single/dummydata/*.feather"
    )
  ),
  
  action(
    name = "process_timevarying",
    run = "r:latest analysis/single/process/process_timevarying.R",
    needs = namelesslst(
      "process_single",
      "extract_timevarying",
      "dummydata_timevarying"
    ),
    highly_sensitive = lst(
      processed = "output/single/process/*.rds"
    )
  ),
  
  action(
    name = "process_stset",
    run = "r:latest analysis/single/process/process_stset.R",
    needs = namelesslst(
      "process_single",
      "process_timevarying"
    ),
    highly_sensitive = lst(
      processed = "output/single/stset/*.rds"
    )
  ),
  
  # model actions
  expand_grid(
    brand=model_brands,
    subgroup=model_subgroups,
    outcome=model_outcomes,
  ) %>%
    pmap(
      function(brand, subgroup, outcome) model_single(brand, subgroup, outcome, 3000, 1000)#150000, 50000)
    ) %>%
    unlist(recursive = FALSE),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "Combine model outputs",
          "# # # # # # # # # # # # # # # # # # #"),
  action(
    name = "msm_combine",
    run = glue("r:latest analysis/model/msm_combine.R"),
    needs = splice(
      as.list(
        glue_data(
          .x=expand_grid(
            brand=model_brands,
            subgroup=model_subgroups,
            outcome=model_outcomes,
          ),
          "msm_postprocess_{brand}_{subgroup}_{outcome}"
        )
      )
    ),
    moderately_sensitive = lst(
      csv = "output/single/combine/*.csv",
      svg = "output/single/combine/*.svg"
    )
  ),
  
  comment("# # # # # # # # # # # # # # # # # # #", 
          "REPORT", 
          "# # # # # # # # # # # # # # # # # # #"),
  
  action(
    name = "flowchart",
    run = glue("r:latest analysis/report/flowchart.R"),
    needs = namelesslst(
      "process_treated",
      "process_controlfinal_pfizer",
      "process_controlfinal_az",
      "process_single"
    ),
    moderately_sensitive = lst(
      flow_matching = "output/report/flowchart/*.csv"
    )
  ),
  
  action(
    name = "brand12counts",
    run = glue("r:latest analysis/report/brand12counts.R"),
    needs = namelesslst(
      "process_stset"
    ),
    moderately_sensitive = lst(
      csv = "output/report/brand12counts/*.csv",
      plots = "output/report/brand12counts/*.png"
    )
  ),
  
  # 
  # comment("# # # # # # # # # # # # # # # # # # #", 
  #         "Move files for release", 
  #         "# # # # # # # # # # # # # # # # # # #"),
  # 
  # action(
  #   name = "release",
  #   run = glue("r:latest analysis/release_objects.R"),
  #   needs = namelesslst(
  #     glue("combine_km_pfizer"),
  #     glue("table1_pfizer"),
  #     glue("combine_km_under12"),
  #     glue("table1_under12"),
  #   ),
  #   highly_sensitive = lst(
  #     txt = glue("output/meta-release/*.txt"),
  #     csv = glue("output/release/*.csv"),
  #   ),
  # ),

  comment("#### End ####")
)

project_list <- splice(
  defaults_list,
  list(actions = actions_list)
)

## convert list to yaml, reformat comments and whitespace ----
thisproject <- as.yaml(project_list, indent=2) %>%
  # convert comment actions to comments
  convert_comment_actions() %>%
  # add one blank line before level 1 and level 2 keys
  str_replace_all("\\\n(\\w)", "\n\n\\1") %>%
  str_replace_all("\\\n\\s\\s(\\w)", "\n\n  \\1")


# if running via opensafely, check that the project on disk is the same as the project created here:
if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("expectations", "tpp")){

  thisprojectsplit <- str_split(thisproject, "\n")
  currentproject <- readLines(here("project.yaml"))

  stopifnot("project.yaml is not up-to-date with create-project.R.  Run create-project.R before running further actions." = identical(thisprojectsplit, currentproject))

# if running manually, output new project as normal
} else if (Sys.getenv("OPENSAFELY_BACKEND") %in% c("")){

  ## output to file ----
  writeLines(thisproject, here("project.yaml"))
  #yaml::write_yaml(project_list, file =here("project.yaml"))
  
  ## grab all action names and send to a txt file
  
  names(actions_list) %>% tibble(action=.) %>%
    mutate(
      model = action==""  & lag(action!="", 1, TRUE),
      model_number = cumsum(model),
    ) %>%
    group_by(model_number) %>%
    summarise(
      sets = str_trim(paste(action, collapse=" "))
    ) %>% pull(sets) %>%
    paste(collapse="\n") %>%
    writeLines(here("actions.txt"))

# fail if backend not recognised
} else {
  stop("Backend not recognised by create.project.R script")
}

