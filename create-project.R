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
action_1matchround <- function(brand, matching_round){
  
  control_extract_date <- study_dates[[brand]][[glue("control_extract_dates")]][matching_round]
  
  if (matching_round == 1) {
    process_controlpotential_hsoutputs <- lst(
      rds = glue("output/sequential/{brand}/matchround{matching_round}/process/*.rds"),
      csv = glue("output/sequential/{brand}/matchround{matching_round}/process/*.csv.gz")
    )
  } else {
    process_controlpotential_hsoutputs <- lst(
      rds = glue("output/sequential/{brand}/matchround{matching_round}/process/*.rds")
    )
  }
  
  splice(
    
    comment(
      "",
      "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
      glue("Matching round {matching_round}:"),
      "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"
      ),
    
    action(
      name = glue("extract_controlpotential_{brand}_{matching_round}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlpotential", 
        " --output-file output/sequential/{brand}/matchround{matching_round}/extract/input_controlpotential.feather", 
        " --param brand={brand}",
        " --param matching_round={matching_round}",
        " --param index_date={control_extract_date}"
      ),
      needs = c(
        "design",
        if(matching_round>1) {
          c(glue("process_controlpotential_{brand}_1"),glue("process_controlactual_{brand}_{matching_round-1}"))
        } else {
            NULL
          }
      ) %>% as.list,
      highly_sensitive = lst(
        cohort = glue("output/sequential/{brand}/matchround{matching_round}/extract/input_controlpotential.feather")
      )
    ),
    
    action(
      name = glue("process_controlpotential_{brand}_{matching_round}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("potential", brand, matching_round),
      needs = namelesslst(
        glue("extract_controlpotential_{brand}_{matching_round}"),
      ),
      highly_sensitive = process_controlpotential_hsoutputs,
      moderately_sensitive = lst(
        input_controlpotential_skim = glue("output/sequential/{brand}/matchround{matching_round}/extract/potential/*.txt"),
        data_processed_skim = glue("output/sequential/{brand}/matchround{matching_round}/potential/*.txt"),
        data_controlpotential_skim = glue("output/sequential/{brand}/matchround{matching_round}/process/*.txt")
      )
    ),
    
    action(
      name = glue("match_potential_{brand}_{matching_round}"),
      run = glue("r:latest analysis/sequential/matching/match_potential.R"),
      arguments = c(brand, matching_round),
      needs = c(
        glue("process_treated"), 
        glue("process_controlpotential_{brand}_{matching_round}"),
        if(matching_round>1) {glue("process_controlactual_{brand}_{matching_round-1}")} else {NULL}
      ) %>% as.list,
      highly_sensitive = lst(
        rds = glue("output/sequential/{brand}/matchround{matching_round}/potential/*.rds"),
        csv = glue("output/sequential/{brand}/matchround{matching_round}/potential/*.csv.gz"),
      )
    ),
    
    action(
      name = glue("extract_controlactual_{brand}_{matching_round}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlactual", 
        " --output-file output/sequential/{brand}/matchround{matching_round}/extract/input_controlactual.feather", 
        " --param brand={brand}",
        " --param matching_round={matching_round}",
      ),
      needs = namelesslst(
        "design",
        glue("match_potential_{brand}_{matching_round}"), 
      ),
      highly_sensitive = lst(
        cohort = glue("output/sequential/{brand}/matchround{matching_round}/extract/input_controlactual.feather")
      )
    ),
    
    
    action(
      name = glue("process_controlactual_{brand}_{matching_round}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("actual", brand, matching_round),
      needs = c(
        glue("process_treated"),
        glue("match_potential_{brand}_{matching_round}"), 
        glue("extract_controlpotential_{brand}_{matching_round}"),  # this is only necessary for the dummy data
        glue("process_controlpotential_{brand}_{matching_round}"), # this is necessary for the vaccine data
        glue("extract_controlactual_{brand}_{matching_round}"),
        if(matching_round>1){glue("process_controlactual_{brand}_{matching_round-1}")} else {NULL}
      ) %>% as.list,
      highly_sensitive = lst(
        rds = glue("output/sequential/{brand}/matchround{matching_round}/actual/*.rds"),
        csv = glue("output/sequential/{brand}/matchround{matching_round}/actual/*.csv.gz"),
      ),
      moderately_sensitive = lst(
        input_controlactual_skim = glue("output/sequential/{brand}/matchround{matching_round}/extract/actual/*.txt"),
        data_actual_skim = glue("output/sequential/{brand}/matchround{matching_round}/actual/*.txt"),
      )
    )

  )
}

# create all necessary actions for n matching rounds
action_extract_and_match <- function(brand, n_matching_rounds) {
  
  allrounds <- map(seq_len(n_matching_rounds), ~action_1matchround(brand, .x)) %>% flatten
  
  splice(
    
    allrounds,
    
    comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            glue("Extract and process data from final controls in the {brand} trials"),
            "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
    
    comment(glue("`extract_controlfinal_{brand}` extracts data from successful matches"),
            "across all matching rounds:"),
    action(
      name = glue("extract_controlfinal_{brand}"),
      run = glue(
        "cohortextractor:latest generate_cohort", 
        " --study-definition study_definition_controlfinal", 
        " --output-file output/sequential/{brand}/extract/input_controlfinal.feather",
        " --param brand={brand}",
        " --param n_matching_rounds={n_matching_rounds}",
      ),
      needs = namelesslst(
        "design",
        glue("process_controlactual_{brand}_{n_matching_rounds}")
      ),
      highly_sensitive = lst(
        extract = glue("output/sequential/{brand}/extract/input_controlfinal.feather")
      )
    ),
    
    comment(glue("`dummydata_controlfinal_{brand}` creates dummy data to represent"),
            glue("the dummy data extraced in `extract_controlfinal_{brand}`(for "),
            "testing only):"),
    action(
      name = glue("dummydata_controlfinal_{brand}"),
      run = glue("r:latest analysis/dummy/dummydata_controlfinal.R"),
      arguments = c(brand),
      needs =map(
        seq_len(n_matching_rounds),
        ~glue("process_controlactual_{brand}_",.x)
      ),
      highly_sensitive = lst(
        dummydata_controlfinal = glue("output/sequential/{brand}/dummydata/dummy_control_final.feather")
      ),
    ),
    
    comment(glue("`process_controlfinal_{brand}` processes the data extracted in"),
            glue("extract_controlfinal_{brand}:")),
    action(
      name = glue("process_controlfinal_{brand}"),
      run = glue("r:latest analysis/process/process_data.R"),
      arguments = c("final", brand),
      needs = c(
        map(
          seq_len(n_matching_rounds),
          ~glue("process_controlactual_{brand}_",.x)
        ),
        glue("extract_controlfinal_{brand}"),
        glue("process_treated"),
        glue("dummydata_controlfinal_{brand}")
      ),
      highly_sensitive = lst(
        extract = glue("output/sequential/{brand}/match/*.rds")
      ),
      moderately_sensitive = lst(
        input_controlfinal_skim = glue("output/sequential/{brand}/extract/*.txt"),
        data_matched_skim = glue("output/sequential/{brand}/match/*.txt")
      )
    )
    
  )
  
}

action_kmcox <- function(brand, subgroup, outcome){
  
  action(
    name = glue("kmcox_{brand}_{subgroup}_{outcome}"),
    run = glue("r:latest analysis/sequential/model/kmcox.R"),
    arguments = c(brand, subgroup, outcome),
    needs = namelesslst(
      glue("process_controlfinal_{brand}"),
    ),
    moderately_sensitive= lst(
      rds= glue("output/sequential/{brand}/model/{subgroup}/{outcome}/*.rds"),
      png= glue("output/sequential/{brand}/model/{subgroup}/{outcome}/*.png"),
    )
  )
  
}

brand_seqtrial <- function(brand) {
  
  if (brand == "pfizer") {
    
    matching_comment <- 
      comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
              "",
              "Due to constraints in the way that data are extracted using the",
              "opensafely cohort extractor (i.e. one-row-per-patient), we",
              glue("conduct the matching over {n_matching_rounds} rounds. Each round (denoted "),
              "{matching_round}) implements the following actions:",
              "",
              "- `extract_controlpotential_pfizer_{matching_round}` extracts data",
              "  from individuals who are potential controls for matching in ",
              "  matching round {matching_round}",
              "",
              "- `process_controlpotential_pfizer_{matching_round}` processes the",
              "  extracted data and applies the eligibility criteria",
              "",
              "- `match_potential_pfizer_{matching_round}` matches the potential",
              "  controls the the treated individuals",
              "",
              "- `extract_controlactual_pfizer_{matching_round}` re-extracts data",
              "  from the individuals who were matched as controls in",
              "  `match_potential_pfizer_{matching_round}`, with data re-defined",
              "  on `trial_date` (the start date fo the sequential trial to which",
              "  they were assigned)",
              "",
              "- `process_controlactual_pfizer_{matching_round}` processes the",
              "  data extracted in `extract_controlactual_pfizer_{matching_round}`",
              "  and checks that the matches made in",
              "  `match_potential_pfizer_{matching_round}` still match bases on",
              "  the re-extracted data",
              "",
              "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #")
  } else {
    
    matching_comment <- 
    comment("",
            "See comment at the start of the pfizer matching round for a",
            "description of the actions in the matching round section.")
    
  }
  
  splice(
    
    comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            glue("Extract control data, match and model for the {brand} trials"),
            "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
    
    matching_comment,
    
    action_extract_and_match(brand, n_matching_rounds),
    
    comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            glue("{brand} trial summary"),
            "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
    
    comment(glue("`coverage_{brand}` summarises the matching coverage on each day of"),
            "the recruitment period:"),
    action(
      name = glue("coverage_{brand}"),
      run = glue("r:latest analysis/sequential/matching/coverage.R"),
      arguments = c(brand),
      needs = namelesslst(
        "process_treated",
        glue("process_controlfinal_{brand}"),
      ),
      moderately_sensitive= lst(
        coverage = glue("output/report/coverage/coverage_{brand}.csv")
      )
    ),
    
    comment(glue("`table1_sequential_{brand}` summarises matching variables and"),
            glue("baseline covariates for individuals included in the {brand} trials:")),
    action(
      name = glue("table1_sequential_{brand}"),
      run = "r:latest analysis/report/table1.R",
      arguments = c("sequential", brand),
      needs = namelesslst(
        "process_treated",
        glue("process_controlfinal_{brand}"),
      ),
      moderately_sensitive= lst(
        table1 = glue("output/report/table1/table1_sequential_{brand}_rounded.csv")
      )
    ),
    
    comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            "Fit models to the sequential trials data",
            "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            "",
            "`kmcox_pfizer_{subgroup}_{outcome}` actions fit models to the",
            paste0(glue("{brand}"), " trial data for a given {subgroup} and {outcome}:")),
    
    expand_grid(
      brand=brand,
      subgroup=model_subgroups,
      outcome=model_outcomes,
    ) %>%
      pmap(
        function(brand, subgroup, outcome) action_kmcox(brand, subgroup, outcome)
      ) %>%
      unlist(recursive = FALSE)
    
  )
  
}

model_single <- function(brand, subgroup, outcome) {
  
  splice(
    
    comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
            glue("Model: brand = {brand}; subgroup = {subgroup}; outcome = {outcome}"),
            "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
    
    action(
      name = glue("preflight_{brand}_{subgroup}_{outcome}_vaccine"),
      run = "r:latest analysis/single/model/preflight.R",
      arguments = c(brand, subgroup, outcome, "vaccine"),
      needs = splice(
        "process_stset",
        "process_data_days"
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/vaccine/*.csv"),
        html = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/vaccine/*.html")
      )
    ),
    
    action(
      name = glue("preflight_{brand}_{subgroup}_{outcome}_outcome"),
      run = "r:latest analysis/single/model/preflight.R",
      arguments = c(brand, subgroup, outcome, "outcome"),
      needs = splice(
        "process_stset",
        "process_data_days"
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/outcome/*.csv"),
        html = glue("output/single/{brand}/{subgroup}/{outcome}/preflight/outcome/*.html")
      )
    ),
    
    action(
      name = glue("ipw_{brand}_{subgroup}_{outcome}"),
      run = "r:latest analysis/single/model/ipw.R",
      arguments = c(brand, subgroup, outcome),
      needs = splice(
        "process_stset",
        "process_data_days"#,
        # glue("preflight_{brand}_{subgroup}_{outcome}_vaccine"),
        # glue("preflight_{brand}_{subgroup}_{outcome}_outcome")
      ),
      highly_sensitive = lst(
        rds = glue("output/single/{brand}/{subgroup}/{outcome}/ipw/*.rds")
      ),
      moderately_sensitive = lst(
        # csv = glue("output/single/{brand}/{subgroup}/{outcome}/ipw/*.csv"),
        svg = glue("output/single/{brand}/{subgroup}/{outcome}/ipw/*.svg"),
        txt = glue("output/single/{brand}/{subgroup}/{outcome}/ipw/*.txt")
      )
    ),
    
    action(
      name = glue("msm_{brand}_{subgroup}_{outcome}"),
      run = "r:latest analysis/single/model/msm.R",
      arguments = c(brand, subgroup, outcome),
      needs = splice(
        glue("ipw_{brand}_{subgroup}_{outcome}")
        # glue("preflight_{brand}_{subgroup}_{outcome}_vaccine"),
        # glue("preflight_{brand}_{subgroup}_{outcome}_outcome"),
        
      ),
      highly_sensitive = lst(
        rds = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.rds")
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.csv")#,
        # svg = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.svg"),
        # txt = glue("output/single/{brand}/{subgroup}/{outcome}/msm/*.txt")
      )
    ),
    
    action(
      name = glue("postprocess_{brand}_{subgroup}_{outcome}"),
      run = "r:latest analysis/single/model/postprocess.R",
      arguments = c(brand, subgroup, outcome),
      needs = namelesslst(
        glue("ipw_{brand}_{subgroup}_{outcome}"),
        glue("msm_{brand}_{subgroup}_{outcome}")
      ),
      highly_sensitive = lst(
        rds = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.rds")
      ),
      moderately_sensitive = lst(
        csv = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.csv"),
        svg = glue("output/single/{brand}/{subgroup}/{outcome}/postprocess/*.svg")
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

  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "DO NOT EDIT project.yaml DIRECTLY",
          "This file is created by create-project.R",
          "Edit and run create-project.R to update the project.yaml",
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
           " "
          ),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
          "PRELIIMINARIES", 
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "`design` defines study metadata:"),
  
  action(
    name = "design",
    run = glue("r:latest analysis/design.R"),
    moderately_sensitive = lst(
      lib = glue("lib/design/study-dates.json")
    ),
  ),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
          "SEQUENTIAL TRIAL APPROACH", 
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "`extract_treated` extracts data from individuals who received a",
          "vaccine dose during the study recruitment period:"),
  
  action(
    name = "extract_treated",
    run = glue(
      "cohortextractor:latest generate_cohort", 
      " --study-definition study_definition_treated", 
      " --output-file output/sequential/treated/extract/input_treated.feather",
    ),
    needs = namelesslst(
      "design",
      "process_single"
    ),
    highly_sensitive = lst(
      extract = "output/sequential/treated/extract/input_treated.feather"
    ),
  ),
  
  comment("`process_treated` processes data and apply eligibility criteria:"),
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

  brand_seqtrial("pfizer"),
  brand_seqtrial("az"),
  
  comment("`combine_kmcox` combines output from all actions that run `kmcox.R`:"),
  
  action(
    name = "combine_kmcox",
    run = glue("r:latest analysis/sequential/model/kmcox_combine.R"),
    needs = splice(
      as.list(
        glue_data(
          .x=expand_grid(
            brand=model_brands,
            subgroup=model_subgroups,
            outcome=model_outcomes,
          ),
          "kmcox_{brand}_{subgroup}_{outcome}"
        )
      )
    ),
    moderately_sensitive = lst(
      rds = "output/sequential/combine/*.csv",
      png = "output/sequential/combine/*.png"
    )
  ),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
          "SINGLE TRIAL APPROACH", 
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
  comment("Extract and process data", 
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "`process_single` processes data extracted in the ", 
          "`extract_controlpotential_pfizer_1` action for the",
          "single trial approach:"),
  
  action(
    name = "process_single",
    run = "r:latest analysis/process/process_data.R",
    arguments = "single",
    needs = namelesslst(
      # "extract_treated",
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
  
  comment("`table1_single_any` calculates summary statistics for the single", 
          "trial cohort:"), 
  
  action(
    name = "table1_single_any",
    run = "r:latest analysis/report/table1.R",
    arguments = c("single", "any"),
    needs = namelesslst(
      "process_single"
    ),
    moderately_sensitive= lst(
      table1 = glue("output/report/table1/table1_single_any_rounded.csv")
    )
  ),
  
  comment("`extract_timevarying` extracts the data needed to derive", 
          "time-varying covariates and outcome variables for the single trial",
          "approach:"), 
  
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
  
  comment("`dummydata_timevarying` creates dummy data represent the data", 
          "extracted in `extract_timevarying` (for testing only):"), 
  
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
  
  comment("`process_timevarying` processes the data extracted in ", 
          "`extract_timevarying` to create time-varying covariates and outcome",
          "variables:"), 
  
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
  
  comment("`process_stset` creates time-to-event datasets necessary to derive", 
          "data_days (one-row-per-patient-per-day):"), 
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
  
  comment("`process_data_days` creates a one-row-per-patient-per-day dataset;",
          "this is split across `process_data_days_n` defined in",
          "due to memory constraints:"), 
  
  action(
    name = "process_data_days",
    run = "r:latest analysis/single/process/process_data_days.R",
    needs = namelesslst(
      "process_stset"
    ),
    highly_sensitive = lst(
      processed = "output/single/stset/data_days_*.rds"
    )
  ),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "Fit models to the single trials data",
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "The actions in this section do the following:",
          "",
          "- `preflight_{brand}_{subgroup}_{outcome}_{stage}",
          "  checks for separation between covariates and outcomes",
          "",
          "- `ipw_{brand}_{subgroup}_{outcome} fits models to calculate the",
          "  probability of vaccination models",
          "",
          "- `msm_{brand}_{subgroup}_{outcome} fits marginal structural models",
          "   to predict the odds of outcome events",
          "",
          "- `postprocess_{brand}_{subgroup}_{outcome}` processes the",
          "  output from the `ipw_preflight_{brand}_{subgroup}_{outcome}`",
          "  and `msm_preflight_{brand}_{subgroup}_{outcome}` actions",
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          ""
  ),
  
  # model actions
  expand_grid(
    brand=model_brands,
    subgroup=model_subgroups,
    outcome=model_outcomes,
  ) %>%
    pmap(
      function(brand, subgroup, outcome) model_single(brand, subgroup, outcome)
    ) %>%
    unlist(recursive = FALSE),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
          "Combine model outputs",
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "`msm_combine` combines the output from all actions that run `msm.R`:"),
  action(
    name = "single_combine",
    run = glue("r:latest analysis/single/model/combine.R"),
    needs = splice(
      as.list(
        glue_data(
          .x=expand_grid(
            brand=model_brands,
            subgroup=model_subgroups,
            outcome=model_outcomes,
          ),
          "postprocess_{brand}_{subgroup}_{outcome}"
        )
      )
    ),
    moderately_sensitive = lst(
      csv = "output/single/combine/*.csv",
      svg = "output/single/combine/*.svg"
    )
  ),
  
  comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
          "REPORT", 
          "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #",
          "",
          "`flowchart` prepares the data to be used in the particiant flow",
          "diagram in the paper (Supplementary Figure xxx):"),
  
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
      check_NAs = "output/report/flowchart/*.txt",
      flow_matching = "output/report/flowchart/*.csv"
    )
  ),
  
  comment("`brand12counts` plots the cumulative incidence of first and second",
  "vaccine doses (Figure xxx):"),
  
  action(
    name = "brand12counts",
    run = glue("r:latest analysis/report/brand12counts.R"),
    needs = namelesslst(
      "process_data_days"
    ),
    moderately_sensitive = lst(
      csv = "output/report/brand12counts/*.csv",
      plots = "output/report/brand12counts/*.png"
    )
  ),
  
  comment("Check vax data mismatch:"),
  action(
    name = "check_vax_data",
    run = glue("r:latest analysis/checks/vax_data.R"),
    needs = namelesslst(
      "extract_treated",
      "extract_controlpotential_pfizer_1"
    ),
    moderately_sensitive = lst(
      vax_data = "output/checks/vax_data.txt"
    )
  ),
  
  # 
  # comment("# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #", 
  #         "Move files for release", 
  #         "# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #"),
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

