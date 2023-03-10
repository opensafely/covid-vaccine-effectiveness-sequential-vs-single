
# combine here() and glue() functionality
ghere <- function(...){
  here::here(glue::glue(..., .sep=.Platform$file.sep))
}

ceiling_any <- function(x, to=1){
  # round to nearest 100 millionth to avoid floating point errors
  ceiling(plyr::round_any(x/to, 1/100000000))*to
}

roundmid_any <- function(x, to=1){
  # like ceiling_any, but centers on (integer) midpoint of the rounding points
  ceiling(x/to)*to - (floor(to/2)*(x!=0))
}


fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

# for relabelling variables
# use like this:
# fct_recoderelevel(variable_coded,  c(`code1`="full name 1", `code2` = "full name 2"))
fct_recoderelevel <- function(x, lookup){
  stopifnot(!is.na(names(lookup)))
  factor(x, levels=lookup, labels=names(lookup))
}

my_skim <- function(
  .data, # dataset to be summarised
  path,
  id_suffix = "_id" # (set to NULL if no id columns)
) {
  
  # specify summary function for each class
  my_skimmers <- list(
    logical = skimr::sfl(
    ),
    # numeric applied to numeric and integer
    numeric = skimr::sfl(
      mean = ~ mean(.x, na.rm=TRUE),
      sd = ~ sd(.x, na.rm=TRUE),
      min = ~ min(.x, na.rm=TRUE),
      p10 = ~ quantile(.x, p=0.1, na.rm=TRUE, type=1),
      p25 = ~ quantile(.x, p=0.25, na.rm=TRUE, type=1),
      p50 = ~ quantile(.x, p=0.5, na.rm=TRUE, type=1),
      p75 = ~ quantile(.x, p=0.75, na.rm=TRUE, type=1),
      p90 = ~ quantile(.x, p=0.9, na.rm=TRUE, type=1),
      max = ~ max(.x, na.rm=TRUE)
    ),
    character = skimr::sfl(),
    factor = skimr::sfl(),
    Date = skimr::sfl(
      # wrap in as.Date to avoid errors when all missing
      min = ~ as.Date(min(.x, na.rm=TRUE)),
      p50 = ~ as.Date(quantile(.x, p=0.5, na.rm=TRUE, type=1)),
      max = ~ as.Date(max(.x, na.rm=TRUE))
    ),
    POSIXct = skimr::sfl(
      # wrap in as.POSIXct to avoid errors when all missing
      min = ~ as.POSIXct(min(.x, na.rm=TRUE)),
      p50 = ~ as.POSIXct(quantile(.x, p=0.5, na.rm=TRUE, type=1)),
      max = ~ as.POSIXct(max(.x, na.rm=TRUE))
    )
  )
  
  my_skim_fun <- skimr::skim_with(
    !!!my_skimmers,
    append = FALSE
  )
  
  # summarise factors as the printing is not very nice or flexible in skim
  summarise_factor <- function(var) {
    
    out <- .data %>%
      group_by(across(all_of(var))) %>%
      count() %>%
      ungroup() %>%
      mutate(across(n, ~roundmid_any(.x, to = 7))) %>%
      mutate(percent = round(100*n/sum(n),2)) %>%
      arrange(!! sym(var)) 
    
    total <- nrow(out)
    
    out %>%
      slice(1:min(total, 10)) %>% 
      knitr::kable(
        format = "pipe",
        caption = glue::glue("{min(total, 10)} of {total} factor levels printed")
        ) %>% 
      print()
    
  }
  
  vars <- .data %>% 
    select(-ends_with(id_suffix)) %>% 
    select(where(~ is.factor(.x) | is.character(.x))) %>%
    names()
  
  options(width = 120)
  capture.output(
    {
      cat("The following id variables are removed from this summary:\n")
      print(.data %>% select(ends_with(id_suffix)) %>% names())
      cat("\n")
      print(my_skim_fun(.data, -ends_with(id_suffix)))
      cat("\n")
      cat("--- counts for factor and character variables ---")
      for (v in vars) {
        summarise_factor(v)
      }
    },
    file = path,
    append = FALSE
  )
  
}

# functions for sampling ----

# function to sample non-outcome patients
sample_nonoutcomes_prop <- function(had_outcome, id, proportion){
  # TRUE if outcome occurs,
  # TRUE with probability of `prop` if outcome does not occur
  # FALSE with probability `prop` if outcome does occur
  # based on `id` to ensure consistency of samples
  
  # `had_outcome` is a boolean indicating if the subject has experienced the outcome or not
  # `id` is a identifier with the following properties:
  # - a) consistent between cohort extracts
  # - b) unique
  # - c) completely randomly assigned (no correlation with practice ID, age, registration date, etc etc) which should be true as based on hash of true IDs
  # - d) is an integer strictly greater than zero
  # `proportion` is the proportion of nonoutcome patients to be sampled
  
  (dplyr::dense_rank(dplyr::if_else(had_outcome, 0L, id)) - 1L) <= ceiling(sum(!had_outcome)*proportion)
  
}

sample_nonoutcomes_n <- function(had_outcome, id, n){
  # TRUE if outcome occurs,
  # TRUE with probability of `prop` if outcome does not occur
  # FALSE with probability `prop` if outcome does occur
  # based on `id` to ensure consistency of samples
  
  # `had_outcome` is a boolean indicating if the subject has experienced the outcome or not
  # `id` is a identifier with the following properties:
  # - a) consistent between cohort extracts
  # - b) unique
  # - c) completely randomly assigned (no correlation with practice ID, age, registration date, etc etc) which should be true as based on hash of true IDs
  # - d) is an integer strictly greater than zero
  # `n` is the number of nonoutcome patients to be sampled
  if(!any(had_outcome)) {dplyr::dense_rank(id) <= n}
  else {  (dplyr::dense_rank(dplyr::if_else(had_outcome, 0L, id)) - 1L) <= n}
  
}


sample_random_prop <- function(id, proportion){
  # TRUE with probability of `prop`
  # FALSE with probability `prop`
  # based on `id` to ensure consistency of samples
  
  # `id` is a identifier with the following properties:
  # - a) consistent between cohort extracts
  # - b) unique
  # - c) completely randomly assigned (no correlation with practice ID, age, registration date, etc etc) which should be true as based on hash of true IDs
  # - d) is an integer strictly greater than zero
  # `proportion` is the proportion patients to be sampled
  
  dplyr::dense_rank(id) <= ceiling(length(id)*proportion)
}

sample_random_n <- function(id, n){
  # select n rows
  # based on `id` to ensure consistency of samples
  
  # `id` is a identifier with the following properties:
  # - a) consistent between cohort extracts
  # - b) unique
  # - c) completely randomly assigned (no correlation with practice ID, age, registration date, etc etc) which should be true as based on hash of true IDs
  # - d) is an integer strictly greater than zero
  # `proportion` is the proportion patients to be sampled
  
  dplyr::dense_rank(id) <= n
}


sample_weights <- function(had_outcome, sampled){
  # `had_outcome` is a boolean indicating if the subject has experienced the outcome or not
  # `sampled` is a boolean indicating if the patient is to be sampled or not
  case_when(
    had_outcome ~ 1,
    !had_outcome & !sampled ~ 0,
    !had_outcome & sampled ~ sum(!had_outcome)/sum((sampled) & !had_outcome),
    TRUE ~ NA_real_
  )
}






tidy_parglm <- function(x, conf.int = FALSE, conf.level = .95,
                        exponentiate = FALSE, ...) {
  
  # nicked from https://github.com/tidymodels/broom/blob/443ebd995760c6674f122d75ceb2e6b82f055439/R/stats-glm-tidiers.R#L12
  # and adapted for parglm glm objects
  
  ret <- as_tibble(summary(x)$coefficients, rownames = "term")
  colnames(ret) <- c("term", "estimate", "std.error", "statistic", "p.value")
  
  # summary(x)$coefficients misses rank deficient rows (i.e. coefs that summary.lm() sets to NA)
  # catch them here and add them back
  
  coefs <- tibble::enframe(stats::coef(x), name = "term", value = "estimate")
  ret <- left_join(coefs, ret, by = c("term", "estimate"))
  
  if (conf.int) {
    ci <- confint.default(x, level = conf.level) # not ideal -- change for more robust conf intervals! - see tidy_plr below
    ci <- as_tibble(ci, rownames = "term")
    names(ci) <- c("term", "conf.low", "conf.high")
    
    ret <- dplyr::left_join(ret, ci, by = "term")
  }
  
  if (exponentiate) {
    ret <- exponentiate(ret)
  }
  
  ret
}

tidy_plr <- function(model, conf.int=TRUE, conf.level=0.95, exponentiate=FALSE, cluster){
  # create tidy dataframe for coefficients of pooled logistic regression
  mod_tidy <- tidy_parglm(model, conf.int=conf.int, conf.level=conf.level, exponentiate=exponentiate)
  robustSEs <- coeftest(model, vcov. = vcovCL(model, cluster = cluster, type = "HC0")) %>% broom::tidy()
  robustCIs <- coefci(model, vcov. = vcovCL(model, cluster = cluster, type = "HC0")) %>% as_tibble(rownames="term")
  robust <- inner_join(robustSEs, robustCIs, by="term")
  
  robust %>%
    rename(
      conf.low=`2.5 %`,
      conf.high=`97.5 %`
    ) %>%
    mutate(
      or = exp(estimate),
      or.ll = exp(conf.low),
      or.ul = exp(conf.high),
    )
  
}


formatpercent100 <- function(x,accuracy){
  formatx <- scales::label_percent(accuracy)(x)
  
  if_else(
    formatx==scales::label_percent(accuracy)(1),
    paste0(">",scales::label_percent(1)((100-accuracy)/100)),
    formatx
  )
}

# define function to add descriptive variables ----
add_descr <- function(.data) {
  
  brand_descr <- brand_lookup %>% filter(brand %in% model_brands) %>% pull(brand_descr)
  subgroup_descr <- names(recoder$subgroups[recoder$subgroups %in% model_subgroups])
  outcome_descr <- events_lookup %>% filter(event %in% model_outcomes) %>% pull(event_descr)
  
  .data %>%
    mutate( 
      brand_descr = factor(brand, levels = model_brands, labels = brand_descr),
      subgroup_descr = factor(subgroup, levels = model_subgroups, labels = subgroup_descr),
      outcome_descr = factor(outcome, levels = model_outcomes, labels = outcome_descr),
    ) %>%
    mutate(across(brand, factor, levels = model_brands)) %>%
    mutate(across(subgroup, factor, levels = model_subgroups)) %>%
    mutate(across(outcome, factor, levels = model_outcomes)) 
  
}
