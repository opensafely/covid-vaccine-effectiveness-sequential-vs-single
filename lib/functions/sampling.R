# functions for sub-sampling ----

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
  # `proportion` is the proportion of nonoutcome patients to be sampled

  (dplyr::dense_rank(dplyr::if_else(had_outcome, 0L, id)) - 1L) <= n

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
