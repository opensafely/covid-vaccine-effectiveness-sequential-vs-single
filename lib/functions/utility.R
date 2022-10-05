
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

