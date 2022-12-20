library(tidyverse)

flowchart_treatedeligible_any_rounded <- read_csv("release20221201/flowchart/flowchart_treatedeligible_any_rounded.csv")

flowchart_treatedeligible_any_rounded %>% 
  transmute(
    criteria,
    n=scales::comma(n, accuracy = 1), 
    n_exclude=scales::comma(n_exclude, accuracy = 1),
    pct_all = round(100*pct_all,1),
    pct_exclude = round(100*pct_exclude,1)
  )

flowchart_matching_any_rounded <- read_csv("release20221201/flowchart/flowchart_matching_any_rounded.csv")

flowchart_matching_any_rounded %>% 
  transmute(
    criteria,
    n=scales::comma(n, accuracy = 1), 
    n_exclude=scales::comma(n_exclude, accuracy = 1),
    pct_all = round(100*pct_all,1),
    pct_exclude = round(100*pct_exclude,1)
  )

total <- flowchart_treatedeligible_any_rounded[1,][["n"]]
eligible <- flowchart_treatedeligible_any_rounded[7,][["n"]]
matched_pfizer <- flowchart_matching_any_rounded[2,][["n"]]
matched_az <- flowchart_matching_any_rounded[5,][["n"]]
unmatched <- total-matched_pfizer-matched_az

100*matched_pfizer/total
100*matched_az/total
100*unmatched/total


