# Import codelists from codelists.py
import codelists

# import json module
import json

from cohortextractor import (
  StudyDefinition,
  patients,
  codelist_from_csv,
  codelist,
  filter_codes_by_category,
  combine_codelists,
  params
)

cohort = params["cohort"]
n_matching_rounds = params["n_matching_rounds"]

# import study dates defined in "./analysis/design.R" script
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)

# change these in design.R if necessary
index_date = study_dates["index_date"]


############################################################
## outcome variables
from variables_outcome import generate_outcome_variables 
outcome_variables = generate_outcome_variables(index_date="trial_date")
############################################################


# Specify study defeinition
study = StudyDefinition(
  
  # Configure the expectations framework
  default_expectations={
    "date": {"earliest": "2020-01-01", "latest": "today"},
    "rate": "uniform",
    "incidence": 0.2,
    "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
    "float": {"distribution": "normal", "mean": 25, "stddev": 5},
  },
  
  # This line defines the study population
  population = patients.which_exist_in_file(f_path=f"output/{cohort}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz"),

  trial_date = patients.with_value_from_file(f_path=f"output/{cohort}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz", returning="trial_date", returning_type="date", date_format='YYYY-MM-DD'),
  
  match_id = patients.with_value_from_file(f_path=f"output/{cohort}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz", returning="match_id", returning_type="int"),

  
  ###############################################################################
  # post variables
  ##############################################################################
  **outcome_variables,
  
  
)
