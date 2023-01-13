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

brand = params["brand"]
n_matching_rounds = params["n_matching_rounds"]

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
  population = patients.which_exist_in_file(
    f_path=f"output/sequential/{brand}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz"
    ),

  trial_date = patients.with_value_from_file(
    f_path=f"output/sequential/{brand}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz", 
    returning="trial_date", 
    returning_type="date", 
    date_format='YYYY-MM-DD'
    ),
  
  match_id = patients.with_value_from_file(
    f_path=f"output/sequential/{brand}/matchround{n_matching_rounds}/actual/cumulative_matchedcontrols.csv.gz", 
    returning="match_id", 
    returning_type="int"
    ),

  ###############################################################################
  # outcome variables
  ##############################################################################
  **outcome_variables,

)
