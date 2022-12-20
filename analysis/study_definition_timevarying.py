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

# import study dates defined in "./analysis/design.R" script
with open("./lib/design/study-dates.json") as f:
  study_dates = json.load(f)
index_date = study_dates['pfizer']['start_date']

############################################################
## outcome variables
from variables_outcome import generate_outcome_variables 
outcome_variables = generate_outcome_variables(index_date="index_date")
############################################################
## timevarying variables
from variables_timevarying import generate_timevarying_variables 
timevarying_variables = generate_timevarying_variables(index_date="index_date", n=6)
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
  population = patients.which_exist_in_file(f_path="output/single/eligible/data_singleeligible.csv.gz"),
  
  index_date = index_date,
  
  ###############################################################################
  # time varying covariates
  ##############################################################################
  **timevarying_variables,

  ###############################################################################
  # outcome variables
  ##############################################################################
  **outcome_variables,
  
)
