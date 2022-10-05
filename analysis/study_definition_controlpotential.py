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

# define params
cohort = params["cohort"]
matching_round = params["matching_round"]
previousmatching_round = int(matching_round)-1
index_date = params["index_date"]


############################################################
## inclusion variables
from variables_vax import generate_vax_variables 
vax_variables = generate_vax_variables(index_date="1900-01-01", n=3)
############################################################
# vax variables
from variables_inclusion import generate_inclusion_variables 
inclusion_variables = generate_inclusion_variables(index_date="index_date")
############################################################
## jcvi variables
from variables_jcvi import generate_jcvi_variables 
jcvi_variables = generate_jcvi_variables(index_date="index_date")
############################################################
## demographic variables
from variables_demo import generate_demo_variables 
demo_variables = generate_demo_variables(index_date="index_date")
############################################################
## pre variables
from variables_pre import generate_pre_variables 
pre_variables = generate_pre_variables(index_date="index_date")



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
  
  index_date = index_date,
  
  # This line defines the study population
  population=patients.satisfying(
    """
    registered
    AND
    age >= 18
    AND
    NOT has_died
    AND 
    covid_vax_disease_2_date
    """,
    
    **inclusion_variables,    

  ),
  
  #################################################################
  ## Covid vaccine dates
  #################################################################
  **vax_variables,
  # all vaccination variables for first three doses to apply selection criteria
    
  ###############################################################################
  # jcvi variables
  ##############################################################################
  **jcvi_variables, 
  
  ###############################################################################
  # demographic variables
  ##############################################################################
  **demo_variables,   

  ###############################################################################
  # pre variables
  ##############################################################################
  **pre_variables,    

)
