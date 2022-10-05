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

############################################################
## inclusion variables
from variables_vax import generate_vax_variables 
vax_variables = generate_vax_variables(index_date="1900-01-01", n=4)
############################################################
# vax variables
from variables_inclusion import generate_inclusion_variables 
inclusion_variables = generate_inclusion_variables(index_date="covid_vax_disease_3_date")
############################################################
## jcvi variables
from variables_jcvi import generate_jcvi_variables 
jcvi_variables = generate_jcvi_variables(index_date="covid_vax_disease_3_date")
############################################################
## demographic variables
from variables_demo import generate_demo_variables 
demo_variables = generate_demo_variables(index_date="covid_vax_disease_3_date")
############################################################
## pre variables
from variables_pre import generate_pre_variables 
pre_variables = generate_pre_variables(index_date="covid_vax_disease_3_date")
############################################################
## outcome variables
from variables_outcome import generate_outcome_variables 
outcome_variables = generate_outcome_variables(index_date="covid_vax_disease_3_date")
############################################################

# Specify study definition
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
  population=patients.satisfying(
    """
    registered
    AND
    age >= 18
    AND
    NOT has_died
    AND 
    covid_vax_disease_3_date
    """,
    
    **inclusion_variables,    

  ),
  
  #################################################################
  ## Covid vaccine dates
  #################################################################
  **vax_variables,
    
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
  
  ###############################################################################
  # posts
  ##############################################################################
  **outcome_variables,      
  
)
