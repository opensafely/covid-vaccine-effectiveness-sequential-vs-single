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
# vax variables
from variables_inclusion import generate_inclusion_variables 
inclusion_variables = generate_inclusion_variables(index_date="vax1_date")
############################################################
## jcvi variables
from variables_jcvi import generate_jcvi_variables 
jcvi_variables = generate_jcvi_variables(index_date="vax1_date")
############################################################
## demographic variables
from variables_demo import generate_demo_variables 
demo_variables = generate_demo_variables(index_date="vax1_date")
############################################################
## pre variables
from variables_pre import generate_pre_variables 
pre_variables = generate_pre_variables(index_date="vax1_date")
############################################################
## outcome variables
from variables_outcome import generate_outcome_variables 
outcome_variables = generate_outcome_variables(index_date="vax1_date")
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
    age31aug2020 >= 70
    AND
    NOT has_died
    AND 
    eligible_single
    AND
    vax1_date
    """,
    
    **inclusion_variables,   

    eligible_single = patients.which_exist_in_file(
      f_path="output/single/eligible/data_singleeligible.csv.gz"
      ),
  ),
  
  #################################################################
  ## Covid vaccine dates
  #################################################################
  vax1_date = patients.with_value_from_file(
    f_path="output/single/eligible/data_singleeligible.csv.gz", 
    returning="vax1_date", 
    returning_type="date", 
    date_format='YYYY-MM-DD'
    ),
  vax2_date = patients.with_value_from_file(
    f_path="output/single/eligible/data_singleeligible.csv.gz", 
    returning="vax2_date", 
    returning_type="date", 
    date_format='YYYY-MM-DD'
    ),
  vax1_type = patients.with_value_from_file(
    f_path="output/single/eligible/data_singleeligible.csv.gz", 
    returning="vax1_type", 
    returning_type="str", 
    ),
  vax2_type = patients.with_value_from_file(
    f_path="output/single/eligible/data_singleeligible.csv.gz", 
    returning="vax2_type", 
    returning_type="str", 
    ),
    
  ###############################################################################
  # jcvi variables
  ###############################################################################
  **jcvi_variables, 
  
  ###############################################################################
  # demographic variables
  ###############################################################################
  **demo_variables,   

  ###############################################################################
  # pre variables
  ###############################################################################
  **pre_variables,      
  
  ###############################################################################
  # posts
  ###############################################################################
  **outcome_variables,      
  
)
