
from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists

############################################################
## functions
from variables_functions import *
############################################################

def generate_vax_variables(index_date, n):

  vax_variables = dict(

# pfizer
  **vaccination_date_X(
    name = "covid_vax_pfizer",
    # use 1900 to capture all possible recorded covid vaccinations, including date errors
    # any vaccines occurring before national rollout are later excluded
    index_date = index_date, 
    n = n,
    product_name_matches="COVID-19 mRNA Vaccine Comirnaty 30micrograms/0.3ml dose conc for susp for inj MDV (Pfizer)"
  ),
  
  # az
  **vaccination_date_X(
    name = "covid_vax_az",
    index_date = index_date,
    n = n,
    product_name_matches="COVID-19 Vac AstraZeneca (ChAdOx1 S recomb) 5x10000000000 viral particles/0.5ml dose sol for inj MDV"
  ),
  
  # moderna
  **vaccination_date_X(
    name = "covid_vax_moderna",
    index_date = index_date,
    n = n,
    product_name_matches="COVID-19 mRNA Vaccine Spikevax (nucleoside modified) 0.1mg/0.5mL dose disp for inj MDV (Moderna)"
  ),
  
  # any covid vaccine
  **vaccination_date_X(
    name = "covid_vax_disease",
    index_date = index_date,
    n = n,
    target_disease_matches="SARS-2 CORONAVIRUS"
  ),
  
  )
  return vax_variables

