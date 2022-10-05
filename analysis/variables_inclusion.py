from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_inclusion_variables(index_date):
  inclusion_variables = dict(
    
    registered = patients.registered_as_of(
        f"{index_date} - 1 day",
    ), 

    has_died = patients.died_from_any_cause(
      on_or_before=f"{index_date} - 1 day",
      returning="binary_flag",
    ),
          
  )
  return inclusion_variables

