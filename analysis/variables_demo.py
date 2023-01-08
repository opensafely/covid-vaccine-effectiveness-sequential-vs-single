from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_demo_variables(index_date):

  demo_variables = dict(

  has_follow_up_previous_year=patients.registered_with_one_practice_between(
    start_date=f"{index_date} - 1 year",
    end_date=f"{index_date} - 1 day",
  ),

  age31aug2020=patients.age_as_of( 
    "2020-08-31",
  ),

  age=patients.age_as_of( 
    f"{index_date} - 1 day",
  ),
    
  sex=patients.sex(
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
      "incidence": 1,
    }
  ),

  # Ethnicity in 6 categories
  ethnicity = patients.with_these_clinical_events(
    codelists.ethnicity,
    returning="category",
    find_last_match_in_period=True,
    include_date_of_match=False,
    return_expectations={
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.75,
    },
  ),
  
  # ethnicity variable that takes data from SUS
  ethnicity_6_sus = patients.with_ethnicity_from_sus(
    returning="group_6",  
    use_most_frequent_code=True,
    return_expectations={
      "category": {"ratios": {"1": 0.2, "2": 0.2, "3": 0.2, "4": 0.2, "5": 0.2}},
      "incidence": 0.8,
    },
  ),
  
  ################################################################################################
  ## Practice and patient ID variables
  ################################################################################################

  # NHS administrative region
  region=patients.registered_practice_as_of(
    f"{index_date} - 1 day",
    returning="nuts1_region_name",
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.2,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.1,
          "South East": 0.1,
          "South West": 0.1
          #"" : 0.01
        },
      },
    },
  ),

  ## IMD - quintile
  imd_Q5=patients.categorised_as(
    {
      "Unknown": "DEFAULT",
      "1 (most deprived)": "imd >= 0 AND imd < 32844*1/5",
      "2": "imd >= 32844*1/5 AND imd < 32844*2/5",
      "3": "imd >= 32844*2/5 AND imd < 32844*3/5",
      "4": "imd >= 32844*3/5 AND imd < 32844*4/5",
      "5 (least deprived)": "imd >= 32844*4/5 AND imd <= 32844",
    },
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"Unknown": 0.02, "1 (most deprived)": 0.18, "2": 0.2, "3": 0.2, "4": 0.2, "5 (least deprived)": 0.2}},
    },
    imd=patients.address_as_of(
      f"{index_date} - 1 day",
      returning="index_of_multiple_deprivation",
      round_to_nearest=100,
      return_expectations={
      "category": {"ratios": {c: 1/320 for c in range(100, 32100, 100)}}
      }
    ),
  ),

  # flu vaccine in flu seasons 16-17, 17-18, 18-19, 19-20 or 20 (only up to 2020-12-08)
  flu_vaccine=patients.satisfying(
    """
    flu_vaccine_tpp_table>0 OR
    flu_vaccine_med>0 OR
    flu_vaccine_clinical>0
    """,
        
    flu_vaccine_tpp_table=patients.with_tpp_vaccination_record(
      target_disease_matches="INFLUENZA",
      between=["2016-07-01", "2020-12-08"], 
      returning="binary_flag",
      ),
        
    flu_vaccine_med=patients.with_these_medications(
        flu_med_codes,
        between=["2016-07-01", "2020-12-08"], 
        returning="binary_flag",
      ),
    flu_vaccine_clinical=patients.with_these_clinical_events(
        flu_clinical_given_codes,
        ignore_days_where_these_codes_occur=flu_clinical_not_given_codes,
        between=["2016-07-01", "2020-12-08"], 
        returning="binary_flag",
        ),
      return_expectations={"incidence": 0.5},
    ),

  )
  return demo_variables

