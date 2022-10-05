from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_matching_variables(index_date):

  matching_variables = dict(

  has_follow_up_previous_6weeks=patients.registered_with_one_practice_between(
    start_date=f"{index_date} - 42 days",
    end_date=f"{index_date} - 1 day",
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
  # practice pseudo id
  practice_id=patients.registered_practice_as_of(
    f"{index_date} - 1 day",
    returning="pseudo_id",
    return_expectations={
      "int": {"distribution": "normal", "mean": 1000, "stddev": 100},
      "incidence": 1,
    },
  ),
  
  # msoa
  msoa=patients.address_as_of(
    f"{index_date} - 1 day",
    returning="msoa",
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {"E02000001": 0.0625, "E02000002": 0.0625, "E02000003": 0.0625, "E02000004": 0.0625,
        "E02000005": 0.0625, "E02000007": 0.0625, "E02000008": 0.0625, "E02000009": 0.0625, 
        "E02000010": 0.0625, "E02000011": 0.0625, "E02000012": 0.0625, "E02000013": 0.0625, 
        "E02000014": 0.0625, "E02000015": 0.0625, "E02000016": 0.0625, "E02000017": 0.0625}},
    },
  ),    

  # stp is an NHS administration region based on geography
  stp=patients.registered_practice_as_of(
    f"{index_date} - 1 day",
    returning="stp_code",
    return_expectations={
      "rate": "universal",
      "category": {
        "ratios": {
          "STP1": 0.1,
          "STP2": 0.1,
          "STP3": 0.1,
          "STP4": 0.1,
          "STP5": 0.1,
          "STP6": 0.1,
          "STP7": 0.1,
          "STP8": 0.1,
          "STP9": 0.1,
          "STP10": 0.1,
        }
      },
    },
  ),

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

  #rurality
  rural_urban=patients.address_as_of(
    f"{index_date} - 1 day",
    returning="rural_urban_classification",
    return_expectations={
      "rate": "universal",
      "category": {"ratios": {1: 0.125, 2: 0.125, 3: 0.125, 4: 0.125, 5: 0.125, 6: 0.125, 7: 0.125, 8: 0.125}},
    },
  ),

  ################################################################################################
  ## Pre-study event dates
  ################################################################################################

  # any covid test 
  covid_test_0_date=patients.with_test_result_in_sgss(
    pathogen="SARS-CoV-2",
    test_result="any",
    on_or_before=f"{index_date} - 1 day",
    returning="date",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
    restrict_to_earliest_specimen_date=False,
  ),

  # positive covid test
  positive_test_0_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="positive",
      returning="date",
      date_format="YYYY-MM-DD",
      on_or_before=f"{index_date} - 1 day",
      # no earliest date set, which assumes any date errors are for tests occurring before study start date
      find_last_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
  ),

  # any emergency attendance for covid
  covidemergency_0_date=patients.attended_emergency_care(
    returning="date_arrived",
    on_or_before=f"{index_date} - 1 day",
    with_these_diagnoses = codelists.covid_emergency,
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),

  # unplanned hospital admission
  admitted_unplanned_0_date=patients.admitted_to_hospital(
    returning="date_admitted",
    on_or_before=f"{index_date} - 1 day",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_patient_classification = ["1"], # ordinary admissions only
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),
  discharged_unplanned_0_date=patients.admitted_to_hospital(
    returning="date_discharged",
    on_or_before=f"{index_date} - 1 day",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_patient_classification = ["1"], # ordinary admissions only
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ), 

  # planned hospital admission
  admitted_planned_0_date=patients.admitted_to_hospital(
    returning="date_admitted",
    on_or_before=f"{index_date} - 1 day",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["11", "12", "13", "81"],
    with_patient_classification = ["1"], # ordinary admissions only 
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),
  discharged_planned_0_date=patients.admitted_to_hospital(
    returning="date_discharged",
    on_or_before=f"{index_date} - 1 day",
    # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["11", "12", "13", "81"],
    with_patient_classification = ["1"], # ordinary admissions only
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True
  ), 
  
  # Positive covid admission prior to study start date
  admitted_covid_0_date=patients.admitted_to_hospital(
    returning="date_admitted",
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
    with_these_diagnoses=codelists.covid_icd10,
    on_or_before=f"{index_date} - 1 day",
    date_format="YYYY-MM-DD",
    find_last_match_in_period=True,
  ),

  # Positive case identification prior to study start date
  primary_care_covid_case_0_date=patients.with_these_clinical_events(
    combine_codelists(
      codelists.covid_primary_care_code,
      codelists.covid_primary_care_positive_test,
      codelists.covid_primary_care_sequelae,
    ),
    returning="date",
    date_format="YYYY-MM-DD",
    on_or_before=f"{index_date} - 1 day",
    find_last_match_in_period=True,
  ),

    # test frequency in six months prior to earliest start date in cohort
    prior_covid_test_frequency=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="any",
      between=[f"{index_date} - 182 days", f"{index_date} - 1 day"], # 182 days = 26 weeks
      returning="number_of_matches_in_period", 
      date_format="YYYY-MM-DD",
      restrict_to_earliest_specimen_date=False,
    ),
  
  )
  return matching_variables

