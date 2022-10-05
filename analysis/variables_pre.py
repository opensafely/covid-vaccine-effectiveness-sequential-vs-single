from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists


def generate_pre_variables(index_date):

  pre_variables = dict(

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
  return pre_variables

