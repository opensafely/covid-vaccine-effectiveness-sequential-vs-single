from cohortextractor import patients, combine_codelists
from codelists import *
import codelists


def generate_outcome_variables(index_date):
  outcome_variables = dict(
  
    # deregistration date
    dereg_date=patients.date_deregistered_from_all_supported_practices(
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
    ),
  
    # Positive case identification after study start date
    primary_care_covid_case_date=patients.with_these_clinical_events(
      combine_codelists(
        codelists.covid_primary_care_code,
        codelists.covid_primary_care_positive_test,
        codelists.covid_primary_care_sequelae,
      ),
      returning="date",
      date_format="YYYY-MM-DD",
      on_or_after=index_date,
      find_first_match_in_period=True,
    ),
    
    
    # covid PCR test dates from SGSS
    covid_test_date=patients.with_test_result_in_sgss(
      pathogen="SARS-CoV-2",
      test_result="any",
      on_or_after=index_date,
      find_first_match_in_period=True,
      restrict_to_earliest_specimen_date=False,
      returning="date",
      date_format="YYYY-MM-DD",
    ),
    
    # positive covid test
    postest_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_after=index_date,
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
    ),
    
    # emergency attendance for covid, as per discharge diagnosis
    covidemergency_date=patients.attended_emergency_care(
      returning="date_arrived",
      date_format="YYYY-MM-DD",
      on_or_after=index_date,
      with_these_diagnoses = codelists.covid_emergency,
      find_first_match_in_period=True,
    ),
    
    # emergency attendance for covid, as per discharge diagnosis, resulting in discharge to hospital
    covidemergencyhosp_date=patients.attended_emergency_care(
      returning="date_arrived",
      date_format="YYYY-MM-DD",
      on_or_after=index_date,
      find_first_match_in_period=True,
      with_these_diagnoses = codelists.covid_emergency,
      discharged_to = codelists.discharged_to_hospital,
    ),
    
    # emergency attendance for respiratory illness
    # FIXME -- need to define codelist
    # respemergency_date=patients.attended_emergency_care(
    #   returning="date_arrived",
    #   date_format="YYYY-MM-DD",
    #   on_or_after=index_date,
    #   with_these_diagnoses = codelists.resp_emergency,
    #   find_first_match_in_period=True,
    # ),
    
    # emergency attendance for respiratory illness, resulting in discharge to hospital
    # FIXME -- need to define codelist
    # respemergencyhosp_date=patients.attended_emergency_care(
    #   returning="date_arrived",
    #   date_format="YYYY-MM-DD",
    #   on_or_after=index_date,
    #   find_first_match_in_period=True,
    #   with_these_diagnoses = codelists.resp_emergency,
    #   discharged_to = codelists.discharged_to_hospital,
    # ),
    
    # any emergency attendance
    emergency_date=patients.attended_emergency_care(
      returning="date_arrived",
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
      find_first_match_in_period=True,
    ),
    
    # emergency attendance resulting in discharge to hospital
    emergencyhosp_date=patients.attended_emergency_care(
      returning="date_arrived",
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
      find_last_match_in_period=True,
      discharged_to = codelists.discharged_to_hospital,
    ),
    
    
    # unplanned hospital admission
    admitted_unplanned_date=patients.admitted_to_hospital(
      returning="date_admitted",
      on_or_after=index_date,
      # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
      # see https://docs.opensafely.org/study-def-variables/#sus for more info
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_patient_classification = ["1"], # ordinary admissions only
      date_format="YYYY-MM-DD",
      find_first_match_in_period=True,
    ),
    
    # # planned hospital admission
    # admitted_planned_date=patients.admitted_to_hospital(
    #   returning="date_admitted",
    #   on_or_after=index_date,
    #   # see https://github.com/opensafely-core/cohort-extractor/pull/497 for codes
    #   # see https://docs.opensafely.org/study-def-variables/#sus for more info
    #   with_admission_method=["11", "12", "13", "81"],
    #   with_patient_classification = ["1"], # ordinary admissions only 
    #   date_format="YYYY-MM-DD",
    #   find_first_match_in_period=True,
    # ),
    
    # Positive covid admission prior to study start date
    covidadmitted_date=patients.admitted_to_hospital(
      returning="date_admitted",
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10,
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
      find_first_match_in_period=True,
    ),
    
    covidcritcare_date=patients.admitted_to_hospital(
      returning="date_admitted",
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10,
      with_at_least_one_day_in_critical_care=True,
      on_or_after=index_date,
      date_format="YYYY-MM-DD",
      find_first_match_in_period=True,
    ),
    
    # Covid-related death
    coviddeath_date=patients.with_these_codes_on_death_certificate(
      codelists.covid_icd10,
      returning="date_of_death",
      date_format="YYYY-MM-DD",
    ),
    
    # All-cause death
    death_date=patients.died_from_any_cause(
      returning="date_of_death",
      date_format="YYYY-MM-DD",
    ),
  )
  
  return outcome_variables