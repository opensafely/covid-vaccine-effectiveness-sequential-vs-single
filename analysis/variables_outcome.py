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
  
    # positive case identification 
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

    # positive covid admission
    covidadmitted_date=patients.admitted_to_hospital(
      returning="date_admitted",
      with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
      with_these_diagnoses=codelists.covid_icd10,
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