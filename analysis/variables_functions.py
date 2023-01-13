from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

####################################################################################################
## function to add days to a string date
# from datetime import datetime, timedelta
# def days(datestring, days):
  
#   dt = datetime.strptime(datestring, "%Y-%m-%d").date()
#   dt_add = dt + timedelta(days)
#   datestring_add = datetime.strftime(dt_add, "%Y-%m-%d")

#   return datestring_add

####################################################################################################
def vaccination_date_X(name, index_date, n, product_name_matches=None, target_disease_matches=None):
  # vaccination date, given product_name
  def var_signature(
    name,
    on_or_after,
    product_name_matches,
    target_disease_matches
  ):
    return {
      name: patients.with_tpp_vaccination_record(
        product_name_matches=product_name_matches,
        target_disease_matches=target_disease_matches,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
  variables = var_signature(f"{name}_1_date", index_date, product_name_matches, target_disease_matches)
  for i in range(2, n+1):
    variables.update(var_signature(
      f"{name}_{i}_date", 
      f"{name}_{i-1}_date + 1 days",
      # pick up subsequent vaccines occurring one day or later -- people with unrealistic dosing intervals are later excluded
      product_name_matches,
      target_disease_matches
    ))
  return variables

####################################################################################################
def covid_test_date_X(name, index_date, n, test_result):
  # covid test date (result can be "any", "positive", or "negative")
  def var_signature(name, on_or_after, test_result):
    return {
      name: patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result=test_result,
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
        returning="date",
        date_format="YYYY-MM-DD"
      ),
    }
  variables = var_signature(f"{name}_1_date", index_date, test_result)
  for i in range(2, n+1):
    variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", test_result))
  return variables

####################################################################################################
# define recurrent event variables
def with_these_clinical_events_date_X(name, codelist, index_date, n):
    
    def var_signature(name, on_or_after, codelist):
        return {
            name: patients.with_these_clinical_events(
                    codelist,
                    returning="date",
                    on_or_after=on_or_after,
                    date_format="YYYY-MM-DD",
                    find_first_match_in_period=True
	    ),
        }
    variables = var_signature(f"{name}_1_date", index_date, codelist)
    for i in range(2, n+1):
        variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", codelist))
    return variables

####################################################################################################
def admitted_date_X(
    name, index_name, index_date, n, returning, 
    with_these_diagnoses=None, 
    with_admission_method=None, 
    with_patient_classification=None, 
    return_expectations=None
):
    def var_signature(
        name, on_or_after, returning, 
        with_these_diagnoses, 
        with_admission_method, 
        with_patient_classification, 
        return_expectations
    ):
        return {
            name: patients.admitted_to_hospital(
                    returning = returning,
                    on_or_after = on_or_after,
                    find_first_match_in_period = True,
                    date_format = "YYYY-MM-DD",
                    with_these_diagnoses = with_these_diagnoses,
                    with_admission_method = with_admission_method,
                    with_patient_classification = with_patient_classification,
                    return_expectations = return_expectations
	        ),
        }
    variables = var_signature(
        f"{name}_1_date", 
        index_date, 
        returning, 
        with_these_diagnoses,
        with_admission_method,
        with_patient_classification,
        return_expectations
    )
    for i in range(2, n+1):
        variables.update(var_signature(
            f"{name}_{i}_date", f"{index_name}_{i-1}_date + 1 day", returning, 
            with_these_diagnoses,
            with_admission_method,
            with_patient_classification,
            return_expectations
        ))
    return variables

####################################################################################################
def emergency_attendance_date_X(
  name, 
  index_date,
  n, 
  with_these_diagnoses,
  return_expectations=None
  ):
    def var_signature(
      name, 
      on_or_after, 
      ):
        return {
            name: patients.attended_emergency_care(
                    returning="date_arrived",
                    on_or_after=on_or_after,
                    with_these_diagnoses = with_these_diagnoses,
                    find_first_match_in_period=True,
                    date_format="YYYY-MM-DD",
                    return_expectations=return_expectations
	        ),
        }
    variables = var_signature(f"{name}_1_date", index_date)
    for i in range(2, n+1):
        variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day"))
    return variables