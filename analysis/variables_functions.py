from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

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
def emergency_attendance_date_X(
  name, index_date, n, with_these_diagnoses=None, discharged_to=None
):
  # emeregency attendance dates
  def var_signature(name, on_or_after, with_these_diagnoses, discharged_to):
    return {
      name: patients.attended_emergency_care(
        returning="date_arrived",
        on_or_after=on_or_after,
        find_first_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=with_these_diagnoses,
        discharged_to=discharged_to
      ),
    }
  variables = var_signature(f"{name}_1_date", index_date, with_these_diagnoses, discharged_to)
  for i in range(2, n+1):
      variables.update(var_signature(f"{name}_{i}_date", f"{name}_{i-1}_date + 1 day", with_these_diagnoses, discharged_to))
  return variables

####################################################################################################
def admitted_date_X(
  # hospital admission and discharge dates, given admission method and patient classification
  # note, it is not easy/possible to pick up sequences of contiguous episodes,
  # because we cannot reliably identify a second admission occurring on the same day as an earlier admission
  # some episodes will therefore be missed
  name, index_date, n,  
  with_these_diagnoses=None, 
  with_admission_method=None, 
  with_patient_classification=None,
):
  def var_signature(
    name, 
    on_or_after, 
    returning,
    with_these_diagnoses, 
    with_admission_method, 
    with_patient_classification
  ):
    return {
      name: patients.admitted_to_hospital(
        returning = returning,
        on_or_after = on_or_after,
        find_first_match_in_period = True,
        date_format = "YYYY-MM-DD",
        with_these_diagnoses = with_these_diagnoses,
        with_admission_method = with_admission_method,
        with_patient_classification = with_patient_classification
	   ),
    }
  variables = var_signature(
    name=f"admitted_{name}_1_date", 
    on_or_after=index_date, 
    returning="date_admitted", 
    with_these_diagnoses=with_these_diagnoses,
    with_admission_method=with_admission_method,
    with_patient_classification=with_patient_classification
  )
  variables.update(var_signature(
    name=f"discharged_{name}_1_date", 
    on_or_after=index_date, 
    returning="date_discharged", 
    with_these_diagnoses=with_these_diagnoses,
    with_admission_method=with_admission_method,
    with_patient_classification=with_patient_classification
  ))
  for i in range(2, n+1):
    variables.update(var_signature(
      name=f"admitted_{name}_{i}_date", 
      on_or_after=f"discharged_{name}_{i-1}_date + 1 day", 
      # we cannot pick up more than one admission per day
      # but "+ 1 day" is necessary to ensure we don't always pick up the same admission
      # some one day admissions will therefore be lost
      returning="date_admitted", 
      with_these_diagnoses=with_these_diagnoses,
      with_admission_method=with_admission_method,
      with_patient_classification=with_patient_classification
    ))
    variables.update(var_signature(
      name=f"discharged_{name}_{i}_date", 
      on_or_after=f"admitted_{name}_{i}_date", 
      returning="date_discharged", 
      with_these_diagnoses=with_these_diagnoses,
      with_admission_method=with_admission_method,
      with_patient_classification=with_patient_classification
    ))
  return variables

####################################################################################################
def admitted_daysincritcare_X(
  # days in critical care for a given admission episode
  name, index_name, index_date, n,  
  with_these_diagnoses=None, 
  with_admission_method=None, 
  with_patient_classification=None
):
  def var_signature(
    name, on_or_after,  
    with_these_diagnoses, 
    with_admission_method, 
    with_patient_classification
  ):
    return {
      name: patients.admitted_to_hospital(
        returning = "days_in_critical_care",
        on_or_after = on_or_after,
        find_first_match_in_period = True,
        date_format = "YYYY-MM-DD",
        with_these_diagnoses = with_these_diagnoses,
        with_admission_method = with_admission_method,
        with_patient_classification = with_patient_classification,
        return_expectations={
        "category": {"ratios": {"0": 0.75, "1": 0.20,  "2": 0.05}},
        "incidence": 0.5,
      },
	   )
    }
  variables = var_signature(
    f"admitted_{name}_ccdays_1", 
    f"admitted_{index_name}_1_date", 
    with_these_diagnoses,
    with_admission_method,
    with_patient_classification
  )
  for i in range(2, n+1):
    variables.update(var_signature(
      f"admitted_{name}_ccdays_{i}", 
      f"admitted_{index_name}_{i}_date", 
      with_these_diagnoses,
      with_admission_method,
      with_patient_classification
    ))
  return variables