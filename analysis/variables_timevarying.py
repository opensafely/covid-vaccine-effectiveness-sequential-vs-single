
from cohortextractor import patients, combine_codelists
from codelists import *
import json
import codelists

############################################################
## functions
from variables_functions import *
############################################################

def generate_timevarying_variables(index_date, n):

  timevarying_variables = dict(

  # for deriving suspected covid
  ## positive test
   positive_test_0_date=patients.with_test_result_in_sgss(
        pathogen="SARS-CoV-2",
        test_result="positive",
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"{index_date} - 1 day",
        find_last_match_in_period=True,
        restrict_to_earliest_specimen_date=False,
    ),
   **covid_test_date_X(
    name="positive_test", 
    index_date=index_date,
     n=n, 
     test_result="positive"
    ),

    ## suspected covid
    primary_care_suspected_covid_0_date=patients.with_these_clinical_events(
        primary_care_suspected_covid_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"{index_date} - 1 day",
        find_last_match_in_period=True,
    ),
    **with_these_clinical_events_date_X(
        name = "primary_care_suspected_covid",
        n = n,
        index_date = index_date,
        codelist = primary_care_suspected_covid_combined,
    ),

    ## probable covid
    primary_care_probable_covid_0_date=patients.with_these_clinical_events(
        covid_primary_care_probable_combined,
        returning="date",
        date_format="YYYY-MM-DD",
        on_or_before=f"{index_date} - 1 day",
        find_last_match_in_period=True,
    ),
    **with_these_clinical_events_date_X(
        name = "primary_care_probable_covid",
        n = n,
        index_date = index_date,
        codelist = covid_primary_care_probable_combined,
    ),

    # for deriving unplaned hospital admissions
    ## upplanned hospital admission
    admitted_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before=f"{index_date} - 1 day",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": "2020-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    **admitted_date_X(
        name = "admitted_unplanned",
        n = n,
        index_name = "admitted_unplanned",
        index_date = "index_date",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        return_expectations={
            "date": {"earliest": "2021-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    discharged_unplanned_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before=f"{index_date} - 1 day",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": "2020-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ), 
    **admitted_date_X(
        name = "discharged_unplanned",
        n = n,
        index_name = "admitted_unplanned",
        index_date = "index_date",
        returning="date_discharged",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        return_expectations={
            "date": {"earliest": "2021-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),


    ## unplanned infectious hospital admission
    admitted_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_admitted",
        on_or_before=f"{index_date} - 1 day",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": "2020-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    **admitted_date_X(
        name = "admitted_unplanned_infectious",
        n = n,
        index_name = "admitted_unplanned_infectious",
        index_date = "index_date",
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        return_expectations={
            "date": {"earliest": "2021-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    discharged_unplanned_infectious_0_date=patients.admitted_to_hospital(
        returning="date_discharged",
        on_or_before=f"{index_date} - 1 day",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        date_format="YYYY-MM-DD",
        find_first_match_in_period=True,
        return_expectations={
            "date": {"earliest": "2020-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ), 
    **admitted_date_X(
        name = "discharged_unplanned_infectious",
        n = n,
        index_name = "admitted_unplanned_infectious",
        index_date = "index_date",
        returning="date_discharged",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        with_these_diagnoses = ICD10_I_codes,
        return_expectations={
            "date": {"earliest": "2021-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),

  )

  return timevarying_variables

