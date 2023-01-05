
from cohortextractor import patients, combine_codelists
from codelists import *
import codelists

############################################################
## functions
from variables_functions import *
############################################################

def generate_timevarying_variables(index_date, n):

  timevarying_variables = dict(

  # for deriving suspected covid
  ## positive test
   # positive_test_0_date defined in variables_pre
#    **covid_test_date_X(
#     name="positive_test", 
#     index_date=index_date,
#      n=n, 
#      test_result="positive"
#     ),

    # for deriving unplaned hospital admissions
    ## upplanned hospital admission
    # admitted_unplanned_0_date defined in variables_pre
    **admitted_date_X(
        name = "admitted_unplanned",
        n = n,
        index_name = "admitted_unplanned",
        index_date = index_date,
        returning="date_admitted",
        with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"],
        with_patient_classification = ["1"],
        return_expectations={
            "date": {"earliest": "2021-05-01", "latest" : "2021-06-01"},
            "rate": "uniform",
            "incidence": 0.05,
        },
    ),
    # discharged_unplanned_0_date defined in variables_pre
    **admitted_date_X(
        name = "discharged_unplanned",
        n = n,
        index_name = "admitted_unplanned",
        index_date = index_date,
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
        index_date = index_date,
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
        index_date = index_date,
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

    # emergency attendance
    # covidemergency_0_date defined in variables_pre
    # **emergency_attendance_date_X(
    #     name = "covidemergency",
    #     index_date = index_date,
    #     n = n,
    #     with_these_diagnoses = codelists.covid_emergency,
    # ),

  )

  return timevarying_variables

