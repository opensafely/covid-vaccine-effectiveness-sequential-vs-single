from cohortextractor import (codelist, codelist_from_csv, combine_codelists)


covid_icd10 = codelist_from_csv(
    "codelists/opensafely-covid-identification.csv",
    system="icd10",
    column="icd10_code",
)

covid_emergency = codelist(
    ["1240751000000100"],
    system="snomed",
)


covid_primary_care_positive_test = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-positive-test.csv",
    system="ctv3",
    column="CTV3ID",
)

covid_primary_care_code = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-clinical-code.csv",
    system="ctv3",
    column="CTV3ID",
)

covid_primary_care_sequelae = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-probable-covid-sequelae.csv",
    system="ctv3",
    column="CTV3ID",
)

covid_primary_care_probable_combined = combine_codelists(
    covid_primary_care_positive_test,
    covid_primary_care_code,
    covid_primary_care_sequelae,
)
covid_primary_care_suspected_covid_advice = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-advice.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_had_test = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-had-test.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_isolation = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-isolation-code.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_nonspecific_clinical_assessment = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-suspected-covid-nonspecific-clinical-assessment.csv",
    system="ctv3",
    column="CTV3ID",
)
covid_primary_care_suspected_covid_exposure = codelist_from_csv(
    "codelists/opensafely-covid-identification-in-primary-care-exposure-to-disease.csv",
    system="ctv3",
    column="CTV3ID",
)
primary_care_suspected_covid_combined = combine_codelists(
    covid_primary_care_suspected_covid_advice,
    covid_primary_care_suspected_covid_had_test,
    covid_primary_care_suspected_covid_isolation,
    covid_primary_care_suspected_covid_exposure,
)



ethnicity = codelist_from_csv(
    "codelists/opensafely-ethnicity.csv",
    system="ctv3",
    column="Code",
    category_column="Grouping_6",
)
ethnicity_16 = codelist_from_csv(
    "codelists/opensafely-ethnicity.csv",
    system="ctv3",
    column="Code",
    category_column="Grouping_16",
)

## PRIMIS
# Patients in long-stay nursing and residential care
carehome = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-longres.csv", 
    system="snomed", 
    column="code",
)


# High Risk from COVID-19 code
shield = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-shield.csv",
    system="snomed",
    column="code",
)

# Lower Risk from COVID-19 codes
nonshield = codelist_from_csv(
    "codelists/primis-covid19-vacc-uptake-nonshield.csv",
    system="snomed",
    column="code",
)



# Asthma Diagnosis code
ast = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ast.csv",
  system="snomed",
  column="code",
)

# Asthma Admission codes
astadm = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-astadm.csv",
  system="snomed",
  column="code",
)

# Asthma systemic steroid prescription codes
astrx = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-astrx.csv",
  system="snomed",
  column="code",
)

# Chronic Respiratory Disease
resp_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-resp_cov.csv",
  system="snomed",
  column="code",
)

# Chronic heart disease codes
chd_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-chd_cov.csv",
  system="snomed",
  column="code",
)

# Chronic kidney disease diagnostic codes
ckd_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd_cov.csv",
  system="snomed",
  column="code",
)

# Chronic kidney disease codes - all stages
ckd15 = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd15.csv",
  system="snomed",
  column="code",
)

# Chronic kidney disease codes-stages 3 - 5
ckd35 = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-ckd35.csv",
  system="snomed",
  column="code",
)

# Chronic Liver disease codes
cld = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cld.csv",
  system="snomed",
  column="code",
)

# Diabetes diagnosis codes
diab = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-diab.csv",
  system="snomed",
  column="code",
)

# Immunosuppression diagnosis codes
immdx_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immdx_cov.csv",
  system="snomed",
  column="code",
)

# Immunosuppression medication codes
immrx = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-immrx.csv",
  system="snomed",
  column="code",
)

# Chronic Neurological Disease including Significant Learning Disorder
cns_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-cns_cov.csv",
  system="snomed",
  column="code",
)

# Asplenia or Dysfunction of the Spleen codes
spln_cov = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-spln_cov.csv",
  system="snomed",
  column="code",
)

# BMI
bmi = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-bmi.csv",
  system="snomed",
  column="code",
)

# All BMI coded terms
bmi_stage = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-bmi_stage.csv",
  system="snomed",
  column="code",
)

# Severe Obesity code recorded
sev_obesity = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-sev_obesity.csv",
  system="snomed",
  column="code",
)

# Diabetes resolved codes
dmres = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-dmres.csv",
  system="snomed",
  column="code",
)

# Severe Mental Illness codes
sev_mental = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-sev_mental.csv",
  system="snomed",
  column="code",
)

# Remission codes relating to Severe Mental Illness
smhres = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-smhres.csv",
  system="snomed",
  column="code",
)


# to represent household contact of shielding individual
hhld_imdef = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-hhld_imdef.csv",
  system="snomed",
  column="code",
)

# Wider Learning Disability
learndis = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-learndis.csv",
  system="snomed",
  column="code",
)



# Carer codes
carer = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-carer.csv",
  system="snomed",
  column="code",
)

# No longer a carer codes
notcarer = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-notcarer.csv",
  system="snomed",
  column="code",
)

# Employed by Care Home codes
carehomeemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-carehome.csv",
  system="snomed",
  column="code",
)

# Employed by nursing home codes
nursehomeemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-nursehome.csv",
  system="snomed",
  column="code",
)

# Employed by domiciliary care provider codes
domcareemployee = codelist_from_csv(
  "codelists/primis-covid19-vacc-uptake-domcare.csv",
  system="snomed",
  column="code",
)


eol = codelist_from_csv(
    "codelists/nhsd-primary-care-domain-refsets-palcare_cod.csv",
    system="snomed",
    column="code",
)

midazolam = codelist_from_csv(
    "codelists/opensafely-midazolam-end-of-life.csv",
    system="snomed",
    column="dmd_id",   
)

housebound = codelist_from_csv(
    "codelists/opensafely-housebound.csv", 
    system="snomed", 
    column="code"
)

no_longer_housebound = codelist_from_csv(
    "codelists/opensafely-no-longer-housebound.csv", 
    system="snomed", 
    column="code"
)

discharged_to_hospital = codelist(
    ["306706006", "1066331000000109", "1066391000000105"],
    system="snomed",
)

ICD10_I_codes = codelist_from_csv(
    "codelists/opensafely-icd-10-chapter-i.csv",
    system="icd10",
    column="code",
)