### description for each procedure


# feature_bucketing_addressing_data.py: it is used to extract features related to SOFA
#     input:  lab_new,  respiratoryCharting_new, nurseCharting_new  intakeOutput_new
#     output: "feature_bucket_with_missing.json"

# feature_bucketing_fill_missing.py: it is used to fill the missing values based on the method of ffill and backfill
#     input:  "feature_bucket_with_missing.json"
#     output: "feature_bucket_filled_missing.json"

# Cardiovascular_score_drug.py: generate each drugs start date, end date, and rate
#     input:  "medication_rate.csv", 'patient.csv'
#     output: "cardiovascular_score_drug.json"

# feature_bucketing_scoring.py: it is used to compute SOFA score based on Sepsis-3 criteria
#     input:  feature_bucket_filled_missing, cardiovascular_score_drug, ventilation_duration
#     output: "feature_bucket_filled_missing_scoring.json"

# sepsis_onset.py:  it is used to compute the onset of sepsis.
#     input:  'feature_bucket_filled_missing_scoring.json';  'list_infection_time.json'
#     output: sepsis_onset_df.csv --->  with two column: ['patientunitstayid','sepsis_onset']