# extract infection itme for all icuStayID patients  from eICU cohort

import json
import pickle
import pandas as pd


def obtain_infection(patient_icd_code_list,each_infe_codes, patient_times):
    # each_infe_codes is list
    infe_0_or_1 = 0

    for index, i in enumerate(patient_icd_code_list):
        if (i in each_infe_codes) and (patient_times[index] < 1440):
            infe_0_or_1 = 1

    return infe_0_or_1


path_paitent_ID = '/Users/suraj/Documents/CBM/Fei/Trial/Data/'
path_diagnose_ICD9_code = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eicu-database-2.0/'
result_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/sepsis_code/data/'

object_patient_ID = pd.read_csv(path_paitent_ID+'sepsis.csv',index_col=False)
patientunitstayid_s = list(object_patient_ID['patientunitstayid'])

DIAGNOSES_ICD = pd.read_csv(path_diagnose_ICD9_code+'diagnosis.csv',index_col=False)
print('before_len(DIAGNOSES_ICD)',len(DIAGNOSES_ICD))
DIAGNOSES_ICD = DIAGNOSES_ICD[DIAGNOSES_ICD['icd9code'].notna()]
DIAGNOSES_ICD = DIAGNOSES_ICD[['patientunitstayid','icd9code', 'diagnosisoffset']]

# DIAGNOSES_ICD.to_csv('/Users/xuzhenxing/Documents/eICU_AKI_Sepsis/testing/DIAGNOSES_ICD_updated.csv',index=False)
# print('after_len(DIAGNOSES_ICD)',len(DIAGNOSES_ICD))

# Get the indexes which are repetative with the split
df = DIAGNOSES_ICD
df['icd9code'] = df['icd9code'].str.split(',')
df = df.explode('icd9code')
# remove blank in each string
df['icd9code'] = df['icd9code'].str.strip()
DIAGNOSES_ICD = df



infection_name = ['cfi_cns','cfi_intabd','cfi_pneumonia','cfi_septicemia_bacteremia','cfi_skin_softtis','cfi_uti']


infection_ICD9_code = {
    'cfi_cns': ["006.5",
                "013.00","013.01","013.02","013.03","013.04","013.05","013.06",
                "013.10","013.11","013.12","013.13","013.14","013.15","013.16",
                "013.20","013.21","013.22","013.23","013.24","013.25","013.26",
                "013.30","013.31","013.32","013.33","013.34","013.35","013.36","013.40",
                "013.41","013.42","013.43","013.44","013.45","013.46","013.50","013.51",
                "013.52","013.53","013.54","013.55","013.56","013.60","013.61","013.62",
                "013.63","013.64","013.65","013.66","013.80","013.81","013.82","013.83",
                "013.84","013.85","013.86","013.90","013.91","013.92","013.93","013.94",
                "013.95","013.96","036.0","036.1","049.8","052.0","052.2","053.0",
                "053.14","054.3","054.72","055.0","056.01","064","072.1","072.2",
                "091.81","094.2","098.82","100.81","112.83","114.2","130.0","320.0",
                "320.1","320.2","320.3","320.7","320.81","320.82","320.89","320.9",
                "321.0","321.1","321.2","321.3","321.8","323.01","323.02","323.1",
                "323.2","323.41","323.42","324.0","324.1","324.9"],

    'cfi_intabd': ["008.00","008.01","008.02","008.03","008.04","008.09","008.1",
                  "008.2","008.3","008.41","008.42","008.43","008.44",
                  "008.45","008.46","008.47","008.49","008.5","009.0",
                  "009.1","009.2","009.3","540.0","540.1","540.9",
                  "541","542","543.9","562.01","562.03","562.11","562.13",
                  "567.0","567.1","567.21","567.22","567.23","567.29","567.31",
                  "567.38","567.39","567.81","567.82","567.89","567.9",
                  "569.5","569.61","569.71","569.83","572.0","572.1",
                  "574.00","574.01","574.10","574.11","574.30","574.31","574.40",
                  "574.41","574.60","574.61","574.70","574.71","574.80",
                  "574.81","575.0","575.10","575.12","575.4","576.1","576.2",
                  "576.3","576.8","576.9","614.0","614.2","614.3",
                  "614.5","614.6","614.8","614.9"],

    'cfi_pneumonia': ["003.22","021.2","039.1","052.1","055.1",
                    "112.4","114.0","115.05","130.4","136.3","480.0",
                    "480.1","480.2","480.3","480.8","480.9","481",
                    "482.0","482.1","482.2","482.30","482.31","482.32",
                    "482.39","482.40","482.41","482.42","482.49","482.81",
                    "482.82","482.83","482.84","482.89","482.9","483.0",
                    "483.1","483.8","484.1","484.3","484.5","484.6",
                    "484.7","484.8","485","486"],

    'cfi_septicemia_bacteremia': ["003.1","020.2","022.3","036.2","038.0","038.10",
                              "038.11","038.12","038.19","038.2","038.3","038.40",
                              "038.41","038.42","038.43","038.44","038.49","038.8",
                              "038.9","054.5","112.5","785.52","790.7","995.91",
                              "995.92"],

    'cfi_skin_softtis':["035","376.01","680.0","680.1","680.2","680.3","680.4",
                          "680.5","680.6","680.7","680.8","680.9","681.00",
                          "681.01","681.02","681.10","681.11","681.9","682.0",
                          "682.1","682.2","682.3","682.4","682.5",
                          "682.6","682.7","682.8","682.9","683","684","685.0",
                          "686.00","686.01","686.09","686.1","686.8","686.9",
                          "728.86"],

    'cfi_uti':["590.10","590.11","590.2","590.3","590.80","590.81",
              "590.9","595.0","595.3","595.4","595.89",
              "595.9","597.0","597.80","597.89","598.00","598.01",
              "599.0"]

}


all_infection_icd9_code_df = pd.DataFrame(columns=['patientunitstayid']+infection_name)

j = 0
# #testing
patient_icd_code_list = DIAGNOSES_ICD.loc[DIAGNOSES_ICD['patientunitstayid'] == 1576563]
patient_icd_code_list.to_csv(result_path+'1576563.csv', index=False)
# print('patient_icd_code_list',patient_icd_code_list)


for p_id in patientunitstayid_s:
    print('j',j)

    person_infection_icd9_code_df = pd.DataFrame(columns=['patientunitstayid']+infection_name)
    person_result_dic = {}

    # print('p_id',p_id)
    person_result_dic['patientunitstayid'] = p_id

    patient_icd_code_list = list(DIAGNOSES_ICD.loc[DIAGNOSES_ICD['patientunitstayid'] == p_id, 'icd9code'])

    patient_times = list(DIAGNOSES_ICD.loc[DIAGNOSES_ICD['patientunitstayid'] == p_id, 'diagnosisoffset'])

    # print('patient_icd_code_list',patient_icd_code_list)

    for infe in infection_ICD9_code:
        infe_name = infe
        each_infe_codes = infection_ICD9_code[infe_name]
        # print('infe_name', infe_name)
        # print('each_infe_codes', each_infe_codes)
        infection_result = obtain_infection(patient_icd_code_list,each_infe_codes, patient_times)

        person_result_dic[infe_name] = infection_result

    person_infection_icd9_code_df = person_infection_icd9_code_df.append(person_result_dic, ignore_index=True, sort=False)

    all_infection_icd9_code_df = all_infection_icd9_code_df.append(person_infection_icd9_code_df)

    j = j + 1
    # if j >1000:
    #     break

all_infection_icd9_code_df.to_csv(result_path+'patient_infection_item_df.csv',index=False)

print('ss')







