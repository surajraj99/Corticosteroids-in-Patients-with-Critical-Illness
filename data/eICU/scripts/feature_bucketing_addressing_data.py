#  to extract the features original information based on lab_new,  respiratoryCharting_new, nurseCharting_new  intakeOutput_new

import pandas as pd
import numpy as  np
import json
import math
import time

file_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/eicu-database-2.0/'
result_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/sepsis_code/data/'

class JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj,'to_dict'):
            return obj.to_dict()
        return json.JSONEncoder.default(self, obj)

def extracting_fea(admit_f,discharge_f,feature_name_f, feature_item_f,interval_f,p_ID_f,lab_f,featureDirection_f):
    feature_name_f = feature_name_f
    num_time_slice = int(math.ceil(discharge_f/(interval_f*1.0)))
    p_ori_infor = lab_f[lab_f['patientunitstayid']==p_ID_f]

    # p_ori_infor.to_csv(result_path+'testing/'+str(p_ID_f)+'ori_feature_information.csv',index=False)
    p_result = []
    try:
        for i in range(num_time_slice):
            # print('slice_i',i)
            if i==0:
                start = -1440 # lab values is also used before 24h,(1440 min) admitted in ICU
            else:
                start = i * interval_f

            end = (i+1) * interval_f
            if start <= discharge_f:
                fea_ori = list(p_ori_infor.loc[(p_ori_infor['feature_item']==feature_item_f)&(p_ori_infor['feature_offset']>=start)&(p_ori_infor['feature_offset']<=end),'feature_result'])

                # print('fea',feature_item_f,fea_ori)
                # print('len()',len(fea_ori))

                if len(fea_ori)>0:
                    fea_ori = list(map(float, fea_ori))
                    if featureDirection_f=='min':
                        if feature_name_f =='Urine':
                            max_min_fea = sum(fea_ori) # for urine, we need to sum
                        else:
                            max_min_fea = min(fea_ori)
                    else:
                        max_min_fea = max(fea_ori)

                    fea_final = max_min_fea
                else:
                    fea_final = np.nan # empty
                # print('fea_final',fea_final)

                p_result.append(fea_final)

            else:
                break
        # print('p_result',p_result)

    except:
        print('there is an exception.')
        # generate empty for abnormal samples
        p_result = [np.nan] * num_time_slice

    return p_result


start_time = time.time()
patient = pd.read_csv(file_path+'patient.csv',index_col=False)
patient_ICU_ID = list(set(patient['patientunitstayid']))
print('len(patient_ICU_ID)',len(patient_ICU_ID))

lab_feature = ['PaO2','Bilirubin','Platelet','Creatinine']
lab_feature_source = ['paO2','total bilirubin','platelets x 1000','creatinine']

respiratoryCharting_feature = ['FiO2']
respiratoryCharting_feature_source = ['FIO2 (%)','FiO2']

nurseCharting_feature = ['MAP','GCS']
nurseCharting_feature_source = ['MAP (mmHg)','Arterial Line MAP (mmHg)','Score (Glasgow Coma Scale)','Glasgow coma score'] # in nurseCharting.csv,nurseCharting_feature_source in column "nursingchartcelltypevallabel"

intakeOutput_feature = ['Urine']
intakeOutput_feature_source = [
'Urine, L neph:',
'Urine Output-Urine Output',
'Urine Output-RIGHT Nephrouretero Stent Urine Output',
'Urine Output-Nephrostomy',
'Urine Output-Foley',
'Urine Output-foley',
'Urine Output-FOLEY',
'Urine Output (mL)-Urethral Catheter',
'Urine Output (mL)-Urethral Catheter',
'Urine Output (mL)-Urethral Catheter',
'URINE CATHETER',
'Urine',
'Urinary Catheter Output: Urethral 16 Fr.',
'Urinary Catheter Output: Urethral',
'Urinary Catheter Output: Suprapubic',
'Urinary Catheter Output: Nephrostomy R',
'Urinary Catheter Output: Nephrostomy L',
'Urinary Catheter Output: Coude Urethral 18 Fr.',
'Urinary Catheter Output: Coude Urethral 16 Fr.',
'Urinary Catheter Output: Coude Urethral 14 Fr.',
'Urinary Catheter Output: Coude Urethral',
'Urinary Catheter Output: Condom Urethral',
'Urinary Catheter Output: Condom External'
]

# read lab_new,  respiratoryCharting_new, nurseCharting_MAP_GCS_new  intakeOutput_new

lab = pd.read_csv(result_path+'lab_new.csv',index_col=False)
# lab = lab.loc[~lab['patientunitstayid'].isin(noise_patient)]

respiratoryCharting = pd.read_csv(result_path+'respiratoryCharting_new.csv',index_col=False)
respiratoryCharting['feature_result']=respiratoryCharting['feature_result'].str.split('%').str[0] # addressing
respiratoryCharting['feature_result'] = pd.to_numeric(respiratoryCharting['feature_result'], errors='coerce')
respiratoryCharting = respiratoryCharting.loc[(respiratoryCharting['feature_result']>0)&(respiratoryCharting['feature_result']<=100)]

# respiratoryCharting.to_csv(result_path+'ssss.csv',index=False)

nurseCharting = pd.read_csv(result_path+'nurseCharting_new.csv',index_col=False)
# address one abnormal "Unable to score due to medication" for GCS  noise_patient = [524364,1048660]
nurseCharting = nurseCharting[nurseCharting['feature_result']!='Unable to score due to medication']

intakeOutput = pd.read_csv(result_path+'intakeOutput_new.csv',index_col=False)
# intakeOutput = intakeOutput.loc[~intakeOutput['patientunitstayid'].isin(noise_patient)]

# read feature list
# feature_info_df = pd.read_csv(result_path+ 'variable_list.csv', sep=',', header=0)
feature_info_df = pd.read_excel(result_path + 'variable_list.xlsx', index_col=False)

feature_bucket = {}
interval = 6*60 # interval 24 h, 1440 mins
count = 0

for i in patient_ICU_ID:
    print('count=',count)
    p_ID = i
    # print('p_ID',p_ID)
    admit = 0 # admit offset is 0 for all patient
    discharge = list(patient.loc[patient['patientunitstayid']==p_ID,'unitdischargeoffset'])[0]

    # print('discharge',discharge)

    time_slice = int(math.ceil(discharge/(interval*1.0)))

    # create column names for DF, which is used to save the all features result for one person
    columns_name = []
    for j in range(time_slice):
        col_name = 'time_' + str(j)
        columns_name.append(col_name)
    columns_name = ['FeatureName'] + columns_name
    # p_feature_bucket_DF is used to save the all features result for one person
    p_feature_bucket_DF = pd.DataFrame(columns=columns_name)

    # print('p_feature_bucket_DF', p_feature_bucket_DF)
    # print('p_ID,discharge',p_ID,discharge,time_slice)

    for idx, row in feature_info_df.iterrows():
        featureName, feature_item, featureDirection = row['FeatureName'], row['name_in_source'],  row['FeatureDirection']
        # print('featureName:',featureName)

        flag_sofa_feature = 0
        if featureName in lab_feature:
            feature_source = lab
            flag_sofa_feature = 1
        elif featureName in respiratoryCharting_feature:
            feature_source = respiratoryCharting
            flag_sofa_feature = 1
        elif featureName in nurseCharting_feature:
            feature_source = nurseCharting
            flag_sofa_feature = 1
        elif featureName in intakeOutput_feature:
            feature_source = intakeOutput
            flag_sofa_feature = 1

        if flag_sofa_feature == 1:
            fea_value = extracting_fea(admit,discharge,featureName,feature_item,interval,p_ID,feature_source,featureDirection)
            p_feature_bucket_DF.loc[len(p_feature_bucket_DF)] = [featureName] + fea_value

    # print('p_feature_bucket_DF',p_feature_bucket_DF)
    feature_bucket[str(p_ID)] = p_feature_bucket_DF

    # if count>1000:
    #     break

    count = count + 1

# save dic file
with open(result_path + "feature_bucket_with_missing.json", "w") as outfile:
    json.dump(feature_bucket, outfile, cls=JSONEncoder, indent=4)

print('done.')
print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

















