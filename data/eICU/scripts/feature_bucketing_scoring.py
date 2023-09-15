# scoring based on the result from feature_bucket_filled_missing, cardiovascular_score_drug, ventilation_duration

import pandas as pd
import numpy as np
import json
import time

file_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/eicu-database-2.0/'
result_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/sepsis_code/data/'

class JSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if hasattr(obj,'to_dict'):
            return obj.to_dict()
        return json.JSONEncoder.default(self, obj)

def add_vent_flag(feature_df,p_ventilation,time_interval):
    interval = time_interval # 24 hour, 1440 mins
    df = feature_df
    p_vent = p_ventilation

    if len(p_vent)==0:
        # print('p_vent',p_vent)
        df = df
    else:
        vent_start = p_vent['start']
        vent_end = p_vent['end']

        # print('len(vent_end)', vent_end)
        # print('len(vent_start)', vent_start)
        for i in range(len(df)):
            i_start = i*interval
            i_end = (i+1)*interval
            # print('000_i_start', i_start)
            # print('000_i_end', i_end)

            for j in range(len(vent_start)):
                j_start = vent_start[j]
                j_end = vent_end[j]

                if (i_start>j_end) or (i_end<j_start):

                    vent_flag = 0
                    # print('there is not vent')
                    # print('000_j_start', j_start)
                    # print('000_j_end', j_end)

                else:
                    vent_flag = 1
                    # assign 1 for i th row, means the patient use vent during i th time duration ;  index is " time_"
                    # print('111_this is vent')
                    # print('111_i_start', i_start)
                    # print('111_i_end', i_end)
                    # print('111_j_start', j_start)
                    # print('111_j_end', j_end)

                    df.loc['time_'+str(i),'vent'] = 1

    return df

def compute_MAP_score(MAP_value):
    score = 0
    map_value = MAP_value

    if (map_value>=70) or (map_value==0): # 0 is used to consider when filling missing values
        score = 0
    else:
        score = 1

    return score

def compute_Dopamine_score(rate):
    score = 0
    if rate <=5:
        score=2
    if (rate>5) and (rate<=15):
        score=3
    if rate >15:
        score=4

    return score

def compute_Dobutamine_score(rate):
    score = 0
    if rate >0:
        score = 2
    return score

def compute_Epinephrine_score(rate):
    score = 0
    if rate <=0.1:
        score=3
    else:
        score = 4

    return score

def compute_Norepinephrine_score(rate):
    score = 0
    if rate <=0.1:
        score=3
    else:
        score = 4

    return score

def compute_Phenylephrine_score(rate):
    score =0
    if (rate>0) and (rate<=200):
        score = 2
    if rate>200:
        score =3
    return score

def compute_Vasopressin_score(rate):
    score =0
    if rate>0:
        score =3
    return score


def compute_drug_score(drug_name, drug_value, i_start, i_end):
    score = 0
    drug_name = drug_name
    i_start = i_start
    i_end = i_end
    drug_value = drug_value
    drug_start = drug_value['drug_start']
    drug_end = drug_value['drug_end']
    drug_rate = drug_value['rate']

    if len(drug_start)==0:
        score = 0
        # print('drug_rate',len(drug_start))
    else:
        for j in range(len(drug_start)):
            j_start = drug_start[j]
            j_end = drug_end[j]
            if (i_start>j_end) or (i_end <j_start):
                score = 0
                # print('there in no overlap two time slices')
            else:

                rate = drug_rate[j]
                if drug_name=='Dopamine':
                    score = compute_Dopamine_score(rate)
                if drug_name=='Dobutamine':
                    score = compute_Dobutamine_score(rate)
                if drug_name=='Epinephrine':
                    score = compute_Epinephrine_score(rate)
                if drug_name=='Norepinephrine':
                    score = compute_Norepinephrine_score(rate)
                if drug_name=='Phenylephrine':
                    score = compute_Phenylephrine_score(rate)

    return score

def compute_drug_Vasopressin_score(drug_name, drug_value, i_start, i_end):
    score = 0
    rate = 0
    drug_name = drug_name
    i_start = i_start
    i_end = i_end
    drug_value = drug_value
    drug_start = drug_value['drug_start']
    drug_end = drug_value['drug_end']
    drug_rate = drug_value['rate']

    if len(drug_start)==0:
        score = 0
        # print('dop_rate',len(drug_start))
    else:
        for j in range(len(drug_start)):
            j_start = drug_start[j]
            j_end = drug_end[j]
            if (i_start>j_end) or (i_end <j_start):
                score = 0
                # print('there in no overlap two time slices')
            else:

                rate = drug_rate[j]
                score = compute_Vasopressin_score(rate)

    return score, rate


def compute_Cardiovascular_score(feature_df,p_drug,time_interval):

    df = feature_df
    for i in range(len(df)):
        i_start = i*time_interval
        i_end = (i+1)*time_interval
# extract MAP value and compute MAP score
        MAP_value = df.loc['time_'+str(i),'MAP']
        MAP_score = compute_MAP_score(MAP_value)
        # print('MAP_value', MAP_value)
        # print('MAP_score', MAP_score)

# compute Dopamine
        Dopamine_value = p_drug['Dopamine']
        Dopamine_name = 'Dopamine'
        Dopamine_score = compute_drug_score(Dopamine_name,Dopamine_value, i_start, i_end)
        # print('Dopamine_score',Dopamine_score)

# compute Dobutamine
        Dobutamine_value = p_drug['Dobutamine']
        Dobutamine_name = 'Dobutamine'
        Dobutamine_score = compute_drug_score(Dobutamine_name,Dobutamine_value, i_start, i_end)
        # print('Dobutamine_score', Dobutamine_score)

# compute Epinephrine
        Epinephrine_value = p_drug['Epinephrine']
        Epinephrine_name = 'Epinephrine'
        Epinephrine_score = compute_drug_score(Epinephrine_name, Epinephrine_value, i_start, i_end)
        # print('Epinephrine_score', Epinephrine_score)

# compute Norepinephrine
        Norepinephrine_value = p_drug['Norepinephrine']
        Norepinephrine_name = 'Norepinephrine'
        Norepinephrine_score = compute_drug_score(Norepinephrine_name, Norepinephrine_value, i_start, i_end)
        # print('Norepinephrine_score', Norepinephrine_score)

# compute Phenylephrine
        Phenylephrine_value = p_drug['Phenylephrine']
        Phenylephrine_name = 'Phenylephrine'
        Phenylephrine_score = compute_drug_score(Phenylephrine_name, Phenylephrine_value, i_start, i_end)
        # print('Phenylephrine_score', Phenylephrine_score)

# compute Vasopressin  Note that, we need obtain  Vasopressin_score, Vasopressin_rate
        Vasopressin_value = p_drug['Vasopressin']
        Vasopressin_name = 'Vasopressin'
        Vasopressin_score, Vasopressin_rate = compute_drug_Vasopressin_score(Vasopressin_name, Vasopressin_value, i_start, i_end)
        # print('Vasopressin_score', Vasopressin_score)
        # print('Vasopressin_rate', Vasopressin_rate)

# summary score based on MAP and drugs
        Cardiovascular_score = max([MAP_score,Dopamine_score,Dobutamine_score,Epinephrine_score,Norepinephrine_score,Phenylephrine_score,Vasopressin_score])

        if (Dopamine_score==3) or (Epinephrine_score==3) or (Norepinephrine_score==3):
            if (Phenylephrine_score==3) or Vasopressin_rate>0.4:
                Cardiovascular_score = 4

        df.loc['time_' + str(i), 'Cardiovascular_score'] = Cardiovascular_score

    return df

start_time = time.time()
# read feature_bucket_with_filled_missing
with open(result_path + 'feature_bucket_filled_missing.json', 'r') as load_f:
    feature_bucket_filled_missing = json.load(load_f)
print('len(feature_bucket_filled_missing)',len(feature_bucket_filled_missing))

# read ventilation duration information
with open(result_path + 'ventilation_duration.json', 'r') as load_f:
    ventilation_duration = json.load(load_f)

print('len(ventilation_duration)',len(ventilation_duration))

# read drug to obtain drug_rate
with open(result_path + 'cardiovascular_score_drug.json', 'r') as load_f:
    cardiovascular_score_drug = json.load(load_f)

print('len(cardiovascular_score_drug)',len(cardiovascular_score_drug))

count = 0
time_interval = 360 # 24h, 1440 min

feature_bucket_filled_missing_scoring = {}

for i in feature_bucket_filled_missing:
    print('count',count)

    p_ID = i
    # print('p_ID',p_ID)
    feature = feature_bucket_filled_missing[p_ID]
    feature_df = pd.DataFrame.from_dict(feature)

    # feature_df.to_csv(result_path+'testing/person_feature_'+p_ID+'_lenth.csv',index=False)

    # set FeatureName as index
    feature_df = feature_df.set_index('FeatureName')
    # df1.T  transpose   column  and row:  column->
    feature_df = feature_df.T
    # Columns: # Bilirubin  # Creatinine # FiO2 # GCS # MAP # PaO2 #   # Urine

    # compute_SOFA score
    # Coagulation score--Platelet
    feature_df['Coagulation_score'] = 0

    feature_df.loc[feature_df['Platelet']>=150,'Coagulation_score'] = 0
    feature_df.loc[(feature_df['Platelet']>=100)&(feature_df['Platelet']< 150), 'Coagulation_score'] = 1
    feature_df.loc[(feature_df['Platelet']>=50)&(feature_df['Platelet']< 100), 'Coagulation_score'] = 2
    feature_df.loc[(feature_df['Platelet']>=20)&(feature_df['Platelet']< 50), 'Coagulation_score'] = 3
    feature_df.loc[(feature_df['Platelet']>0)&(feature_df['Platelet']<=20), 'Coagulation_score'] = 4
    feature_df.loc[feature_df['Platelet']==0, 'Coagulation_score'] = 0 # missing value is filled with 0 and scoring 0

    # Liver score--Bilirubin
    feature_df['Liver_score'] = 0

    feature_df.loc[feature_df['Bilirubin']<1.2,'Liver_score'] = 0 # (including # missing value is filled with 0 and scoring 0
    feature_df.loc[(feature_df['Bilirubin']>=1.2)&(feature_df['Bilirubin']<=1.9), 'Liver_score'] = 1
    feature_df.loc[(feature_df['Bilirubin']>1.9)&(feature_df['Bilirubin']<=5.9), 'Liver_score'] = 2
    feature_df.loc[(feature_df['Bilirubin']>5.9)&(feature_df['Bilirubin']<=11.9), 'Liver_score'] = 3
    feature_df.loc[feature_df['Bilirubin'] >11.9, 'Liver_score'] = 4

    # CNS score--GCS
    feature_df['CNS_score'] = 0

    feature_df.loc[feature_df['GCS']==15,'CNS_score'] = 0
    feature_df.loc[(feature_df['GCS']>=13)&(feature_df['GCS']<=14), 'CNS_score'] = 1
    feature_df.loc[(feature_df['GCS']>=10)&(feature_df['GCS']<=12), 'CNS_score'] = 2
    feature_df.loc[(feature_df['GCS']>=6)&(feature_df['GCS']<=9 ), 'CNS_score'] = 3
    feature_df.loc[(feature_df['GCS'] >0) & (feature_df['GCS'] <=6), 'CNS_score'] = 4
    feature_df.loc[feature_df['GCS']==0, 'CNS_score'] = 0 # missing value is filled with 0 and scoring 0

    # Creatinine_score--Creatinine
    feature_df['Creatinine_score'] = 0

    feature_df.loc[feature_df['Creatinine']<1.2,'Creatinine_score'] = 0  #(including # missing value is filled with 0 and scoring 0)
    feature_df.loc[(feature_df['Creatinine']>=1.2)&(feature_df['Creatinine']<=1.9), 'Creatinine_score'] = 1
    feature_df.loc[(feature_df['Creatinine']>1.9)&(feature_df['Creatinine']<=3.4), 'Creatinine_score'] = 2
    feature_df.loc[(feature_df['Creatinine']>3.4)&(feature_df['Creatinine']<=4.9), 'Creatinine_score'] = 3
    feature_df.loc[feature_df['Creatinine'] >4.9, 'Creatinine_score'] = 4

    # Urine_score, for urine, only 0, 3, 4 score
    feature_df['Urine_score'] = 0

    feature_df.loc[feature_df['Urine']>=500,'Urine_score'] = 0
    feature_df.loc[(feature_df['Urine']>=200)&(feature_df['Urine']<500), 'Urine_score'] = 3
    feature_df.loc[(feature_df['Urine']>0)&(feature_df['Urine']<=200), 'Urine_score'] = 4
    feature_df.loc[feature_df['Urine']==0, 'Urine_score'] = 0 # missing value is filled with 0 and scoring 0

    # Renal_score based on Creatinine_score and Urine_score
    feature_df['Renal_score'] = feature_df[['Creatinine_score', 'Urine_score']].max(axis=1)

    # Respiration score
    feature_df['vent'] = 0 # add one column that is used for vent flag  1 means the patient uses, and 0 means the patient does not use ventilation
    p_ventilation = ventilation_duration[p_ID]# extracting ventilation information
    # print('p_ventilation',p_ventilation)
    feature_df = add_vent_flag(feature_df,p_ventilation,time_interval)

    # compute respiration score based on PaO2 and FiO2
    feature_df['PaO2_FiO2'] = feature_df['PaO2']/((feature_df['FiO2']+1)*0.01) # FiO2 add 1 to avoid  division; FiO2 is fraction
    feature_df['PaO2_FiO2'] = feature_df['PaO2_FiO2'].round(4)

    feature_df['Respiration_score'] = 0
    feature_df.loc[(feature_df['PaO2_FiO2'] >= 400)&(feature_df['vent'] == 0), 'Respiration_score'] = 0
    feature_df.loc[(feature_df['PaO2_FiO2'] >= 400) & (feature_df['vent'] == 1), 'Respiration_score'] = 3

    feature_df.loc[(feature_df['PaO2_FiO2'] <400) & (feature_df['PaO2_FiO2'] >=300)& (feature_df['vent'] == 0), 'Respiration_score'] = 1
    feature_df.loc[(feature_df['PaO2_FiO2'] <400) & (feature_df['PaO2_FiO2'] >=300)& (feature_df['vent'] == 1), 'Respiration_score'] = 3

    feature_df.loc[(feature_df['PaO2_FiO2'] <300) & (feature_df['PaO2_FiO2'] >0)& (feature_df['vent'] == 0), 'Respiration_score'] = 2
    feature_df.loc[(feature_df['PaO2_FiO2'] <300) & (feature_df['PaO2_FiO2'] >0)& (feature_df['vent'] == 1), 'Respiration_score'] = 3

    feature_df.loc[(feature_df['PaO2_FiO2'] <=200) & (feature_df['PaO2_FiO2'] >=100)& (feature_df['vent'] == 1), 'Respiration_score'] = 3
    feature_df.loc[(feature_df['PaO2_FiO2'] <100) & (feature_df['PaO2_FiO2'] >0)& (feature_df['vent'] == 1), 'Respiration_score'] = 4

    feature_df.loc[feature_df['PaO2_FiO2'] == 0, 'Respiration_score'] = 0  # missing value is filled with 0 and scoring 0

   # compute Cardiovascular score based on MAP and  drugs Dopamine # Dobutamine, Phenylephrine  Vasopressin # Epinephrine # Norepinephrine
    feature_df['Cardiovascular_score'] = 0
    p_drug = cardiovascular_score_drug[p_ID] # extracting drug information, to obtain drug_rate
    feature_df = compute_Cardiovascular_score(feature_df,p_drug,time_interval)

    # Compute SOFA score based on six subscore

    feature_df['SOFA_score'] = feature_df['Respiration_score'] + feature_df['Cardiovascular_score'] + feature_df['CNS_score'] + \
                               feature_df['Liver_score'] + feature_df['Coagulation_score'] + feature_df['Renal_score']

    # feature_df.to_csv(result_path + 'testing/person_feature_' + p_ID + '_finished.csv', index=True)
    feature_bucket_filled_missing_scoring[p_ID] = feature_df

    # feature_df.to_csv(result_path+'testing/person_feature_'+p_ID+'_scoring.csv',index=True)

    # if count>1000:
    #     break
    count = count + 1

# # save results feature_bucket_filled_missing_scoring
with open(result_path + "feature_bucket_filled_missing_scoring.json", "w") as outfile:
    json.dump(feature_bucket_filled_missing_scoring, outfile, cls=JSONEncoder, indent=4)

print('done.')
print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
