# build a json file to save ventilation duration information (start, end) for each patient
# this reference code from paper: A Machine-Learning Approach for Dynamic Prediction of Sepsis-Induced Coagulopathy in Critically Ill Patients With Sepsis

import pandas as pd
import numpy as  np
import json

file_path = '/Users/xuzhenxing/Documents/eicu-database-2.0/'
result_path = '/Users/xuzhenxing/Documents/eICU_AKI_Sepsis/'

# read  patients table to obtain discharged offset (time)
patient = pd.read_csv(file_path+'patient.csv',index_col=False)
patient_ICU_ID = list(set(patient['patientunitstayid']))
patient = patient[['patientunitstayid','unitdischargeoffset']]
patient['on_off'] = 0 # # add one column to be used for indicating the use of discharged and they do not use ventilation
patient = patient.rename(columns={'unitdischargeoffset': 'ventoffset'})

# patient.to_csv(result_path+'testing/patient_testing.csv',index=False)


## read respiratoryCare
ori_respiratoryCare = pd.read_csv(file_path+'respiratoryCare.csv',index_col=False)

respiratoryCare = ori_respiratoryCare
respiratoryCare = respiratoryCare.loc[respiratoryCare['ventstartoffset']!=0,['patientunitstayid','ventstartoffset']]
respiratoryCare['on_off']=1 # add one column to be used for indicating the use of ventilation
# change name of column "ventstartoffset"  to "ventoffset"
respiratoryCare = respiratoryCare.rename(columns={'ventstartoffset': 'ventoffset'})

# respiratoryCare.to_csv(result_path+'testing/respiratoryCare_testing.csv',index=False)

# read privous ventilation start information
prior_start_respiratoryCare = ori_respiratoryCare.loc[(ori_respiratoryCare['priorventstartoffset']!=0)|(ori_respiratoryCare['priorventendoffset']!=0),['patientunitstayid','priorventstartoffset']]
prior_start_respiratoryCare['on_off'] = 1
prior_start_respiratoryCare = prior_start_respiratoryCare.rename(columns={'priorventstartoffset': 'ventoffset'})

prior_start_respiratoryCare.to_csv(result_path+'testing/prior_start_respiratoryCare_testing.csv',index=False)

# read privous ventilation end information
prior_end_respiratoryCare = ori_respiratoryCare.loc[(ori_respiratoryCare['priorventstartoffset']!=0)|(ori_respiratoryCare['priorventendoffset']!=0),['patientunitstayid','priorventendoffset']]
prior_end_respiratoryCare['on_off'] = 0
prior_end_respiratoryCare = prior_end_respiratoryCare.rename(columns={'priorventendoffset': 'ventoffset'})

prior_end_respiratoryCare.to_csv(result_path+'testing/prior_end_respiratoryCare_testing.csv',index=False)


## read respiratoryCharting
respiratoryCharting = pd.read_csv(file_path+'respiratoryCharting.csv',index_col=False)
respchartvaluelabel = [
    'NIV Setting Set_RR'
    , 'Backup RR (Set)'
    , 'Plateau Pressure'
    , 'Total RR'
    , 'FIO2 (%)'
    , 'Oxygen Flow Rate'
    , 'Tidal Volume Observed (VT)'
    , 'VS RESP RATE'
    , 'PEEP'
    , 'RR Spont'
    , 'RSBI'
    , 'Set Fraction of Inspired Oxygen (FIO2)'
    , '4. RSBI (RR/Vt)'
    , 'Total RSBI'
    , 'Mean Airway Pressure'
    , 'Tidal Volume, Delivered'
    , 'Resp Rate Total'
    , 'RT Vent On/Off'
    , 'Spontaneous Respiratory Rate'
    , 'Adult Con Pt/Vent Spont Rate'
    , 'Minute Ventilation Set(L/min)'
    , 'Minute Volume, Spontaneous'
    , 'Mechanical Ventilator Mode'
    , 'RR (patient)'
    , 'FiO2'
    , 'Tidal Volume (set)'
    , 'Adult Con Alarms Backup RR'
    , 'Mechanical Ventilator High Tidal Volume Alarm'
    , 'Adult Con Setting Set RR'
]
respiratoryCharting = respiratoryCharting.loc[respiratoryCharting['respchartvaluelabel'].isin(respchartvaluelabel),['patientunitstayid','respchartoffset','respchartvaluelabel','respchartvalue']]
respiratoryCharting['on_off'] = 2 # 2 is "continued" for ventilation, 1 start, 0 end and Suspended
respiratoryCharting.loc[(respiratoryCharting['respchartvaluelabel']=='RT Vent On/Off')&(respiratoryCharting['respchartvalue']=='Start'),'on_off'] = 1
respiratoryCharting.loc[(respiratoryCharting['respchartvaluelabel']=='RT Vent On/Off')&(respiratoryCharting['respchartvalue']=='Suspended'),'on_off'] = 0
respiratoryCharting.loc[(respiratoryCharting['respchartvaluelabel']=='RT Vent On/Off')&(respiratoryCharting['respchartvalue'].isin(['Off','off'])),'on_off'] = 0
respiratoryCharting.loc[(respiratoryCharting['respchartvaluelabel']=='RT Vent On/Off')&(respiratoryCharting['respchartvalue']=='Continued'),'on_off'] = 2

respiratoryCharting = respiratoryCharting[['patientunitstayid','respchartoffset','on_off']]
respiratoryCharting = respiratoryCharting.rename(columns={'respchartoffset': 'ventoffset'})

respiratoryCharting.to_csv(result_path+'testing/respiratoryCharting_testing.csv',index=False)


# combine 5 tables: patients, respiratoryCare, prior_start_respiratoryCare, prior_end_respiratoryCare, respiratoryCharting

frames = [respiratoryCare, patient, prior_start_respiratoryCare, prior_end_respiratoryCare, respiratoryCharting]
respiratoryCare = pd.concat(frames,ignore_index=True)

respiratoryCare.to_csv(result_path+'testing/respiratoryCare_total_ori.csv',index=False)
print('len',len(respiratoryCare))

# delete replicated rows
respiratoryCare = respiratoryCare.drop_duplicates()

print('len(respiratoryCare)',len(respiratoryCare))
respiratoryCare.to_csv(result_path+'respiratoryCare_total.csv',index=False)


