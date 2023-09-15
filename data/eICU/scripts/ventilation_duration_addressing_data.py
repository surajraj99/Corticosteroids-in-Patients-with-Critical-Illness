# build a json file to save ventilation duration information (start, end) for each patient
# this reference code from paper: A Machine-Learning Approach for Dynamic Prediction of Sepsis-Induced Coagulopathy in Critically Ill Patients With Sepsis

import pandas as pd
import numpy as  np
import json

file_path = '/Users/xuzhenxing/Documents/eicu-database-2.0/'
result_path = '/Users/xuzhenxing/Documents/eICU_AKI_Sepsis/'

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super(NpEncoder, self).default(obj)

# read respiratoryCare_total table
respiratoryCare = pd.read_csv(result_path+'respiratoryCare_total.csv',index_col=False)
patient_ICU_ID = list(set(respiratoryCare['patientunitstayid']))
print('len',len(patient_ICU_ID))

data = respiratoryCare
data['next_offset'] = data.sort_values(['patientunitstayid','ventoffset'],ascending=True).groupby('patientunitstayid')['ventoffset'].shift(-1)# lead()
data['vent_lag'] = data.sort_values(['patientunitstayid','ventoffset'],ascending=True).groupby('patientunitstayid')['on_off'].shift(1)# lag()
data.loc[data['on_off']==0,'next_offset'] = data.loc[data['on_off']==0,'ventoffset']

# patient_ID = 141584# 141853,141436,     3145787,1572929,524351
# patient_sample = data[data['patientunitstayid']==patient_ID]
# patient_sample.to_csv(result_path+'testing/'+ str(patient_ID) +'_before_new_vent.csv')
#

# fill missing for next_offset and vent_lag
data.loc[~data['next_offset'].notna(),'next_offset'] = data.loc[~data['next_offset'].notna(),'ventoffset']
data.loc[~data['vent_lag'].notna(),'vent_lag'] = 0


# new_vent
data['new_vent'] = 0
data.loc[(data['vent_lag']==0)&(data['on_off']>0),'new_vent'] = 1
data.loc[(data['vent_lag']>0)&(data['on_off']==0),'new_vent'] = 0
data.loc[(data['vent_lag']>0)&(data['on_off']==1),'new_vent'] = 1
data.loc[(data['vent_lag']>0)&(data['on_off']==2),'new_vent'] = 0

data['num_vent'] = data.sort_values(['patientunitstayid','ventoffset'],ascending=True).groupby('patientunitstayid')['new_vent'].cumsum()

# delete rows if 'vent_lag'==0 and 'on_off'==0
data = data.loc[~((data['vent_lag']==0)&(data['on_off']==0))]

patient_ven_num = len(set(data['patientunitstayid']))
print('patient_ven_num',patient_ven_num)
data_ventilation = data

# compute ventilation duration
# patient_ID = 2100636# 141853,141436,3145787
# patient_sample = data[data['patientunitstayid']==patient_ID]
# patient_sample.to_csv(result_path+'testing/'+ str(patient_ID) +'num_vent.csv')

ventilation_duration = {} # dic is used to save vetilation duration for each patient

count = 0
for i in patient_ICU_ID:

    print('count',count)
    p_ven_duration = {}

    p_ID = i # i
    # print('p_ID',p_ID)
    p_ven = data_ventilation[data_ventilation['patientunitstayid']==p_ID]

    if len(p_ven)>0: # there is ventilation information for this patient
        num_vent = list(set(p_ven['num_vent']))
        # start and end are used to save start time and end time for each vent
        start = []
        end = []
        for j in num_vent:
            # print('j',j)
            each_start = (p_ven.loc[p_ven['num_vent']==j,'ventoffset']).to_numpy()
            each_start_min = min(each_start)

            each_end = (p_ven.loc[p_ven['num_vent']==j,'next_offset']).to_numpy()
            each_end_max = max(each_end)

            start.append(each_start_min)
            end.append(each_end_max)
            # print('each_start_min',each_start_min)
            # print('each_end_max', each_end_max)

        p_ven_duration['start'] = start
        p_ven_duration['end'] = end

# save ventilation duration for one patient
    ventilation_duration[str(p_ID)] = p_ven_duration

    # if count > 500:
    #     break
    # p_ven.to_csv(result_path + 'testing/' + str(p_ID) + '.csv', index=False)


    count = count + 1

# save dic file
with open(result_path + "ventilation_duration.json", "w") as outfile:
    json.dump(ventilation_duration, outfile,cls=NpEncoder,indent=4)

print('it is over.')


















