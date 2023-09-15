# obtaining sepsis onset based on infection and SOFA socre > 2
import pandas as pd
import numpy as  np
import json
import time

file_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/eicu-database-2.0/'
result_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eICU/sepsis_code/data/'

start_time = time.time()
# read feature_bucket_with_filled_missing
with open(result_path + 'feature_bucket_filled_missing_scoring.json', 'r') as load_f:
    SOFA_score = json.load(load_f)
print('len(SOFA_score)',len(SOFA_score))

with open(result_path + 'list_infection_time.json', 'r') as load_f:
    infection_time = json.load(load_f)

print('len(infection_time)',len(infection_time))


count = 0
interval = 360 # 24hour, 1440 min
sepsis_onset_df = pd.DataFrame(columns=['patientunitstayid','sepsis_onset'])

for i in SOFA_score:
    count = count + 1
    print('count',count)
    p_ID = i
    # print('p_ID',p_ID)
    p_sofa_score = SOFA_score[p_ID]
    p_sofa_score_df = pd.DataFrame.from_dict(p_sofa_score)
    # reset index, start from 0,1,...
    p_sofa_score_df = p_sofa_score_df.reset_index(drop=False)

# obtain SOFA >= 2 time index
    p_SOFA_time_index = p_sofa_score_df.index[p_sofa_score_df['SOFA_score'] >= 2].tolist()
    p_infection_time = infection_time[p_ID]

    if (len(p_SOFA_time_index)>0) & (len(p_infection_time)>0): # the patient needs to have infection time and SOFA
        p_SOFA_time_index = np.array(p_SOFA_time_index)
        p_SOFA_time = list((p_SOFA_time_index+1)*interval) # convert time to min, * interval

        # print('p_SOFA_time',p_SOFA_time)
        # print('p_infection_time',p_infection_time)

        SOFA_before24h = []
        SOFA_after12h  = []

# the rules for obtaining sepsis onset time based on infection time and SOFA score time
# reference: Early Prediction of Sepsis from Clinical Data the PhysioNet/Computing in Cardiology Challenge 2019
# tsepsis: Onset of sepsis identi ed as the earlier of tsuspicion and  tSOFA as long as tSOFA occurred no more than 24 hours before or 12 hours after tsuspicion.

        for i in p_infection_time:
            for j in p_SOFA_time:
                if (i-j>0) & (i-j<=1440):# before 24h of infection time
                    SOFA_before24h.append(j)# save SOFA time which is before 24 hour of infection time
                    # print('SOFA_before24h_ i,j',i,j)

                if (j-i>0) & (j-i<720):# after 12h of infection time
                    SOFA_after12h.append(i) # save infection time
                    # print('SOFA_after12h i,j', i, j)
        p_sepsis_candidate_time = SOFA_before24h + SOFA_after12h
        if len(p_sepsis_candidate_time)>0:
            p_sepsis_onset = min(p_sepsis_candidate_time)
            sepsis_onset_df.loc[len(sepsis_onset_df)] = [int(p_ID),p_sepsis_onset]
            # print('p_sepsis_onset',p_sepsis_onset)
    # p_sofa_score_df.to_csv(result_path+'testing/person_feature_'+p_ID+'_sofa_score_df.csv',index=False)

    # if count>1000:
    #     break

sepsis_onset_df.to_csv(result_path+'sepsis_onset_df.csv',index=False)

print('done.')
print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
