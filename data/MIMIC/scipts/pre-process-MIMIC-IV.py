#  MIMIC-IV
import pandas as pd
import numpy as np
import time
import math
from collections import defaultdict
import pickle as pkl
import ast
import datetime
from tslearn.metrics import cdist_dtw
import scipy.spatial as sp, scipy.cluster.hierarchy as hc
import seaborn as sns;
import matplotlib.pyplot as plt
sns.set(color_codes=True)
import matplotlib
matplotlib.use('TkAgg')
from sklearn.cluster import AgglomerativeClustering
import json



DURATION = 24  # hour #changed this to 24
TOTAL_DURATION = 24  # hour #changed this to 24
CHART_WINDOW_LEFT = -24  # hour

infection_time_threshold = 24  # hour
left_infection_time = -1 * 24  # hour  # -2 * 24, 0
right_infection_time = 1 * 24  # hour 7 * 24


class TimestampEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, pd.Timestamp):
            return obj.isoformat()
        return super(TimestampEncoder, self).default(obj)
    
def timestamp_hook(obj):
    if 'timestamp' in obj:
        obj['timestamp'] = pd.Timestamp(obj['timestamp'])
    return obj

def func_renal_score(creatinine, urine):
    renal_score = 0
    creatinine_score = 0
    urine_score = 0
    if pd.notna(creatinine):
        if creatinine > 0:
            if creatinine < 1.2:
                creatinine_score = 0
            elif creatinine >= 1.2 and creatinine < 2.0:
                creatinine_score = 1
            elif creatinine >= 2.0 and creatinine < 3.5:
                creatinine_score = 2
            elif creatinine >= 3.5 and creatinine < 5:
                creatinine_score = 3
            else:
                creatinine_score = 4

    if pd.notna(urine):
        if urine > 0:
            if urine <= 125:  # (125 = 500/4), 24 hours--500; 6 hours--125;
                urine_score = 3
            if urine <= 50:  # (50 = 200/4), 24 hours--200; 6 hours--50;
                urine_score = 4

    renal_score = max([creatinine_score, urine_score])

    return renal_score


def func_CNS_score(gcsEye, gcsVerbal, gcsMotor):
    GCS_total = 0
    CNS_score = 0
    gcsEye_score = 0
    gcsVerbal_score = 0
    gcsMotor_score = 0
    if pd.notna(gcsEye):
        gcsEye_score = gcsEye
    if pd.notna(gcsVerbal):
        gcsVerbal_score = gcsVerbal
    if pd.notna(gcsMotor):
        gcsMotor_score = gcsMotor

    GCS_total = gcsEye_score + gcsVerbal_score + gcsMotor_score
    if GCS_total == 15:
        CNS_score = 0
    elif (GCS_total >= 13) and (GCS_total <= 14):
        CNS_score = 1
    elif (GCS_total >= 10) and (GCS_total <= 12):
        CNS_score = 2
    elif (GCS_total >= 6) and (GCS_total <= 9):
        CNS_score = 3
    else:
        CNS_score = 4

    return CNS_score


def func_MAP_score(MAP):
    # MAP
    MAP_score = 0
    if pd.notna(MAP):
        if MAP > 0:
            if MAP >= 70:
                MAP_score = 0
            else:
                MAP_score = 1
    return MAP_score


def func_dopamine_score(dopamine):
    # dopamine
    dopamine_score = 0
    if pd.notna(dopamine):
        if dopamine > 0:
            if dopamine < 5:
                dopamine_score = 2
            elif (dopamine >= 5) and (dopamine <= 15):
                dopamine_score = 3
            else:
                dopamine_score = 4
    return dopamine_score


def func_dobutamine_score(dobutamine):
    # dobutamine
    dobutamine_score = 0
    if pd.notna(dobutamine):
        if dobutamine > 0:
            dobutamine_score = 2
    return dobutamine_score


def func_epinephrine_score(epinephrine):
    # epinephrine
    epinephrine_score = 0
    if pd.notna(epinephrine):
        if epinephrine > 0:
            if epinephrine <= 0.1:
                epinephrine_score = 3
            else:
                epinephrine_score = 4

    return epinephrine_score


def func_vasopressin_score(vasopressin):
    vasopressin_score = 0
    flag_vasopressin_4 = 0
    if pd.notna(vasopressin):
        if vasopressin > 0:
            vasopressin_score = 3
        if vasopressin > 0.4:
            flag_vasopressin_4 = 1
    return vasopressin_score, flag_vasopressin_4


def func_phenylephrine_score(phenylephrine):
    phenylephrine_score = 0
    if pd.notna(phenylephrine):
        if phenylephrine > 0:
            if phenylephrine < 200:
                phenylephrine_score = 2
            else:
                phenylephrine_score = 3
    return phenylephrine_score


def func_cardiovascular_score(MAP, dopamine, dobutamine, epinephrine, norepinephrine, vasopressin, phenylephrine):
    card_score = 0
    MAP_score = func_MAP_score(MAP)
    dopamine_score = func_dopamine_score(dopamine)
    dobutamine_score = func_dobutamine_score(dobutamine)
    epinephrine_score = func_epinephrine_score(epinephrine)
    norepinephrine_score = func_epinephrine_score(
        norepinephrine)  # the criteria of epinephrine and norepinephrine  are same.
    vasopressin_score, flag_vas_4 = func_vasopressin_score(vasopressin)
    phenylephrine_score = func_phenylephrine_score(phenylephrine)

    card_score = max(
        [MAP_score, dopamine_score, dobutamine_score, epinephrine_score, norepinephrine_score, vasopressin_score,
         phenylephrine_score])
    # card_score==4:  (Norepinephrine between 0-0.1 ug/kg/min or Epinephrine between 0-0.1 ug/kg/min or Dopamine between 5-15 ug/kg/min) and (Vasopressin > 0.4 units/min or Phenylephrine > 200 ug/min)
    if flag_vas_4 == 1 or phenylephrine_score == 3:
        if max(dopamine_score, epinephrine_score, norepinephrine_score) == 3:
            card_score = 4

    return card_score


def func_liver_score(bilirubin):
    liver_score = 0
    if pd.notna(bilirubin):
        if bilirubin > 0:
            if bilirubin < 1.2:
                liver_score = 0
            if (bilirubin >= 1.2) and (bilirubin < 2.0):
                liver_score = 1
            if (bilirubin >= 2.0) and (bilirubin < 6.0):
                liver_score = 2
            if (bilirubin >= 6.0) and (bilirubin < 12.0):
                liver_score = 3
            if bilirubin >= 12.0:
                liver_score = 4

    return liver_score


def func_coagulation_score(platelet):
    coagulation_score = 0
    if pd.notna(platelet):
        if platelet > 0:
            if platelet >= 150:
                coagulation_score = 0
            if platelet < 150:
                coagulation_score = 1
            if platelet < 100:
                coagulation_score = 2
            if platelet < 50:
                coagulation_score = 3
            if platelet < 20:
                coagulation_score = 4

    return coagulation_score


def func_respiration_score(pao2, fio2, ventilation):
    respiration_score = 0
    if pd.isna(ventilation):
        ventilation = 0
    if pd.notna(pao2) and pd.notna(fio2):
        if fio2 > 0:
            fio2 = fio2 / 100  # Note that, the unit of fio2 is %, we need to convert FiO2 into fraction
            pao2_fio2 = pao2 / fio2
            if pao2_fio2 >= 400:
                respiration_score = 0
            if pao2_fio2 < 400:
                respiration_score = 1
            if pao2_fio2 < 300:
                respiration_score = 2
            if (pao2_fio2 < 200) and (ventilation == 1):
                respiration_score = 3
            if (pao2_fio2 < 100) and (ventilation == 1):
                respiration_score = 4

    return respiration_score


def func_compute_subscore(fea_sofa_df):
    df = fea_sofa_df
    row_name_s = list(df.index.values)
    subscore_s = ['Respiration_score', 'Coagulation_score', 'Liver_score', 'Cardiovascular_score', 'CNS_score',
                  'Renal_score', 'SOFA_score']
    for subscore in subscore_s:
        row_index_name = subscore
        # compute sofa subscore every DURATION = 6 hours
        for i in range(int(TOTAL_DURATION / DURATION)):
            column_index_name = 'hour_' + str((i + 1) * DURATION)
            if subscore == 'Respiration_score':
                pao2 = df.loc['PaO2', column_index_name]
                fio2 = df.loc['FiO2', column_index_name]
                ventilation = df.loc['Ventilation', column_index_name]
                resp_score = func_respiration_score(pao2, fio2, ventilation)
                df.loc[row_index_name, column_index_name] = resp_score

            elif subscore == 'Coagulation_score':
                platelet = df.loc['Platelet', column_index_name]
                coag_score = func_coagulation_score(platelet)
                df.loc[row_index_name, column_index_name] = coag_score

            elif subscore == 'Liver_score':
                bilirubin = df.loc['Bilirubin', column_index_name]
                liver_score = func_liver_score(bilirubin)
                df.loc[row_index_name, column_index_name] = liver_score

            elif subscore == 'Cardiovascular_score':
                MAP = df.loc['MAP', column_index_name]
                dopamine = df.loc['Dopamine', column_index_name]
                dobutamine = df.loc['Dobutamine', column_index_name]
                epinephrine = df.loc['Epinephrine', column_index_name]
                norepinephrine = df.loc['Norepinephrine', column_index_name]
                vasopressin = df.loc['Vasopressin', column_index_name]
                phenylephrine = df.loc['Phenylephrine', column_index_name]
                card_score = func_cardiovascular_score(MAP, dopamine, dobutamine, epinephrine, norepinephrine,
                                                       vasopressin, phenylephrine)
                df.loc[row_index_name, column_index_name] = card_score

            elif subscore == 'CNS_score':
                gcsEye = df.loc['gcsEye', column_index_name]
                gcsVerbal = df.loc['gcsVerbal', column_index_name]
                gcsMotor = df.loc['gcsMotor', column_index_name]
                CNS_score = func_CNS_score(gcsEye, gcsVerbal, gcsMotor)
                df.loc[row_index_name, column_index_name] = CNS_score

            elif subscore == 'Renal_score':
                creatinine = df.loc['Creatinine', column_index_name]
                urine = df.loc['Urine', column_index_name]
                renal_score = func_renal_score(creatinine, urine)
                df.loc[row_index_name, column_index_name] = renal_score

            else:  # SOFA score
                df.loc[row_index_name, column_index_name] = df.loc['Respiration_score', column_index_name] + df.loc[
                    'Coagulation_score', column_index_name] + \
                                                            df.loc['Liver_score', column_index_name] + df.loc[
                                                                'Cardiovascular_score', column_index_name] + \
                                                            df.loc['CNS_score', column_index_name] + df.loc[
                                                                'Renal_score', column_index_name]

    sofa_list = list(df.loc['SOFA_score'])  # select sofa
    subscore_df = df
    return sofa_list, subscore_df


def compute_sofa(input_file, output_file):
    id_icustays_f = open(input_file + "icustays_features.pkl", "rb")
    id_icustays = pkl.load(
        id_icustays_f)  # id_icustays is dictionary, and icustay_id is key; value is a list [(),(without_fill_miss, with_fill_miss)] with each item: ()-> (intime, outtime, subject_id, hadm_id, first_careunit, last_careunit, los)
    print('len(id_icustays)', len(id_icustays))

    # build df for saving all patient's sofa score during 72 hours
    # column_names = []
    # for i in range(int(TOTAL_DURATION / DURATION)):
    #     column_names.append('hour_' + str((i + 1) * DURATION))
    #
    # sofa_result_df = pd.DataFrame(columns=['icustay_id']+column_names)

    pat_no = 0
    for key, value in id_icustays.items():
        print('pat_no = ', pat_no)
        # if pat_no > 100:
        #     break
        pat_no = pat_no + 1
        icustay_id = key
        feature_df = value[1][1]

        fea_add_row_index = ['Respiration_score', 'Coagulation_score', 'Liver_score', 'Cardiovascular_score',
                             'CNS_score', 'Renal_score', 'SOFA_score']  # it is used to save SOFA subscore
        row_index = list(feature_df.index.values) + fea_add_row_index
        feature_df = feature_df.reindex(row_index)  # add several rows with fea_add_row_index
        # choose features for compute SOFA score
        obj_fea_row_index = ['PaO2', 'FiO2', 'Ventilation', 'Respiration_score',
                             'Platelet', 'Coagulation_score',
                             'Bilirubin', 'Liver_score',
                             'MAP', 'Dopamine', 'Dobutamine', 'Epinephrine', 'Norepinephrine', 'Vasopressin',
                             'Phenylephrine', 'Cardiovascular_score',
                             'GCS', 'gcsEye', 'gcsVerbal', 'gcsMotor', 'CNS_score',
                             'Urine', 'Creatinine', 'Renal_score',
                             'SOFA_score']

        feature_df = feature_df.loc[obj_fea_row_index]
        fea_sofa_df = feature_df.copy()
        sofa_score, sofa_subscore_df = func_compute_subscore(fea_sofa_df)
        id_icustays[icustay_id].append((sofa_subscore_df))

        # sofa_score_list = [icustay_id] + sofa_score
        # sofa_result_df.loc[len(sofa_result_df)] = sofa_score_list

    # save id_icustays, which contains all subsofa scores
    # dump
    icustays_features_sofa = id_icustays
    pkl.dump(icustays_features_sofa, open(output_file + 'icustays_features_sofa.pkl', 'wb'))
    print('len(icustays_features_sofa):', len(icustays_features_sofa))

    # # read pkl
    # testing_f = open(output_file + "icustays_features_sofa.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # icustay_id = '35479615'
    # p_df_withmiss = testing[icustay_id][1][0]
    # p_df_withmiss.to_csv(output_file + 'testing/' + icustay_id + '_nonfill_testing.csv', index=True)
    # p_df_withoutmiss = testing[icustay_id][1][1]
    # p_df_withoutmiss.to_csv(output_file + 'testing/' + icustay_id + '_filled_testing.csv', index=True)
    #
    # p_df_sofa = testing[icustay_id][2]
    # p_df_sofa.to_csv(output_file + 'testing/' + icustay_id + '_sofa.csv', index=True)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
    return icustays_features_sofa


def func_fill_missing(df_nonfill):
    # variables need to be fillfed
    # Note that, do not fill ---- 'Dobutamine', 'Dopamine', 'Epinephrine','Norepinephrine', 'Phenylephrine', 'Vasopressin', 'Ventilation', 'Urine',
    obj_row = ['Albumin', 'ALT', 'AST', 'Bands', 'Bicarbonate', 'Bilirubin', 'BMI', 'BUN', 'Chloride', 'Creatinine',
               'CRP',
               'FiO2', 'GCS', 'gcsEye', 'gcsVerbal', 'gcsMotor', 'Glucose', 'Heart_rate', 'Hemoglobin', 'INR',
               'Lactate', 'Lymphocyte_count', 'Lymphocyte_percent',
               'MAP', 'PaO2', 'Platelet', 'RDW', 'Respiratory_rate', 'SO2', 'Sodium', 'Systolic_ABP', 'Temperature',
               'Troponin I', 'Troponin T', 'WBC']
    df_nonfill.loc[obj_row, :] = df_nonfill.loc[obj_row, :].fillna(method='pad', axis=1)
    df_nonfill.loc[obj_row, :] = df_nonfill.loc[obj_row, :].fillna(method='bfill', axis=1)
    df_filled = df_nonfill
    return df_filled


def func_fill_procedureevents_feature(pat_result_df, pat_intime, pat_outtime, pat_procedures_list, fea_name,
                                      item_id_list):
    # pat_procedures_list is list, each item is tuple: (starttime, endtime, itemid, value, valueuom)
    pat_intime = int(pat_intime.timestamp())  # convert to integer, seconds
    pat_outtime = int(pat_outtime.timestamp())

    for procedures in pat_procedures_list:
        starttime = procedures[0]
        starttime = int(starttime.timestamp())
        endtime = procedures[1]
        endtime = int(endtime.timestamp())
        itemid = procedures[2]

        if (itemid in item_id_list):
            # fill pat_result_df based on fea_name (row_index name)
            for i in range(int(TOTAL_DURATION / DURATION)):
                column_index_name = 'hour_' + str((i + 1) * DURATION)
                row_index_name = fea_name
                left_time = pat_intime + i * DURATION * 60 * 60
                right_time = pat_intime + (i + 1) * DURATION * 60 * 60

                pre_fea_value = pat_result_df.loc[row_index_name, column_index_name]
                if pd.isna(pre_fea_value) or (pre_fea_value == 0):
                    pat_result_df.loc[row_index_name, column_index_name] = 0
                    # check overlap between [left_time,right_time] and [starttime, endtime]
                    if left_time < pat_outtime:  # consider early death or discharged
                        if starttime >= right_time or endtime <= left_time:
                            no_overlap = 0
                        else:
                            pat_result_df.loc[row_index_name, column_index_name] = 1

    return pat_result_df


def func_fill_inputevents_feature(pat_result_df, pat_intime, pat_outtime, pat_inputs_list, fea_name, item_id_list,
                                  fea_direction):
    # pat_inputs is list, each item is tuple: (starttime, endtime, itemid, amount, amountuom, rate, rateuom, patientweight)
    pat_intime = int(pat_intime.timestamp())  # convert to integer, seconds
    pat_outtime = int(pat_outtime.timestamp())

    for inputs in pat_inputs_list:
        starttime = inputs[0]
        starttime = int(starttime.timestamp())
        endtime = inputs[1]
        endtime = int(endtime.timestamp())
        itemid = inputs[2]
        rate = inputs[5]

        dif_start_end = endtime - starttime
        flag_inclusion = True
        if (itemid in ['221662', '221289', '229617', '221906']) and (
                dif_start_end < 60 * 60):  # Note that, Catecholamine (Dopamine, Epinephrine, Norepinephrine) doses are given as ug/kg/min for at least 1 hour.
            flag_inclusion = False

        if (itemid in item_id_list) and flag_inclusion:
            # fill pat_result_df based on fea_name (row_index name)
            for i in range(int(TOTAL_DURATION / DURATION)):
                column_index_name = 'hour_' + str((i + 1) * DURATION)
                row_index_name = fea_name
                left_time = pat_intime + i * DURATION * 60 * 60
                right_time = pat_intime + (i + 1) * DURATION * 60 * 60
                # check overlap between [left_time,right_time] and [starttime, endtime]
                if left_time < pat_outtime:  # consider early death or discharged
                    pre_fea_value = pat_result_df.loc[row_index_name, column_index_name]
                    if starttime >= right_time or endtime <= left_time:
                        no_overlap = 0
                    else:
                        if pd.notna(pre_fea_value):
                            if fea_direction == 'max':
                                if pre_fea_value < float(rate):
                                    pat_result_df.loc[row_index_name, column_index_name] = float(rate)
                            if fea_direction == 'min':
                                if pre_fea_value > float(rate):
                                    pat_result_df.loc[row_index_name, column_index_name] = float(rate)
                        else:
                            pat_result_df.loc[row_index_name, column_index_name] = float(rate)
    return pat_result_df


def func_build_column_index(start_date, pat_intime, pat_outtime):
    column_index_name = ''
    left_sec = (start_date - pat_intime).total_seconds()
    if (left_sec >= CHART_WINDOW_LEFT * 60 * 60) and (left_sec <= TOTAL_DURATION * 60 * 60) and (
            start_date <= pat_outtime):  # lab [-24 hour, 72 hour]
        column_index_tem = left_sec / (60 * 60 * DURATION)
        if column_index_tem <= 0:
            column_index = 0
        elif column_index_tem > 0:
            column_index = 0
        # elif column_index_tem > 0 and column_index_tem < 12:
        #     column_index = int(column_index_tem)
        # else:
        #     column_index = 11  # column_index_tem == 12
        column_index_name = 'hour_' + str((column_index + 1) * DURATION)
    return column_index_name


def func_fill_charts_feature(flag_source, pat_result_df, pat_intime, pat_outtime, pat_charts, pat_outputs, fea_name,
                             item_id_list, fea_direction):
    # pat_charts_list[], each item is (start_date, itemid, value, valueuom)
    if flag_source == 'chartevents':
        pat_charts_list = pat_charts
    else:
        pat_charts_list = pat_outputs

    for chart in pat_charts_list:
        start_date = chart[0]
        start_date = pd.Timestamp(start_date)
        itemid = chart[1]
        value = chart[2]

        if itemid in item_id_list:
            column_index_name = func_build_column_index(start_date, pat_intime, pat_outtime)
            if len(column_index_name) > 0:  # generate column name successfully
                row_index_name = fea_name
                pre_fea_value = pat_result_df.loc[row_index_name, column_index_name]
                if pd.notna(pre_fea_value):
                    if flag_source == 'chartevents':  # 'chartevents', # if there are multiple values in duration, we need to choose one based on their direction
                        if fea_direction == 'max':
                            if pre_fea_value < float(value):
                                pat_result_df.loc[row_index_name, column_index_name] = float(value)
                        if fea_direction == 'min':
                            if pre_fea_value > float(value):
                                pat_result_df.loc[row_index_name, column_index_name] = float(value)
                    else:  # flag_source = 'outputevents' # if there are multiple values in duration, we need to Sum Urine
                        pat_result_df.loc[row_index_name, column_index_name] = float(value) + pre_fea_value
                else:
                    pat_result_df.loc[row_index_name, column_index_name] = float(value)

    return pat_result_df


def func_create_df_icustay(feature_df):
    ori_df = feature_df
    index_name = list(ori_df['FeatureName'])
    column_name = []
    for i in range(int(TOTAL_DURATION / DURATION)):
        name = 'hour_' + str((i + 1) * DURATION)
        column_name.append(name)
    new_df = pd.DataFrame(columns=column_name, index=index_name)
    return new_df


def pre_features(input_file, output_file):
    # read original information
    id_icustays, id_inputevents, id_outputevents, id_chartevents, id_procedureevents, feature_df = func_read_ori_data(
        input_file)
    pat_no = 0
    print('len(id_icustays) = ', len(id_icustays))
    for key, value in id_icustays.items():
        print('pat_no = ', pat_no)
        # if pat_no > 20:
        #     break
        pat_no = pat_no + 1
        icustay_id = key
        intime = value[0][0]  # value is list where each item is tuple
        outtime = value[0][1]

        # extract all information of an icustay
        p_charts = id_chartevents.get(icustay_id, [])
        p_inputs = id_inputevents.get(icustay_id, [])
        p_outputs = id_outputevents.get(icustay_id, [])
        p_procedures = id_procedureevents.get(icustay_id, [])

        # create empty dataframe p_df for an icustay
        p_df = func_create_df_icustay(feature_df)
        # fill p_df
        for i, row in feature_df.iterrows():
            # print('i = ', i)
            featureName = row['FeatureName']
            source = row['Source']
            itemid_in_source = row['Itemid_in_source']
            featureDirection = row['FeatureDirection']
            unit = row['Unit']

            if pd.notna(itemid_in_source):
                itemid_str = ast.literal_eval(itemid_in_source)  # convert str to list
                itemid_list = list(map(str, itemid_str))  # convert each int item into str
            else:
                itemid_list = []

            if (source in ['chartevents', 'outputevents']) and itemid_list:
                p_df = func_fill_charts_feature(source, p_df, intime, outtime, p_charts, p_outputs, featureName,
                                                itemid_list, featureDirection)
            if (source in ['inputevents']) and itemid_list:
                p_df = func_fill_inputevents_feature(p_df, intime, outtime, p_inputs, featureName, itemid_list,
                                                     featureDirection)
            if (source in ['procedureevents']) and itemid_list:
                p_df = func_fill_procedureevents_feature(p_df, intime, outtime, p_procedures, featureName, itemid_list)

        p_df_withmiss = p_df.copy()
        p_df_withoutmiss = func_fill_missing(p_df)
        # p_df_withmiss.to_csv(output_file + 'testing/' + icustay_id + '_nonfill.csv', index=True)
        # p_df_withoutmiss.to_csv(output_file + 'testing/' + icustay_id + '_filled.csv', index=True)

        id_icustays[icustay_id].append((p_df_withmiss, p_df_withoutmiss))

    # dump
    icustays_features = id_icustays
    pkl.dump(icustays_features, open(output_file + 'icustays_features.pkl', 'wb'))
    print('len(icustays_features):', len(icustays_features))

    # # read pkl
    # testing_f = open(output_file + "icustays_features.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # icustay_id = '39553978'
    # p_df_withmiss = testing[icustay_id][1][0]
    # p_df_withmiss.to_csv(output_file + 'testing/' + icustay_id + '_nonfill_testing.csv', index=True)
    # p_df_withoutmiss = testing[icustay_id][1][1]
    # p_df_withoutmiss.to_csv(output_file + 'testing/' + icustay_id + '_filled_testing.csv', index=True)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
    return icustays_features


def pre_inputevents(input_file, output_file):
    input = pd.read_csv(input_file + 'ICU/' + 'inputevents.csv.gz', compression='gzip', index_col=False, dtype=str)
    # input = pd.read_csv(input_file + 'MIMIC_IV_others/' + 'Med_DDENVP.csv',  index_col=False, dtype=str)
    print('len(input)', len(input))
    icustayid = set(list(input['stay_id']))
    patid = set(list(input['subject_id']))
    print('len(patid)', len(patid))
    print('len(icustayid)', len(icustayid))

    id_indexrecord = input
    id_inputevents = defaultdict(list)  # save result
    n_no_stay_id = 0
    n_no_starttime = 0
    n_no_endtime = 0
    n_no_itemid = 0
    n_no_amountuom = 0
    n_no_amount = 0
    negative_no_amount = 0
    n_no_patientweight = 0

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>100000:
        #     break
        stay_id = row['stay_id']
        starttime = pd.to_datetime(row['starttime'])  # there is no nan value for starttime and endtime
        endtime = pd.to_datetime(row['endtime'])

        itemid = row['itemid']
        amount = row['amount']
        amountuom = row['amountuom']

        rate = row['rate']
        rateuom = row['rateuom']
        patientweight = row['patientweight']

        if pd.isna(patientweight):
            n_no_patientweight += 1

        if pd.isna(stay_id):
            n_no_stay_id += 1
        elif pd.isna(starttime):
            n_no_starttime += 1
        elif pd.isna(endtime):
            n_no_endtime += 1
        elif pd.isna(itemid):
            n_no_itemid += 1
        elif pd.isna(amountuom):
            n_no_amountuom += 1
        elif pd.isna(amount):
            n_no_amount += 1
        elif float(amount) <= 0:
            negative_no_amount += 1
        else:
            # address Epinephrine,  itemid = [221289, 229617], itemid= 229617, there is no rate, need to compute it. Rate = amount / ((endtime - starttime)*patientweight), note that, rate unit: mcg/kg/min. The unit of amount is mg or mcg.
            if itemid == '229617':
                amount = float(amount)
                if amountuom == 'mg':
                    amount = 1000 * amount
                time_seconds = (endtime - starttime).total_seconds()
                time_mins = time_seconds / 60
                rate = amount / (time_mins * float(patientweight))

                # update amount, amountuom, rate, rateuom
                amount = str(amount)
                amountuom = 'mcg'
                rate = str(rate)
                rateuom = 'mcg/kg/min'

            # address Norepinephrine, itemid = [221906] # Including rate unit: mg/kg/min and mcg/kg/min. Neet to convert rate unit.
            if itemid == '221906':
                if rateuom == 'mg/kg/min':
                    rate = 1000 * float(rate)
                    rate = str(rate)
                    rateuom = 'mcg/kg/min'

            # address Vasopressin: itemid = [222315] # Including rate unit: units/hour and units/min. Need to convert rate unit.
            if itemid == '222315':
                if rateuom == 'units/hour':
                    rate = float(rate) / 60
                    rate = str(rate)
                    rateuom = 'units/min'

            # address Phenylephrine: itemid = [221749, 229630, 229632]  #
            # Items = 221749, including rate unit: mcg/kg/min and mcg/min; for mcg/kg/min, it need to be converted to mcg/min by considering weight.
            if itemid in ['221749', '229630', '229632']:
                if rateuom == 'mcg/kg/min':
                    rate = float(rate) * float(patientweight)
                    rate = str(rate)
                    rateuom = 'mcg/min'

            id_inputevents[stay_id].append((starttime, endtime, itemid, amount, amountuom, rate, rateuom,
                                            patientweight))  # note that, we used stay_id as key

    print('n_no_stay_id', n_no_stay_id, 'n_no_starttime', n_no_starttime, 'n_no_endtime', n_no_endtime, 'n_no_itemid',
          n_no_itemid, 'n_no_amountuom', n_no_amountuom, 'n_no_amount', n_no_amount, 'negative_no_amount',
          negative_no_amount, 'n_no_patientweight', n_no_patientweight)

    # sort
    print('Sort id_inputevents lists in id_inputevents by starttime')
    for patid, item_list in id_inputevents.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_inputevents[patid] = item_list_sorted

    # dump
    pkl.dump(id_inputevents, open(output_file + 'inputevents.pkl', 'wb'))
    print('len(id_inputevents):', len(id_inputevents))

    # # read pkl
    # testing_f = open(output_file + "inputevents.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['39765666'] # 28662225, 39553978
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_inputevents


def pre_outputevents(input_file, output_file):
    output = pd.read_csv(input_file + 'ICU/' + 'outputevents.csv.gz', compression='gzip', index_col=False, dtype=str)
    print('len(output)', len(output))
    icustayid = set(list(output['stay_id']))
    patid = set(list(output['subject_id']))
    print('len(patid)', len(patid))
    print('len(icustayid)', len(icustayid))

    id_indexrecord = output
    id_outputevents = defaultdict(list)  # save result

    n_no_stay_id = 0
    n_no_start_date = 0
    n_no_itemid = 0
    n_no_value = 0

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>100000:
        #     break
        stay_id = row['stay_id']
        charttime = pd.to_datetime(row['charttime'])  # there is no nan value for starttime and storetime
        storetime = pd.to_datetime(row['storetime'])

        itemid = row['itemid']
        value = row['value']
        valueuom = row['valueuom']

        # start_date
        if pd.notna(charttime):
            start_date = charttime
        elif pd.notna(storetime):
            start_date = storetime
        else:
            start_date = np.nan

        if pd.isna(stay_id):
            n_no_stay_id += 1
        elif pd.isna(start_date):
            n_no_start_date += 1
        elif pd.isna(itemid):
            n_no_itemid += 1
        elif pd.isna(value):
            n_no_value += 1
        else:
            id_outputevents[stay_id].append((start_date, itemid, value, valueuom))  # note that, we used stay_id as key

    print('n_no_stay_id', n_no_stay_id, 'n_no_start_date', n_no_start_date, 'n_no_itemid', n_no_itemid, 'n_no_value',
          n_no_value)

    # sort
    print('Sort id_outputevents lists in id_outputevents by starttime')
    for patid, item_list in id_outputevents.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_outputevents[patid] = item_list_sorted

    # dump
    pkl.dump(id_outputevents, open(output_file + 'outputevents.pkl', 'wb'))
    print('len(id_outputevents):', len(id_outputevents))

    # # read pkl
    # testing_f = open(output_file + "outputevents.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['39765666'] # 28662225, 39553978
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_outputevents


def pre_procedureevents(input_file, output_file):
    proc = pd.read_csv(input_file + 'ICU/' + 'procedureevents.csv.gz', compression='gzip', index_col=False, dtype=str)
    print('len(chart)', len(proc))
    icustayid = set(list(proc['stay_id']))
    patid = set(list(proc['subject_id']))
    print('len(patid)', len(patid))
    print('len(icustayid)', len(icustayid))

    id_indexrecord = proc
    id_procedureevents = defaultdict(list)  # save result

    n_no_stay_id = 0
    n_no_starttime = 0
    n_no_endtime = 0
    n_no_itemid = 0
    n_no_value = 0

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>1000:
        #     break
        stay_id = row['stay_id']
        starttime = pd.to_datetime(row['starttime'])  # there is no nan value for starttime and endtime
        endtime = pd.to_datetime(row['endtime'])

        itemid = row['itemid']
        value = row['value']
        valueuom = row['valueuom']

        if pd.isna(stay_id):
            n_no_stay_id += 1
        elif pd.isna(starttime):
            n_no_starttime += 1
        elif pd.isna(endtime):
            n_no_endtime += 1
        elif pd.isna(itemid):
            n_no_itemid += 1
        elif pd.isna(value):
            n_no_value += 1
        else:
            id_procedureevents[stay_id].append(
                (starttime, endtime, itemid, value, valueuom))  # note that, we used stay_id as key

    print('n_no_stay_id', n_no_stay_id, 'n_no_starttime', n_no_starttime, 'n_no_endtime', n_no_endtime, 'n_no_itemid',
          n_no_itemid, 'n_no_value', n_no_value)

    # sort
    print('Sort id_procedureevents lists in id_procedureevents by starttime')
    for patid, item_list in id_procedureevents.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_procedureevents[patid] = item_list_sorted

    # dump
    pkl.dump(id_procedureevents, open(output_file + 'procedureevents.pkl', 'wb'))
    print('len(id_procedureevents):', len(id_procedureevents))

    # # read pkl
    # testing_f = open(output_file + "procedureevents.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['31205490'] # 28662225
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_procedureevents


def pre_chartevents(input_file, output_file):
    chart = pd.read_csv(input_file + 'ICU/' + 'chartevents.csv.gz', compression='gzip', index_col=False, dtype=str)
    # chart = pd.read_csv(input_file + 'MIMIC_IV_others/'+'platelets.csv',index_col=False, dtype=str)
    print('len(chart)', len(chart))
    icustayid = set(list(chart['stay_id']))
    patid = set(list(chart['subject_id']))
    print('len(patid)', len(patid))
    print('len(icustayid)', len(icustayid))

    # chart['charttime'] = pd.to_datetime(chart['charttime'])
    # chart['storetime'] = pd.to_datetime(chart['storetime'])

    id_indexrecord = chart
    id_chartevents = defaultdict(list)  # save result
    n_no_stay_id = 0
    n_no_time = 0
    n_no_itemid = 0
    n_no_value = 0

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>1000:
        #     break
        stay_id = row['stay_id']

        # charttime = row['charttime']
        # storetime = row['storetime']

        charttime = pd.to_datetime(row['charttime'])
        storetime = pd.to_datetime(row['storetime'])
        itemid = row['itemid']
        value = row['value']
        valueuom = row['valueuom']

        # address GCS value, using valuenum
        if int(itemid) in [220739, 223900, 223901, 226756, 226757, 226758]:
            value = row['valuenum']

        # start_date
        if pd.notna(charttime):
            start_date = charttime
        elif pd.notna(storetime):
            start_date = storetime
        else:
            start_date = np.nan

        if pd.isna(stay_id):
            n_no_stay_id += 1
        elif pd.isna(start_date):
            n_no_time += 1
        elif pd.isna(itemid):
            n_no_itemid += 1
        elif pd.isna(value):
            n_no_value += 1
        else:
            id_chartevents[stay_id].append((start_date, itemid, value, valueuom))  # note that, we used stay_id as key

    print('Readlines:', i, 'n_no_stay_id:', n_no_stay_id, 'n_no_time:', n_no_time, 'n_no_itemid:', n_no_itemid,
          'n_no_value:', n_no_value, 'len(id_chartevents):', len(id_chartevents))

    # sort
    print('Sort id_chartevents lists in id_chartevents by start_date')
    for patid, item_list in id_chartevents.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_chartevents[patid] = item_list_sorted

    # dump
    # pkl.dump(id_chartevents, open(output_file + 'chartevents.pkl', 'wb'), protocol=pkl.HIGHEST_PROTOCOL)
    # print('len(id_chartevents):', len(id_chartevents))

    # # Save the large dictionary to a JSON file
    with open(output_file + 'chartevents.json', 'w') as f:
        json.dump(id_chartevents, f, cls=TimestampEncoder)
    print('len(id_chartevents):', len(id_chartevents))

    # # read pkl
    # testing_f = open(output_file + "chartevents.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['31205490'] # 28662225
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_chartevents


def pre_icustays(input_file, output_file):
    ori_icustays = pd.read_csv(input_file + 'ICU/' + 'icustays.csv.gz', compression='gzip', index_col=False, dtype=str)
    print('len(ori_icustays)', len(ori_icustays))
    patid = set(list(ori_icustays['subject_id']))

    print('patid', len(patid))
    icustays_id = set(list(ori_icustays['stay_id']))
    print('icustays_id', len(icustays_id))

    ori_icustays['intime'] = pd.to_datetime(ori_icustays['intime'])
    ori_icustays['outtime'] = pd.to_datetime(ori_icustays['outtime'])

    id_icustays = defaultdict(list)
    id_icustays_final = defaultdict(list)

    id_indexrecord = ori_icustays

    n_no_intime = 0

    for i, row in id_indexrecord.iterrows():
        if i % 100 == 0:
            print('i = ', i)
        # if i>100:
        #     break
        subject_id = row['subject_id']
        hadm_id = row['hadm_id']
        stay_id = row['stay_id']

        first_careunit = row['first_careunit']
        last_careunit = row['last_careunit']
        intime = row['intime']
        outtime = row['outtime']
        los = row['los']

        if pd.isna(intime):
            n_no_intime += 1
        else:
            # id_icustays[hadm_id].append((intime, outtime, subject_id, stay_id, first_careunit, last_careunit, los)) # note that, we used hadm_id as key
            id_icustays[subject_id].append((intime, outtime, hadm_id, stay_id, first_careunit, last_careunit,
                                            los))  # note that, we used subject_id as key

    print('Readlines:', i, 'n_no_intime:', n_no_intime, 'len(id_icustays):', len(id_icustays))

    # sort
    print('Sort id_icustays lists in id_icustays by intime')
    for patid, item_list in id_icustays.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        first_icustay = item_list_sorted[0]  # only select the first icustay
        first_icustay = list(first_icustay)
        pat_icustay_id = first_icustay[3]  # icustay_id
        pat_hadm_id = first_icustay[2]  # hadm_id

        first_icustay[2] = patid  # subject_id
        first_icustay[3] = pat_hadm_id  # hadm_id
        id_icustays_final[pat_icustay_id] = [tuple(
            first_icustay)]  # for id_icustays_final: (intime, outtime, subject_id, hadm_id, first_careunit, last_careunit, los)

    # dump
    pkl.dump(id_icustays_final, open(output_file + 'icustays.pkl',
                                     'wb'))  # for id_icustays_final: (intime, outtime, subject_id, hadm_id, first_careunit, last_careunit, los)

    print('len(id_icustays_final):', len(id_icustays_final))

    # # read pkl
    # testing_f = open(output_file + "icustays.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['35479615'] # 28662225
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_icustays


def pre_patients(input_file, output_file):
    ori_patients = pd.read_csv(input_file + 'Hosp/' + 'patients.csv.gz', compression='gzip', index_col=False, dtype=str)
    print('len(ori_patients)', len(ori_patients))
    patid = set(list(ori_patients['subject_id']))
    print('len(patid)', len(patid))

    id_patients = defaultdict(list)
    id_indexrecord = ori_patients

    n_no_subject_id = 0

    for i, row in id_indexrecord.iterrows():
        if i % 100 == 0:
            print('i = ', i)
        # if i>100:
        #     break
        subject_id = row['subject_id']
        gender = row['gender']
        age = row['anchor_age']
        anchor_year = row['anchor_year']
        anchor_year_group = row['anchor_year_group']
        dod = row['dod']
        if pd.isna(dod):
            dod = np.nan
        else:
            dod = pd.Timestamp(dod)

        if pd.isna(subject_id):
            n_no_subject_id += 1
        else:
            id_patients[subject_id].append((gender, age, anchor_year, anchor_year_group, dod))

    print('Readlines:', i, 'n_no_subject_id:', n_no_subject_id, 'len(ori_patients):', len(ori_patients))

    # dump
    pkl.dump(id_patients, open(output_file + 'patients.pkl', 'wb'))

    print('len(id_patients):', len(id_patients))

    # # read pkl
    # testing_f = open(output_file + "patients.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['35479615'] # 28662225
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_patients


def pre_services(input_file, output_file):
    ori_services = pd.read_csv(input_file + 'Hosp/' + 'services.csv.gz', compression='gzip', index_col=False, dtype=str)
    print('len(ori_services)', len(ori_services))

    id_services = defaultdict(list)
    id_indexrecord = ori_services

    n_no_subject_id = 0

    for i, row in id_indexrecord.iterrows():
        if i % 100 == 0:
            print('i = ', i)
        # if i>100:
        #     break
        subject_id = row['subject_id']
        hadm_id = row['hadm_id']
        transfertime = row['transfertime']
        if pd.isna(transfertime):
            transfertime = np.nan
        else:
            transfertime = pd.Timestamp(transfertime)

        prev_service = row['prev_service']
        curr_service = row['curr_service']

        if pd.isna(subject_id):
            n_no_subject_id += 1
        else:
            id_services[subject_id].append((transfertime, hadm_id, prev_service, curr_service))

    print('Readlines:', i, 'n_no_subject_id:', n_no_subject_id, 'len(id_services):', len(id_services))

    # sort
    print('Sort id_services lists in id_services by transfertime')
    for patid, item_list in id_services.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_services[patid] = item_list_sorted

    # dump
    pkl.dump(id_services, open(output_file + 'services.pkl', 'wb'))
    print('len(id_services):', len(id_services))

    # # read pkl
    # testing_f = open(output_file + "services.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['31205490'] # 28662225
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_services


def pre_culture(input_file, output_file):
    ori_culture = pd.read_csv(input_file + 'Hosp/' + 'microbiologyevents.csv.gz', compression='gzip', index_col=False,
                              dtype=str)
    print('len(ori_culture)', len(ori_culture))
    # patid = set(list(ori_culture['subject_id']))
    # print('patid',len(patid))

    ori_culture['chartdate'] = pd.to_datetime(ori_culture['chartdate'])
    ori_culture['charttime'] = pd.to_datetime(ori_culture['charttime'])

    ori_culture['storedate'] = pd.to_datetime(ori_culture['storedate'])
    ori_culture['storetime'] = pd.to_datetime(ori_culture['storetime'])

    id_culture = defaultdict(list)

    id_indexrecord = ori_culture
    n_no_charttime = 0

    spec_type_desc_no = 0
    spec_type_desc_exception = []

    test_name_no = 0
    test_name_exception = []

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>1000:
        #     break
        subject_id = row['subject_id']
        hadm_id = row['hadm_id']

        charttime = row['charttime']
        chartdate = row['chartdate']

        storetime = row['storetime']
        storedate = row['storedate']

        spec_type_desc = row['spec_type_desc']
        try:
            spec_type_desc = spec_type_desc.strip()
        except:
            spec_type_desc_no += 1
            spec_type_desc_exception.append(spec_type_desc)

        test_name = row['test_name']
        try:
            test_name = test_name.strip()
        except:
            test_name_no += 1
            test_name_exception.append(test_name)

        # start_date
        if pd.notna(charttime):
            start_date = charttime
        elif pd.notna(chartdate):
            start_date = chartdate
        elif pd.notna(storetime):
            start_date = storetime
        elif pd.notna(storedate):
            start_date = storedate
        else:
            start_date = np.nan

        if pd.isna(start_date):
            n_no_charttime += 1
        else:
            id_culture[subject_id].append((start_date, hadm_id, spec_type_desc, test_name))

    print('Readlines:', i, 'n_no_charttime:', n_no_charttime, 'spec_type_desc_no:', spec_type_desc_no, 'test_name_no:',
          test_name_no, 'len(id_culture):', len(id_culture))

    print('len(spec_type_desc_exception)', len(spec_type_desc_exception))
    print('spec_type_desc_exception', spec_type_desc_exception)

    print('len(test_name_exception)', len(test_name_exception))
    print('test_name_exception', test_name_exception)

    # sort
    print('Sort id_culture lists in id_culture by charttime')
    for patid, item_list in id_culture.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_culture[patid] = item_list_sorted

    # dump
    pkl.dump(id_culture, open(output_file + 'microbiologyevents.pkl', 'wb'))

    print('len(id_culture):', len(id_culture))

    # # read pkl
    # testing_f = open(output_file + "microbiologyevents.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['10000032']
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_culture


def pre_prescription(input_file, output_file):
    ori_prescription = pd.read_csv(input_file + 'Hosp/' + 'prescriptions.csv.gz', compression='gzip', index_col=False,
                                   dtype=str)
    print('len(ori_prescription)', len(ori_prescription))

    # patid = set(list(ori_prescription['subject_id']))
    # print('patid',len(patid))

    ori_prescription['starttime'] = pd.to_datetime(ori_prescription['starttime'])
    ori_prescription['stoptime'] = pd.to_datetime(ori_prescription['stoptime'])

    id_prescription = defaultdict(list)

    id_indexrecord = ori_prescription
    n_no_starttime = 0

    drug_type_no = 0
    drug_type_exception = []

    drug_no = 0
    drug_exception = []

    route_no = 0
    route_exception = []

    for i, row in id_indexrecord.iterrows():
        if i % 1000 == 0:
            print('i = ', i)
        # if i>1000:
        #     break
        subject_id = row['subject_id']
        hadm_id = row['hadm_id']

        starttime = row['starttime']
        stoptime = row['stoptime']

        drug_type = row['drug_type']
        try:
            drug_type = drug_type.strip()
        except:
            drug_type_no += 1
            drug_type_exception.append(drug_type)

        drug = row['drug']
        try:
            drug = drug.strip()
        except:
            drug_no += 1
            drug_exception.append(drug)

        route = row['route']
        try:
            route = route.strip()
        except:
            route_no += 1
            route_exception.append(route)

        if pd.isna(starttime):
            n_no_starttime += 1
        else:
            id_prescription[subject_id].append((starttime, stoptime, hadm_id, drug_type, drug, route))

    print('Readlines:', i, 'n_no_starttime:', n_no_starttime, 'drug_type_no:', drug_type_no, 'drug_no:', drug_no,
          'route_no:', route_no, 'len(id_prescription):', len(id_prescription))

    print('len(drug_type_exception)', len(drug_type_exception))
    print('drug_type_exception', drug_type_exception)

    print('len(drug_exception)', len(drug_exception))
    print('drug_exception', drug_exception)

    print('len(route_exception)', len(route_exception))
    print('route_exception', route_exception)

    # sort
    print('Sort id_prescription lists in id_prescription by starttime')
    for patid, item_list in id_prescription.items():
        # add a set operation to reduce duplicates
        # sorted returns a sorted list
        # print('item_list',item_list)
        item_list_sorted = sorted(set(item_list), key=lambda x: x[0])
        # print('item_list_sorted', item_list_sorted)

        id_prescription[patid] = item_list_sorted

        # if patid =='10000032':
        #     break

    # dump
    pkl.dump(id_prescription, open(output_file + 'prescription.pkl', 'wb'))

    print('len(id_prescription):', len(id_prescription))

    # # read pkl
    # testing_f = open(output_file + "prescription.pkl", "rb")
    # testing = pkl.load(testing_f)
    # print('len',len(testing))
    # p = testing['10000032']
    # print('p',p)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return id_prescription


def extract_antibiotic_time(pat_prescription, icu_intime):
    antibiotic_time = []  # save antiobic time (starttime,stoptime)
    antibiotic_drug = ['bicillin l-a', 'rifadin', 'minocycline', 'vibra-tabs', 'clarithromycin', 'zinacef', 'eryc',
                       'monurol', 'cefepime', 'raniclor', 'trimethoprim', 'ceftibuten', 'factive',
                       'vibramycin', 'avidoxy', 'septra', 'proquin', 'sulfisoxazole', 'ery-tab', 'biaxin',
                       'chloramphenicol', 'myrac', 'rocephin', 'zmax', 'sulfamethoxazole', 'periostat',
                       'sulfatrim', 'cubicin', 'cefditoren', 'morgidox', 'ocudox', 'moxatag', 'macrobid', 'amikacin',
                       'bactocill', 'oxacillin', 'cefaclor', 'nafcillin sodium', 'ciprofloxacin',
                       'cefotetan', 'sulfadiazine', 'cedax', 'cefazolin', 'bactrim', 'panixine', 'cefadroxil',
                       'tazobactam', 'spectracef', 'noroxin', 'duricef', 'flagyl', 'streptomycin sulfate',
                       'pfizerpen', 'ceftin', 'cefotaxime', 'cefprozil', 'ofloxacin', 'solodyn', 'cleocin',
                       'macrodantin', 'amikin', 'cipro', 'omnicef', 'adoxa', 'dicloxacillin', 'azithromycin',
                       'vancomycin', 'dynacin', 'oracea', 'amoxicillin/clavulanate', 'rifampin', 'doxycycline',
                       'oraxyl', 'septra ds', 'timentin', 'ampicillin', 'vibativ', 'ala-tet', 'cefpodoxime',
                       'unasyn', 'azactam', 'metronidazole', 'amoxicillin', 'vantin', 'clavulanate', 'cayston',
                       'piperacillin', 'zosyn', 'cefdinir', 'erythrocin', 'nicazel doxy 30', 'minocin',
                       'zyvox', 'alodox', 'axetil', 'doryx', 'tobramycin', 'smz-tmp', 'nitrofurantoin', 'vancocin',
                       'cefoxitin', 'cefuroxime', 'augmentin', 'moxifloxacin', 'primsol', 'suprax',
                       'eryped', 'tetracycline', 'ceftazidime', 'penicillin', 'zithromax', 'cephalexin', 'monodox',
                       'pce dispertab', 'clindamycin', 'aztreonam', 'synercid', 'pediazole', 'tobi',
                       'mefoxin', 'avelox', 'tazicef', 'erythromycin', 'bethkis', 'maxipime']

    # print('len(antibiotic_drug)',len(antibiotic_drug)) # 125 antibiotic drug

    for item in pat_prescription:
        starttime = item[0]
        stoptime = item[1]
        drug_type = item[3]
        drug = item[4]
        route = item[5]

        # print('item',item)
        if pd.notna(drug_type):
            if drug_type in ['MAIN', 'ADDITIVE']:
                if pd.notna(route):
                    if route not in ['OU', 'OS', 'OD', 'AU', 'AS', 'AD', 'TP']:
                        if ('ear' not in route) and ('eye' not in route):
                            if pd.notna(drug):
                                if ('cream' not in drug.lower()) and ('desensitization' not in drug.lower()) and (
                                        'ophth oint' not in drug.lower()) and ('gel' not in drug.lower()):
                                    for antibiotic in antibiotic_drug:
                                        if antibiotic in drug.lower():
                                            if ((
                                                        starttime - icu_intime).total_seconds() >= left_infection_time * 60 * 60) and (
                                                    (
                                                            starttime - icu_intime).total_seconds() <= right_infection_time * 60 * 60):
                                                antibiotic_time.append(starttime)
                                                break
    antibiotic_time = list(set(antibiotic_time))

    return antibiotic_time


def extract_culture_time(pat_culture, icu_intime):
    culture_time = []
    for item in pat_culture:
        pat_start_date = item[0]

        if ((pat_start_date - icu_intime).total_seconds() >= left_infection_time * 60 * 60) and (
                (pat_start_date - icu_intime).total_seconds() <= right_infection_time * 60 * 60):
            culture_time.append(pat_start_date)

    culture_time = list(set(culture_time))

    return culture_time


def identify_infection_time(output_file):
    # create df to save result: icustay_id, flag_infection
    infec_time_df = pd.DataFrame(columns=['icustay_id', 'flag_infection'])

    id_icustays_f = open(output_file + "icustays" + ".pkl", "rb")
    # id_icustays is dictionary, and icustay_id is key; value is a list with each item: such assample: (intime, outtime, subject_id, hadm_id, first_careunit, last_careunit, los),
    id_icustays = pkl.load(id_icustays_f)
    print('len(id_icustays)', len(id_icustays))

    id_prescription_f = open(output_file + "prescription" + ".pkl", "rb")
    # id_prescription is dictionary, and subject_id is key; value is a list with each item: (starttime, stoptime, hadm_id, drug_type, drug, route)
    id_prescription = pkl.load(id_prescription_f)
    print('len(id_prescription)', len(id_prescription))

    id_culture_f = open(output_file + "microbiologyevents" + ".pkl", "rb")
    # id_culture is dictionary, and subject_id is key; value is a list with each item: (start_date, hadm_id, spec_type_desc, test_name)
    id_culture = pkl.load(id_culture_f)
    print('len(id_culture)', len(id_culture))

    j = 0
    for key, value in id_icustays.items():
        if j % 10 == 0:
            print('j = ', j)
            # break
        j = j + 1

        icustay_id = key
        icu_intime = value[0][0]
        icu_outtime = value[0][1]
        subject_id = value[0][2]
        hadm_id = value[0][3]

        pat_prescription = id_prescription.get(subject_id, [])
        antibiotics_time = extract_antibiotic_time(pat_prescription, icu_intime)  # antibiotics_time

        pat_culture = id_culture.get(subject_id, [])
        culture_time = extract_culture_time(pat_culture, icu_intime)

        flag_infection = 0
        flag_infection_middle = 0
        flag_infection_early = 0
        flag_infection_late = 0
        if antibiotics_time and culture_time:  # have antiobiotics and culture at the same time
            anti_cul_time = antibiotics_time + culture_time
            for i in anti_cul_time:
                if (i - icu_intime).total_seconds() < 0:
                    flag_infection_early = 1
                elif ((i - icu_intime).total_seconds() <= infection_time_threshold * 60 * 60) and (
                        (i - icu_intime).total_seconds() > 0):  # 6 hours [-24,6]:>= -24*60*60;
                    flag_infection_middle = 1
                else:
                    flag_infection_late = 1
                    # print('i = ', i)
                    # print('icu_intime',icu_intime)
                    # print('flag_infection===', flag_infection)
                    # print('(i - icu_intime).seconds', (i - icu_intime).total_seconds())
                    # break
        if flag_infection_early == 0:
            flag_infection = flag_infection_middle

        pat_infection = [icustay_id, flag_infection]
        infec_time_df.loc[len(infec_time_df)] = pat_infection  # save
    infec_time_df.to_csv(output_file + 'infection.csv', index=False)

    print('done.')
    print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    return infec_time_df

def func_read_ori_data(candidate_folder_path):
    # read icustays, inputevents, outputevents, chartevents files;
    print("Beginning Loading Data")
    id_icustays_f = open(candidate_folder_path + "icustays.pkl", "rb")
    id_icustays = pkl.load(
        id_icustays_f)  # id_icustays is dictionary, and icustay_id is key; value is a list [] with each item: such as (intime, outtime, subject_id, hadm_id, first_careunit, last_careunit, los)
    print("Loaded ICU Stays")
    id_inputevents_f = open(candidate_folder_path + "inputevents.pkl", "rb")
    id_inputevents = pkl.load(
        id_inputevents_f)  # id_inputevents is dictionary, and icustay_id is key; value is a list with each item: such as (starttime, endtime, itemid, amount, amountuom, rate, rateuom, patientweight)
    print("Loaded Input Events")
    id_outputevents_f = open(candidate_folder_path + "outputevents.pkl", "rb")
    id_outputevents = pkl.load(
        id_outputevents_f)  # id_outputevents is dictionary, and icustay_id is key; value is a list with each item: such as (start_date, itemid, value, valueuom)
    print("Loaded Output Events")
    # id_chartevents_f = open(candidate_folder_path + "chartevents.pkl", "rb")
    # id_chartevents = pkl.load(
    #     id_chartevents_f)  # id_chartevents is dictionary, and icustay_id is key; value is a list with each item: such as (start_date, itemid, value, valueuom)
    with open(candidate_folder_path + 'chartevents.json', 'r') as f:
        id_chartevents = json.load(f, object_hook=timestamp_hook)
    print("Loaded Chart Events")
    id_procedureevents_f = open(candidate_folder_path + "procedureevents.pkl", "rb")
    id_procedureevents = pkl.load(
        id_procedureevents_f)  # id_procedureevents is dictionary, and icustay_id is key; value is a list with each item: such as (starttime, endtime, itemid, value, valueuom)
    print("Loaded Procedure Events")
    feature_file = candidate_folder_path + 'variable_list_MIMIC_IV.xlsx'
    feature_df = pd.read_excel(feature_file, sheet_name='variable_list')

    return id_icustays, id_inputevents, id_outputevents, id_chartevents, id_procedureevents, feature_df

def func_combine_sofa_infection(input_file, output_file):
    infection_df = pd.read_csv(input_file + 'infection.csv', index_col=False, dtype={'icustay_id': str})
    print('len(infection_df)', len(infection_df))
    sofa_feature_f = open(output_file + "icustays_features_sofa.pkl", "rb")
    sofa_feature = pkl.load(sofa_feature_f)
    print('len(sofa_feature)', len(sofa_feature))

    patients_f = open(output_file + "patients.pkl", "rb")
    patients = pkl.load(patients_f)
    print('len(patient)', len(patients))

    icustays_f = open(output_file + "icustays.pkl", "rb")
    icustays = pkl.load(icustays_f)
    print('len(icustays)', len(icustays))

    services_f = open(output_file + "services.pkl", "rb")
    services = pkl.load(services_f)
    print('len(services)', len(services))

    # create df for save all patients sofa and infection,
    column_name = []
    for i in range(int(TOTAL_DURATION / DURATION)):
        name = 'hour_' + str((i + 1) * DURATION)
        column_name.append(name)
    result_df = pd.DataFrame(columns=['icustay_id'] + column_name)

    j = 0
    for key, values in sofa_feature.items():
        print('j = ', j)
        # if j>20:
        #     break
        j = j + 1
        patid = key
        sofa_df = values[-1]
        sofa_value_list = list(sofa_df.loc['SOFA_score', :])
        # save sofa value list
        result_df.loc[len(result_df)] = [patid] + sofa_value_list
    # merge sofa and infection
    result_df = result_df.merge(infection_df, how='left', on='icustay_id')

    # add age, service to result_df
    result_df['age'] = np.nan
    result_df['service'] = np.nan
    result_df['los'] = np.nan
    for i, row in result_df.iterrows():
        icustay_id = row['icustay_id']
        subject_id = icustays[icustay_id][0][2]
        age = patients[subject_id][0][1]
        icustay_intime = icustays[icustay_id][0][0]
        icustay_outtime = icustays[icustay_id][0][1]
        icustay_los = icustays[icustay_id][0][6]
        p_service_type = np.nan
        p_service_list = services[subject_id]
        flag_server = 0
        for p_service in p_service_list:
            p_service_time = p_service[0]
            if ((p_service_time - icustay_intime).days > -1) and ((p_service_time - icustay_outtime).days < 1):
                if flag_server == 0:
                    p_service_type = p_service[3]
                else:
                    flag_server = 1

        result_df.loc[i, 'age'] = int(age)
        result_df.loc[i, 'los'] = float(icustay_los)
        result_df.loc[i, 'service'] = p_service_type

    # sepsis label
    result_df['sepsis'] = 0
    surgical_service = ['CSURG', 'VSURG', 'TSURG']
    surgical_service = surgical_service + ['DENT','ENT','EYE','NB','NBB','NMED','NSURG','OBS','ORTHO','OMED','PSURG','PSYCH','TRAUM']
    result_df.loc[(result_df['los'] >= 0.25) & (result_df['hour_24'] >= 2) & (result_df['flag_infection'] == 1) & (
                result_df['age'] >= 18) & (~result_df['service'].isin(surgical_service)), 'sepsis'] = 1

    # save all ori_sofa_infection
    result_df.to_csv(output_file + 'ori_sofa_infection.csv', index=False)
    ori_sofa_infection = result_df
    # save sepsis result
    sepsis_sofa = result_df.loc[result_df['sepsis'] == 1]
    sepsis_sofa = sepsis_sofa[['icustay_id'] + column_name]

    # delete outliers , no any change in terms of sofa trajectory
    print('before_sepsis_sofa_trajectory_change:',len(sepsis_sofa))
    sepsis_sofa = sepsis_sofa.set_index('icustay_id')
    sepsis_sofa = sepsis_sofa.loc[sepsis_sofa.std(axis=1) != 0]
    print('after_sepsis_sofa_trajectory_change:', len(sepsis_sofa))

    sepsis_sofa.to_csv(output_file + 'sepsis_sofa.csv', index=True)

    return ori_sofa_infection, sepsis_sofa


if __name__ == '__main__':
    start_time = time.time()

    # input_file = '/Users/xuzhenxing/Documents/MIMIC-IV/'
    # output_file = input_file

    # on server:
    input_file = '/Users/suraj/Documents/CBM/Fei/Trial/Data/MIMIC/mimic-iv-2.2/'
    output_file = input_file

    services = pre_services(input_file, output_file)
    patients = pre_patients(input_file, output_file)
    prescription = pre_prescription(input_file, output_file) # patid: 166758
    culture = pre_culture(input_file, output_file) # patid: 190695
    chartevents = pre_chartevents(input_file,output_file) # did it complete?
    procedureevents = pre_procedureevents(input_file, output_file)
    inputevents = pre_inputevents(input_file, output_file)
    outputevents = pre_outputevents(input_file, output_file)
    icustays = pre_icustays(input_file, output_file)

    infection = identify_infection_time(output_file)
    features = pre_features(input_file, output_file)
    sofa = compute_sofa(input_file, output_file)
    ori_sofa_infection, sepsis_sofa = func_combine_sofa_infection(input_file, output_file)


