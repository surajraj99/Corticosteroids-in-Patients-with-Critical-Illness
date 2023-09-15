import json
import pickle
import pandas as pd

def obtain_comorbidity(patient_icd_code_list,each_com_codes):
    # each_com_codes is dictionary
    five_icd9_code = each_com_codes['five_icd9_code']
    four_icd9_code = each_com_codes['four_icd9_code']
    three_icd9_code = each_com_codes['three_icd9_code']
    com_0_or_1 = 0

    for k in patient_icd_code_list:
        # print('k',k)
        i = k.split('.')
        # print('i', i)

        if len(i) == 2:
            i = i[0] + i[1]
        else:
            i = i

        i_5 = i
        i_4 = i[:4]
        i_3 = i[:3]
        # print('i_4',i_4)
        # print('i_3', i_3)
        if i_5 in five_icd9_code:
            com_0_or_1 = 1
        if i_4 in four_icd9_code:
            com_0_or_1 = 1
        if i_3 in three_icd9_code:
            com_0_or_1 = 1

    # print('ss')

    return com_0_or_1

path_diagnose_ICD9_code = '/Users/suraj/Documents/CBM/Fei/Trial/Data/eicu-database-2.0/'
result_path = '/Users/suraj/Documents/CBM/Fei/Trial/Data/'
path_paitent_ID = result_path

object_patient_ID = pd.read_csv(path_paitent_ID+'sepsis.csv',index_col=False)
patientunitstayid_s = list(object_patient_ID['patientunitstayid'])

DIAGNOSES_ICD = pd.read_csv(path_diagnose_ICD9_code+'diagnosis.csv',index_col=False)
print('before_len(DIAGNOSES_ICD)',len(DIAGNOSES_ICD))
DIAGNOSES_ICD = DIAGNOSES_ICD[DIAGNOSES_ICD['icd9code'].notna()]
DIAGNOSES_ICD = DIAGNOSES_ICD[['patientunitstayid','icd9code']]

# DIAGNOSES_ICD.to_csv('/Users/xuzhenxing/Documents/eICU_AKI_Sepsis/testing/DIAGNOSES_ICD_updated.csv',index=False)
# print('after_len(DIAGNOSES_ICD)',len(DIAGNOSES_ICD))

# Get the indexes which are repetative with the split
df = DIAGNOSES_ICD
df['icd9code'] = df['icd9code'].str.split(',')
df = df.explode('icd9code')
df['icd9code'] = df['icd9code'].str.strip()
DIAGNOSES_ICD = df
# DIAGNOSES_ICD.to_csv('/Users/xuzhenxing/Documents/eICU_AKI_Sepsis/testing/DIAGNOSES_ICD_updated_2.csv',index=False)
# print('after_len(DIAGNOSES_ICD)',len(DIAGNOSES_ICD))

# DIAGNOSES_ICD['icd9code'] = DIAGNOSES_ICD['icd9code'].str.split(',',expand=False)

comorbidity_name = ['CHF','VALVE','PULMCIRC','PERIVASC','HTN','PARA','NEURO','CHRNLUNG','DM','DMCX',
                    'HYPOTHY','RENLFAIL','LIVER','ULCER','AIDS','LYMPH','METS','TUMOR','ARTH','COAG',
                    'OBESE','WGHTLOSS','LYTES','BLDLOSS','ANEMDEF','ALCOHOL','DRUG','PSYCH','DEPRESS'
                    ]

comorbidity_ICD9_code = {

    'CHF': {
        'five_icd9_code':['39891','40201','40211','40291','40401','40403','40411','40413','40491','40493'],
        'four_icd9_code': ['4254','4255','4257','4258','4259'],
        'three_icd9_code': ['428'],
        'Full_name':'Congestive heart failure'
    },

    'VALVE': {
        'five_icd9_code': [],
        'four_icd9_code': ['0932','7463','7464','7465','7466','V422','V433'],
        'three_icd9_code': ['394','395','396','397','424'],
        'Full_name': 'Valvular disease'
    },

    'PULMCIRC': {
        'five_icd9_code': [],
        'four_icd9_code': ['4150','4151','4170','4178','4179'],
        'three_icd9_code': ['416'],
        'Full_name': 'Pulmonary circulation disorder'
    },

    'PERIVASC': {
        'five_icd9_code': [],
        'four_icd9_code': ['0930','4373','4431','4432','4438','4439','4471','5571','5579','V434'],
        'three_icd9_code': ['440','441'],
        'Full_name': 'Peripheral vascular disorder'
    },

    'HTN': {
        'five_icd9_code': [],
        'four_icd9_code': [],
        'three_icd9_code': ['401','402','403','404','405'],
        'Full_name': 'Hypertension'
    },

    'PARA': {
        'five_icd9_code': [],
        'four_icd9_code': ['3341','3440','3441','3442','3443','3444','3445','3446','3449'],
        'three_icd9_code': ['342','343'],
        'Full_name': 'Paralysis'
    },

    'NEURO': {
        'five_icd9_code': ['33392'],
        'four_icd9_code': ['3319','3320','3321','3334','3335','3362','3481','3483','7803','7843'],
        'three_icd9_code': ['334','335','340','341','345'],
        'Full_name': 'Other neurological'
    },

    'CHRNLUNG': {
        'five_icd9_code': [],
        'four_icd9_code': ['4168','4169','5064','5081','5088'],
        'three_icd9_code': ['490','491','492','493','494','495','496','500','501','502','503','504','505'],
        'Full_name': 'Chronic pulmonary disease'
    },

    'DM': {
        'five_icd9_code': [],
        'four_icd9_code': ['2500','2501','2502','2503'],
        'three_icd9_code': [],
        'Full_name': 'Diabetes w/o chronic complications'
    },

    'DMCX': {
        'five_icd9_code': [],
        'four_icd9_code': ['2504','2505','2506','2507','2508','2509'],
        'three_icd9_code': [],
        'Full_name': 'Diabetes w/ chronic complications'
    },

    'HYPOTHY': {
        'five_icd9_code': [],
        'four_icd9_code': ['2409','2461','2468'],
        'three_icd9_code': ['243','244'],
        'Full_name': 'Hypothyroidism'
    },

    'RENLFAIL': {
        'five_icd9_code': ['40301','40311','40391','40402','40403','40412','40413','40492','40493'],
        'four_icd9_code': ['5880','V420','V451'],
        'three_icd9_code': ['585','586','V56'],
        'Full_name': 'Renal failure'
    },

    'LIVER': {
        'five_icd9_code': ['07022','07023','07032','07033','07044','07054'],
        'four_icd9_code': ['0706','0709','4560','4561','4562','5722','5723','5724','5728','5733','5734','5738','5739','V427'],
        'three_icd9_code': ['570','571'],
        'Full_name': 'Liver disease'
    },

    'ULCER': {
        'five_icd9_code': [],
        'four_icd9_code': ['5317','5319','5327','5329','5337','5339','5347','5349'],
        'three_icd9_code': [],
        'Full_name': 'Chronic Peptic ulcer disease'
    },

    'AIDS': {
        'five_icd9_code': [],
        'four_icd9_code': [],
        'three_icd9_code': ['042','043','044'],
        'Full_name': 'HIV_AIDS'
    },

    'LYMPH': {
        'five_icd9_code': [],
        'four_icd9_code': ['2030','2386'],
        'three_icd9_code': ['200','201','202'],
        'Full_name': 'Lymphoma'
    },

    'METS': {
        'five_icd9_code': [],
        'four_icd9_code': [],
        'three_icd9_code': ['196','197','198','199'],
        'Full_name': 'Metastatic cancer'
    },

    'TUMOR': {
        'five_icd9_code': [],
        'four_icd9_code': [],
        'three_icd9_code': ['140','141','142','143','144','145','146','147','148','149','150','151','152'
    ,'153','154','155','156','157','158','159','160','161','162','163','164','165'
    ,'166','167','168','169','170','171','172','174','175','176','177','178','179'
    ,'180','181','182','183','184','185','186','187','188','189','190','191','192'
    ,'193','194','195'],
        'Full_name': 'Solid tumor without metastasis'
    },

    'ARTH': {
        'five_icd9_code': ['72889','72930'],
        'four_icd9_code': ['7010','7100','7101','7102','7103','7104','7108','7109','7112','7193','7285'],
        'three_icd9_code': ['446','714','720','725'],
        'Full_name': 'Rheumatoid arthritis/collagen vascular diseases'
    },

    'COAG': {
        'five_icd9_code': [],
        'four_icd9_code': ['2871','2873','2874','2875'],
        'three_icd9_code': ['286'],
        'Full_name': 'Coagulation deficiency'
    },

    'OBESE': {
        'five_icd9_code': [],
        'four_icd9_code': ['2780'],
        'three_icd9_code': [],
        'Full_name': 'Obesity'
    },

    'WGHTLOSS': {
        'five_icd9_code': [],
        'four_icd9_code': ['7832','7994'],
        'three_icd9_code': ['260','261','262','263'],
        'Full_name': 'Weight loss'
    },

    'LYTES': {
        'five_icd9_code': [],
        'four_icd9_code': ['2536'],
        'three_icd9_code': ['276'],
        'Full_name': 'Fluid and electrolyte disorders'
    },

    'BLDLOSS': {
        'five_icd9_code': [],
        'four_icd9_code': ['2800'],
        'three_icd9_code': [],
        'Full_name': 'Blood loss anemia'
    },

    'ANEMDEF': {
        'five_icd9_code': [],
        'four_icd9_code': ['2801','2808','2809'],
        'three_icd9_code': ['281'],
        'Full_name': 'Deficiency anemias'
    },

    'ALCOHOL': {
        'five_icd9_code': [],
        'four_icd9_code': ['2652','2911','2912','2913','2915','2918','2919','3030','3039','3050','3575','4255','5353','5710','5711','5712','5713','V113'],
        'three_icd9_code': ['980'],
        'Full_name': 'Alcohol abuse'
    },

    'DRUG': {
        'five_icd9_code': ['V6542'],
        'four_icd9_code': ['3052','3053','3054','3055','3056','3057','3058','3059'],
        'three_icd9_code': ['292','304'],
        'Full_name': 'Drug abuse'
    },

    'PSYCH': {
        'five_icd9_code': ['29604','29614','29644','29654'],
        'four_icd9_code': ['2938'],
        'three_icd9_code': ['295','297','298'],
        'Full_name': 'Psychoses'
    },

    'DEPRESS': {
        'five_icd9_code': [],
        'four_icd9_code': ['2962','2963','2965','3004'],
        'three_icd9_code': ['309','311'],
        'Full_name': 'Depression'
    },

}

all_comorbidity_icd9_code_df = pd.DataFrame(columns=['patientunitstayid']+comorbidity_name)

j = 0

for p_id in patientunitstayid_s:
    print('j',j)
    # create dataframe for each patient with columns ['ICU_stay_id','Subject_id']+comorbidity_name
    person_comorbidity_icd9_code_df = pd.DataFrame(columns=['patientunitstayid']+comorbidity_name)
    person_result_dic = {}

    # print('p_id',p_id)
    person_result_dic['patientunitstayid'] = int(p_id)

    patient_icd_code_list = list(DIAGNOSES_ICD.loc[DIAGNOSES_ICD['patientunitstayid'] == int(p_id), 'icd9code'])

    # print('patient_icd_code_list',patient_icd_code_list)

    for com in comorbidity_ICD9_code:
        com_name = com
        each_com_codes = comorbidity_ICD9_code[com_name]
        comorbidity_result = obtain_comorbidity(patient_icd_code_list,each_com_codes)
        # print('comorbidity_result',comorbidity_result)
        person_result_dic[com_name] = comorbidity_result

    person_comorbidity_icd9_code_df = person_comorbidity_icd9_code_df.append(person_result_dic, ignore_index=True, sort=False)

    all_comorbidity_icd9_code_df = all_comorbidity_icd9_code_df.append(person_comorbidity_icd9_code_df)

    j = j + 1
    # if j >10:
    #     break

com_df = all_comorbidity_icd9_code_df
# compute the comorbidity score
com_df['comorbidity_score'] = 9*com_df['CHF'] + 0*com_df['VALVE']+5*com_df['PULMCIRC'] + 4*com_df['PERIVASC']+\
                              (-2)*com_df['HTN']+4*com_df['PARA']+5*com_df['NEURO']+3*com_df['CHRNLUNG']+\
                              0*com_df['DM']+1*com_df['DMCX']+0*com_df['HYPOTHY']+7*com_df['RENLFAIL']+7*com_df['LIVER']+\
                              0*com_df['ULCER']+0*com_df['AIDS']+8*com_df['LYMPH']+17*com_df['METS']+10*com_df['TUMOR']+\
                              0*com_df['ARTH']+12*com_df['COAG']+(-5)*com_df['OBESE']+10*com_df['WGHTLOSS']+11*com_df['LYTES']+\
                              (-3)*com_df['BLDLOSS']+0*com_df['ANEMDEF']+0*com_df['ALCOHOL']+(-11)*com_df['DRUG']+(-6)*com_df['PSYCH']+(-5)*com_df['DEPRESS']

com_df.to_csv(result_path+'patient_comorbidity_score_df.csv',index=False)

print('Done!')


