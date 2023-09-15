# fill missing based on the result from feature_bucketing_addressing_data

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

start_time = time.time()
# read feature_bucket_with_missing
with open(result_path + 'feature_bucket_with_missing.json', 'r') as load_f:
    feature_bucket = json.load(load_f)

count = 0

feature_bucket_filled_missing = {}

for i in feature_bucket:
    print('count',count)

    p_ID = i
    feature = feature_bucket[p_ID]
    feature_df = pd.DataFrame.from_dict(feature)

    # feature_df.to_csv(result_path+'testing/person_feature_'+p_ID+'_not_filen.csv',index=False)

    # set FeatureName as index and fill missing value
    feature_df = feature_df.set_index('FeatureName')

    # impute missing value with ffill and backfill
    feature_df = feature_df.fillna(method='ffill', axis=1)
    feature_df = feature_df.fillna(method='backfill', axis=1)
    # filling 0 for variables if there is not any values
    feature_df = feature_df.fillna(0)

    feature_df = feature_df.reset_index(drop=False)
    feature_bucket_filled_missing[p_ID] = feature_df

    # feature_df.to_csv(result_path+'testing/person_feature_'+p_ID+'_fill.csv',index=True)
    # if count>1000:
    #     break
    count = count + 1


# save results filled missing file
with open(result_path + "feature_bucket_filled_missing.json", "w") as outfile:
    json.dump(feature_bucket_filled_missing, outfile, cls=JSONEncoder, indent=4)

print('done.')
print('Time used:', time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
