import os
import re
import pandas as pd
from tqdm import tqdm


df = pd.read_csv('N:/results/rlt_combine_entraplib_G3/rawdata_list.tsv')
rawdata_ls = df.FileName.tolist()

file_ls = [rd.split('\\')[-1].strip('.d') for rd in rawdata_ls]
datetime_ls = []
for rd in tqdm(rawdata_ls):
    search_ls = [re.search(r'^([0-9\-_]+).+separation$', f).group(1).strip('_') for f in os.listdir(rd) if f.endswith('separation')]
    if len(search_ls)!= 1:
        print(rd)
        break
    else:
        datetime_ls.append(search_ls[0])



dfout = pd.DataFrame({'RawDataPath': rawdata_ls,
              'FileName': file_ls,
              'DateTime': datetime_ls})
dfout.to_excel('N:/results/rlt_combine_entraplib_G3/rawdata_datetime.xlsx', index=None)



