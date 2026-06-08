#!/usr/bin/env python
import re
from collections import defaultdict
import pandas as pd
import numpy as np
orbcoef = defaultdict(list)
pattern = r'(\d+\.)\s*([A-Z]+\s*\d+)\s*(\([Spxyzd123456]+\))\s*([\s\-0-9.]+)'
gaulog = open("33TS.log","r")
readnao=10
atomindex=[26,27]
orbrange=[-6,6]
for l in gaulog:
    if "MOs in the NAO basis:" in l:
        readnao=3
        continue
    if readnao==3 and "NAO" in l:
        readnao = 2 
        continue
    if readnao==2 and "------" in l:
        readnao = 1
        continue
    if readnao == 1 and l.strip()=="":
        readnao= 3
    if "NATURAL BOND ORBITAL" in l:
        readnao= 10
    if readnao == 1:
        match= re.search(pattern,l)
        if match:
            orbcoef[(match.group(1),match.group(2),match.group(3))] += [float(i) for i in match.group(4).strip().split()]
        else:
            print(l)
gaulog.close()
df = pd.DataFrame.from_dict(orbcoef, orient='index')
df.index = pd.MultiIndex.from_tuples(df.index)
df_sum = df.groupby(level=1).sum()
print(df_sum)
gaulog = open("33TS.log","r")
ao, av, bo, bv = [], [], [], []
for l in gaulog:
    l = l.rstrip()
    if 'Alpha  occ. eigenvalues' in l:
        ao = ao + [float(l[i:i+10].strip()) for i in range(28, len(l), 10)]
    if 'Beta  occ. eigenvalues' in l:
        bo = bo + [float(l[i:i+10].strip()) for i in range(28, len(l), 10)]
    if 'Alpha virt. eigenvalues' in l:
        av = av + [float(l[i:i+10].strip()) for i in range(28, len(l), 10)]
    if 'Beta virt. eigenvalues' in l:
        bv = bv + [float(l[i:i+10].strip()) for i in range(28, len(l), 10)]
orb = ao+av
fermi_level = (ao[-1]+av[0])/2
print(f"fermi level is {fermi_level}")
hidx=len(ao)
orb_range=np.arange(orbrange[0],orbrange[1])+hidx #h-2 到 lumo+2
print(f"HOMO is {hidx}")
def extract_number(s):
    match = re.search(r'\d+', s)
    return int(match.group()) if match else None

# 应用这个函数到索引上，并创建一个布尔索引
bool_idx = [i in atomindex for i in df_sum.index.map(extract_number)]

# 使用布尔索引来选取行
selected_rows = df_sum[bool_idx]
row_multiply = selected_rows.iloc[0] * selected_rows.iloc[1]
result = pd.Series([value / (fermi_level-orb[index]) for index, value in enumerate(row_multiply)])
df_transposed = selected_rows.transpose()
df_transposed.index = range(1, len(df_transposed) + 1)
df_transposed['orb'] = orb
df_transposed['c1c2'] = df_transposed.iloc[:,0]*df_transposed.iloc[:,1]
df_transposed['delta_orb'] = fermi_level - np.array(orb)
df_transposed['c1c2/dE'] = df_transposed['c1c2']/df_transposed['delta_orb'] 
print(df_transposed.iloc[orb_range,:])
#print(result)
#print(sum(result))
#print(max(result))
#print(result[orb_range])
#print(row_multiply[orb_range])

