#!/bin/env python
import numpy as np
import argparse
import os
import re

parser = argparse.ArgumentParser(description="需要计算SOC、重组能、ΔE(st)、freq、H-R-index/S")
parser.add_argument("NAC",type=str, help="SOC/NAC：cm-1")
parser.add_argument("delta",type=str, help="ΔEst/ΔG：eV")
parser.add_argument("T",default=298.15,type=str, help="温度：K")

args = parser.parse_args()

# 将文件名分别赋值给不同的变量
NAC_cm = args.delta
ΔG = args.delta
T = args.delta

# 读取文件

HRindexfile = "HuangRhys.dat"
Fccfile = "fcclasses.out"

def ReadHRindex(HRindexfile):
    print(f"reading ReadHRindex:{HRindexfile}")
    HRindex = []
    slag = 0
    # snumber = 0
    # mnumber = 0
    # RH = open(path + filename, 'r', encoding='utf8')
    RH = open(HRindexfile, 'r', encoding='utf8')
    for line in RH:
        if 'Huang-Rhys' in line:
            slag = 1
            continue
        if slag == 1 and len(line.strip().split()) ==2:
            # snumber = int(line.strip().split()[0])
            # print(snumber)
            # #line_save = line.split(' ')
            # #line_save = [string.replace('\n', '') for string in line_save]
            # # if snumber > mnumber:
            #     print(line.strip().split())
            HRindex.append(line.strip().split()[-1])
            #     mnumber = snumber
            # else:
            #     flag = 0 
    RH.close()
    return HRindex
if not os.path.isfile(HRindexfile):
    print(f"File {HRindexfile} does not exist.")
else:
    arrHRindex = np.array(ReadHRindex(HRindexfile))

def Readfreq(Fccfile):
    print(f"reading Readfcc.out:{Fccfile}")
    pattern = r'------------------------\s+FREQUENCIES \(cm-1\)\s+------------------------\s+((?:\s*\d+\s+\d+\.\d+\s*)+)\s+------------------------'
    FR = open(Fccfile, 'r', encoding='utf8')
    content = FR.read()
    matches = re.findall(pattern, content)
    all_frequencies = []
    for match in matches:
        frequencies = re.findall(r'\d+\s+(\d+\.\d+)', match)
        all_frequencies.append(frequencies)
    return all_frequencies

arrfrequencies = np.array(Readfreq(Fccfile)).astype(np.float32)

# 计算重组能/λ

Harrf = arrfrequencies[0,arrfrequencies[0]>=400]
HarrHR = arrHRindex[arrfrequencies[0]>=400].astype(np.float32)
Larrf = arrfrequencies[0,arrfrequencies[0]<400] 
LarrHR = arrHRindex[arrfrequencies[0]<400].astype(np.float32) 

λM_sum = np.sum(LarrHR*(Larrf*29979245800*6.6260696E-34)*219474.6363/2625500*6.02214179E+23*1.98644586E-23/1.602176634E-19)

# 计算 w_eff

w_eff = np.sum(Harrf*HarrHR)/np.sum(HarrHR)

# 计算 H-R-Index

S = np.sum(HarrHR)
