#!/bin/env python
import numpy as np
import argparse
import os
import math
import re

parser = argparse.ArgumentParser(description="需要计算SOC、ΔE(st)、T 默认为 298.15 K、fcclasses.out和HuangRhys.dat")
parser.add_argument("NAC",type=float, help="SOC/NAC：cm-1")
parser.add_argument("delta",type=float, help="ΔEst/ΔG：eV")
parser.add_argument("T",nargs='?',default=298.15,type=float, help="温度：K")

args = parser.parse_args()

# 将文件名分别赋值给不同的变量
NAC_cm = args.NAC   # 变量 1 SOC
ΔG = args.delta     # 变量 2 ΔEst
T = args.T          # 变量 3 温度

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

# 计算重组能/λM

arrf = arrfrequencies[0]
arrHR = arrHRindex.astype(np.float32)

λM = np.sum(arrHR*(arrf*29979245800*6.6260696E-34)*219474.6363/2625500*6.02214179E+23*1.98644586E-23/1.602176634E-19)



def marcus_rate(NAC,NAC_cm, λM, ΔG, T):
    """
    计算 Marcus 速率常数
    参数:
    NAC (eV): 非绝热耦合
    λM (eV): 重组能量
    ΔG (eV): 自由能变化
    ω_eff (cm^-1): 有效振动频率
    S (无单位): Huang-Rhys 因子
    T (K): 温度

    返回:
    kT_Marcus (s^-1): Marcus 速率常数
    """
    # 常量
    h = 4.135667696e-15  # eV·s (普朗克常数)
    kb = 8.617333262145e-5  # eV/K (玻尔兹曼常数)
    ħ = h / (2 * np.pi)  # eV·s (约化普朗克常数)
    
    # 计算 Marcus 速率常数
    factor1 = (np.pi / ħ) * (NAC**2)
    factor2 = 1 / np.sqrt(np.pi * λM * kb * T)
    sum_term = 0
    exponent = -((ΔG + λM) ** 2) / (4 * λM * kb * T)
    sum_term =  np.exp(exponent)
    
    kT_Marcus = factor1 * factor2 * sum_term
    
    return kT_Marcus

# 示例参数
NAC = NAC_cm*1.239841E-4  # SOC / eV

print(f'SOC_eV: {NAC}')
print(f"λM_CM-1: {λM}")
print(f"ΔG_eV: {ΔG}")
print(f"K: {T}")

# 计算 Marcus 速率常数
rate_constant = marcus_rate(NAC,NAC_cm, λM, ΔG, T)
print(f"Marcus 速率常数: {rate_constant:.4e} s^-1")
