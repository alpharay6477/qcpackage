#!/usr/bin/env python
"""
orcamonitor.py - 监控ORCA计算输出文件的工具

此脚本用于分析ORCA量子化学计算的输出文件，提取优化循环信息、
虚频模式、收敛情况等关键数据，并以友好的格式展示。

版本: 1.0
改进: 添加类型注解、文档字符串、异常处理等
"""

import re
import sys
import os
import numpy as np
import argparse
from typing import List, Tuple, Dict, Any, Optional, Union

# Regular expression for matching the:
# Redundant Internal Coordinates (Angstroem and degrees)
# Definition  Value    dE/dq     Step     New-Value  [comp.(TS mode)]
pattern = r'^\s*(\d+)\.\s+([BADL]\(([A-Za-z]+\s+\d+,?)+(?:\s+\d+)?\))\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)(?:\s+(-?\d+\.\d+))?'

# extract the atom and the number: 'B(C   1,C   0)' --> [('C', '1'), ('C', '0')]
pattern2 = r'([A-Za-z])+\s+(\d+),?'

# Regular expression for matching the optimization cycle beginning
pattern3 = r'^\s*\*\s*GEOMETRY OPTIMIZATION CYCLE\s*(\d+)\s*\*'

# Regular expression for matching the imaginary mode originated from Hessian Calculation
pattern4 = r'^\s+(\d+):\s+(-?\d+\.\d+)\s+cm\*\*-1\s+(\*\*\*imaginary mode\*\*\*)?'

# Regular expression for matching charge and multiplicity
pattern5 = r'^\|\s+\d+>\s+\*\s*(?:XYZ|XYZFILE)\s+(\d+)\s+(\d+)'

# Regular expression for matching the optimization convergence
pattern6_dict = {
    "E_change": r'^\s+Energy\schange\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)\s+(?:YES|NO)',
    "RMS_grad": r'^\s+RMS\sgradient\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)\s+(?:YES|NO)',
    "MAX_grad": r'^\s+MAX gradient\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)\s+(?:YES|NO)',
    "RMS_step": r'^\s+RMS\sstep\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)\s+(?:YES|NO)',
    "MAX_step": r'^\s+MAX\sstep\s+(-?\d+\.\d+)\s+(?:-?\d+\.\d+)\s+(?:YES|NO)',
}

# Regular expression for matching the num of atoms
# Number of atoms                             ...     26
pattern7 = r'^Number\sof\satoms\s+\.\.\.\s+(\d+)'

# Extract minimum eigen value of the Hessian matrix
pattern8 = r'Lowest eigenvalues of (?:augmented|the) Hessian:\s+(-?\d*\.\d+)\s+(-?\d*\.\d+)\s+(-?\d*\.\d+)'

opt_flag_txt = "GEOMETRY OPTIMIZATION CYCLE"
ts_flag_txt = "TS mode is"
freq_flag_txt = "THERMOCHEMISTRY AT"


def get_read_internal_coords(lines: List[str], ts_flag: bool = True) -> Tuple[List[str], List[str], List[Tuple[int, int]], np.ndarray]:
    """
    读取内部坐标信息
    
    Args:
        lines: ORCA输出文件的行列表
        ts_flag: 是否为过渡态计算
        
    Returns:
        Tuple[List[str], List[str], List[Tuple[int, int]], np.ndarray]:
            - 内部坐标名称列表
            - 列名列表
            - 每行内部坐标的开始和结束行号
            - 内部坐标数据的3D数组
    """
    all_matches: List[List[Any]] = []
    for line_num, line in enumerate(lines, start=1):
        # Search for matches in the current line
        match = re.search(pattern, line)
        if match:
            # If a match is found, print the match and the line number
            match_lis = list(match.groups())
            all_matches.append([line_num, match_lis])
    
    def progress_match_item(match_item: List[str]) -> List[Union[int, str, float]]:
        """处理匹配项"""
        internal_item = match_item[1]  # string
        assert internal_item[0] in ('B', 'A', 'D', 'L'), "Internal coordinate type not recognized"
        
        temp_lis = re.findall(pattern2, internal_item)
        temp_lis = [internal_item[0]] + [item[1] for item in temp_lis]
        new_internal_item = "-".join(temp_lis)
        
        res = [int(match_item[0]), new_internal_item, float(match_item[3]), 
               float(match_item[4]), float(match_item[5]), float(match_item[6])]
        
        if ts_flag:
            if match_item[7]:
                res.append(float(match_item[7]))
            else:
                res.append(0.)
        return res
    
    all_matches = list(map(lambda a: [a[0], progress_match_item(a[1])], all_matches))
    all_matches_line_idx = [item[0] for item in all_matches]
    all_matches_cont = [item[1] for item in all_matches]
    
    internal_idx_lis = [item[0] for item in all_matches_cont]
    internal_coordName_lis = [item[1] for item in all_matches_cont]
    
    internal_size = max(internal_idx_lis)
    internal_coord_list = internal_coordName_lis[:internal_size]
    
    number_mat_lis = [item[2:] for item in all_matches_cont]
    number_mat_arr = np.array(number_mat_lis, dtype=float)
    
    mat_columns = ['Value', 'dE/dq', 'Step', 'New-Value'] if not ts_flag else \
                  ['Value', 'dE/dq', 'Step', 'New-Value', 'comp.(TS mode)']
    number_mat_arr = number_mat_arr.reshape(-1, internal_size, len(mat_columns))
    
    line_idx_begin_end = [(all_matches_line_idx[i], all_matches_line_idx[i+internal_size-1]) 
                          for i in range(0, len(all_matches_line_idx), internal_size)]
    
    return internal_coord_list, mat_columns, line_idx_begin_end, number_mat_arr


def get_ts_mode(lines: List[str], n_modes: int = 2) -> List[List[Union[int, str]]]:
    """
    获取过渡态模式
    
    Args:
        lines: ORCA输出文件的行列表
        n_modes: 要获取的模式数量
        
    Returns:
        List[List[Union[int, str]]]: 过渡态模式信息列表
    """
    full_text = ''.join(lines)
    
    if opt_flag_txt in full_text:
        if ts_flag_txt in full_text:
            ts_flag = True
        else:
            ts_flag = False
    else:
        raise ValueError("Not Optimization Output")
    
    internal_coord_list, mat_columns, line_idx_begin_end, number_mat_arr = \
        get_read_internal_coords(lines, ts_flag=ts_flag)
    assert number_mat_arr.shape == (len(line_idx_begin_end), len(internal_coord_list), len(mat_columns)), "Shape mismatch"
    
    ts_mode_col = number_mat_arr[:, :, -1]
    
    max_n_ts_mode_idx = np.argsort(ts_mode_col, axis=1)[:, -n_modes:][:, ::-1]
    
    max_ts_mode = [", ".join([internal_coord_list[i] for i in sublis]) for sublis in max_n_ts_mode_idx]
    
    result = [[line_idx_begin_end[i][0], max_ts_mode[i], "redundant_internal"] for i in range(len(line_idx_begin_end))]
    return result

def get_opt_cycle_begin(lines: List[str]) -> List[List[Union[int, str]]]:
    """
    获取优化循环开始的行号
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        List[List[Union[int, str]]]: 优化循环开始信息列表
    """
    all_matches: List[List[Union[int, str]]] = []
    for line_num, line in enumerate(lines, start=1):
        # Search for matches in the current line
        match = re.search(pattern3, line)
        if match:
            all_matches.append([line_num, int(match.groups()[0]), "opt_cycle_begin"])
            
    return all_matches

def get_cartesian_coords_begin(lines: List[str]) -> List[int]:
    """
    获取笛卡尔坐标开始的行号
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        List[int]: 笛卡尔坐标开始行号列表
    """
    all_matches: List[int] = []
    for line_num, line in enumerate(lines, start=1):
        # Search for matches in the current line
        match = re.search(r'^CARTESIAN COORDINATES \(ANGSTROEM\)', line)
        if match:
            all_matches.append(line_num + 1)
            
    return all_matches
    
def get_imaginary_freq(lines: List[str]) -> Tuple[List[Tuple[int, int]], np.ndarray, np.ndarray]:
    """
    获取虚频信息
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        Tuple[List[Tuple[int, int]], np.ndarray, np.ndarray]:
            - 每行频率信息的开始和结束行号
            - 频率值数组
            - 虚频标志数组
    """
    all_matches: List[List[Any]] = []
    for line_num, line in enumerate(lines, start=1):
        # Search for matches in the current line
        match = re.search(pattern4, line)
        if match:
            match_lis = list(match.groups())
            all_matches.append([line_num, match_lis])
            
    
    def progress_match_item(match_item: List[str]) -> List[Union[int, float, bool]]:
        """处理匹配项"""
        frequency = float(match_item[1])
        if match_item[2]:
            if match_item[2] == "***imaginary mode***":
                flag = True
            else:
                raise ValueError(f"Unknown imaginary mode information:{match_item[2]}")
        else:
            flag = False  # 0 for normal mode
            
        res = [int(match_item[0]), frequency, flag]
        return res
    
    all_matches = list(map(lambda a: [a[0], progress_match_item(a[1])], all_matches))
    all_matches_line_idx = [item[0] for item in all_matches]
    all_matches_cont = [item[1] for item in all_matches]
    
    mode_idx_lis = [item[0] for item in all_matches_cont]
    
    mode_size = max(mode_idx_lis) + 1
    
    all_matches_freq_arr = np.array([item[1] for item in all_matches_cont], dtype=float)
    all_matches_freq_arr = all_matches_freq_arr.reshape(-1, mode_size)
    
    all_matches_flag_arr = np.array([item[2] for item in all_matches_cont], dtype=bool)
    all_matches_flag_arr = all_matches_flag_arr.reshape(-1, mode_size)
    
    line_idx_begin_end = [(all_matches_line_idx[i], all_matches_line_idx[i+mode_size-1]) 
                          for i in range(0, len(all_matches_line_idx), mode_size)]
    return line_idx_begin_end, all_matches_freq_arr, all_matches_flag_arr

def get_imaginary_mode(lines: List[str]) -> List[List[Union[int, str]]]:
    """
    获取虚频模式
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        List[List[Union[int, str]]]: 虚频模式信息列表
    """
    line_idx_begin_end, all_matches_freq_arr, all_matches_flag_arr = \
        get_imaginary_freq(lines)
    
    imaginary_freq = [all_matches_freq_arr[i][all_matches_flag_arr[i]] 
                     for i in range(all_matches_freq_arr.shape[0])]
    imaginary_freq = list(map(str, imaginary_freq))
    
    result = [[line_idx_begin_end[i][0], imaginary_freq[i], "imaginary_mode"] 
              for i in range(len(line_idx_begin_end))]
    return result

def get_opt_convergence(lines: List[str], get_ratio: bool = True) -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """
    获取优化收敛信息
    
    Args:
        lines: ORCA输出文件的行列表
        get_ratio: 是否返回收敛阈值的比值
        
    Returns:
        Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray]]:
            - 收敛信息字典
            - 收敛阈值比值字典
    """
    pattern6_match_res = {
        key: [] for key in pattern6_dict
    }
    for line_num, line in enumerate(lines, start=1):
        # Search for matches in the current line
        match = re.search(pattern6_dict["RMS_grad"], line)
        if match:
            pattern6_match_res["RMS_grad"].append(float(match.groups()[0]))
            match_E = re.search(pattern6_dict["E_change"], lines[line_num-2])
            pattern6_match_res["E_change"].append(float(match_E.groups()[0]) if match_E else np.nan)
            match_MAX_grad = re.search(pattern6_dict["MAX_grad"], lines[line_num])
            pattern6_match_res["MAX_grad"].append(float(match_MAX_grad.groups()[0]))
            match_RMS_step = re.search(pattern6_dict["RMS_step"], lines[line_num+1])
            pattern6_match_res["RMS_step"].append(float(match_RMS_step.groups()[0]))
            match_MAX_step = re.search(pattern6_dict["MAX_step"], lines[line_num+2])
            pattern6_match_res["MAX_step"].append(float(match_MAX_step.groups()[0]))
    
    opt_conver_normal = {
        "E_change": 0.0000050000,
        "RMS_grad": 0.0001000000,
        "MAX_grad": 0.0003000000,
        "RMS_step": 0.0020000000,
        "MAX_step": 0.0040000000,
    }

    opt_conver_ratio = {
        key: np.array(pattern6_match_res[key]) / value
            for key, value in opt_conver_normal.items()
    }
    
    if get_ratio:
        return opt_conver_ratio, opt_conver_ratio
    else:
        return pattern6_match_res, opt_conver_ratio


def output_pop(lines: List[str], output_filename: str = '', output_ratio: bool = True, 
               format_out: List[int] = [3, 6, 6, 6, 6, 6, 22, 24, 17, 5]) -> str:
    """
    输出优化信息
    
    Args:
        lines: ORCA输出文件的行列表
        output_filename: 输出文件名
        output_ratio: 是否输出收敛阈值的比值
        format_out: 输出格式
        
    Returns:
        str: 输出文本
    """
    opt_cycle_begin = get_opt_cycle_begin(lines)
    opt_conver, opt_conver_rat = get_opt_convergence(lines, get_ratio=output_ratio)
    min_eigen = get_hessian_min_eigen(lines)
    job_type = get_task_type(lines)
    del_col = []
    
    if job_type != "TS":
        del_col.append(6)  # delete the "im_mode"
        ts_mode = []
    else:
        ts_mode = get_ts_mode(lines)  # not this if it's not ts job
    
    if freq_flag_txt in "".join(lines):
        imaginary_mode = get_imaginary_mode(lines)  # not this if we do not have frequency calculation during optimization
    else:
        del_col.append(7)  # delete the "im_freq"
        imaginary_mode = []

    def get_res_cache(res_cache: List[str]) -> List[str]:
        """处理结果缓存"""
        temp = res_cache[1:]
        if 7 in del_col:
            temp.pop(1)
        if 6 in del_col:
            temp.pop(0)
        return temp

    temp_list = opt_cycle_begin + imaginary_mode + ts_mode
    temp_list = sorted(temp_list, key=lambda x: x[0])[::-1]
    
    result = [["n", "dE", "RMS_g", "MAX_g", "RMS_s", "MAX_s", "im_mode", "im_freq", "EigenX100", "MonConv"]]
    
    # delete the columns in del_col
    format_out = [format_out[i] for i in range(len(format_out)) if i not in del_col]
    result[0] = [result[0][i] for i in range(len(result[0])) if i not in del_col]
    
    res_cache = ["", "", ""]
    
    converge_judge = [False] * len(opt_cycle_begin)
    
    for count_cycle, idx in enumerate(range(len(opt_cycle_begin)), start=1):
        thresh = [1, 1, 1, 1, 1] 
        if idx >= len(opt_conver_rat['E_change']):
            converge_judge[idx] = False
        elif (int(abs(opt_conver_rat['E_change'][idx]) < thresh[0]) + 
              int(opt_conver_rat['RMS_grad'][idx] < thresh[1]) + 
              int(opt_conver_rat['MAX_grad'][idx] < thresh[2]) + 
              int(opt_conver_rat['RMS_step'][idx] < thresh[3]) + 
              int(opt_conver_rat['MAX_step'][idx] < thresh[4])) >= 4:
            converge_judge[idx] = True
        elif (int(abs(opt_conver_rat['E_change'][idx]) < thresh[0]) + 
              int(opt_conver_rat['RMS_grad'][idx] < thresh[1]/2) + 
              int(opt_conver_rat['MAX_grad'][idx] < thresh[2]/2)) == 3:
            converge_judge[idx] = True

    if not output_ratio:
        for key in opt_conver:
            opt_conver[key] = list(map(lambda a: a * (10**4), opt_conver[key]))

    while len(temp_list) > 0:
        temp = temp_list.pop()
        if temp[2] == "opt_cycle_begin":
            if res_cache[0] != "":
                auxi_idx = int(res_cache[0]) - 1
                
                auxil_res = [f"{opt_conver['E_change'][auxi_idx]:.1f}", 
                             f"{opt_conver['RMS_grad'][auxi_idx]:.1f}", 
                             f"{opt_conver['MAX_grad'][auxi_idx]:.1f}", 
                             f"{opt_conver['RMS_step'][auxi_idx]:.1f}", 
                             f"{opt_conver['MAX_step'][auxi_idx]:.1f}"]

                if res_cache[2] == "":
                    converge_judge[auxi_idx] = False
                                    
                res_cache = [res_cache[0]] + auxil_res + get_res_cache(res_cache) + \
                    [min_eigen[auxi_idx]] + ["YES" if converge_judge[auxi_idx] else "NO"]

                result.append(res_cache)
                res_cache = ["", "", ""]
            else:
                pass
            res_cache[0] = str(temp[1])
        elif temp[2] == "redundant_internal":
            res_cache[1] = temp[1] if res_cache[1] == "" else res_cache[1] + " | " + temp[1]
        elif temp[2] == "imaginary_mode":
            res_cache[2] = temp[1] if res_cache[2] == "" else res_cache[2] + " | " + temp[1]
            auxi_idx = int(res_cache[0]) - 1
            
            task_type = get_task_type(lines)
            judgement_ = 1 if task_type == "TS" else 0
            if temp[1].count("-") != judgement_:  # more than 1 imaginary mode
                converge_judge[auxi_idx] = False

    auxi_idx = int(res_cache[0]) - 1
    if res_cache[2] == "":
        converge_judge[auxi_idx] = False
    if auxi_idx >= len(opt_conver['E_change']):
        auxil_res = ["None"] * 5
        res_cache = [res_cache[0]] + auxil_res + get_res_cache(res_cache) + \
                    ["None"] + ["NO"]
    else:
        auxil_res = [f"{opt_conver['E_change'][auxi_idx]:.1f}", 
                     f"{opt_conver['RMS_grad'][auxi_idx]:.1f}", 
                     f"{opt_conver['MAX_grad'][auxi_idx]:.1f}", 
                     f"{opt_conver['RMS_step'][auxi_idx]:.1f}", 
                     f"{opt_conver['MAX_step'][auxi_idx]:.1f}"]
        
        res_cache = [res_cache[0]] + auxil_res + get_res_cache(res_cache) + \
            [min_eigen[auxi_idx]] + ["YES" if converge_judge[auxi_idx] else "NO"]

    result.append(res_cache)
    
    output_text = ""
    # use the number in format_out to format the output
    for item in result:
        output_text += '\t'.join([f"{item[i]:<{format_out[i]}}" for i in range(len(item))]) + '\n'
        
    if output_filename:
        with open(output_filename, 'w') as f:
            f.write(output_text)
    else:
        print(output_text)
    return output_text
    
def get_charge_mul(lines: List[str]) -> Tuple[int, int]:
    """
    获取电荷和多重度
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        Tuple[int, int]: 电荷和多重度
    """
    text_all = "".join(lines)
    match = re.search(pattern5, text_all, re.MULTILINE | re.IGNORECASE)
    if match:
        charge, mul = match.groups()
        return int(charge), int(mul)
    else:
        raise ValueError("Charge and Multiplicity not found")
   
def get_num_of_atoms(lines: List[str]) -> int:
    """
    获取原子数量
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        int: 原子数量
    """
    text_all = "".join(lines)
    match = re.search(pattern7, text_all, re.MULTILINE)
    if match:
        return int(match.groups()[0])
    else:
        raise ValueError("Number of Atoms not found")
   
    
def get_task_type(lines: List[str]) -> str:
    """
    获取任务类型
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        str: 任务类型 (Opt, TS, Freq, Unknown)
    """
    full_text = "".join(lines)
    if opt_flag_txt in full_text:
        if ts_flag_txt in full_text:
            return "TS"
        else:
            return "Opt"
    elif freq_flag_txt in full_text:
        return "Freq"
    else:
        return "Unknown"

def get_opt_xyz_frame(lines: List[str], n_frame: int = -1, output_filename: str = '') -> str:
    """
    获取优化的XYZ坐标
    
    Args:
        lines: ORCA输出文件的行列表
        n_frame: 帧编号
        output_filename: 输出文件名
        
    Returns:
        str: XYZ坐标字符串
    """
    opt_cycle_begin = get_cartesian_coords_begin(lines)
    charge, mul = get_charge_mul(lines)
    num_of_atoms = get_num_of_atoms(lines)
    xyz_begin_idx = [int(item) for item in opt_cycle_begin]
    xyz_eng_idx = [item + num_of_atoms for item in xyz_begin_idx]
    
    begin_idx_ = xyz_begin_idx[n_frame]
    end_idx_ = xyz_eng_idx[n_frame]
    
    xyz_string = f"{num_of_atoms}\nCharge: {charge} ,Mul: {mul} ,frame: {n_frame} ,type: {get_task_type(lines).lower()}\n"
    xyz_string += '\n'.join(list(map(lambda a: a.strip(), lines[begin_idx_:end_idx_])))
    
    if output_filename:
        with open(output_filename, 'w') as f:
            f.write(xyz_string)
    else:
        return xyz_string


def get_hessian_min_eigen(lines: List[str]) -> List[str]:
    """
    获取Hessian矩阵的最小特征值
    
    Args:
        lines: ORCA输出文件的行列表
        
    Returns:
        List[str]: 最小特征值列表
    """
    text_all = "".join(lines)
    matches = re.findall(pattern8, text_all, re.MULTILINE)
    if matches:
        matches = list(map(lambda a: list(map(lambda b: f"{float(b)*100:.2f}", a)), matches))
        matches_str = list(map(lambda a: " ".join(a), matches))
        return matches_str
    else:
        raise ValueError("Lowest eigenvalues not found")

def output_freq(lines: List[str], freq_type: str = 'opt') -> None:
    """
    输出频率信息
    
    Args:
        lines: ORCA输出文件的行列表
        freq_type: 频率类型 (opt, ts)
    """
    imaginary_mode = get_imaginary_mode(lines)[-1][1]
    
    num_if = imaginary_mode.count("-")
    valid_flag = None
    
    if freq_type == 'opt':
        valid_flag = True if num_if == 0 else False
    elif freq_type == 'ts':
        valid_flag = True if num_if == 1 else False
    else:
        raise ValueError(f"Unknown freq type: {freq_type}")
    
    print(f"frquency type:{freq_type}, imaginary_freq:  {imaginary_mode}")
    print(f"Validity:{"YES" if valid_flag else "NO"}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python orcamonitor.py <filename> [-o] [-i] [-x xyz_frame] [-f freq_type]")
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Process ORCA output files.')
    parser.add_argument('filename', type=str, help='The ORCA output file to process')
    parser.add_argument('-i', '--interactive', action='store_true', help='Enable interactive mode')
    parser.add_argument('-x', '--xyz', type=int, help='Extract the xyz frame, int.')
    parser.add_argument('-o', '--ongoing', action='store_true', help='never raise error if encoutered')
    parser.add_argument('-f', '--freq_type', type=str, default='opt', 
                        help='Works for Freq-Only task, to judge if if qualifies. opt or ts (default: opt)')

    args = parser.parse_args()
    filename = args.filename
    
    if not os.path.exists(filename):
        print(f"错误: {filename} 文件不存在")
        sys.exit(1)
    
    freq_type = args.freq_type.lower()
    if freq_type not in ['opt', 'ts']:
        print(f"错误: freq_type 必须是 'opt' 或 'ts'，当前值: {freq_type}")
        sys.exit(1)
    
    interactive = args.interactive
    
    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
            text_all = ''.join(lines)
    except Exception as e:
        print(f"错误: 读取文件失败: {e}")
        sys.exit(1)
    
    def work() -> None:
        """主要工作函数"""
        if args.xyz:
            try:
                print(get_opt_xyz_frame(lines, n_frame=args.xyz - 1))
            except Exception as e:
                print(f"警告: 提取指定帧失败: {e}")
                print("尝试提取最后一帧...")
                try:
                    print(get_opt_xyz_frame(lines, n_frame=-1))
                except Exception as e2:
                    print(f"错误: 提取最后一帧也失败: {e2}")
                    sys.exit(1)
            sys.exit(0)
        
        print(f"正在监控 ORCA 计算: {filename}")
        try:
            num_of_atoms = get_num_of_atoms(lines)
            charge, mul = get_charge_mul(lines)
            task_type = get_task_type(lines)
            print(f"任务类型: {task_type.lower()}")
            print(f"原子数量: {num_of_atoms}  电荷: {charge}  多重度: {mul}")
        except Exception as e:
            print(f"错误: 提取基本信息失败: {e}")
            sys.exit(1)
        
        if task_type == "Freq":
            try:
                output_freq(lines, freq_type)
            except Exception as e:
                print(f"错误: 处理频率信息失败: {e}")
                sys.exit(1)
            sys.exit(0)
        elif task_type == "Unknown":
            print("错误: 未知任务类型")
            sys.exit(1)
        elif task_type in ["Opt", "TS"]:
            pass
        
        try:
            opt_cycle_begin = get_opt_cycle_begin(lines)
            print(f"优化循环次数: {len(opt_cycle_begin)}")
        except Exception as e:
            print(f"错误: 提取优化循环信息失败: {e}")
            sys.exit(1)
        
        output_ratio = True
        input_frame = None
        try:
            output_pop(lines, output_ratio=output_ratio, format_out=[2, 6, 6, 6, 6, 6, 24, 24, 18, 5])
        except Exception as e:
            print(f"错误: 输出优化信息失败: {e}")
            sys.exit(1)
        
        while True and interactive:
            if output_ratio:
                print("优化收敛信息以比值形式显示 (值/收敛阈值: 标准)")
            else:
                print("优化收敛信息以原始值形式显示, X 10^4")
                
            input_frame = input("请输入要提取的XYZ帧编号，输入't'切换收敛信息显示方式，或输入'q'退出: ")
            if input_frame == 'q':
                break
            elif input_frame == 't':
                output_ratio = not output_ratio
                try:
                    output_pop(lines, output_ratio=output_ratio, format_out=[2, 6, 6, 6, 6, 6, 24, 24, 18, 5])
                except Exception as e:
                    print(f"错误: 输出优化信息失败: {e}")
                    continue
                continue
            else:
                try:
                    input_frame = int(input_frame)
                except ValueError:
                    print("错误: 无效的输入，请输入数字")
                    continue
                try:
                    res = get_opt_xyz_frame(lines, n_frame=input_frame)
                    print(f"帧 {input_frame} 已提取")
                    print(res)
                except Exception as e:
                    print(f"错误: 提取XYZ帧失败: {e}")
                    continue
        
        sys.exit(0)
        
    if args.ongoing:
        try:
            work()
        except Exception as e:
            print(0)
            sys.exit(1)
    else:
        work()
