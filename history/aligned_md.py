import MDAnalysis as mda
from MDAnalysis.analysis import align
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="分子轨迹.xyz")
parser.add_argument('-f',"--file",required=True,type=str, help="xyz文件")
parser.add_argument('-p',"--part",required=True,type=str, help="原子序号")

args = parser.parse_args()

# 你的轨迹文件和拓扑文件
xyz_file = args.file
group_atoms_str = args.part

print(f'part:{group_atoms_str}')
# 加载轨迹文件
u = mda.Universe(xyz_file)

def parse_atom_indices(atom_indices_str):
    indices = []
    for part in atom_indices_str.split(','):
        if '-' in part:
            start, end = part.split('-')
            indices.extend(range(int(start), int(end) + 1))
        else:
            indices.append(int(part))
    indices = np.array(indices) - 1   # MDAnalysis 与 Gaussian 中编号不同，需要修改
    return indices

group_atoms = parse_atom_indices(group_atoms_str)

# Define the atom selection string for MDAnalysis
selection_string = 'index ' + ' '.join(map(str, group_atoms))

# Select the reference atoms for alignment
ref_atoms = u.select_atoms(selection_string)

# Create a reference Universe from the first frame of the trajectory
ref = u.select_atoms('all').copy()

# Align the trajectory to the reference atoms
aligner = align.AlignTraj(u, ref, select=selection_string, in_memory=True)
aligner.run()

# Save the aligned trajectory
aligned_trajectory_filename = "aligned_trajectory.xyz"
with mda.Writer(aligned_trajectory_filename, ref.n_atoms) as W:
    for ts in u.trajectory:
        W.write(u)

print(f'Aligned trajectory saved as {aligned_trajectory_filename}')