import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 常数定义
H_if_cm = 1  # 单位：cm^-1
h = 6.62607015e-34  # 普朗克常数，单位：J·s
c = 2.998e10  # 光速，单位：cm/s
H_if = H_if_cm * h * c  # 转换为J

T = 298  # 温度，单位：开尔文
h_bar = 1.0545718e-34  # 约化普朗克常数，单位：J·秒
k_B = 1.380649e-23  # 玻尔兹曼常数，单位：J/K

eV_to_J = 1.60218e-19  # eV到J的转换因子

# 创建ΔG和λM的网格，单位为eV
ΔG_range_eV = np.linspace(-0.3, 0.3, 6000)  # ΔG范围，从-2到2 eV
λM_range_eV = np.linspace(0, 0.25, 25000)  # λM范围，从0.01到2 eV
ΔG_eV, λM_eV = np.meshgrid(ΔG_range_eV, λM_range_eV)

# 转换为J
ΔG = ΔG_eV * eV_to_J
λM = λM_eV * eV_to_J

# 根据给定公式计算kT
prefactor = (np.pi / h_bar) * np.abs(H_if)**2
exponential = np.exp(-( (ΔG + λM)**2 ) / (4 * λM * k_B * T))
kT = prefactor * (1 / np.sqrt(np.pi * λM * k_B * T)) * exponential

# 对kT取10的对数
log_kT = np.log10(kT)

# 将超出范围的数据设为NaN
# log_kT = np.where((log_kT < 0) | (log_kT > 10), np.nan, log_kT)

# 绘图
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(ΔG_eV, λM_eV, log_kT, cmap='viridis', edgecolor='none')

ax.set_xlabel('ΔEst (eV)', fontsize=13)
ax.set_ylabel('λM (eV)', fontsize=13)
ax.set_zlabel('log10(krisc) (s^-1)', fontsize=12)
ax.tick_params(axis='both', which='major', labelsize=12)

# ax.set_zlim(0, 10)

plt.title('3D Surface Plot of log10(krisc) as a Function of ΔEst and λM')
plt.colorbar(surf, shrink=0.5, aspect=5)
plt.show()