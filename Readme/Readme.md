## 计算方法一览

### 目前可计算的体系：TADF、（双）自由基、ESIPT (Gaussian、Orca)；单分子动力学 (cp2k)；

### 1. 基态

&emsp;&emsp;分子或原子在最低能量时的电子结构和几何构型。

&emsp;&emsp;可计算：限制性优化、分子轨道模拟、静电势模拟 (ESP)、偶极矩计算、振动分析、分子体积（pdb 盒子）、热力学函数计算 (内能 U、焓 H、功函 A、吉布斯自由能 G)。

### 2. 垂直激发态

&emsp;&emsp;垂直激发态指分子在基态几何结构不变的情况下，通过吸收光子直接跃迁到激发态。该跃迁遵循 Franck–Condon 原理，因此激发过程发生得非常快，分子结构来不及改变。结果是在激发态下，分子通常处于一个非平衡的、高能量的几何构型。

&emsp;&emsp;可计算：激发态能级（LE 态可选择 Orca 使用双杂化计算）、电子-空穴分布、振子强度、电子-空穴距离/重叠程度、CT/LE 指数、偶极矩计算。空穴-电子重叠等值面（即 Sm 或 Sr 函数）、Cele 和 Chole 等值面、跃迁密度分布、跃迁偶极矩（VMD 箭头绘图）、跃迁磁偶极矩密度、电荷密度分布、空穴-电子库仑吸引能。

### 3. 绝热激发态

&emsp;&emsp;绝热激发态指分子在吸收光子后，沿着能量最低路径从基态几何结构缓慢调整到激发态的最低能量几何构型。这个过程保持电子状态不变，但允许分子结构在激发态中重新优化。绝热激发态通常用于描述分子在稳定的激发态几何构型下的能量和性质。

&emsp;&emsp;可计算：激发态能级、电子-空穴分布、振子强度、电子-空穴距离/重叠程度、CT/LE 指数、偶极矩计算、激发态重组能、分子构象差异分析 (RMSD)、振动分析。

### 4. 调控密度泛函方法计算分子（Tuned Density Functional Theory，TDFT）

&emsp;&emsp;介绍：[最优化“调控”区间分离密度泛函理论的研究进展](https://doi.org/10.3866/PKU.WHXB201605301)

&emsp;&emsp;这种方法特别适用于那些标准泛函难以准确描述的系统，比如含有电荷转移激发或较大能带间隙的材料。通过调整密度泛函理论中的参数，如交换-相关泛函（Exchange-Correlation Functional）中的参数，来更精确地描述电子结构，尤其是激发态和光学性质。调控密度泛函的关键在于找到最佳的参数设置，使得计算结果与实验数据或其他高精度计算结果相吻合。这通常涉及到对泛函中的长程校正参数（如ω参数）进行优化，以使得如 HOMO 能量更接近电离能，从而提高激发能、极化率等性质的计算精度。

&emsp;&emsp;这种方法虽然可以提高计算的准确性，但也可能带来计算成本的增加，因为需要对每个体系单独优化参数。调控密度泛函的关键在于正确选择泛函（例如：lc-blyp、lc-whpbe、wb97xd）和调整参数 (ω)。

### 5. 自然跃迁轨道分析 (Natural transition orbital, NTO)

&emsp;&emsp;介绍：[计算化学公社：NTO 计算](http://sobereva.com/91)

&emsp;&emsp;虽然电子跃迁问题本质是体系电子态与电子态之间的跃迁，但是人们往往喜欢通过轨道模型来描述，以使跃迁模式更易于考察和理解。NTO 的思路是将多组态波函数密度矩阵对角化，得到的本征向量就是自然轨道，通常只需要本征值最大的 $n_{occ}$ 个轨道就能很充分地描述体系密度矩阵，冗余的信息就被去掉了，使一大堆轨道跃迁模式转化为一个“紧凑”的轨道跃迁模式。

&emsp;&emsp;需要说明的是电子-空穴分析比 NTO 强大得多、普适性强得多，而且使用很方便，因此强烈建议用空穴-电子分析而非 NTO！NTO 分析唯一的好处仅在于可以给出相位信息，你不需要专门考察相位的话一律应当用空穴-电子分析代替 NTO。

### 6. 片段间电荷转移分析 (Interfragment charge transfer analysis, IFCT)

&emsp;&emsp;介绍：[思想家公社：IFCT 计算](http://sobereva.com/433)

&emsp;&emsp;用于研究分子中不同部分之间电子密度变化的计算化学方法，帮助理解分子激发态内电荷的分布和转移情况。在 IFCT 分析输出中会直接给出电荷转移百分比 CT (%) 和局域激发百分比 LE (%)，在指认电子激发类型的时候极其方便和严格！

### 7. 空间电荷转移分析 (Space charge transfer analysis, TSCT)

&emsp;&emsp;介绍：[计算化学公社：TSCT 计算](http://bbs.keinsci.com/thread-18597-1-1.html)

&emsp;&emsp;在片段间电荷转移分析的基础上研究分子中电子给体（D）与电子受体（A）之间通过空间相互作用实现电子转移的现象。
电荷转移类型：空间电荷转移 (Through-space charge transfer, TSCT) 与化学键电荷转移 (Through-bond charge transfer, TBCT)

### 8. 轨道离域指数分析 (Orbital delocalization index analysis, ODI)

&emsp;&emsp;介绍：[思想家公社：ODI 计算](http://sobereva.com/525)

&emsp;&emsp;轨道离域指数是一个衡量轨道离域性的指标，它通过计算分子中各个原子对特定轨道的贡献来定量描述轨道的离域程度。

### 9. 片段相互作用分析 (Independent gradient model based on Hirshfeld partition, IGMH)

&emsp;&emsp;介绍：[思想家公社：IGMH 计算](http://sobereva.com/621)

&emsp;&emsp;片段相互作用分析（IGMH）是一种基于 Hirshfeld 分子密度分区的独立梯度模型（Independent Gradient Model），用于可视化研究化学系统中的相互作用，如：氢键、空间位阻效应、笼效应、范德华相互作用等。

### 10. 势能面扫描 (Potential energy surface scanning, PES)

介绍：[思想家公社：PES 计算](http://sobereva.com/474)

&emsp;&emsp;它用于研究分子或一组分子在反应过程中 (构象变化过程中) 的能量变化情况。

#### &emsp;10.1 刚性扫描（Rigid scan）

&emsp;&emsp;在这种扫描中，逐渐改变选定的内坐标，并对分子体系进行单点能计算，而分子内其他变量保持不变。

#### &emsp;10.2 柔性扫描（Relaxed scan）

&emsp;&emsp;在柔性扫描中，改变被扫描的变量，在每个结构处固定当前变量，并优化分子结构（实际是一个限制性优化）

### 11. 黄里因子计算 (Huang-Rhys factor)

&emsp;&emsp;介绍：[思想家公社：分解重组能和 Huang-Rhys 因子计算](http://sobereva.com/330)

&emsp;&emsp;分子的电子跃迁过程中，分子的振动模式与电子态之间的相互作用。这种耦合效应对分子光谱以及能量转移等特性有重要的影响。使用 FCclasses 3 软件包计算。

### 12. 计算电子跃迁速率 (Marcus equation)

&emsp;&emsp;介绍：[Libretexts：Marcus Theory for Electron Transfer（我也看不懂）](<https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Time_Dependent_Quantum_Mechanics_and_Spectroscopy_(Tokmakoff)/15%3A_Energy_and_Charge_Transfer/15.05%3A_Marcus_Theory_for_Electron_Transfer>)

&emsp;&emsp;公式 1：$k_{T, \text { Marcus }}=\frac{\pi}{\hbar}\left|\hat{H}_{i f}\right|^{2} \frac{1}{\sqrt{\pi \lambda_{M} k_{\mathrm{B}} T}} \exp \left(\frac{-\left(\Delta G^{\circ}+\lambda_{M}\right)^{2}}{4 \lambda_{M} k_{\mathrm{B}} T}\right)$

&emsp;&emsp;公式 2：$k_{T, M \mathrm{JL}, \mathrm{eff}}=\frac{\pi}{\hbar}|\mathrm{NAC}|^{2} \frac{1}{\sqrt{\pi \lambda_{M} k_{\mathrm{B}} T}} \sum_{n=0}^{\infty} \exp (-S) \frac{S^{n}}{n!} \exp \left(\frac{-\left(\Delta G^{\circ}+\lambda_{M}+n \hbar \omega_{\text {eff }}\right)^{2}}{4 \lambda_{M} k_{\mathrm{B}} T}\right)$

&emsp;&emsp;Marcus方程是用于计算电子跃迁速率的一个经典理论模型，由 R. A. Marcus 在1956年首次提出。这个方程主要用于描述外层电子转移过程的速率常数，后来被扩展和改进以适用于不同类型的电子转移过程。以下是Marcus方程的核心内容和计算电子跃迁速率的关键参数：$SOC(H_{if})$，$ΔE_{st}(ΔG)$， $λ_M$ 有效重组能，$ω_{eff}$ 有效振动频率；$S$ 黄里因子，$\hbar$，$k_B$，$T$ 均为常数。

### 13. 单分子动力学 (Monomolecule dynamics method)

&emsp;&emsp;一种用于研究单个分子运动和动态行为的模拟方法。

### 14. 键解离能 (Bond disscociation energy, BDE)

&emsp;&emsp;介绍：[墨灵格的博客：BDE 计算](http://blog.molcalx.com.cn/2017/11/05/gaussian-tutorial-bond-dissociation-energy.html)

&emsp;&emsp;计算分子中化学键断裂过程的反应焓变，它反映了键断裂过程所需要的能量。

### 15. 单分子反应分析

&emsp;&emsp;介绍：[思想家公社：TS、IRC 计算](http://sobereva.com/44)

&emsp;&emsp;过渡态分析 (Transition state analysis) 、反应路径分析（IRC）

&emsp;&emsp;渡态分析 和 IRC分析 是理解化学反应机理，可以精确地定位过渡态、计算反应的能量障碍，并探讨反应的路径。例如：激发态分子内质子转移 (ESIPT) 、螺手性分子反转能 ...


### 16. 芳香性计算

&emsp;&emsp;芳香性计算方法一览：[思想家公社：芳香性计算](http://sobereva.com/176)

#### &emsp;16.1 核独立化学位移（Nucleus-independent chemical shifts, NICS）

&emsp;&emsp;介绍：[墨灵格的博客：NICS 计算](http://blog.molcalx.com.cn/2021/01/09/calculate-nics-to-evaluate-aromaticity.html)

&emsp;&emsp;NICS是一个在人为设定的不在原子核位置上的磁屏蔽值的负值，负值越大（即对磁场屏蔽越强）则芳香性越强
。最初NICS(0)是将这个位置取在共轭环的几何中心。

&emsp;&emsp;后来，为了更准确地反映 π 电子的贡献，提出了在环中心上方或下方 1 埃的位置计算 NICS，称为 NICS(1)。进一步地，有研究者提出只考察磁屏蔽张量的 zz 值（z 轴垂直于环平面），以更好地体现 π 芳香性，这被称为 NICS(1)_zz。


#### &emsp;16.2 磁感应电流模拟

&emsp;&emsp;介绍：[思想家公社：AICD 计算磁感应电流](http://sobereva.com/294)

&emsp;&emsp;软件部署：[王喆的博客：芳香性可视化 AICD](https://wongzit.github.io/visualization-of-aromaticity-aicd/)

&emsp;&emsp;用于衡量芳香性的方法，它基于分子在外加磁场下的电子响应。当分子体系中的电子在整体或某个部分有很强的离域性时，例如苯环的π电子，在外加磁场作用下，会在相应区域产生一圈明显的感应环形电流。这种电流密度的分布可以用来考察体系的电子离域性，进而研究芳香性。

&emsp;&emsp;环电流的方向与芳香性有直接关系。如果环电流方向与左手规则相同时（称为 diamagnetic 或 diatropic 环电流），电流越大，则芳香性越强。相反，如果电流方向与左手规则相反（paratropic 环电流），电流越大则表现出越强的反芳香性。如果没有产生明显净电流，则为非芳香性。

#### &emsp;16.3 化学屏蔽表面 (Iso-chemical shielding surfaces, ICSS)

&emsp; &emsp;介绍：[思想家公社：ICSS 计算](http://sobereva.com/216)

&emsp;&emsp;ICSS 是指在分子周围某个特定等值面上，该面上的磁屏蔽值相等。这个等值面能够直观地展示分子中电子密度的分布情况，尤其是 π 电子的分布，从而反映分子的芳香性。

&emsp;&emsp;在芳香性区域，由于π电子的离域，磁屏蔽通常较强，因此在 ICSS 图中，芳香性区域通常会显示为屏蔽较强的区域（即低化学位移值的区域）。相反，非芳香性区域或反芳香性区域则显示出较弱的磁屏蔽。

### 17. 分子中的原子极化率 (Atomic polarizabilities in molecules)

&emsp;&emsp;介绍：[思想家公社：原子极化率计算](http://sobereva.com/600)

&emsp;&emsp;原子在孤立状态下的极化率是可以实验测量的，也很容易理论计算，在 [ctcp.massey](http://ctcp.massey.ac.nz/index.php?menu=dipole&page=dipole) 有全周期表各种元素的实验和高精度理论计算数据。分子的极化率可以视为是其中各个原子的有效 (effective) 极化率的总和，但是分子环境中原子极化率通常不是实验可观测的（除非是所有原子等价，比如C60，就用分子极化率除以60就是各个碳的极化率），但可以通过计算说明。

&emsp;&emsp;可以输出的信息：原子有效体积、自由体积和极化率。

<p align="right">Last update: 25/01/10</p>
<p align="right"> —— By Alpha Ray</p>

