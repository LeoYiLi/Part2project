# Section 1 README: MD Quality Check

这个文件夹用于整理 **Result and Discussion 第 1 节：MD Quality Check** 所需的原始数据、分析任务、图表计划和后续写作依据。

本 README 的作用不是写正文，而是明确：

- 第 1 节要回答什么问题；
- 文件夹里需要放哪些数据；
- 后续要做哪些分析；
- main text 和 supplementary 分别放哪些图和表；
- 后续如何衔接 `nature-figure` 作图和正文写作。

## 1. 研究目的

第 1 节的核心目的，是证明用于 **VACF/VDOS** 分析的 MD 轨迹具有可靠的有限温度采样窗口。

这里要证明的是：

1. 选定的 production window 已经热平衡，不再主要受初始结构弛豫影响。
2. 缺陷在 VDOS 分析窗口内的状态是明确的：保持、迁移、复合，或需要分段分析。
3. 主要 VDOS 特征具有数值稳定性，不是由温度漂移、能量漂移、COM drift 或单个时间块噪声造成的。

需要注意的表述边界：

- 本节**不是**证明缺陷具有长期热力学稳定性。
- 不建议写 “the defect is stable”；应写 “metastable over the simulated timescale”。
- 不建议写 “no recombination occurs”；应写 “no recombination was observed within the simulated trajectory”。
- 如果没有多条 independent runs，不要写 recombination rate、defect lifetime 或 diffusion coefficient。

## 2. 文件夹中存放的数据

| 文件类型 | 用途 |
|---|---|
| `initial cif` | MD 前的初始结构 |
| `final cif` | MD 后的最终结构 |
| `log` | 提取 temperature、potential energy、total energy 等随时间变化的数据 |
| `traj` | 提取 RMSD、RDF、COM drift、defect descriptor、VACF/VDOS 等 |
| `md.py` | 记录 MD 参数，例如 timestep、thermostat、friction、目标温度、运行步数、输出间隔 |
| optional VDOS data | 已经算好的 VACF/VDOS 或 block VDOS 数据 |



## 3. 文件配对表

把文件放进文件夹后，在这里补全每个体系对应的文件。

| 体系 | Initial CIF | Final CIF | LOG | TRAJ | MD script | 备注 |
|---|---|---|---|---|---|---|
| pristine |  |  |  |  |  |  |
| vacancy |  |  |  |  |  |  |
| interstitial |  |  |  |  |  |  |
| Frenkel 1 |  |  |  |  |  |  |
| Frenkel 2 |  |  |  |  |  |  |
| Frenkel 3 |  |  |  |  |  |  |
| Frenkel 4 |  |  |  |  |  |  |

## 4. Frenkel defect 说明

在计算 `Oi-Ov distance-t` 之前，需要先手动补全 Frenkel defect 的位置信息。

| Frenkel case | Vacancy site / 被移除的 O | Interstitial O | 初始 Oi-Ov distance | 局部位置说明 | 备注 |
|---|---|---|---|---|---|
| Frenkel 1 |  |  |  |  |  |
| Frenkel 2 |  |  |  |  |  |
| Frenkel 3 |  |  |  |  |  |
| Frenkel 4 |  |  |  |  |  |

这里需要说明：

- 哪个 O 被移除，形成 vacancy；
- 哪个 O 是 interstitial O；
- vacancy site 是如何定义的；
- Frenkel pair 在 MD 过程中是保持、迁移，还是发生 recombination。

## 5. 需要做的主要分析

| 分析内容 | 数据来源 | 目的 |
|---|---|---|
| Temperature vs time, `T-t` | `log` | 判断轨迹是否达到热平衡，并确定 production window |
| Potential energy vs time, `Epot-t` | `log` | 判断初始结构弛豫是否结束 |
| Total energy trend | `log` | 检查是否存在异常数值漂移；如果是 NVT，需要谨慎解释 |
| Production window selection | `log` + `traj` | 去掉初始 heating / relaxation 段，确定用于 VDOS 的时间窗口 |
| COM drift | `traj` | 检查或去除整体平移，避免低频 VDOS 伪信号 |
| RMSD / displacement | `traj` + `initial cif` | 判断结构是否在 production window 内保持合理 |
| RDF / bond-length distribution | `traj` | 检查局部结构是否合理，是否存在非物理短接触 |
| Coordination number | `traj` | 跟踪缺陷附近局部配位环境 |
| Vacancy / interstitial descriptor | `traj` | 判断 vacancy 或 interstitial 在 VDOS 窗口中的状态 |
| Frenkel `Oi-Ov distance-t` | `traj` + Frenkel defect 说明 | 判断 Frenkel pair 是否保持、迁移或复合 |
| VACF decay | `traj` 或 VDOS data | 判断 VACF 是否在 Fourier transform 前衰减到合理基线 |
| VDOS block convergence | `traj` 或 VDOS data | 判断主要 VDOS 峰是否在不同时间块中稳定 |

## 6. Methods 中需要交代的内容

Methods 不需要放太多质量检查结果，但必须交代 MD 和 VDOS 的计算设置，保证别人能复现。

### 6.1 MD simulation 设置

需要从 `md.py`、`log` 和结构文件中整理：

| 内容 | 数据来源 | 写入位置 |
|---|---|---|
| 体系信息：pristine / vacancy / interstitial / Frenkel | 文件名 + `initial cif` | Methods |
| supercell size 和原子数 | `initial cif` | Methods 或 Methods table |
| calculator / potential | `md.py` | Methods |
| ensemble / thermostat | `md.py` | Methods |
| target temperature | `md.py` + `log` | Methods |
| timestep | `md.py` | Methods |
| friction / thermostat parameter | `md.py` | Methods |
| MD 总步数和总时间 | `md.py` + `log` | Methods |
| trajectory / log 输出间隔 | `md.py` | Methods |
| initial relaxation / production window 选择原则 | `log` + 分析结果 | Methods |

### 6.2 VACF/VDOS 设置

需要说明：

- VDOS 来自 velocity autocorrelation function, VACF；
- `traj` 中需要使用 velocities；
- 是否去除了 frame-wise mean velocity / COM drift；
- 使用的 production window；
- VACF 截断长度；
- 是否使用 window function；
- 是否进行 smoothing / normalization；
- 是否做 3-block VDOS convergence。

### 6.3 Methods 中建议做的表

| 表格 | 内容 | 备注 |
|---|---|---|
| Methods Table 1 | MD metadata：体系、原子数、缺陷类型、timestep、temperature、thermostat、friction、run length、output interval | 可以放 Methods；如果太长则放 Appendix/Supplementary |
| Methods Table 2 | Production window：丢弃时间、production start/end、用于 VDOS 的帧数 | 可以放 Methods 或 Appendix/Supplementary |

## 7. Result and Discussion Section 1 中建议放的内容

Result and Discussion 第 1 节只放最关键的质量检查结果。目标是让读者相信：后续 VDOS 分析使用的轨迹窗口是合理的。

### 7.1 Main text 推荐图

main text 建议做一个紧凑的四面板图。

推荐图题：

**Trajectory reliability for finite-temperature VDOS analysis**

| Panel | 图 | 体系 | 支持的结论 |
|---|---|---|---|
| a | `T-t` | 代表性的 pristine / vacancy / interstitial / Frenkel | 轨迹达到有限温度采样窗口 |
| b | `Epot-t` | 同上 | 初始弛豫已经结束 |
| c | Frenkel `Oi-Ov distance-t` | 代表性 Frenkel cases | Frenkel defect 在 VDOS 窗口内状态明确 |
| d | VDOS block convergence | 代表性体系 | 主要 VDOS 特征具有数值稳定性 |

### 7.2 Main text 中可以写的分析点

正文建议围绕这几句话展开：

- `T-t` 显示所选 production window 已经达到目标温度附近，并且没有持续升温或降温漂移。
- `Epot-t` 显示初始结构弛豫已经结束，production window 中能量围绕稳定均值波动。
- Frenkel `Oi-Ov distance-t` 用于判断 Frenkel pair 在 VDOS 窗口中是保持、迁移还是复合。
- VDOS block convergence 显示主要峰位或主要谱特征在不同时间块中一致。
- 因此，这些轨迹可以用于有限时间 VACF/VDOS 分析。

### 7.3 Result and Discussion 中建议做的表

| 表格 | 内容 | 备注 |
|---|---|---|
| Results Table 1 | Defect-state summary：vacancy / interstitial / Frenkel 在 production window 中的状态 | 如果文字能说清楚，可以不放表 |
| Results Table 2 | 代表性 VDOS block convergence summary：主要峰位和 block-to-block variation | 如果图已经足够清楚，可以放 Appendix/Supplementary |

## 8. Appendix / Supplementary 中建议放的内容

Appendix 或 Supplementary 用来放完整质量检查，不占用 main text 篇幅。

### 8.1 Supplementary figures

| 图 | 内容 |
|---|---|
| Fig. S1 | 所有体系的 Temperature vs time |
| Fig. S2 | 所有体系的 Potential energy vs time |
| Fig. S3 | 所有体系的 Total energy trend |
| Fig. S4 | 所有体系的 RMSD 或 mean displacement |
| Fig. S5 | COM drift 或 mean velocity check |
| Fig. S6 | Vacancy 和 interstitial descriptors |
| Fig. S7 | 所有 Frenkel cases 的 `Oi-Ov distance-t` |
| Fig. S8 | RDF / bond-length distributions |
| Fig. S9 | 所有体系的 VACF decay |
| Fig. S10 | 所有体系的 VDOS block convergence |

### 8.2 Supplementary tables

| 表格 | 内容 | 建议位置 |
|---|---|---|
| Table S1 | 完整 MD metadata：体系、原子数、缺陷类型、timestep、temperature、thermostat、friction、run length、output interval | Appendix / Supplementary |
| Table S2 | 完整 production window：丢弃时间、production start/end、用于 VDOS 的帧数 | Appendix / Supplementary |
| Table S3 | Production window 内的 temperature 和 energy 统计量 | Appendix / Supplementary |
| Table S4 | 所有 defect-state summary：保持、迁移、复合，或需要分段分析 | Appendix / Supplementary |
| Table S5 | 所有 VDOS block convergence summary：主要峰位和 block-to-block variation | Appendix / Supplementary |
| Table S6 | 如果有 independent runs，列出 seed、event status、event time | Appendix / Supplementary |

## 9. 写作前最低证据清单

开始写 Section 1 正文前，至少需要完成：

- 文件配对表已补全；
- Frenkel defect 位置信息已补全；
- 每个体系的 production window 已确定；
- 每个体系的 `T-t` 和 `Epot-t` 已检查；
- vacancy、interstitial 和 Frenkel 的 defect descriptors 已检查；
- VACF decay 和 VDOS block convergence 已检查；
- COM drift 已量化，或在 VDOS 后处理中已去除。

## 10. 后续工作流程

1. 把所有 raw files 放入本文件夹。
2. 补全文件配对表和 Frenkel defect 说明。
3. 使用 `nature-figure` 画 main text 和 supplementary figures。
4. 使用 `academic-research-suite` 或 `nature-writing` 写 Section 1 正文。

使用 `nature-figure` 时，开始前需要先确定作图后端：Python or R。
