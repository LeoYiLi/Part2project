# Section 2 README: Frenkel Pair Geometry and Recombination Behavior

这个文件夹用于整理 **Result and Discussion 第 2 节：Frenkel Pair Geometry and Recombination Behavior** 所需的原始数据、分析任务、图表计划和后续写作依据。

本 README 的作用不是写正文，而是明确：

- 第 2 节要回答什么问题；
- Frenkel cases 应如何分类和比较；
- 每个 case 必须记录哪些结构和轨迹信息；
- 哪些数据可以从 `log` 和 `traj` 中提取；
- 已有的数据提取代码放在哪里；
- `Oi-Ov distance` 在 initial、relaxed、production-start 和 MD time series 中分别是什么意思；
- recombination 应如何定义，而不是只凭肉眼判断；
- main text 和 supplementary 分别放哪些图和表；
- 后续如何衔接 `nature-figure` 作图和正文写作。

## 1. 研究目的

第 2 节的核心目的，是把多组 Frenkel pair 模拟组织成一个受控比较，而不是在主文中逐个展示所有 Frenkel cases。

本节要支持的核心主张是：

> Frenkel pair behavior in YBa2Cu3O7 depends on vacancy-interstitial geometry, especially initial separation and local oxygen-site environment. Close pairs may recombine or strongly rearrange, while more separated pairs may remain metastable over the simulated MD window.

这里需要注意：

- 这只是建议主张，最终必须由 `Oi-Ov distance`、coordination、snapshots 和轨迹结果支持。
- 主文不应逐个展示所有 Frenkel cases；主文应展示代表性构型和对比逻辑。
- Supplementary 应放完整 case inventory，证明分析不是 cherry-picking。
- 如果观察到 interstitial 和 vacancy 在同一个 crystallographic unit cell 时发生 recombination，应将其作为一个可检验的几何分类结果，而不是单独依赖这个描述下结论。

建议写作边界：

- 可以写 “This Frenkel configuration remains separated over the 40 ps simulated window.”
- 可以写 “No recombination was observed within the sampled trajectory.”
- 可以写 “The configuration is metastable on the simulated timescale.”
- 不建议写 “stable defect”。
- 不建议写 “will not recombine”。
- 不建议写 “recombination rate”，除非有多个 independent seeds 和足够事件数。
- 不建议写 “diffusion coefficient”，除非有长时间、统计充分的迁移轨迹。
- 不建议写 “thermodynamic stability”，除非有自由能或更系统的采样。

## 2. 文件夹中存放的数据

| 文件类型 | 用途 |
|---|---|
| `initial cif` | Frenkel pair 初始结构，记录 vacancy 和 interstitial 的初始几何关系 |
| `relaxed cif` | MD 前结构弛豫后的 Frenkel pair 几何 |
| `production-start cif` | VDOS 或 production window 起点结构 |
| `final cif` | MD 后最终结构，用于判断有限时间演化结果 |
| `log` | 提取 temperature、potential energy、total energy 和事件时间 |
| `traj` | 提取 `Oi-Ov distance-t`、coordination、site occupancy、snapshots 和 defect migration |
| `md.py` | 记录 MD 参数，例如 timestep、thermostat、friction、目标温度、运行步数、输出间隔 |
| optional analysis data | 已经提取好的 distance、coordination、occupancy 或 event summary 数据 |
| optional snapshots | 代表性 initial / relaxed / production-start / final snapshots |

### 2.1 与 Section 1 数据的关系

Section 2 使用的 `traj` 和 `log` 与 Section 1 的 MD quality check 是同一批原始 MD 输出，不需要重新生成一套独立轨迹。

两节的区别是：

| 数据 | Section 1 用途 | Section 2 用途 |
|---|---|---|
| `log` | 检查 temperature、energy、production window | 提取每个 Frenkel case 的 MD 时长、温度、能量和事件时间参考 |
| `traj` | 检查 trajectory quality、VACF / VDOS 稳定性 | 提取 `Oi-Ov distance-t`、vacancy site occupancy、local coordination 和 recombination / migration event summary |
| `cif` | 记录初始、弛豫和最终结构 | 定义 vacancy site、interstitial site、same-cell / different-cell 关系和初始几何 |

因此，后续只需要把所有 Frenkel cases 对应的 `cif`、`log` 和 `traj` 文件放入项目中，就可以从同一套 MD 输出中提取 Section 2 所需的 distance、coordination、occupancy 和 event summary。

### 2.2 已有提取代码

本文件夹中已经放置了一个可复用的数据提取脚本：

| 脚本 | 作用 |
|---|---|
| `extract_frenkel_metrics.py` | 从 Frenkel `cif`、`log` 和 `traj` 中提取 case inventory、MD log summary、`Oi-Ov distance-t`、vacancy occupancy、local Cu-O coordination 和 event summary |

运行后会在本文件夹下生成：

| 输出 | 内容 |
|---|---|
| `analysis_outputs/frenkel_static_case_inventory.csv` | 每个 Frenkel case 的 vacancy/interstitial site、same-cell 分类、cell offset 和 initial `Oi-Ov distance` |
| `analysis_outputs/frenkel_log_summary.csv` | 每个 case 的 MD 时长、温度和能量摘要 |
| `analysis_outputs/frenkel_event_summary.csv` | 每个有 `traj` 的 case 的 recombination / metastable / missing-trajectory 判断摘要 |
| `analysis_outputs/*_timeseries.csv` | 每条 trajectory 的逐帧 `Oi-Ov distance`、vacancy occupancy 和 coordination time series |

如果某个 case 只有 `log` 但缺少 `traj`，脚本只能提取 temperature / energy / MD duration，不能可靠判断 recombination。recombination、coordination recovery 和 vacancy occupancy 必须依赖包含原子坐标的 `traj`。

## 3. 文件配对表

把文件放进文件夹后，在这里补全每个 Frenkel case 对应的文件。

| Case ID | Initial CIF | Relaxed CIF | Production-start CIF | Final CIF | LOG | TRAJ | MD script | 备注 |
|---|---|---|---|---|---|---|---|---|
| F1 |  |  |  |  |  |  |  |  |
| F2 |  |  |  |  |  |  |  |  |
| F3 |  |  |  |  |  |  |  |  |
| F4 |  |  |  |  |  |  |  |  |
| F5 |  |  |  |  |  |  |  |  |

## 4. Frenkel case inventory

每个 Frenkel case 必须有统一索引。后续所有图、轨迹、VDOS 文件和正文讨论都使用同一个 Case ID。

| Case ID | Vacancy site | Vacancy type | Interstitial site | Oi type | Same parent unit cell? | Cell offset | Initial Oi-Ov distance (A) | Relaxed Oi-Ov distance (A) | Production-start distance (A) | Minimum MD distance (A) | Final MD distance (A) | Recombined? | Event time (ps) | Outcome class | Notes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| F1 |  | chain / plane / apical |  | interstitial channel / other | yes / no |  |  |  |  |  |  | yes / no / unclear |  | recombined / metastable / migrated / unclear |  |
| F2 |  | chain / plane / apical |  | interstitial channel / other | yes / no |  |  |  |  |  |  | yes / no / unclear |  | recombined / metastable / migrated / unclear |  |
| F3 |  | chain / plane / apical |  | interstitial channel / other | yes / no |  |  |  |  |  |  | yes / no / unclear |  | recombined / metastable / migrated / unclear |  |
| F4 |  | chain / plane / apical |  | interstitial channel / other | yes / no |  |  |  |  |  |  | yes / no / unclear |  | recombined / metastable / migrated / unclear |  |
| F5 |  | chain / plane / apical |  | interstitial channel / other | yes / no |  |  |  |  |  |  | yes / no / unclear |  | recombined / metastable / migrated / unclear |  |

这里需要说明：

- 哪个 O 被移除，形成 vacancy；
- 哪个 O 是 interstitial O；
- vacancy site 属于 chain、plane、apical 还是其他局部环境；
- interstitial 的初始位置属于哪个局部通道或邻近环境；
- vacancy 和 interstitial 是否来自同一个 crystallographic unit cell；
- 如果不是同一个 unit cell，需要记录 cell offset，例如 `(0,0,0)`、`(1,0,0)`、`(0,1,0)`；
- Frenkel pair 在 MD 过程中是 recombined、metastable、migrated，还是 unclear。

## 5. `Oi-Ov distance` 的定义和作图逻辑

本节的核心图不是静态结构中某一个距离，而是 **MD 时间内 `Oi-Ov distance` 如何随时间变化**。

推荐主图：

> `Oi-Ov distance vs time`

其中：

- x-axis: MD time, ps；
- y-axis: distance to original vacancy site, A；
- 不同曲线表示 close / intermediate / separated Frenkel cases；
- 用水平虚线标出 recombination / vacancy-occupancy cutoff；
- 如果发生 recombination，用竖线标出 event time；
- 如果空间允许，可在下方面板加入 local Cu-O coordination vs time。

### 5.1 `Oi` 和 `Ov` 分别是什么

| 符号 | 含义 | 注意事项 |
|---|---|---|
| `Oi` | interstitial oxygen，即插入的间隙氧，或轨迹中被识别为 interstitial-like 的氧 | 如果发生 oxygen exchange，不应只依赖最初插入的 atom ID |
| `Ov` | original vacancy site，即最初被移除的氧晶格位置 | vacancy 是一个 site coordinate，不是真实原子 |

因此，`Oi-Ov distance` 不是“两个不同 oxygen sites 之间的距离”，而是 **interstitial oxygen 到原始 vacancy site 的距离**。

因为 `Ov` 是一个晶格 site，recombination 后距离不一定严格等于 `0 A`。MD 中氧原子会热振动，所以更合理的判断是：距离下降到预先定义的 occupancy cutoff 以下，并且持续一段时间。

### 5.2 静态距离和动态距离的区别

| 距离 | 含义 | 主要用途 |
|---|---|---|
| `Initial Oi-Ov distance` | 刚构造 Frenkel pair 时，interstitial O 到 original vacancy site 的距离 | 定义 close / intermediate / separated 初始几何 |
| `Relaxed Oi-Ov distance` | 静态结构弛豫后，interstitial O 到 original vacancy site 的距离 | 判断 relaxation 是否已经改变 Frenkel pair geometry |
| `Production-start Oi-Ov distance` | MD production window 起点时，interstitial O 到 original vacancy site 的距离 | 确定用于 VDOS 和后续分析的真实起始构型 |
| `Minimum MD distance` | MD 过程中该距离达到的最小值 | 判断是否接近 vacancy site 或可能发生 recombination |
| `Final MD distance` | MD 最后一帧中该距离 | 总结有限时间演化结果 |
| `Oi-Ov distance vs time` | MD 中逐帧距离随时间变化 | 主文最关键的动态证据 |

更清楚的表述方式是把 `Relaxed Oi-Ov distance` 写成：

> Relaxed interstitial-to-original-vacancy-site distance

### 5.3 建议同时跟踪两种距离

为了避免 oxygen exchange 造成误判，建议同时提取并比较两种距离：

| 距离 | 定义 | 用途 |
|---|---|---|
| `tracked interstitial O -> original vacancy site distance` | 最初插入的那个 O 原子到原 vacancy site 的距离 | 判断最初的 interstitial O 是否向 vacancy site 移动 |
| `nearest O -> original vacancy site distance` | 任意 O 原子到原 vacancy site 的最短距离 | 判断 vacancy site 是否被任何氧重新占据，更适合作为 recombination 判据 |

如果两条曲线都下降并持续保持在 cutoff 以下，且 local Cu-O coordination 恢复，则 recombination 判断更强。

如果 tracked interstitial O 没有回到 vacancy site，但 nearest O distance 下降到 cutoff 以下，说明可能发生了 oxygen exchange；这种 case 不应简单写成“最初的 interstitial 回去了”，而应写成 vacancy-site recovery 或 oxygen-mediated recombination。

### 5.4 `distance vs time` 的典型解读

| 行为 | `distance vs time` 表现 | 结构解释 | 建议 outcome |
|---|---|---|---|
| Recombined | 距离快速下降到 cutoff 以下，并持续保持 | 原 vacancy site 被重新占据，局部 Cu-O coordination 恢复 | recombined |
| Not recombined / separated | 距离围绕一个有限值上下振动 | Frenkel pair 在 MD 窗口内保持分离 | metastable over sampled window |
| Metastable but flexible | 平均距离保持有限值，但振动幅度较大 | interstitial 周围局部环境较软，热运动更强 | metastable / flexible |
| Migration / rearrangement | 距离出现突然跳变，或从一个平台跳到另一个平台 | O 迁移、oxygen exchange 或局部结构重排 | migrated / rearranged |
| Unclear | 短暂接近 cutoff 后又离开 | 可能只是瞬时热振动接近 | unclear |

对未复合的 cases，不应只说 “no recombination”。更有信息量的比较是：

- 平均 `Oi-Ov distance` 是否保持在有限值；
- 振动幅度是否随 initial separation 或 site environment 改变；
- 是否有突然跳变，提示 migration 或 oxygen exchange；
- close、intermediate、separated cases 的 distance fluctuation 是否存在系统差异。

可以记录以下统计量：

| 统计量 | 用途 |
|---|---|
| mean distance during production window | 表示 metastable separation 的中心值 |
| standard deviation of distance | 表示热振动幅度 |
| peak-to-peak distance range | 表示最大结构波动 |
| time below cutoff | 判断是否为持续接近 |
| event time | 标记 recombination、migration 或跳变发生时间 |

## 6. Recombination 判据

Frenkel defect 的核心问题之一是 vacancy 与 interstitial 是否保持分离、迁移，或复合。因此本节必须分析 recombination。

不要只追踪某一个氧原子 ID。高温或缺陷附近可能发生 oxygen exchange，因此更稳妥的定义应基于 site occupancy 和局部配位恢复。

推荐定义：

> Recombination is assigned when the original vacancy site becomes re-occupied according to a site-occupancy criterion and the local Cu-O coordination around the original vacancy is restored for a sustained time window.

至少需要同时满足两个条件：

| 条件 | 目的 |
|---|---|
| `Oi-Ov distance` 低于预先定义的 cutoff，并持续一段时间 | 避免把瞬时热振动接近误判为 recombination |
| 原 vacancy site 的局部配位环境恢复，或 interstitial oxygen 占据 / 接近原 vacancy lattice site | 确认结构上确实发生 vacancy recovery |

需要记录：

| 信息 | 为什么需要 |
|---|---|
| Minimum Oi-Ov distance during MD | 判断是否发生接近或可能复合 |
| Time below recombination cutoff | 判断是否为持续事件，而不是瞬时接近 |
| Vacancy site occupancy | 判断原 vacancy site 是否被重新占据 |
| Cu-O coordination recovery | 判断局部结构是否恢复 |
| Event time | 记录 recombination、migration 或明显结构重排发生的时间 |
| Notes | 记录 oxygen exchange、异常迁移、局部重构等情况 |

## 7. Close vs separated 分组

不要先人为假设“近的一定复合，远的一定稳定”。应先按照初始 `Oi-Ov distance` 和局部几何环境分组，再比较 MD 演化结果。

建议分组方式：

| 分组 | 定义方式 | 要比较什么 |
|---|---|---|
| Close | 小于第一近邻或局部结构 cutoff | 是否快速复合、是否局部结构恢复 |
| Intermediate | 第一与第二邻近壳层之间 | 是否迁移、是否出现不确定行为 |
| Separated | 大于第二邻近壳层或明显远离 | 是否在 40 ps 内保持分离 |

cutoff 不要凭空设。应从 pristine / relaxed YBa2Cu3O7 的 O-O、Cu-O 距离或 vacancy-site neighbor geometry 推出来。

额外建议记录 same-cell 分类：

| 分组 | 定义方式 | 用途 |
|---|---|---|
| Same-cell Frenkel pair | vacancy 和 interstitial 位于同一个 crystallographic unit cell，cell offset 为 `(0,0,0)` | 检验你观察到的 same-unit-cell recombination |
| Different-cell Frenkel pair | vacancy 和 interstitial 位于不同 crystallographic unit cell | 与 same-cell case 对比有限时间亚稳性 |

如果当前结果显示 same-cell cases 都发生 recombination，可以在正文中谨慎写成：

> In the sampled trajectories, Frenkel pairs initialized within the same crystallographic unit cell, corresponding to short Oi-Ov separation and strong local overlap of the vacancy-interstitial environment, recombined during MD.

如果 separated cases 没有复合，可以谨慎写成：

> By contrast, Frenkel pairs initialized with larger spatial separation remained separated over the 40 ps simulated window, indicating metastability on the sampled timescale.

## 8. 需要做的主要分析

| 分析内容 | 数据来源 | 目的 |
|---|---|---|
| Initial Oi-Ov distance | `initial cif` + defect site labels | 定义 close / intermediate / separated |
| Relaxed Oi-Ov distance | `relaxed cif` | 判断结构弛豫是否已经改变 pair geometry |
| Production-start Oi-Ov distance | `production-start cif` 或 `traj` | 确定 VDOS / production window 的真实起点 |
| Tracked interstitial O to vacancy-site distance vs time | `traj` | 判断最初的 interstitial O 是否向 vacancy site 移动 |
| Nearest O to vacancy-site distance vs time | `traj` | 判断原 vacancy site 是否被任何氧重新占据 |
| Minimum and final Oi-Ov distance | `traj` | 总结有限时间演化结果 |
| Mean distance and fluctuation amplitude | `traj` | 比较未复合 cases 的热振动幅度和局部柔性 |
| Vacancy site occupancy | `traj` + vacancy site definition | 判断原 vacancy site 是否恢复 |
| Local Cu-O coordination | `traj` | 判断 vacancy 周围配位是否恢复 |
| Same-cell vs different-cell classification | `initial cif` + unit-cell indices | 检验 same-unit-cell recombination 观察 |
| Snapshots before / after event | `traj` + selected frames | 支持结构性判断 |
| Event time extraction | `traj` + distance / occupancy / coordination | 记录 recombination 或 migration 发生时间 |
| Outcome classification | all analysis outputs | 形成完整 case summary，避免选择性展示 |

## 9. Methods 中需要交代的内容

Methods 不需要放太多 Frenkel case 结果，但必须交代 defect construction、几何分类和 recombination 判据。

### 9.1 Frenkel pair construction

需要从结构文件和 case inventory 中整理：

| 内容 | 数据来源 | 写入位置 |
|---|---|---|
| YBa2Cu3O7 supercell size，例如 5x5x3 | `initial cif` | Methods |
| Frenkel pair 构造方法 | defect construction notes | Methods |
| vacancy site label 和 site type | case inventory | Methods 或 Supplementary |
| interstitial site label 和 Oi type | case inventory | Methods 或 Supplementary |
| same-cell / different-cell 定义 | unit-cell indices | Methods |
| initial separation 的计算方法 | `initial cif` + analysis script | Methods |
| relaxed 和 production-start geometry 的定义 | `relaxed cif` / `traj` | Methods |

### 9.2 Recombination and outcome classification

需要说明：

- recombination 不只由单个氧原子 ID 判断；
- 使用 site occupancy / vacancy recovery 作为主要结构判据；
- `Oi-Ov distance` 使用 sustained cutoff，而不是瞬时最小距离；
- 同时区分 tracked interstitial-to-vacancy-site distance 和 nearest-O-to-vacancy-site distance；
- local Cu-O coordination 用于辅助确认 vacancy-site recovery；
- 没有复合的 case 只能称为在 simulated window 内 separated 或 metastable；
- 如果存在 oxygen exchange，需要说明如何处理 atom identity 和 site occupancy 的区别。

### 9.3 Methods 中建议做的表

| 表格 | 内容 | 备注 |
|---|---|---|
| Methods Table 1 | Frenkel construction metadata：Case ID、vacancy site、interstitial site、site type、same-cell / different-cell、initial distance | 可以放 Methods；如果太长则放 Supplementary |
| Methods Table 2 | Recombination criteria：distance cutoff、sustained time window、site occupancy criterion、coordination criterion | 可以放 Methods 或 Appendix/Supplementary |

## 10. Result and Discussion Section 2 中建议放的内容

Result and Discussion 第 2 节只放最关键的 Frenkel pair 对比结果。目标是让读者看到：Frenkel pair behavior 由初始几何和局部氧环境控制，而不是随机列举 defect cases。

### 10.1 Main text 推荐图和表

main text 建议选 3-5 个代表性 cases，而不是全部。

| 类别 | 选择目的 | 主文角色 |
|---|---|---|
| Close Frenkel pair | 测试短距离复合 / 局部恢复 | 关键对照 |
| Intermediate pair | 看是否存在过渡行为 | 可选 |
| Far / separated pair | 测试有限时间亚稳性 | 关键对照 |
| 不同 vacancy site | 比较 chain / plane / apical 等局部环境影响 | 若差异明显则主文 |
| 异常 case | 展示意外迁移 / 复合 / 结构重排 | 若结果重要，可作为讨论亮点 |

推荐 main text 内容：

| 图 / 表 | 作用 |
|---|---|
| Fig. 1: representative pristine / vacancy / interstitial / Frenkel structures | 让读者理解缺陷构型 |
| Fig. 2: selected Oi-Ov distance vs time during MD | 比较 close / intermediate / far behavior；展示 recombination、metastability 和 distance fluctuation |
| Table 1: Frenkel case summary | 展示所有 case 的分类和 outcome |
| Fig. 3: before / after snapshots for recombined or metastable representative cases | 支持结构性判断 |

Fig. 2 建议至少包含：

| Panel | 内容 | 目的 |
|---|---|---|
| a | tracked interstitial O to original vacancy site distance vs time | 看最初 interstitial 是否靠近 vacancy |
| b | nearest O to original vacancy site distance vs time | 判断 vacancy site 是否重新被氧占据 |
| c | local Cu-O coordination around original vacancy site vs time | 支持 coordination recovery 判据 |
| d | representative snapshots before / after event | 把 distance 变化和真实结构联系起来 |

如果版面有限，main text 可以只放 panel a/b 的代表性曲线，把所有 coordination time series 放 Supplementary。

### 10.2 Main text 中可以写的分析点

正文建议围绕这几句话展开：

- Frenkel cases 不是随机列举，而是按 `Oi-Ov` 初始距离、same-cell / different-cell 关系和局部 oxygen-site environment 系统比较。
- 完整 case summary 记录了每个 case 的初始条件、弛豫后几何、production-start 几何和 MD 后 outcome。
- same-cell Frenkel pairs 如果都发生 recombination，应表现为 vacancy-site distance 快速下降到 cutoff 以下，并伴随 sustained site occupancy 和 coordination recovery。
- separated Frenkel pairs 如果没有复合，应表现为 distance 围绕有限值热振动；这只能说明它们在 simulated window 内保持分离或呈现 metastability。
- 未复合 cases 也值得比较：distance fluctuation amplitude 的差异可以反映不同局部 oxygen-site environment 下 Frenkel pair 的有限温度柔性。
- 如果曲线出现突然跳变或平台切换，应检查是否发生 migration、oxygen exchange 或局部重构。
- 如果 chain / plane / apical vacancy site 表现不同，应作为 local oxygen-site environment effect 讨论。
- 如果存在 abnormal migration、oxygen exchange 或 local reconstruction，应单独说明，并放入 Supplementary 的完整 case notes。

可直接用于正文的句子模板：

> Recombined cases show a rapid decrease in the vacancy-site distance, followed by sustained occupation of the original vacancy site and recovery of the local Cu-O coordination.

> In contrast, non-recombined cases retain a finite Oi-Ov separation and exhibit thermal oscillations around a metastable geometry over the simulated MD window.

> The oscillation amplitude differs between Frenkel configurations, suggesting that the local oxygen-site environment influences the finite-temperature flexibility of the defect pair.

### 10.3 Result and Discussion 中建议做的表

| 表格 | 内容 | 备注 |
|---|---|---|
| Results Table 1 | 所有 Frenkel case summary：Case ID、site type、same-cell / different-cell、initial distance、mean MD distance、fluctuation amplitude、minimum distance、final distance、recombined、event time、outcome | 建议 main text 放精简版 |
| Results Table 2 | 代表性 cases 的 coordination recovery summary | 如果 Fig. 3 已经足够清楚，可以放 Supplementary |
| Results Table 3 | close / intermediate / separated 分组统计 | 如果 cases 数量足够，可以放 main text 或 Supplementary |

## 11. Appendix / Supplementary 中建议放的内容

Appendix 或 Supplementary 用来放完整 Frenkel case inventory，不占用 main text 篇幅。

### 11.1 Supplementary figures

| 图 | 内容 |
|---|---|
| Fig. S1 | 所有 Frenkel cases 的 initial / relaxed / final structure snapshots |
| Fig. S2 | 所有 Frenkel cases 的 tracked interstitial O to vacancy-site distance vs time |
| Fig. S3 | 所有 Frenkel cases 的 nearest O to vacancy-site distance vs time |
| Fig. S4 | 所有 Frenkel cases 的 vacancy site occupancy time series |
| Fig. S5 | 所有 Frenkel cases 的 local Cu-O coordination time series |
| Fig. S6 | 所有 same-cell Frenkel cases 的 recombination snapshots |
| Fig. S7 | 所有 separated Frenkel cases 的 finite-time metastable snapshots |
| Fig. S8 | 异常 case 的迁移、oxygen exchange 或局部重构说明 |

### 11.2 Supplementary tables

| 表格 | 内容 | 建议位置 |
|---|---|---|
| Table S1 | 完整 Frenkel case inventory：site labels、site type、same-cell / different-cell、cell offset、所有距离和 outcome | Supplementary |
| Table S2 | Recombination criteria 和 cutoff 来源：O-O、Cu-O 或 neighbor-shell analysis | Supplementary |
| Table S3 | 所有 case 的 initial / relaxed / production-start / mean / standard deviation / minimum / final `Oi-Ov` distances | Supplementary |
| Table S4 | 所有 case 的 coordination recovery summary | Supplementary |
| Table S5 | 所有 case 的 event time、event type 和 notes | Supplementary |
| Table S6 | 如果有 independent runs，列出 seed、event status、event time | Supplementary |

## 12. 写作前最低证据清单

开始写 Section 2 正文前，至少需要完成：

- 所有 Frenkel cases 的文件配对表已补全；
- vacancy site label、interstitial site label 和 site type 已补全；
- same-cell / different-cell 分类和 cell offset 已补全；
- initial、relaxed、production-start、minimum 和 final `Oi-Ov distance` 已计算；
- tracked interstitial O to vacancy-site distance vs time 已提取；
- nearest O to vacancy-site distance vs time 已提取；
- 未复合 cases 的 mean distance 和 fluctuation amplitude 已统计；
- close / intermediate / separated 的 cutoff 来源已经说明；
- 每个 case 的 `Oi-Ov distance-t` 已检查；
- 每个 case 的 vacancy site occupancy 已检查；
- 每个 case 的 local Cu-O coordination recovery 已检查；
- recombination、migration、metastable 或 unclear 的 outcome 已分类；
- same-cell recombination 的观察已经由 site occupancy 和 coordination recovery 支持；
- main text 代表性 cases 已选定；
- Supplementary 完整 case inventory 已准备好。

## 13. 后续工作流程

1. 把所有 Frenkel raw files 放入本文件夹。
2. 补全文件配对表和 Frenkel case inventory。
3. 从 pristine / relaxed YBa2Cu3O7 的 O-O、Cu-O 或 neighbor-shell geometry 推导 close / intermediate / separated cutoff。
4. 运行 `extract_frenkel_metrics.py`，提取每个 case 的 `Oi-Ov distance-t`、site occupancy 和 coordination time series。
5. 判断每个 case 的 outcome：recombined、metastable、migrated 或 unclear。
6. 选择 3-5 个 main text representative cases。
7. 使用 `nature-figure` 画 main text 和 supplementary figures。
8. 使用 `academic-research-suite` 或 `nature-writing` 写 Section 2 正文。

使用 `nature-figure` 时，开始前需要先确定作图后端：Python or R。
