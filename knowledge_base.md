---
name: Part2 Project - YBCO Defect Study
description: Comprehensive knowledge base for the YBCO superconductor defect simulation project using ASE+MACE
type: project
originSessionId: bd8023f0-7957-4564-ab1c-00a46dd930bd
---
# Part 2 Project: YBCO Defect Simulation Knowledge Base

---

## Current Direction Update (2026-05-24)

The previous four-stage plan based on NEB, quenching, MSD diffusion extraction,
and VDOS-derived formation free energy is no longer the main project route.

The current supervisor-guided direction is:

- run higher-temperature MD, especially around 600 K and 800 K;
- compare oxygen interstitials placed at different candidate sites;
- focus on whether oxygen interstitials produce c-axis elongation in REBCO;
- use c-axis-flexible MD as the main workflow, with fixed-cell NVT only as a
  control;
- do not prioritize quenching, diffusion coefficients, vibrational entropy, or
  defect formation free energy at this stage.

The active workflow is documented in:

- `high_temperature_interstitial_README.md`

Older VDOS and free-energy notes are retained as background only, not as the
current execution plan.

---

## 1. Research Context & Motivation

### 1.1 One-sentence Summary
**This MSc project uses a MACE machine-learning interatomic potential, benchmarked against DFT, to study the room-temperature stability of irradiation-induced defects in REBCO through molecular dynamics simulations.**

### 1.2 Background
In controlled nuclear fusion applications, high-Tc superconductor REBCO (specifically YBa₂Cu₃O₇) is exposed to neutron/ion irradiation, causing performance changes in Tc and Jc. The project addresses:
- What defects form after irradiation in REBCO
- Whether these defects are stable at room temperature (after removal from cryogenic environment)
- How defect presence, migration, and recombination relate to local bonding, vibration, and structure
- Whether MACE can replace expensive DFT for large-scale MD studies

**All user questions are about this MSc project — no other topics expected. All files under `/home/leo/Part2project/Part2project/`.**

### 1.3 Key Methodological Positioning
- **DFT is NOT the main method — it serves as benchmark only**
- **MACE ML potential is the core tool**, enabling larger systems and longer timescales
- The project evolved from static structure optimization → dynamic stability via MD → vibrational analysis (VDOS)

### 1.4 Literature Context
- At low neutron fluences, Tc decreases nearly linearly
- Jc increases correlate with flux pinning regions
- Nicholls et al. used He⁺ irradiation as a controllable proxy for neutron irradiation
- Combined HERFD-XAS, XRD, STEM with DFT spectral simulation to study irradiation effects

---

## 2. Crystal Structure: YBa₂Cu₃O₇

- Orthorhombic, space group Pmmm (#47)
- Lattice parameters (unit cell): a=3.819 Å, b=3.885 Å, c=11.68 Å
- Key oxygen sites:
  - **O1** (chain oxygen): (0, 0.5, 0) in Cu1 chains along b-axis, occupancy ~0.971
  - **O2** (apical oxygen): (x, 0, z) bridging CuO₂ planes and Cu1 chains
  - **O3** (plane oxygen): (0, y, z) in CuO₂ planes
  - **O4** (chain oxygen): same as O1 in some notations
- Cu1: chain copper at (0, 0, 0)
- Cu2: plane copper at (0, 0, ~0.118)

### Why 5×5×3 Supercell
- Single unit cell defects under PBC create artificially high defect concentrations
- 5×5×3 (975 atoms) is large enough to reduce defect-image interactions
- Better represents low-concentration irradiation defects
- Balances computational cost vs physical realism
- Cell: a=19.095 Å, b=19.425 Å, c=35.04 Å
- Formula: Y75Ba150Cu225O525 (perfect)

---

## 3. Defect Structures

### 3.1 Defect Types Studied
- **O vacancy** — removed from chain O1 site
- **O interstitial** — inserted at specific crystallographic site
- **Frenkel defect** — vacancy + interstitial pair (net atom count unchanged)
- Different Frenkel pair separations to study recombination vs stability

### 3.2 Key Concept: Frenkel Pair Distance & Recombination Results
A central question in the project:
- If vacancy and interstitial are far apart, do they still constitute a Frenkel pair?
- Close pairs more easily recombine; distant pairs behave like independent defects

**MD Results (current as of 2026-04-22):**
- **Frenkel 1** (0.4, 0.4, 0.5): **RECOMBINES** during 298K MD
- **Frenkel 2** (0.4, 0.4, 0.167): **does NOT recombine** — remains separated
- **Frenkel 3** (0.2, 0.4, 0.5): **does NOT recombine** — remains separated

### 3.3 CIF Files

#### Perfect Structures
| File | Description |
|------|-------------|
| `YBa2Cu3O7.cif` | Pristine unit cell |
| `YBa2Cu3O7_553.cif` | Perfect 5×5×3 supercell (Y75Ba150Cu225O525) |
| `Vacancy_553.cif` | VESTA export with full symmetry (619 symops), O4 occupancy=0.971 |

#### Relaxed Structures (MACE)
| File | Description |
|------|-------------|
| `YBa2Cu3O7_relaxed_MACE.cif` | Pristine unit cell relaxed with MACE |
| `YBa2Cu3O7_relaxed_sym_MACE.cif` | Relaxed with symmetry constraints |
| `YBa2Cu3O7_Int1_MACE.cif` | O-interstitial structure relaxed |
| `YBa2Cu3O7_Ovac1_relaxed_MACE.cif` | O-vacancy structure relaxed |

#### Oxygen Vacancy
| File | Description |
|------|-------------|
| `YBa2Cu3O7_Ovac1.cif` | Unit cell with O1 vacancy |
| `YBa2Cu3O7_Ovac1_md.cif` | After MD simulation |

#### Oxygen Interstitial
| File | Description |
|------|-------------|
| `YBa2Cu3O7_Int1.cif` | Unit cell with O interstitial |
| `Interstitial_553.cif` | 553 supercell, O added at (0.4, 0.4, 0.167) between Cu2 pair, **O526 total** |

#### Frenkel Defects (Vacancy + Interstitial pair)
| File | Formula | Interstitial O Position | Description |
|------|---------|------------------------|-------------|
| `YBa2Cu3O7_553_Frenkel.cif` | Y75Ba150Cu225O**525** | **(0.4, 0.4, 0.5)** | Frenkel 1: interstitial at c=0.5 (BaO spacer region, midpoint of two Cu2) |
| `YBa2Cu3O7_553_Frenkel2.cif` | Y75Ba150Cu225O**525** | **(0.4, 0.4, 0.167)** | Frenkel 2: interstitial at c~1/6 (Y layer, between CuO₂ planes) |
| `YBa2Cu3O7_553_Frenkel3.cif` | Y75Ba150Cu225O**525** | **(0.2, 0.4, 0.5)** | Frenkel 3: interstitial at c=0.5, different ab-plane position |

All three have net O count unchanged (525): one O removed from chain, one O inserted at interstitial site.

#### MD-Relaxed Structures
| File | Description |
|------|-------------|
| `Frankel_md.cif` | Frenkel defect after 298K MD (all atoms displaced from ideal) |
| `YBa2Cu3O7_553_md.cif` | Perfect supercell after MD |
| `YBa2Cu3O7_md.cif` | Unit cell after MD |

---

## 4. DFT vs MACE Benchmark

### 4.1 Purpose
Validate MACE reliability before using it for production MD. DFT is the accuracy reference.

### 4.2 What Was Compared
- Geometric structures after relaxation
- Bond lengths (most important metric — sensitive to local defect environment)
- Both pristine and defect configurations

### 4.3 Bond Length Benchmark Approach
- x-axis: different bond types (Cu-O, Ba-O, Y-O, etc.)
- y-axis: bond length
- Two colors: DFT (blue) vs MACE (orange)
- Sorted by length for visual clarity
- Conclusion: MACE and DFT are close enough → MACE is trustworthy for MD

### 4.4 Observed Deviations
- Some deviations exist near defect sites
- Errors concentrate around defect neighbours — physically meaningful, not random
- This validates that MACE captures local structural distortions

### 4.5 DFT Parameters (for reference)
- Pseudo potential, k-point sampling, cutoff energy
- These determine DFT accuracy and explain why DFT is expensive → motivates MACE

---

## 5. Molecular Dynamics

### 5.1 Research Rationale for MD
- REBCO has been irradiated and contains defects
- Material is removed from cryogenic environment, now at room temperature
- Question: are defects stable at 298K? Do they migrate or recombine?
- NVT ensemble appropriate: fixed volume, studying thermal equilibrium behavior

### 5.2 Why NVT + Langevin + 298K
- **NVT**: fixed particles, volume, temperature — studying thermal stability, not pressure effects
- **298K**: room temperature — the post-irradiation storage/testing condition
- **Langevin thermostat**: adds damping + random force, maintains temperature, simulates heat bath coupling

### 5.3 Primary Production Run (MD.py)
- Thermostat: Langevin
- Timestep: 0.5 fs
- Temperature: 298 K
- Friction: 0.0005 / fs (light coupling for good dynamics)
- Steps: 10,000 (= 5 ps)
- Trajectory write interval: every step
- Supercell: 5×5×3
- Cell: fixed (NVT)

### 5.4 Annealing Run (MD_anealing.py)
- Timestep: 1.0 fs
- Temperatures: 20K, 77K, 300K
- Steps per temperature: 50,000 (= 50 ps)
- Friction: 0.01 / fs (stronger coupling for equilibration)
- Trajectory write interval: every 100 steps

### 5.5 Trajectory Files

| File | System | Temperature | Notes |
|------|--------|-------------|-------|
| `YBa2Cu3O7_md_298K.traj` | Perfect 553 | 298K | Main pristine reference |
| `Frankel_md_298K.traj` | Frenkel 1 | 298K | Frenkel defect (z=0.5), **recombines** |
| `Frankel_md_77K.traj` | Frenkel 1 | 77K | Frenkel defect at liquid N₂ temp |
| `Frankel2_md_298K.traj` | Frenkel 2 | 298K | Frenkel defect (z=0.167), **does not recombine** |
| `YBa2Cu3O7_553_Frenkel3.traj` | Frenkel 3 | 298K | Frenkel defect (0.2,0.4,0.5), **does not recombine** |
| `Vacancy_md_298K.traj` | Vacancy | 298K | O-vacancy 553 supercell |
| `Interstitial_md_298K.traj` | Interstitial | 298K | O-interstitial 553 supercell |
| `yba553_md_30K.traj` | Perfect 553 | 30K | Low temperature reference |
| `yba_md_77K.traj` | Perfect 553 | 77K | Intermediate temperature |
| `yba_md_Ovac1_30K.traj` | O-vacancy | 30K | Low temp vacancy |

---

## 6. Analysis Methods

### 6.1 VDOS (Vibrational Density of States)

#### Methodology
1. Read MD trajectory (skip equilibration frames, typically first 100)
2. Extract velocities from all frames
3. Optionally remove center-of-mass drift (V -= V.mean per frame)
4. Compute velocity autocorrelation function (VACF):
   - **Direct method**: time-origin averaging `VACF(t) = <v(0)·v(t)>`
   - **FFT method** (preferred, faster): `VACF = IFFT(|FFT(v)|²)` with zero-padding to power of 2
5. Apply Hanning window (optional, reduces spectral leakage)
6. FFT of VACF → power spectrum = VDOS
7. Convert frequency axis to THz
8. Element-resolved: repeat for subsets of atoms (O, Cu, Ba, Y)
9. Verify: sum of element VDOS ~ total VDOS

#### Physical Interpretation
- Peak positions correspond to specific vibrational modes
- Peak broadening → increased local disorder
- Smoothing of spectrum → more dispersed vibrational modes
- Defects modify local bonding → change vibrational spectrum
- If Frenkel pair recombines, VDOS changes compared to separated pair
- Comparing pristine vs defect VDOS reveals defect-induced phonon changes

### 6.2 Bond Length Analysis
- Cu-O: 1.8–2.5 Å (or up to 3.0 Å for extended analysis)
- Ba-O: 2.5–3.2 Å (or up to 3.3 Å)
- Y-O: 2.2–2.8 Å (or up to 3.0 Å)
- O-O: up to 3.2 Å
- Used for DFT/MACE benchmark and defect local structure characterization

### 6.3 OVITO Visualization
- Used to observe defect movement in trajectories
- Challenges with large supercells: internal defects hard to see
- Solution: cut out layers containing defects (`cut layer.py`)
- Boundary atom flickering/teleporting is a PBC display artifact, not physical
- Key insight: "stable" does not mean "motionless" — atoms vibrate at 298K; stability means defect identity and spatial relationships are preserved

### 6.4 Temperature & Energy Monitoring
- MD logs track T and E over time
- Used to verify equilibration and reasonable fluctuations
- NVT + Langevin should show T fluctuating around target with reasonable amplitude

---

## 7. Python Scripts

### Defect Creation
| Script | Function |
|--------|----------|
| `make_O1_vac.py` | Creates O1 vacancy; finds nearest O to target fractional coord, deletes it |
| `add_O_interstitial.py` | Inserts O at specified fractional position |
| `# make_frenkel.py` | Creates Frenkel pair: interstitial at midpoint of two Cu2 atoms; vacancy creation commented out. Outputs `Interstitial_553.cif` |

### Structure Relaxation
| Script | Function |
|--------|----------|
| `ASE_0.py` | BFGS + UnitCellFilter (fixed angles) relaxation with MACE, fmax=0.005 eV/Å. Also determines space group with spglib. |
| `ASE/ase_mace_symmetry.py` | Relaxation with FixSymmetry constraint |

### Molecular Dynamics
| Script | Function |
|--------|----------|
| `MD.py` | Main MD: Langevin 298K, 0.5fs step, friction=0.0005/fs, 10000 steps, CUDA |
| `MD_anealing.py` | Multi-temp annealing: 20K→77K→300K, 1.0fs, 50000 steps/temp |
| `plot_md.py` | Plots temperature and energy vs time from MD logs |

### VDOS Analysis
| Script | Function |
|--------|----------|
| `VDOS.py` | Basic VDOS: direct VACF with time-origin averaging → FFT |
| `VDOS_Oxygen.py` | O-specific VDOS from Frenkel trajectory |
| `import numpy as np.py` | **Advanced VDOS**: FFT-based autocorrelation, element-resolved, drift removal, zero-padding |
| `VDOS_tet.py` | Plots external DOS data (cm⁻¹ → THz conversion) |

### Bond Analysis
| Script | Function |
|--------|----------|
| `bond length.py` | Extracts bonds preserving CIF labels, pair cutoffs |
| `simple_bond_analysis.py` | YBCO-specific bond statistics |
| `extract_bond_lengths.py` | General utility with CLI |

### Utilities
| Script | Function |
|--------|----------|
| `cut layer.py` | Extracts supercell sub-layers from trajectory (XYZ output) for OVITO |
| `check.py` | Verifies trajectory has momenta/velocity data |

---

## 8. MACE Calculator Setup

- Model: `MACE-MP-0a.model` (foundation model, pre-trained on Materials Project)
- Linux/CUDA: `/home/leo/Part2project/Part2project/MACE-MP-0a.model` (device="cuda")
- macOS/CPU: `/Users/liyi/Desktop/Part 2 project/MACE-MP-0a.model` (device="cpu")
- Originally ran on CPU (slow), migrated to GPU machine for production runs

---

## 9. Current Open Issue: MACE-MD VDOS vs DFT Perturbation DOS

A collaborator (same group) computed the phonon DOS of pristine YBa₂Cu₃O₇ using DFT perturbation theory (DFPT / finite displacement). This is saved as `graph/VDOS_perturbation.png`. The MACE-MD VDOS (from velocity autocorrelation) shows **significant differences** compared to the DFT perturbation result. Possible causes to investigate:
- MD VDOS is inherently broadened (finite temperature, anharmonicity, finite trajectory length)
- DFT perturbation DOS is 0K harmonic — sharper peaks expected
- MACE potential accuracy limitations for YBCO phonons
- Convergence issues (trajectory length, equilibration, sampling)
- Different normalization or units

This discrepancy is a key issue in the current stage of the project.

---

## 10. Output Plots (graph/ folder)

| File | Description |
|------|-------------|
| `VDOS_perturbation.png` | DFT perturbation phonon DOS of pristine YBa₂Cu₃O₇ (from collaborator) |
| `O only frankel2.png` | O-element VDOS for Frenkel 2 |
| `Ba only frankel 2.png` | Ba-element VDOS for Frenkel 2 |
| `Cu only frankel 2.png` | Cu-element VDOS for Frenkel 2 |
| `Y only frankel 2.png` | Y-element VDOS for Frenkel 2 |
| `O vs total frankel.png` | O vs total VDOS for Frenkel 1 |
| `O vs total frankel 2.png` | O vs total VDOS for Frenkel 2 |
| `O vs total pristine.png` | O vs total VDOS for pristine |
| `O vs total Interstitial.png` | O vs total VDOS for interstitial |
| `total vs sum frankel2.png` | Total vs sum-of-elements VDOS for Frenkel 2 |
| `strict sum frankel 2.png` | Strict sum of element VDOS for Frenkel 2 |

### Frenkel 3 additional files
| File | Description |
|------|-------------|
| `YBa2Cu3O7_553_Frenkel3.traj` | MD trajectory |
| `YBa2Cu3O7_553_Frenkel3.log` | MD log |
| `YBa2Cu3O7_553_Frenkel3.xyz` | Full supercell XYZ |
| `YBa2Cu3O7_553_Frenkel3_113.traj/.xyz` | Cut layer 1×1×3 around defect |
| `YBa2Cu3O7_553_Frenkel3_331.traj/.xyz` | Cut layer 3×3×1 around defect |

---

## 11. Key Open Physical Questions

1. How do different O defects (vacancy vs interstitial vs Frenkel) modify the phonon spectrum?
2. Where in the VDOS do defect-induced modes appear?
3. Which element's vibrations are most affected by each defect type?
4. How does the Frenkel interstitial position (z=0.5 vs z=0.167) change the vibrational signature?
5. Does Frenkel pair separation affect recombination probability? Is there a critical distance?
6. Temperature dependence of defect dynamics (20K, 77K, 298K, 300K)
7. Are defects stable at room temperature, or do they migrate/recombine?

---

## 10. Project Evolution & Milestones

### Phase 1: Structure Building & DFT Benchmark
- Built pristine REBCO structure and supercells
- Introduced vacancy, interstitial, Frenkel defects
- Compared DFT vs MACE relaxation results via bond lengths

### Phase 2: MD Simulations
- Established NVT + 298K + Langevin protocol
- Migrated from CPU to GPU for performance
- Ran MD for pristine, vacancy, interstitial, and two Frenkel configurations
- Observed: Frenkel 1 (z=0.5) can recombine; Frenkel 2 (z=0.167) stays separated

### Phase 3: Trajectory Analysis & VDOS
- Developed VACF → VDOS pipeline
- Element-resolved VDOS analysis
- Used OVITO for visual defect tracking
- Cut layers from large supercell for focused visualization
- Compared MACE-MD VDOS with collaborator's DFT phonon DOS (`VDOS_perturbation.png`)
- **Key finding: significant discrepancy between MACE-MD VDOS and DFT perturbation phonon DOS** — this is an open issue under investigation

### Phase 4: Writing & Presentation (current)
- PPT and presentation scripts refined
- Entering thesis/paper writing stage

---

## 11. Planned/Proposed Paper Structure

1. **Introduction**: REBCO in fusion, irradiation defects, DFT limitations, ML potential motivation
2. **Methodology**: Structure building, defect construction, DFT benchmark setup, MACE, MD parameters (NVT/298K/Langevin/553), analysis methods (bond length, OVITO, VACF, VDOS)
3. **Results**: MACE vs DFT benchmark, defect structure stability, Frenkel pair separation effects, MD trajectory observations, VDOS/VACF vibrational analysis
4. **Discussion**: Room-temperature defect stability, recombination conditions, local disorder and vibrational signatures, ML-MD advantages and limitations
5. **Conclusion**: Main findings, method value, future work

---

## 12. Reusable Text Snippets

### English
> The simulations were performed in the NVT ensemble at 298 K using a Langevin thermostat. A 5×5×3 REBCO supercell was employed to better represent low-concentration irradiation defects while reducing artificial interactions between periodic images.

### Chinese
> 模拟采用 NVT 系综，在 298 K 下使用 Langevin thermostat。结构模型选用 5×5×3 的 REBCO supercell，以更真实地表征低浓度辐照缺陷，并尽可能减小周期镜像之间的相互作用。

---

## 13. File Naming Conventions

- `*_553.cif`: 5×5×3 supercell
- `*_MACE.cif`: MACE-relaxed
- `*_Ovac1*`: O1 vacancy
- `*_Int1*`: O interstitial
- `*_Frenkel*` / `*Frankel*`: Frenkel defect (both spellings used)
- `*_md.cif`: structure after MD
- `*_md_XXXK.traj`: MD trajectory at XXX Kelvin

---

## 14. Future Directions

1. **Systematic Frenkel pair separation study**: vary vacancy-interstitial distance, find critical recombination distance
2. **Standardize MD analysis per case**: trajectory snapshots, T/E evolution, bond-length evolution, VACF, VDOS, final defect state
3. **Systematize benchmark**: pristine + vacancy + interstitial + Frenkel, error statistics focused on defect-neighbour bonds
4. **Unified narrative for thesis**: REBCO → irradiation defects → room-temp stability → MACE (not DFT) → structure + dynamics + vibrational analysis
