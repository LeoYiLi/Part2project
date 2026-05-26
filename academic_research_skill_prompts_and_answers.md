# Academic Research Skill Prompts and Answer Framework

Date: 2026-05-25

Project context: YBa2Cu3O7 5x5x3 supercell, oxygen vacancy, oxygen interstitial, Frenkel defects, MD, MACE/DFT-related relaxation, VDOS, quantum-corrected VDOS, and high-temperature MD.

This file has two uses:

1. Copy the prompt blocks into a new Codex conversation if you want `academic-research-suite` to answer each layer separately.
2. Use the answer framework below as a planning checklist for your manuscript.

Important limitation: Codex cannot open a new conversation window and collect results across windows automatically. The prompt blocks below are ready to paste into a new window.

---

## 0. Recommended Overall Prompt

Use this first if you want one integrated planning answer.

```text
ars-plan

I am planning a computational materials manuscript on YBa2Cu3O7 using a 5x5x3 supercell. I study oxygen vacancy, oxygen interstitial, and Frenkel defect configurations using MD, MACE/DFT-related structural relaxation, VDOS, quantum-corrected VDOS, and high-temperature MD.

Please use academic-research-suite plan mode. Do not draft the manuscript yet. Help me plan the article structure and analysis logic.

My current questions are:
1. What MD quality checks are necessary before using trajectories for VDOS?
2. How should I organize multiple Frenkel defect configurations, including different O interstitial and O vacancy combinations?
3. Which bond length or local structure analyses should I do, and what is their scientific purpose?
4. What conclusions can VDOS and quantum-corrected VDOS support?
5. What is the purpose of high-temperature MD? Should it be used to compare with low-temperature MD, detect recombination, or analyze defect stability?

Please output:
1. Recommended central research question.
2. Recommended article storyline.
3. Results section structure.
4. For each section: claim, required evidence, figure/table, analysis method, reviewer risk.
5. Which analyses are main-text essential and which should go to Supplementary Information.
6. Which claims must be phrased cautiously.
7. A data organization checklist for all 5x5x3 Frenkel defect cases.
```

---

## 1. MD Quality Check

### Corresponding ARS Functions

- `deep-research`
- `literature_strategist_agent`
- `methodology_reviewer_agent`
- `source_verification_agent`
- `argument_builder_agent`

### Copy-Paste Prompt

```text
Please use academic-research-suite deep-research and methodology-review style.

My manuscript uses MD trajectories of YBa2Cu3O7 oxygen vacancy, oxygen interstitial, and Frenkel defect systems to compute VDOS. I need to justify that the MD trajectories are reliable before discussing VDOS.

Please answer:
1. What MD quality checks are standard or defensible in computational materials papers before extracting dynamical observables such as VDOS?
2. Are temperature vs time, potential/total energy vs time, RMSD/displacement vs time, and center-of-mass drift checks appropriate?
3. What additional checks should I consider: equilibration window, velocity autocorrelation convergence, timestep stability, thermostat effect, RDF, bond-length stability, defect stability?
4. Which checks should be in Methods, Results opening, or Supplementary Information?
5. What figures and quantitative criteria should I report?
6. What kind of literature should I cite to support these checks?
7. Please flag overclaims and give cautious wording suitable for a manuscript.
```

### Answer Framework

The purpose of MD quality checks is to establish that the trajectory is usable for structural and vibrational analysis. It does not prove long-term thermodynamic stability.

Recommended checks:

| Check | Purpose | Output |
|---|---|---|
| Temperature vs time | Verify thermostat reaches and maintains target temperature | T-t plot |
| Potential/total energy vs time | Check equilibration and absence of abrupt structural collapse | E-t plot |
| RMSD or displacement vs time | Detect global structural drift or instability | RMSD-t plot |
| Center-of-mass drift | Prevent translational drift contaminating VACF/VDOS | Text or small table |
| Defect descriptor vs time | Confirm whether defect remains separated, migrates, or recombines | Oi-Ov distance / coordination |
| VACF block convergence | Check VDOS peaks are not dominated by one trajectory segment | VDOS from 3 blocks |
| RDF or bond-length distribution | Confirm local structure remains physically plausible | RDF / histogram |

Suggested manuscript role:

- Methods: timestep, thermostat, equilibration protocol, production window, COM removal, VACF/VDOS calculation.
- Results opening: representative T-t and E-t stability, defect descriptor vs time.
- Supplementary: full RMSD, COM drift, all trajectories, block convergence.

Cautious wording:

> The trajectories remain thermally equilibrated and structurally stable over the simulated production window, supporting their use for finite-temperature VDOS analysis.

Avoid:

> The defects are permanently stable.

Use instead:

> The defect configuration remains metastable on the simulated timescale.

---

## 2. Frenkel Defect Dynamics and 5x5x3 Supercell Organization

### Corresponding ARS Functions

- `structure_architect_agent`
- `argument_builder_agent`
- `experiment-agent`
- `state_tracker_agent`

### Copy-Paste Prompt

```text
Please use academic-research-suite structure_architect_agent and argument_builder_agent.

My supercell is 5x5x3. I simulated multiple Frenkel defect configurations in YBa2Cu3O7, with different oxygen interstitial positions and oxygen vacancy positions. I will organize all raw data later.

Please help me design how to present these configurations in the manuscript.

Questions:
1. Should all Frenkel cases be shown in the main text, or should I select representative cases?
2. What information should be recorded for each case?
3. Should I analyze recombination? If yes, how should recombination be defined?
4. What figures should be main text and what should be supplementary?
5. How should I compare close Frenkel pairs and more separated Frenkel pairs?
6. What cautious claims can be made from a 40 ps or similar MD window?
7. Please design a table template for all Frenkel cases.
```

### Answer Framework

The goal is to turn many defect simulations into a controlled comparison rather than a collection of unrelated cases.

Minimum data table:

| Case ID | Vacancy site | Interstitial site | Initial Oi-Ov distance | Relaxed Oi-Ov distance | MD final Oi-Ov distance | Recombined? | Notes |
|---|---|---|---:|---:|---:|---|---|
| F1 | O site label | Oi site label |  |  |  | yes/no/unclear |  |
| F2 | O site label | Oi site label |  |  |  | yes/no/unclear |  |
| F3 | O site label | Oi site label |  |  |  | yes/no/unclear |  |

Main figure recommendation:

1. Representative structures: pristine, vacancy, interstitial, Frenkel.
2. Vacancy-interstitial distance vs time for Frenkel 1/2/3.
3. A summary table of recombination/metastability outcomes.

Possible central claim:

> Close Frenkel pairs recombine rapidly, whereas more separated configurations remain metastable on the simulated timescale.

Important wording:

- Say "metastable on the simulated timescale", not "stable".
- Say "no recombination was observed within 40 ps", not "will not recombine".
- Do not estimate recombination rates unless you have multiple independent seeds and enough events.

Definition of recombination:

Use a structural criterion, for example:

- O interstitial occupies or approaches the vacancy site within a defined cutoff.
- Original vacancy coordination is restored.
- Oi-Ov distance falls below a threshold and remains low for a sustained time window.

The cutoff should be justified using nearest-neighbor O-O or Cu-O distances in pristine/relaxed structures.

---

## 3. Bond Length and Local Structure Analysis

### Corresponding ARS Functions

- `methodology_reviewer_agent`
- `visualization_agent`
- `argument_builder_agent`

### Copy-Paste Prompt

```text
Please use academic-research-suite methodology_reviewer_agent and visualization_agent.

I do not fully understand which bond length analyses are meaningful for YBa2Cu3O7 oxygen defects. I want bond-length analysis to support the connection between defect structure and VDOS.

Please answer:
1. For oxygen vacancy, which local bonds should I track?
2. For oxygen interstitial, which local bonds should I track?
3. For Frenkel defects, should I track vacancy-neighbor bonds, interstitial-neighbor bonds, and Oi-Ov distance?
4. Should I distinguish Cu-O plane, Cu-O chain, apical O, Ba-O, and Y-O environments?
5. Should I plot time series, histogram, box plot, RDF, or before/after comparison?
6. What is the scientific purpose of each bond-length figure?
7. Which local-structure analyses are main-text essential and which are supplementary?
8. Please give a minimal but convincing local-structure analysis package.
```

### Answer Framework

Bond-length analysis is not decorative. Its purpose is to show that VDOS changes arise from local bonding distortion and coordination changes.

Recommended analyses:

| Analysis | Purpose | Suggested figure |
|---|---|---|
| Cu-O bond-length histogram | Main local bonding distortion in Cu-O network | Pristine vs vacancy vs interstitial vs Frenkel |
| Defect-neighbor Cu-O time series | Check local structural stability around vacancy/interstitial | Time series or box plot |
| Oi-neighbor distance | Identify whether interstitial forms a stable local environment | Time series + histogram |
| Oi-Ov distance | Track Frenkel recombination or separation | Time series |
| RDF: Cu-O and O-O | Broader local structure comparison | Supplementary or main if strong |
| Local displacement vs distance from defect | Show spatial extent of defect perturbation | Main or Supplementary |

Most useful main-text local-structure figure:

> Cu-O bond-length distribution for pristine, vacancy, interstitial, and representative Frenkel systems.

Why it matters:

- Vacancy can under-coordinate nearby Cu or alter Cu-O network.
- Interstitial oxygen can create new local bonding environments.
- Frenkel defects can combine vacancy-like and interstitial-like distortions.
- These local distortions can rationalize shifts, broadening, or new features in oxygen-dominated VDOS.

Avoid trying to plot all possible bonds. Select bonds by physical role:

1. Nearest neighbors to the vacancy.
2. Nearest neighbors to the interstitial.
3. Cu-O chain/plane/apical distinctions if your structure labels allow this.
4. Bulk-like reference atoms far from the defect.

---

## 4. VDOS and Quantum Correction

### Corresponding ARS Functions

- `argument_builder_agent`
- `visualization_agent`
- `methodology_reviewer_agent`
- `citation_compliance_agent`

### Copy-Paste Prompt

```text
Please use academic-research-suite argument_builder_agent and visualization_agent.

I have computed VDOS for pristine YBa2Cu3O7 and defect systems, including vacancy, interstitial, and Frenkel defects. I also applied quantum correction to the VDOS.

Please help me design the VDOS Results section.

Questions:
1. What should classical VDOS and quantum-corrected VDOS each be used to argue?
2. What conclusions can be supported by total VDOS, element-resolved VDOS, oxygen-only VDOS, difference VDOS, and local defect-neighbor VDOS?
3. How should I interpret peak shift, broadening, low-frequency enhancement, high-frequency oxygen-mode changes, or spectral weight redistribution?
4. How should I connect VDOS features to local bond-length distortion?
5. What figures should be in the main text?
6. What should be treated cautiously or left as supplementary?
7. Please give a claim-by-figure map for this section.
```

### Answer Framework

VDOS should support a mechanistic claim:

> Oxygen defects modify the local bonding environment of YBa2Cu3O7, producing oxygen-dominated changes in vibrational density of states.

Recommended figure order:

| Figure | Claim |
|---|---|
| Total VDOS: pristine vs defect | Defects modify the overall vibrational spectrum |
| O-only VDOS overlay | Main spectral changes are associated with oxygen sublattice |
| O/Cu/Ba/Y projected VDOS | Element-specific contribution to spectral changes |
| Defect - pristine difference VDOS | Direct visualization of defect-induced spectral redistribution |
| Local oxygen VDOS | Interstitial O, vacancy-neighbor O, and bulk-like O have distinct vibrational fingerprints |
| Block convergence VDOS | Spectral features are robust across trajectory segments |

If possible, split oxygen sites:

- chain O
- plane O
- apical O
- interstitial O
- vacancy-neighbor O
- bulk-like O

Potential interpretations:

| Observation | Possible interpretation |
|---|---|
| Peak broadening | Larger structural disorder or anharmonic sampling |
| Peak shift to lower frequency | Weaker/longer local bonds or softened local modes |
| Peak shift to higher frequency | Shorter/stronger local bonds or constrained local environment |
| Low-frequency enhancement | Defect-induced soft modes, local rattling, or increased disorder |
| Oxygen-mode changes | Defect perturbation mainly affects oxygen sublattice |

Quantum correction:

- Use it to compare classical MD-derived spectra with quantum occupation-weighted vibrational response.
- Do not overclaim that quantum correction alone proves quantum phonon physics.
- Report the correction formula and temperature.
- Compare whether qualitative defect trends remain after correction.

Cautious wording:

> The quantum-corrected spectra preserve the main defect-induced trends, indicating that the observed differences are not solely artifacts of classical spectral weighting.

Avoid:

> Quantum correction proves the exact experimental phonon spectrum.

---

## 5. High-Temperature MD and c-axis Analysis

### Corresponding ARS Functions

- `experiment-agent`
- `research_question_agent`
- `argument_builder_agent`
- `methodology_reviewer_agent`

### Copy-Paste Prompt

```text
Please use academic-research-suite experiment-agent and argument_builder_agent.

I am running high-temperature MD for YBa2Cu3O7 oxygen interstitial, vacancy, and Frenkel defect systems. I need to define the scientific purpose of this part.

Please answer:
1. Should high-temperature MD be framed as defect stability, recombination, oxygen migration, temperature-dependent VDOS, or c-axis/lattice response?
2. How should it connect to my lower-temperature MD and VDOS results?
3. What temperature points and trajectory lengths are defensible?
4. What analyses should I perform: Oi-Ov distance, interstitial z-coordinate, c-axis length, volume/strain, MSD, coordination number, bond breaking/forming, high-T VDOS?
5. How should I define recombination or migration?
6. What claims are safe if recombination is observed?
7. What claims are safe if no recombination is observed?
8. Which results should be main text and which should be supplementary?
```

### Answer Framework

High-temperature MD can serve three different manuscript purposes. Pick one main purpose.

| Purpose | What it proves | Analyses |
|---|---|---|
| Defect stability | Whether a defect remains intact over the simulated high-T window | Oi-Ov distance, coordination, snapshots |
| Recombination/migration | Whether interstitial oxygen moves toward vacancy or along c-axis | Oi-Ov distance, Oi z-coordinate, trajectory snapshots |
| Temperature-dependent vibrational response | Whether temperature softens/broadens oxygen modes | high-T VDOS vs 298 K VDOS |
| Lattice response | Whether interstitial defects couple to c-axis or volume strain | c-axis length, volume/strain vs time |

Recommended storyline:

> Low-temperature MD and VDOS establish defect-induced vibrational fingerprints. High-temperature MD then tests whether these defect configurations remain metastable or undergo recombination/migration under stronger thermal perturbation.

Useful analyses:

- c-axis length vs time
- average c-axis change by interstitial site
- volume/strain vs time
- interstitial z-coordinate vs time
- vacancy-interstitial distance vs time
- high-temperature VDOS vs 298 K VDOS
- snapshots before and after any recombination event

Safe claims if recombination occurs:

> The close Frenkel configuration recombines within the simulated high-temperature trajectory, suggesting a thermally accessible local recovery pathway for this geometry.

Safe claims if recombination does not occur:

> No recombination is observed within the simulated high-temperature window, suggesting metastability over the sampled timescale.

Avoid:

- "recombination rate" without multiple seeds/events.
- "diffusion coefficient" from short trajectories unless migration is extensive and sampling is adequate.
- "long-term stability" from tens of ps.

---

## 6. Mapping Your Existing Six-Layer Plan to Manuscript Sections

Recommended Results structure:

1. MACE validation by structural relaxation
   - Claim: MACE-relaxed structures reproduce reasonable local bonding and are suitable for MD setup.
   - Evidence: DFT vs MACE bond lengths or relaxed structure comparison, if available.

2. MD equilibration and defect stability
   - Claim: Production trajectories are thermally equilibrated and structurally usable for VDOS.
   - Evidence: T-t, E-t, RMSD, COM drift, Oi-Ov distance.

3. Defect dynamics in the 5x5x3 supercell
   - Claim: Frenkel behavior depends on initial vacancy-interstitial geometry.
   - Evidence: Oi-Ov distance vs time for Frenkel 1/2/3, recombination snapshots.

4. Local structural distortion around oxygen defects
   - Claim: Oxygen defects perturb local Cu-O and neighboring coordination environments.
   - Evidence: Cu-O bond histograms, Oi-neighbor distances, RDF, local displacement.

5. Defect-induced VDOS changes
   - Claim: Defect configurations alter VDOS, primarily through oxygen-related modes.
   - Evidence: total VDOS, O-only VDOS, element-resolved VDOS, difference spectrum.

6. Local oxygen vibrational fingerprints
   - Claim: Interstitial O, vacancy-neighbor O, and bulk-like O show distinct vibrational signatures.
   - Evidence: local projected VDOS, site-resolved oxygen VDOS.

7. Temperature and site dependence, if high-T data are ready
   - Claim: High temperature reveals whether defect configurations are metastable, recombine, migrate, or alter c-axis response.
   - Evidence: high-T Oi-Ov distance, Oi z-coordinate, c-axis length, high-T VDOS.

---

## 7. Main-Text vs Supplementary Recommendation

### Main Text

- Representative structure images.
- T-t and E-t, possibly compact.
- Frenkel Oi-Ov distance vs time for key cases.
- Cu-O bond-length histogram.
- Total VDOS overlay.
- O-only or element-resolved VDOS.
- Difference VDOS.
- Local oxygen VDOS.
- High-T c-axis or recombination figure only if results are clean and central.

### Supplementary

- Full RMSD and COM drift.
- All individual trajectory stability plots.
- All Frenkel case tables.
- All RDFs.
- Full bond-length distributions for Ba-O, Y-O, Cu-O separated by site.
- Block convergence VDOS.
- Additional temperature points.
- Raw or extended projected VDOS.

---

## 8. Claims You Can and Cannot Make Yet

### Stronger, defensible claims

- The simulated trajectories are stable over the production window.
- Oxygen defects induce local structural distortion.
- VDOS changes are concentrated in oxygen-related modes.
- Different Frenkel geometries show different finite-temperature behavior.
- Some Frenkel pairs recombine or remain separated within the simulated window, depending on geometry.

### Cautious claims only

- Defect lifetime.
- Recombination rate.
- Diffusion coefficient.
- Long-term thermal stability.
- Free energy stability.
- Experimental phonon agreement.

### Better wording

Use:

> finite-temperature behavior over the simulated timescale

Instead of:

> long-term stability

Use:

> no recombination was observed within the sampled trajectory

Instead of:

> the defect does not recombine

Use:

> defect-induced vibrational fingerprints

Instead of:

> definitive phonon signatures

---

## 9. Next Data Organization Checklist

For each system:

- System name: pristine / vacancy / interstitial / Frenkel ID.
- Supercell: 5x5x3.
- Defect site labels.
- Initial structure file.
- Relaxed structure file.
- MD trajectory file.
- Temperature.
- Timestep.
- Thermostat.
- Total MD length.
- Equilibration window removed.
- Production window used for VDOS.
- Temperature mean and fluctuation.
- Energy drift or plateau description.
- RMSD/displacement.
- COM drift result.
- Oi-Ov distance vs time, if relevant.
- Recombined? yes/no/unclear.
- Local bond statistics.
- RDF files.
- VDOS files.
- Quantum correction formula and temperature.
- Block convergence status.

For each figure:

- Figure ID.
- Data source.
- Claim supported.
- Main text or Supplementary.
- Missing analysis.
- Risk or reviewer concern.

---

## 10. Suggested Immediate Next Prompt

After organizing your data table, use:

```text
ars-outline

I have organized my YBa2Cu3O7 5x5x3 defect MD data using the table below: [paste table].

Please generate a detailed manuscript outline and evidence map. For each section, provide:
1. section title;
2. main claim;
3. required figure or table;
4. exact analysis needed;
5. expected result pattern;
6. cautious wording;
7. likely reviewer concern;
8. whether the item belongs in main text or Supplementary Information.
```

