# High-Temperature O-Interstitial MD Plan

## Updated project direction

The previous staged plan is no longer the main route. The useful next step is to
run molecular dynamics at temperatures above room temperature and compare how
oxygen interstitials at different sites affect the YBCO lattice, especially the
c-axis length.

This follows the new supervisor guidance:

- Existing MD is mostly at room temperature, so the temperature dependence is
  not yet well sampled.
- Run higher-temperature MD, for example 600 K and 800 K.
- REBCO c-axis elongation may be linked to oxygen interstitials.
- Oxygen interstitials can occupy different candidate sites, so the site
  dependence should be tested directly.
- There is no need at this stage to calculate quenching energies, diffusion
  coefficients, vibrational entropy, or defect formation free energies.

## Main scientific question

Which oxygen interstitial site gives the strongest and most persistent c-axis
expansion in YBa2Cu3O7, and how does this response depend on temperature?

## Proposed workflow

1. Start from the 5x5x3 pristine supercell, `YBa2Cu3O7_553.cif`.
2. Build several O-interstitial structures by adding one oxygen atom at
   different candidate fractional coordinates.
3. For each interstitial site, run MD at several temperatures:
   - 298 K for comparison with the previous room-temperature work
   - 600 K
   - 800 K
   - optionally 1000 K if the structure remains physically meaningful
4. Use an NPT-style run that allows the c axis to change while keeping a and b
   fixed, because the target observable is c-axis elongation.
5. Track:
   - final and time-averaged cell lengths a, b, c
   - percent change in c relative to pristine
   - volume change
   - total energy and temperature stability
   - whether the interstitial remains near its starting site or migrates

## Why not fixed-volume NVT only?

The earlier 298 K simulations used fixed-cell NVT, which is useful for comparing
defect stability under identical volume. However, fixed-volume NVT cannot show a
real c-axis elongation, because the cell is not allowed to relax. For the new
question, use c-axis NPT as the main calculation and keep NVT only as a control.

## Suggested interstitial sites

The exact sites should be refined based on structural visualization and
literature, but the current project already uses several candidate positions:

| Label | Fractional coordinate | Notes |
| --- | --- | --- |
| site_y_layer | (0.4, 0.4, 0.167) | Existing interstitial/Frenkel2-type site |
| site_bao_mid | (0.4, 0.4, 0.5) | Existing Frenkel1-type site |
| site_shifted_ab | (0.2, 0.4, 0.5) | Existing Frenkel3-type site |
| site_origin_mid | (0.0, 0.0, 0.5) | Used by the older `add_O_interstitial.py` helper |

The new script allows these coordinates to be edited in one dictionary.

## Outputs to keep

For each site and temperature, keep lightweight outputs in Git:

- input CIF for the interstitial structure
- MD log file
- final CIF
- cell-vs-time CSV
- summary CSV

Large trajectory files should stay local or on the GPU server and should not be
uploaded to GitHub.

## Large trajectory data

Some previous trajectories are not in the GitHub repository because they are too
large. This does not necessarily mean the calculations were not run. It means
the repository is incomplete as a data archive.

The `.gitignore` file intentionally excludes large trajectory and model files.
If a trajectory is needed later, copy it directly from the GPU machine or keep a
separate local data folder outside GitHub.

## Program

Use the existing project MD program for these runs. The notes above define the
scientific direction and the outputs to track; they are not tied to a newly
written script.
