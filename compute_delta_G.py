#!/usr/bin/env python3
"""
Compute defect formation free energy ΔG for REBCO Frenkel defects.

Pipeline:
  1. Read MD trajectory (with velocities) for pristine & defect systems
  2. Compute VDOS via mass-weighted VACF + FFT
  3. Apply quantum statistical mechanics (NOT classical equipartition)
     → This is the "quantum correction": high-freq modes (>6 THz O vibrations)
       are suppressed because they are not fully excited at 298K
  4. Integrate over frequency to get S_vib and F_vib
  5. Compute ΔG = (E_defect - E_pristine) + (F_vib_defect - F_vib_pristine)
     For Frenkel defect: no chemical potential correction (atom count unchanged)
  6. If ΔG > 0 → defect is thermodynamically unstable (will heal)
     If ΔG < 0 → defect is thermodynamically stable

Usage:
  python3 compute_delta_G.py

Requirements:
  pip install ase numpy scipy matplotlib

Author: Andy (generated for kedaya's Part 2 Project)
"""

import numpy as np
from scipy.integrate import trapezoid
import matplotlib.pyplot as plt
import os
import sys

# ============================================================
# CONFIGURATION — edit these paths and values
# ============================================================

# Trajectory files (need velocities stored)
# Supports ASE .traj or extended .xyz with momenta
PRISTINE_TRAJ = "pristine_298K.traj"     # TODO: copy from GPU machine
DEFECT_TRAJ   = "frenkel2_298K.traj"     # TODO: copy from GPU machine

# Static (0K) energies from relaxation or MD t=0 (eV)
# From log files:
#   YBa2Cu3O7_md_298k.log  →  Epot(t=0) = -6049.321 eV
#   Frankel2_md_298k.log   →  Epot(t=0) = -6044.843 eV
E_PRISTINE = -6049.321  # eV
E_DEFECT   = -6044.843  # eV (Frenkel 2, did not recombine)

# MD parameters
MD_DT_FS       = 0.5    # MD timestep in femtoseconds
WRITE_INTERVAL  = 1      # trajectory written every N steps
SKIP_FRAMES     = 500    # skip initial equilibration frames
TEMPERATURE     = 298    # K

# Defect type: 'frenkel' | 'vacancy' | 'interstitial'
DEFECT_TYPE = 'frenkel'

# Output
OUTPUT_DIR = "delta_G_results"

# ============================================================
# Physical constants
# ============================================================
kB_eV   = 8.617333e-5       # Boltzmann constant (eV/K)
h_eV_s  = 4.135668e-15      # Planck constant (eV·s)
hbar_eV = 6.582119e-16      # ħ (eV·s)

# ============================================================
# Step 1: Compute VDOS from MD trajectory
# ============================================================

def compute_vdos(traj_path, md_dt_fs=0.5, write_interval=1,
                 skip_frames=0, max_THz=20.0):
    """
    Compute mass-weighted VDOS from MD trajectory.

    Method:
      1. Read velocities from trajectory
      2. Remove center-of-mass drift
      3. Mass-weight: v_i → √(m_i) × v_i
      4. Compute VACF via Wiener-Khinchin theorem (FFT autocorrelation)
      5. FFT of VACF → VDOS g(ν)
      6. Normalize so ∫g(ν)dν = 3N (total degrees of freedom)

    Returns: (freq_THz, vdos, n_atoms)
    """
    from ase.io import read

    print(f"  Reading: {traj_path}")
    frames = read(traj_path, index=":")

    # Skip equilibration
    if skip_frames > 0:
        frames = frames[skip_frames:]
    print(f"  Using {len(frames)} frames (skipped {skip_frames})")

    V = np.array([a.get_velocities() for a in frames], dtype=float)
    masses = frames[0].get_masses()
    K, N, _ = V.shape  # K=frames, N=atoms

    dt_sample_fs = md_dt_fs * write_interval
    total_ps = K * dt_sample_fs / 1000
    print(f"  {N} atoms, dt={dt_sample_fs} fs, total={total_ps:.2f} ps")

    # Remove COM drift
    V = V - V.mean(axis=1, keepdims=True)

    # Mass weighting: v_i → √(m_i) · v_i
    sqrt_m = np.sqrt(masses)[np.newaxis, :, np.newaxis]  # (1, N, 1)
    V = V * sqrt_m

    # VACF via Wiener-Khinchin: pad to next power of 2
    nfft = 1
    while nfft < 2 * K:
        nfft *= 2

    X = V.reshape(K, 3 * N)
    F = np.fft.rfft(X, n=nfft, axis=0)
    ac = np.fft.irfft(F * np.conjugate(F), n=nfft, axis=0)[:K]
    vacf = ac.sum(axis=1) / (K - np.arange(K))

    # VDOS = FFT of VACF
    Y = np.fft.rfft(vacf)
    vdos = np.real(Y)
    vdos[vdos < 0] = 0

    # Frequency axis
    freq_THz = np.fft.rfftfreq(K, d=dt_sample_fs * 1e-15) / 1e12

    # Frequency cutoff
    mask = freq_THz <= max_THz
    freq_THz = freq_THz[mask]
    vdos = vdos[mask]

    # Normalize: ∫g(ν)dν = 3N
    integral = trapezoid(vdos, freq_THz)
    if integral > 0:
        vdos *= (3 * N) / integral

    print(f"  VDOS integral = {trapezoid(vdos, freq_THz):.1f} (target: {3*N})")

    return freq_THz, vdos, N


# ============================================================
# Step 2: Quantum thermodynamics from VDOS
# ============================================================

def compute_vibrational_thermodynamics(freq_THz, vdos, temperature):
    """
    Compute vibrational thermodynamic quantities using QUANTUM harmonic
    oscillator statistics (this is the 'quantum correction').

    Classical MD gives all modes kBT energy (equipartition).
    Quantum mechanics: high-freq modes with hν >> kBT are NOT fully excited.

    The Bose-Einstein occupation n(ν) = 1/(exp(hν/kBT) - 1) naturally
    suppresses high-frequency contributions.

    Formulas (integrated over g(ν)):
      ZPVE  = ∫ (hν/2) g(ν) dν
      E_vib = ∫ hν [1/2 + n_BE(ν)] g(ν) dν
      F_vib = ZPVE + kBT ∫ ln(1 - e^{-x}) g(ν) dν     where x = hν/(kBT)
      S_vib = (E_vib - F_vib) / T

    Returns dict with all quantities in eV (entropy in eV/K and kB units).
    """
    # Skip ν=0 (acoustic modes at Γ point → divergence)
    mask = freq_THz > 0.01  # THz
    nu = freq_THz[mask]
    g = vdos[mask]

    # Convert frequency to energy
    h_nu = h_eV_s * nu * 1e12  # eV (h × ν)
    x = h_nu / (kB_eV * temperature)  # dimensionless hν/kBT

    # Clamp x to avoid overflow
    x_safe = np.minimum(x, 500)
    exp_x = np.exp(x_safe)

    # Zero-point vibrational energy
    zpve = trapezoid(0.5 * h_nu * g, nu)

    # Bose-Einstein occupation
    n_bose = 1.0 / (exp_x - 1)

    # Total vibrational energy (ZPE + thermal)
    e_vib = trapezoid(h_nu * (0.5 + n_bose) * g, nu)

    # Vibrational Helmholtz free energy
    log_term = np.log(1.0 - np.exp(-x_safe))
    f_vib = zpve + kB_eV * temperature * trapezoid(log_term * g, nu)

    # Vibrational entropy: S = (E - F) / T
    s_vib = (e_vib - f_vib) / temperature

    results = {
        'ZPVE_eV':    zpve,
        'E_vib_eV':   e_vib,
        'F_vib_eV':   f_vib,
        'S_vib_eV_K': s_vib,         # eV/K
        'S_vib_kB':   s_vib / kB_eV,  # in units of kB
        'T_K':        temperature,
    }

    print(f"  T = {temperature} K:")
    print(f"    ZPVE     = {zpve:.4f} eV")
    print(f"    E_vib    = {e_vib:.4f} eV")
    print(f"    F_vib    = {f_vib:.4f} eV")
    print(f"    S_vib    = {s_vib/kB_eV:.2f} kB  ({s_vib*1e3:.4f} meV/K)")

    return results


# ============================================================
# Step 3: ΔG calculation
# ============================================================

def compute_delta_G(thermo_pristine, thermo_defect,
                    E_pristine, E_defect, temperature,
                    defect_type='frenkel'):
    """
    Defect formation Gibbs free energy:

      ΔG = (E_defect - E_pristine) + Δμ + (F_vib_defect - F_vib_pristine)

    Equivalently (since F = E_vib - TS):
      ΔG = ΔE_static + ΔF_vib
         = ΔE_static + ΔZPVE - TΔS_vib + Δ(thermal E)

    For Frenkel defect: Δμ = 0 (no atoms added/removed)
    """
    delta_E = E_defect - E_pristine
    delta_F_vib = thermo_defect['F_vib_eV'] - thermo_pristine['F_vib_eV']
    delta_S_vib = thermo_defect['S_vib_eV_K'] - thermo_pristine['S_vib_eV_K']
    delta_ZPVE = thermo_defect['ZPVE_eV'] - thermo_pristine['ZPVE_eV']

    # Chemical potential correction
    delta_mu = 0.0  # Frenkel: no atoms added/removed
    # For vacancy: delta_mu = +mu_O
    # For interstitial: delta_mu = -mu_O

    delta_G = delta_E + delta_mu + delta_F_vib

    print(f"\n{'='*55}")
    print(f"  ΔG Calculation ({defect_type}, T={temperature} K)")
    print(f"{'='*55}")
    print(f"  ΔE (static 0K)    = {delta_E:+.4f} eV")
    print(f"  ΔZPVE             = {delta_ZPVE:+.4f} eV")
    print(f"  ΔF_vib            = {delta_F_vib:+.4f} eV")
    print(f"  ΔS_vib            = {delta_S_vib/kB_eV:+.2f} kB")
    print(f"  −TΔS_vib          = {-temperature*delta_S_vib:+.4f} eV")
    print(f"  ─────────────────────────────────")
    print(f"  ΔG({temperature}K)         = {delta_G:+.4f} eV")
    print(f"{'='*55}")

    if delta_G > 0:
        print(f"  → ΔG > 0: Defect is UNSTABLE, will tend to recombine/heal")
    else:
        print(f"  → ΔG < 0: Defect is STABLE at {temperature} K")

    return {
        'delta_G_eV':       delta_G,
        'delta_E_eV':       delta_E,
        'delta_F_vib_eV':   delta_F_vib,
        'delta_S_vib_eV_K': delta_S_vib,
        'delta_S_vib_kB':   delta_S_vib / kB_eV,
        'delta_ZPVE_eV':    delta_ZPVE,
        'T_K':              temperature,
    }


# ============================================================
# Step 4: Temperature scan
# ============================================================

def temperature_scan(freq_p, vdos_p, freq_d, vdos_d,
                     E_p, E_d, T_range, defect_type='frenkel'):
    """Compute ΔG vs temperature."""
    results = []
    for T in T_range:
        tp = compute_vibrational_thermodynamics(freq_p, vdos_p, T)
        td = compute_vibrational_thermodynamics(freq_d, vdos_d, T)
        delta_F = td['F_vib_eV'] - tp['F_vib_eV']
        delta_G = (E_d - E_p) + delta_F
        results.append(delta_G)
    return np.array(results)


# ============================================================
# Plotting
# ============================================================

def plot_vdos_comparison(freq_p, vdos_p, freq_d, vdos_d, output_dir):
    """Plot pristine vs defect VDOS."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

    # Normalize for visual comparison
    vp_norm = vdos_p / np.max(vdos_p)
    vd_norm = vdos_d / np.max(vdos_d)

    ax1.plot(freq_p, vp_norm, 'b-', lw=1.5, label='Pristine')
    ax1.plot(freq_d, vd_norm, 'r-', lw=1.5, alpha=0.8, label='Defect')
    ax1.set_xlabel('Frequency (THz)', fontsize=12)
    ax1.set_ylabel('Normalized VDOS', fontsize=12)
    ax1.set_title('VDOS: Pristine vs Frenkel Defect', fontsize=13)
    ax1.legend(fontsize=11)
    ax1.set_xlim(0, 20)
    ax1.grid(True, alpha=0.3)

    # Quantum weight visualization at 298K
    nu_plot = np.linspace(0.1, 20, 500)
    h_nu = h_eV_s * nu_plot * 1e12
    x = h_nu / (kB_eV * 298)
    # Entropy weight: x/(e^x - 1) - ln(1 - e^{-x})
    exp_x = np.exp(np.minimum(x, 500))
    weight = x / (exp_x - 1) - np.log(1 - np.exp(-np.minimum(x, 500)))
    weight_norm = weight / weight.max()

    ax2.fill_between(nu_plot, weight_norm, alpha=0.3, color='green',
                     label='Quantum weight (298K)')
    ax2.plot(nu_plot, weight_norm, 'g-', lw=1.5)
    ax2.axvline(x=6, color='red', linestyle='--', alpha=0.7,
                label='~6 THz (O modes suppressed)')
    ax2.set_xlabel('Frequency (THz)', fontsize=12)
    ax2.set_ylabel('Weight (normalized)', fontsize=12)
    ax2.set_title('Quantum Correction Weight at 298K', fontsize=13)
    ax2.legend(fontsize=11)
    ax2.set_xlim(0, 20)
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(f"{output_dir}/vdos_comparison.png", dpi=300)
    plt.close()
    print(f"  Saved: {output_dir}/vdos_comparison.png")


def plot_delta_G_vs_T(T_range, delta_Gs, temperature, output_dir):
    """Plot ΔG vs temperature."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(T_range, delta_Gs, 'b-o', lw=2, markersize=4)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=temperature, color='red', linestyle=':',
               alpha=0.7, label=f'{temperature} K (target)')

    # Mark the target temperature value
    from scipy.interpolate import interp1d
    f_interp = interp1d(T_range, delta_Gs)
    dG_at_T = f_interp(temperature)
    ax.plot(temperature, dG_at_T, 'r*', markersize=15,
            label=f'ΔG({temperature}K) = {dG_at_T:.3f} eV')

    ax.set_xlabel('Temperature (K)', fontsize=12)
    ax.set_ylabel('ΔG (eV)', fontsize=12)
    ax.set_title('Defect Formation Free Energy vs Temperature', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/delta_G_vs_T.png", dpi=300)
    plt.close()
    print(f"  Saved: {output_dir}/delta_G_vs_T.png")


# ============================================================
# MAIN
# ============================================================

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("  ΔG Calculation for REBCO Frenkel Defect")
    print("=" * 60)

    # --- Check files exist ---
    for f, label in [(PRISTINE_TRAJ, "Pristine"), (DEFECT_TRAJ, "Defect")]:
        if not os.path.exists(f):
            print(f"\n❌ {label} trajectory not found: {f}")
            print(f"   Copy the 298K trajectory file from the GPU machine.")
            print(f"   Expected: 5×5×3 supercell (975 atoms), 298K, ≥5000 steps")
            sys.exit(1)

    # --- Step 1: Compute VDOS ---
    print("\n--- Step 1a: Pristine VDOS ---")
    freq_p, vdos_p, N_p = compute_vdos(
        PRISTINE_TRAJ, md_dt_fs=MD_DT_FS,
        write_interval=WRITE_INTERVAL, skip_frames=SKIP_FRAMES)

    print("\n--- Step 1b: Defect VDOS ---")
    freq_d, vdos_d, N_d = compute_vdos(
        DEFECT_TRAJ, md_dt_fs=MD_DT_FS,
        write_interval=WRITE_INTERVAL, skip_frames=SKIP_FRAMES)

    # --- Step 2: Vibrational thermodynamics (quantum) ---
    print("\n--- Step 2a: Pristine thermodynamics ---")
    thermo_p = compute_vibrational_thermodynamics(freq_p, vdos_p, TEMPERATURE)

    print("\n--- Step 2b: Defect thermodynamics ---")
    thermo_d = compute_vibrational_thermodynamics(freq_d, vdos_d, TEMPERATURE)

    # --- Step 3: ΔG ---
    print("\n--- Step 3: ΔG ---")
    result = compute_delta_G(thermo_p, thermo_d,
                             E_PRISTINE, E_DEFECT,
                             TEMPERATURE, DEFECT_TYPE)

    # --- Step 4: Temperature scan ---
    print("\n--- Step 4: Temperature scan ---")
    T_range = np.arange(50, 601, 25)
    delta_Gs = temperature_scan(freq_p, vdos_p, freq_d, vdos_d,
                                E_PRISTINE, E_DEFECT, T_range,
                                DEFECT_TYPE)

    # --- Plots ---
    print("\n--- Generating plots ---")
    plot_vdos_comparison(freq_p, vdos_p, freq_d, vdos_d, OUTPUT_DIR)
    plot_delta_G_vs_T(T_range, delta_Gs, TEMPERATURE, OUTPUT_DIR)

    # --- Save results ---
    np.savez(f"{OUTPUT_DIR}/results.npz",
             freq_pristine=freq_p, vdos_pristine=vdos_p,
             freq_defect=freq_d, vdos_defect=vdos_d,
             temperatures=T_range, delta_G_vs_T=delta_Gs,
             E_pristine=E_PRISTINE, E_defect=E_DEFECT,
             delta_G_298K=result['delta_G_eV'],
             delta_S_vib_kB=result['delta_S_vib_kB'],
             N_pristine=N_p, N_defect=N_d)

    # Summary text file
    with open(f"{OUTPUT_DIR}/summary.txt", 'w') as f:
        f.write("REBCO Frenkel Defect ΔG Calculation Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Pristine traj: {PRISTINE_TRAJ} ({N_p} atoms)\n")
        f.write(f"Defect traj:   {DEFECT_TRAJ} ({N_d} atoms)\n")
        f.write(f"Temperature:   {TEMPERATURE} K\n")
        f.write(f"Defect type:   {DEFECT_TYPE}\n\n")
        f.write(f"E_pristine (static) = {E_PRISTINE:.3f} eV\n")
        f.write(f"E_defect (static)   = {E_DEFECT:.3f} eV\n")
        f.write(f"ΔE (static)         = {result['delta_E_eV']:+.4f} eV\n\n")
        f.write(f"ΔZPVE               = {result['delta_ZPVE_eV']:+.4f} eV\n")
        f.write(f"ΔF_vib              = {result['delta_F_vib_eV']:+.4f} eV\n")
        f.write(f"ΔS_vib              = {result['delta_S_vib_kB']:+.2f} kB\n")
        f.write(f"−TΔS_vib            = {-TEMPERATURE*result['delta_S_vib_eV_K']:+.4f} eV\n\n")
        f.write(f"ΔG({TEMPERATURE}K)           = {result['delta_G_eV']:+.4f} eV\n\n")
        if result['delta_G_eV'] > 0:
            f.write("→ Defect is UNSTABLE at room temperature\n")
        else:
            f.write("→ Defect is STABLE at room temperature\n")

    print(f"\n{'='*60}")
    print(f"  All results saved to {OUTPUT_DIR}/")
    print(f"  Key result: ΔG({TEMPERATURE}K) = {result['delta_G_eV']:+.4f} eV")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
