import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

# ========= settings =========
traj_path = "Frankel_md_298K.traj"
md_dt_fs = 0.5
write_interval = 1      # !!! 改成你实际写 traj 的间隔
skip_frames = 100
use_last_frames = 1000
remove_drift = True
use_window = True
smooth_pts = 11          # 1 表示不平滑；建议 11 或 21（必须是奇数）
xmax_THz = 20          # 例如 20
# ===========================

def moving_average(y, k):
    if k is None or k <= 1:
        return y
    if k % 2 == 0:
        k += 1
    kernel = np.ones(k) / k
    return np.convolve(y, kernel, mode="same")

def spectrum_from_vel(V, dt_fs, remove_drift=True, use_window=True):
    # V: (n_steps, n_atoms, 3)
    V = V.astype(float)

    if remove_drift:
        V = V - V.mean(axis=1, keepdims=True)

    n_steps, n_atoms, _ = V.shape

    # VACF with time-origin averaging
    vacf = np.zeros(n_steps, dtype=float)
    for t in range(n_steps):
        prod = (V[:n_steps - t] * V[t:]).sum(axis=(1, 2))  # dot sum per origin
        vacf[t] = prod.mean()

    if use_window:
        vacf = vacf * np.hanning(n_steps)

    dt = dt_fs * 1e-15
    freq = np.fft.rfftfreq(n_steps, d=dt)    # Hz
    # Power-spectrum-like (non-negative)
    spec = np.real(np.fft.rfft(vacf))
    #spec = np.real(np.fft.rfft(vacf)) ** 2
    freq_THz = freq * 1e-12
    return freq_THz, spec

# ---- load & trim
traj = read(traj_path, index=":")
traj = traj[skip_frames:]
if use_last_frames is not None and len(traj) > use_last_frames:
    traj = traj[-use_last_frames:]

n_steps = len(traj)
if n_steps < 10:
    raise ValueError("Too few frames after trimming.")

dt_fs = md_dt_fs * write_interval
symbols = traj[0].get_chemical_symbols()
O_idx = [i for i, s in enumerate(symbols) if s == "O"]

print("total atoms =", len(symbols))
print("O atoms     =", len(O_idx))
print("frames used =", n_steps)
print(f"dt_sampling = {dt_fs} fs, f_max ~ {1/(2*dt_fs*1e-15)*1e-12:.2f} THz")

V_all = np.array([a.get_velocities() for a in traj], dtype=float)
V_O = V_all[:, O_idx, :]

f, spec_all = spectrum_from_vel(V_all, dt_fs, remove_drift=remove_drift, use_window=use_window)
_, spec_O = spectrum_from_vel(V_O, dt_fs, remove_drift=remove_drift, use_window=use_window)

# smooth (optional)
spec_all_s = moving_average(spec_all, smooth_pts)
spec_O_s = moving_average(spec_O, smooth_pts)

# normalize for shape comparison
# spec_all_n = spec_all_s / spec_all_s.max() if spec_all_s.max() != 0 else spec_all_s
# spec_O_n = spec_O_s / spec_O_s.max() if spec_O_s.max() != 0 else spec_O_s

# diff = spec_O_n - spec_all_n

# ---- plot overlay
plt.figure()
plt.plot(f, spec_all_s, label="Total (norm, |FFT|)")
plt.plot(f, spec_O_s, label="O-only (norm, |FFT|)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (normalized a.u.)")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# ============================
# 在你现有代码基础上，新增：Cu/Ba/Y 以及 8 张图
# （不修改你已有函数和计算逻辑）
# ============================

# ---- element indices
Cu_idx = [i for i, s in enumerate(symbols) if s == "Cu"]
Ba_idx = [i for i, s in enumerate(symbols) if s == "Ba"]
Y_idx  = [i for i, s in enumerate(symbols) if s == "Y"]

print("Cu atoms    =", len(Cu_idx))
print("Ba atoms    =", len(Ba_idx))
print("Y atoms     =", len(Y_idx))

# ---- slice velocities
V_Cu = V_all[:, Cu_idx, :]
V_Ba = V_all[:, Ba_idx, :]
V_Y  = V_all[:, Y_idx, :]

# ---- compute spectra (same function, same settings)
_, spec_Cu = spectrum_from_vel(V_Cu, dt_fs, remove_drift=remove_drift, use_window=use_window)
_, spec_Ba = spectrum_from_vel(V_Ba, dt_fs, remove_drift=remove_drift, use_window=use_window)
_, spec_Y  = spectrum_from_vel(V_Y,  dt_fs, remove_drift=remove_drift, use_window=use_window)

# ---- smooth (same smoothing)
spec_Cu_s = moving_average(spec_Cu, smooth_pts)
spec_Ba_s = moving_average(spec_Ba, smooth_pts)
spec_Y_s  = moving_average(spec_Y,  smooth_pts)

# ---- sum of 4 elements (for visual comparison)
spec_sum4_s = spec_O_s + spec_Cu_s + spec_Ba_s + spec_Y_s

# ============================
# 8 figures
# ============================

# 1) O only
plt.figure()
plt.plot(f, spec_O_s, label="O-only (|FFT| of VACF)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("O-only")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 2) Cu only
plt.figure()
plt.plot(f, spec_Cu_s, label="Cu-only (|FFT| of VACF)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Cu-only")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 3) Ba only
plt.figure()
plt.plot(f, spec_Ba_s, label="Ba-only (|FFT| of VACF)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Ba-only")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 4) Y only
plt.figure()
plt.plot(f, spec_Y_s, label="Y-only (|FFT| of VACF)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Y-only")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 5) total
plt.figure()
plt.plot(f, spec_all_s, label="Total (|FFT| of VACF)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Total")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 6) total + O
plt.figure()
plt.plot(f, spec_all_s, label="Total")
plt.plot(f, spec_O_s, label="O-only")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Total vs O-only")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 7) sum4 vs total (visual check)
# 注意：由于你用的是 abs(FFT)，此图不保证严格重合（这是方法本身的数学性质，不是你数据的问题）
plt.figure()
plt.plot(f, spec_all_s, label="Total (direct)")
plt.plot(f, spec_sum4_s, "--", label="O+Cu+Ba+Y (summed)")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Direct Total vs Sum of Element Contributions (visual check)")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

# 8) O/Cu/Ba/Y overlay
plt.figure()
plt.plot(f, spec_O_s,  label="O")
plt.plot(f, spec_Cu_s, label="Cu")
plt.plot(f, spec_Ba_s, label="Ba")
plt.plot(f, spec_Y_s,  label="Y")
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS-like spectrum (a.u.)")
plt.title("Element-resolved spectra overlay")
if xmax_THz is not None:
    plt.xlim(0, xmax_THz)
plt.legend()
plt.tight_layout()

plt.show()