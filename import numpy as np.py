import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

# ========= settings =========
traj_path = "YBa2Cu3O7_md_298K.traj"
md_dt_fs = 0.5
write_interval = 1      # !!! 改成你实际写 traj 的间隔
skip_frames = 100
use_last_frames = None
remove_drift = True
use_window = True
smooth_pts = 1          # 1 表示不平滑；建议 11 或 21（必须是奇数）
max_THz = 20          # 例如 20
# ===========================


# ---------- 读取轨迹 ----------
frames = read(traj_path, index=":")
V = np.array([a.get_velocities() for a in frames], dtype=float)
symbols = np.array(frames[0].get_chemical_symbols())
masses = frames[0].get_masses()

# ---------- 去除漂移 ----------
if remove_drift:
    V = V - V.mean(axis=1, keepdims=True)

# ---------- 质量加权: v_i -> sqrt(m_i) * v_i ----------
sqrt_m = np.sqrt(masses)[np.newaxis, :, np.newaxis]
V = V * sqrt_m

K, N, _ = V.shape
dt_sample_fs = md_dt_fs * write_interval
dt_sample_ps = dt_sample_fs / 1000.0

# ---------- total VACF ----------
X = V.reshape(K, 3 * N)

nfft = 1
while nfft < 2 * K:
    nfft *= 2

F = np.fft.rfft(X, n=nfft, axis=0)
ac = np.fft.irfft(F * np.conjugate(F), n=nfft, axis=0)[:K]
vacf_total = ac.sum(axis=1) / (K - np.arange(K))

tau_ps = np.arange(K) * dt_sample_ps

# ---------- total VDOS ----------
Y = np.fft.rfft(vacf_total)
vdos_total = np.real(Y)
vdos_total[vdos_total < 0] = 0

f_THz = np.fft.rfftfreq(K, d=dt_sample_fs * 1e-15) / 1e12
mask = f_THz <= max_THz

# ---------- 按元素计算 ----------
elements = np.unique(symbols)
vdos_sum = np.zeros_like(vdos_total)

for el in elements:

    idx = np.where(symbols == el)[0]

    V_el = V[:, idx, :]
    X = V_el.reshape(K, 3 * len(idx))

    F = np.fft.rfft(X, n=nfft, axis=0)
    ac = np.fft.irfft(F * np.conjugate(F), n=nfft, axis=0)[:K]

    vacf = ac.sum(axis=1) / (K - np.arange(K))

    Y = np.fft.rfft(vacf)
    vdos = np.real(Y)
    vdos[vdos < 0] = 0

    vdos_sum += vdos

    plt.figure()
    plt.plot(f_THz[mask], vdos[mask])
    plt.xlabel("frequency (THz)")
    plt.ylabel("VDOS")
    plt.title(el)
    plt.tight_layout()

# ---------- 单独画 total ----------
vdos_total = vdos_total / vdos_total.max()
plt.figure()
plt.plot(f_THz[mask], vdos_total[mask])
plt.xlabel("frequency (THz)")
plt.ylabel("VDOS")
plt.title("Total VDOS")
plt.tight_layout()

# ---------- 元素总和 vs total ----------

plt.figure()
plt.plot(f_THz[mask], vdos_total[mask], label="total")
plt.plot(f_THz[mask], vdos_sum[mask], label="sum(elements)")
plt.xlabel("frequency (THz)")
plt.ylabel("VDOS")
plt.legend()
plt.tight_layout()

plt.show()