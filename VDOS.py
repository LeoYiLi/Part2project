import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

traj_path = "yba553_md_30K.traj"   # 改成你的
dt_fs = 0.5             # 你说的是 0.5 fs
skip = 100              # 丢掉前200帧(按需改：比如 5-20% 的轨迹)

traj = read(traj_path, index=":")

# 取平衡后的片段
traj = traj[skip:]
n_steps = len(traj)
n_atoms = len(traj[0])

# 提取 velocities: (n_steps, n_atoms, 3)
V = np.array([a.get_velocities() for a in traj], dtype=float)

# VACF(t) = <v(0)·v(t)> 做 time-origin 平均
vacf = np.zeros(n_steps)
for t in range(n_steps):
    prod = (V[:n_steps - t] * V[t:]).sum(axis=(1,2))  # 每个 time-origin 的 dot sum
    vacf[t] = prod.mean() / n_atoms

# 加一个简单窗函数减少频谱泄漏（可选但推荐）
#window = np.hanning(n_steps)
#vacf_w = vacf * window

# FFT -> VDOS
dt = dt_fs * 1e-15
freq = np.fft.rfftfreq(n_steps, d=dt)          # Hz
vdos = np.real(np.fft.rfft(vacf))

# 转成 THz 更直观
freq_THz = freq * 1e-12

plt.figure()
plt.plot(freq_THz, vdos)
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS (a.u.)")
plt.xlim(0, freq_THz.max())
plt.show()
