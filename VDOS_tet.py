import numpy as np
import matplotlib.pyplot as plt

# ========= settings =========
dos_path = "YBa2Cu3O7shuju2"  # 你现在的输入文件
x_unit = "THz"         # "cm^-1" 或 "THz"
xmax = 20              # 如果 x_unit="THz"，这里是 THz；如果 "cm^-1"，这里是 cm^-1
normalize = False      # True: 归一化到最大值=1（只看形状）
smooth_pts = 1         # 1 表示不平滑；可用 11/21（奇数）
# ===========================

def moving_average(y, k):
    if k is None or k <= 1:
        return y
    if k % 2 == 0:
        k += 1
    kernel = np.ones(k) / k
    return np.convolve(y, kernel, mode="same")

# --- read data (skip comment lines beginning with '#')
data = np.loadtxt(dos_path, comments="#")
if data.shape[1] < 2:
    raise ValueError("文件列数不足：需要至少两列 (Energy[cm^-1], TOTAL-DOS)")

energy_cm = data[:, 0]
dos_total = data[:, 1]

# --- optional smoothing
dos_total_s = moving_average(dos_total, smooth_pts)

# --- optional normalize
if normalize and dos_total_s.max() != 0:
    dos_total_s = dos_total_s / dos_total_s.max()

# --- unit conversion
# 1 THz ≈ 33.3564095 cm^-1  => THz = (cm^-1) / 33.3564095
CM_PER_THz = 33.3564095

if x_unit.lower() in ["thz", "t", "freq"]:
    x = energy_cm / CM_PER_THz
    xlabel = "Frequency (THz)"
else:
    x = energy_cm
    xlabel = r"Energy / Frequency (cm$^{-1}$)"

# --- plot
plt.figure()
plt.plot(x, dos_total_s, label="TOTAL-DOS")
plt.xlabel(xlabel)
plt.ylabel("VDOS / phonon DOS (a.u.)" if normalize else "VDOS / phonon DOS")
plt.title("YBa2Cu3O7 DOS (from file)")
if xmax is not None:
    plt.xlim(0, xmax)
plt.legend()
plt.tight_layout()
plt.show()