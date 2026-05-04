import numpy as np
import matplotlib
matplotlib.use("Agg")  # 在无图形界面环境下渲染
import matplotlib.pyplot as plt

log_file = "md.log"

# 自动根据列名读取
data = np.genfromtxt(log_file, names=True)

# genfromtxt 会把列名中的 [] 或 / 替换成下划线
time = data['Timeps']
T = data['TK']
Etot = data['EtoteV']

# -------------------------
# Plot 1: Temperature vs Time
# -------------------------
plt.figure(figsize=(10,4))
plt.plot(time, T, lw=1)
plt.xlabel("Time (ps)", fontsize=12)
plt.ylabel("Temperature (K)", fontsize=12)
plt.title("Temperature vs Time", fontsize=14)
plt.tight_layout()
plt.savefig("temperature_vs_time.png", dpi=200)
plt.close()

# -------------------------
# Plot 2: Total Energy vs Time
# -------------------------
plt.figure(figsize=(10,4))
plt.plot(time, Etot, lw=1)
plt.xlabel("Time (ps)", fontsize=12)
plt.ylabel("Total Energy (eV)", fontsize=12)
plt.title("Total Energy vs Time", fontsize=14)
plt.tight_layout()
plt.savefig("energy_vs_time.png", dpi=200)
plt.close()
