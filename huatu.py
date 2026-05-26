import matplotlib.pyplot as plt
import numpy as np

# =========================
# Pristine DFT data: defines x-axis labels and ordering
# =========================
dft = {
    "Cu1--O2": 1.87545,
    "Cu1--O1": 1.87545,
    "Cu3--O4": 1.94252,
    "Cu2--O3": 1.94252,
    "Cu1--O7": 1.96188,
    "Cu2--O5": 1.98220,
    "Cu3--O6": 1.98220,
    "Cu2--O1": 2.31937,
    "Cu3--O2": 2.31937,
    "Y1--O6": 2.39493,
    "Y1--O5": 2.39493,
    "Y1--O3": 2.41963,
    "Y1--O4": 2.41963,
    "O1--O7": 2.71409,
    "O2--O7": 2.71409,
    "O4--O6": 2.74464,
    "O3--O5": 2.74464,
    "Ba2--O2": 2.75785,
    "Ba1--O1": 2.75785,
    "O3--O4": 2.83241,
    "O5--O6": 2.86489,
    "Ba2--O7": 2.87880,
    "Ba1--O7": 2.87880,
}

# =========================
# Pristine ASE data: taken from screenshot in ascending order
# then mapped onto DFT labels by position
# =========================
ase_sorted_values = [
    1.8786, 1.8786, 1.9496, 1.9496, 1.9606, 1.9830, 1.9830,
    2.3067, 2.3067, 2.4015, 2.4015, 2.4268, 2.4268,
    2.7154, 2.7154, 2.7485, 2.7485, 2.7609, 2.7609,
    2.8604, 2.8685, 2.8792, 2.8792
]

# 按 DFT 数值升序排列
dft_sorted = sorted(dft.items(), key=lambda kv: kv[1])
bond_labels = [k for k, v in dft_sorted]
dft_y = [v for k, v in dft_sorted]
ase_y = ase_sorted_values

x = np.arange(len(bond_labels))

plt.figure(figsize=(16, 6))
plt.plot(x, dft_y, marker='o', linewidth=1.8, label='DFT')
plt.plot(x, ase_y, marker='s', linewidth=1.8, label='ASE')

plt.xticks(x, bond_labels, rotation=60, ha='right')
plt.xlabel("Bond type ")
plt.ylabel("Bond length (Å)")
plt.title("Pristine bond length comparison ")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()