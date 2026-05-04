# make_frenkel_from_coords.py
import numpy as np
from ase.io import read, write
from ase import Atom

IN_CIF  = "YBa2Cu3O7_553.cif"
OUT_CIF = "Interstitial_553.cif"

# ====== 你给定的坐标（fractional） ======
# 图1：要删除的 O2
# O_del_frac = np.array([0.50000, 0.40000, 0.45944])

# 图2/图3：两个 Cu2，中点放 interstitial O
Cu_low_frac  = np.array([0.40000, 0.40000, 0.21496])
Cu_high_frac = np.array([0.40000, 0.40000, 0.11838])
O_ins_frac   = 0.5 * (Cu_low_frac + Cu_high_frac)  # -> [0.4, 0.4, 0.5]

# ====== 工具：PBC 下 fractional 距离 ======
def frac_delta(a, b):
    d = a - b
    d -= np.round(d)  # minimum image in fractional space
    return d

atoms = read(IN_CIF)
cell = atoms.get_cell().array
frac = atoms.get_scaled_positions(wrap=True)
sym  = np.array(atoms.get_chemical_symbols())


# ====== 1) 删除最接近指定坐标的那个 O 原子 ======
# O_idx = np.where(sym == "O")[0]
# if len(O_idx) == 0:
#     raise RuntimeError("No O atoms found in structure.")

# d2 = np.array([np.dot(frac_delta(frac[i], O_del_frac), frac_delta(frac[i], O_del_frac)) for i in O_idx])
# del_idx = int(O_idx[np.argmin(d2)])

# print(f"Deleting O index: {del_idx}")
# print(f"Deleting O frac (found): {frac[del_idx]}")
# print(f"Target O frac (given):  {O_del_frac}")

# del atoms[del_idx]

# ====== 2) 添加 interstitial O（在两个 Cu2 中点） ======
# frac -> cart
O_ins_cart = O_ins_frac @ cell
atoms.append(Atom("O", position=O_ins_cart))

print(f"Inserted O frac (given): {O_ins_frac}")

# ====== 3) 写出 ======
write(OUT_CIF, atoms)
print("Written:", OUT_CIF)
print("Total atoms:", len(atoms))
