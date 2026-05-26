# rebco_bonds.py
from ase.io import read
from ase.neighborlist import neighbor_list

# === 1) 读 CIF ===
atoms = read("YBa2Cu3O7_relaxed_MACE.cif")

# === 2) 取标签：优先用 CIF 的 _atom_site_label ===
labels = None
for k in ("labels", "atom_site_label", "cif_labels"):
    if k in atoms.arrays:
        labels = atoms.get_array(k)
        break

if labels is None:  # 兜底：按元素计数生成 Y1、Ba1…
    counts = {}
    labels = []
    for s in atoms.get_chemical_symbols():
        counts[s] = counts.get(s, 0) + 1
        labels.append(f"{s}{counts[s]}")

sym = atoms.get_chemical_symbols()

# === 3) 只保留 M–O 键，并给每类一个合理的 cutoff（Å）===
pair_cut = {
    tuple(sorted(("Cu", "O"))): 2.60,
    tuple(sorted(("Y",  "O"))): 2.50,
    tuple(sorted(("Ba", "O"))): 3.30,
    tuple(sorted(("O",  "O"))): 3.20,
}
Rmax = max(pair_cut.values())

# === 4) 找邻居（最小影像），按类型与 cutoff 过滤，去重 ===
i_idx, j_idx, dists = neighbor_list("ijd", atoms, Rmax)
seen = set()
rows = []
for i, j, d in zip(i_idx, j_idx, dists):
    a, b = sym[i], sym[j]
    key = tuple(sorted((a, b)))
    if key not in pair_cut:
        continue
    if d > pair_cut[key]:
        continue
    pair = (min(i, j), max(i, j))
    if pair in seen:
        continue
    seen.add(pair)
    rows.append((labels[pair[0]], labels[pair[1]], d))

# === 5) 排序并输出 ===
rows.sort(key=lambda x: x[2])

with open("rebco_bonds.txt", "w") as f:
    for L1, L2, d in rows:
        line = f"{L1:>4} -- {L2:<4}  {d:7.4f} Å"
        print(line)
        f.write(line + "\n")