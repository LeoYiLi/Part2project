# python make_O1_vac.py input.cif
# This Python script creates an oxygen-vacancy structure at the O1 site of YBa₂Cu₃O₇ (YBCO). 
# It reads a pristine .cif file, identifies the oxygen atom nearest to the O1 fractional coordinate (0, 0, 0.15919)
# deletes it, and writes a new .cif file containing the defective structure.
import sys
import os
import numpy as np
from ase.io import read, write

if len(sys.argv) < 2:
    print("Usage: python make_O1_vac.py input.cif")
    sys.exit(1)

infile = sys.argv[1]
root, _ = os.path.splitext(infile)
outfile = f"{root}_Ovac1.cif"


def extract_labels_from_cif(path):
    """返回 CIF 中 _atom_site_label 顺序列表"""
    labels = []
    with open(path, "r") as f:
        lines = [ln.rstrip() for ln in f]

    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        if not ln.startswith("loop_"):
            i += 1
            continue

        headers = []
        rows = []
        i += 1
        while i < len(lines) and lines[i].strip().startswith("_"):
            headers.append(lines[i].strip())
            i += 1

        while i < len(lines):
            s = lines[i].strip()
            if not s or s.startswith("loop_") or s.startswith("_"):
                break
            rows.append(s.split())
            i += 1

        if (
            "_atom_site_label" in headers
            and "_atom_site_type_symbol" in headers
        ):
            label_idx = headers.index("_atom_site_label")
            for row in rows:
                if len(row) > label_idx:
                    labels.append(row[label_idx])
            break

    return labels


atoms = read(infile)
labels = extract_labels_from_cif(infile)

target_label = "O1"
del_idx = None
if len(labels) == len(atoms):
    for idx, label in enumerate(labels):
        if label == target_label:
            del_idx = idx
            break

frac = atoms.get_scaled_positions(wrap=True)
symbols = atoms.get_chemical_symbols()

if del_idx is None:
    # 如果无法从 label 匹配，就退回坐标查找
    target_frac = np.array([0.0, 0.5, 0.0])
    tolerance = 5e-3
    o_idx = [i for i, s in enumerate(symbols) if s == "O"]
    if not o_idx:
        raise RuntimeError("结构中没有找到任何氧原子，无法删除 O1。")

    diff = ((frac[o_idx] - target_frac + 0.5) % 1.0) - 0.5
    dist = np.linalg.norm(diff, axis=1)
    closest = int(np.argmin(dist))
    if dist[closest] > tolerance:
        raise RuntimeError(
            f"未找到接近坐标 {target_frac.tolist()} 的氧原子，"
            "请检查输入 CIF 或调高 tolerance。"
        )
    del_idx = o_idx[closest]

print(f"Deleting atom #{del_idx} (label={labels[del_idx] if del_idx < len(labels) else 'N/A'}, "
      f"frac={frac[del_idx]})")

del atoms[del_idx]
write(outfile, atoms)

print(f"Wrote: {outfile}")
