"""在给定的 CIF 结构中，于分数坐标 (0, 0, 0.5) 处插入一个氧间隙原子。

用法:
    python add_O_interstitial.py input.cif [output.cif]
"""

import sys
import os
import numpy as np
from ase import Atom
from ase.io import read, write

TARGET_FRAC = np.array([0.0, 0.0, 0.5])
NEW_LABEL = "O_int"


def extract_labels(path):
    """提取 CIF 中 atom_site_label 列，保持原顺序。"""
    labels = []
    with open(path, "r") as handle:
        lines = [ln.rstrip() for ln in handle]

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("loop_"):
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

        if "_atom_site_label" in headers and "_atom_site_type_symbol" in headers:
            label_idx = headers.index("_atom_site_label")
            for row in rows:
                if len(row) > label_idx:
                    labels.append(row[label_idx])
            break

    return labels


def rewrite_labels(path, new_labels):
    """将输出 CIF 中的 label 列替换为新列表。"""
    with open(path, "r") as handle:
        lines = handle.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line.startswith("loop_"):
            i += 1
            continue

        headers = []
        j = i + 1
        while j < len(lines) and lines[j].strip().startswith("_"):
            headers.append(lines[j].strip())
            j += 1
        data_start = j
        while j < len(lines):
            s = lines[j].strip()
            if not s or s.startswith("loop_") or s.startswith("_"):
                break
            j += 1
        data_end = j

        if "_atom_site_label" in headers and "_atom_site_type_symbol" in headers:
            rows = lines[data_start:data_end]
            if len(rows) != len(new_labels):
                raise RuntimeError("label 数量与原子数不一致，无法重写 CIF。")
            label_idx = headers.index("_atom_site_label")
            new_rows = []
            for row_line, label in zip(rows, new_labels):
                parts = row_line.split()
                if len(parts) <= label_idx:
                    new_rows.append(row_line)
                    continue
                parts[label_idx] = label
                new_rows.append("  " + "  ".join(parts) + "\n")
            lines[data_start:data_end] = new_rows
            break

        i = j

    with open(path, "w") as handle:
        handle.writelines(lines)


def main():
    if len(sys.argv) < 2:
        print("用法: python add_O_interstitial.py input.cif [output.cif]")
        sys.exit(1)

    infile = sys.argv[1]
    outfile = sys.argv[2] if len(sys.argv) >= 3 else os.path.splitext(infile)[0] + "_Oint.cif"

    atoms = read(infile)
    labels = extract_labels(infile)
    if len(labels) != len(atoms):
        labels = atoms.get_chemical_symbols()

    cart = np.dot(TARGET_FRAC, atoms.get_cell())
    atoms.append(Atom("O", position=cart))
    labels.append(NEW_LABEL)

    write(outfile, atoms)
    rewrite_labels(outfile, labels)

    print(f"已在分数坐标 {TARGET_FRAC.tolist()} 处添加 O 间隙原子，写入 {outfile}")


if __name__ == "__main__":
    main()
