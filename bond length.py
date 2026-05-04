from ase.io import read
import itertools

cif_path = "YBa2Cu3O7_Int1_MACE.cif"

# 1) 读结构
atoms = read(cif_path)
symbols = atoms.get_chemical_symbols()

# 2) 在 CIF 的 loop_ 区块里定位列名并抽取 _atom_site_label
def extract_labels_from_cif(path):
    with open(path, "r") as f:
        lines = [ln.rstrip() for ln in f]

    labels = []
    in_loop = False
    headers = []
    rows = []

    i = 0
    while i < len(lines):
        ln = lines[i].strip()
        # 进入一个 loop_，重置
        if ln.startswith("loop_"):
            in_loop = True
            headers, rows = [], []
            i += 1
            # 收集本 loop_ 的所有 header 行
            while i < len(lines) and lines[i].strip().startswith("_"):
                headers.append(lines[i].strip())
                i += 1
            # 收集本 loop_ 的数据行（到下一段空行或新 loop_/新 header 停）
            while i < len(lines):
                s = lines[i].strip()
                if not s or s.startswith("loop_") or s.startswith("_"):
                    break
                rows.append(s.split())
                i += 1

            # 如果这个 loop_ 同时有 type_symbol 和 label 两列，就抽取
            if ("_atom_site_type_symbol" in headers) and ("_atom_site_label" in headers):
                ti = headers.index("_atom_site_type_symbol")
                li = headers.index("_atom_site_label")
                # 有些 CIF 某些行会被换行/缺列，做下保护
                for r in rows:
                    if len(r) > max(ti, li):
                        labels.append(r[li])
            # 继续下一轮 while，不要 i += 1（上面已在内部推进或停在分隔行）
            continue

        i += 1

    return labels

labels = extract_labels_from_cif(cif_path)
print(labels)

# 3) 简单的键型与 cutoff（Å）——和你前面一致
pair_cut = {
    tuple(sorted(("Cu", "O"))): 3.0,
    tuple(sorted(("Y", "O"))): 3.0,
    tuple(sorted(("Ba", "O"))): 3.3,
    tuple(sorted(("O", "O"))): 3.2,
}

# 4) 遍历所有原子对并输出（使用 MIC 以处理周期性）
bond_entries = []
for i, j in itertools.combinations(range(len(atoms)), 2):
    a, b = symbols[i], symbols[j]
    key = tuple(sorted((a, b)))
    if key not in pair_cut:
        continue
    d = atoms.get_distance(i, j, mic=True)
    if d <= pair_cut[key]:
        bond_entries.append((d, labels[i], labels[j]))

# 按距离升序排序并格式化文本
bond_entries.sort(key=lambda item: item[0])
out_lines = [
    f"{lab_i:>4} -- {lab_j:<4}   {dist:7.4f} Å"
    for dist, lab_i, lab_j in bond_entries
]

# 打印并写文件
for ln in out_lines:
    print(ln)
with open("rebco_bonds.txt", "w") as f:
    f.write("\n".join(out_lines))
