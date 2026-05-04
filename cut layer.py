import numpy as np
from ase.io import read, write

INFILE = "YBa2Cu3O7_553_Frenkel3.traj"
OUT    = "YBa2Cu3O7_553_Frenkel3_331.traj"

Na, Nb, Nc = 5, 5, 3
ia, ib, ic = 2, 3, 2
r_ab, r_c  = 1, 0   # 改这里：331 -> (1,0); 113 -> (0,1); 333 -> (1,1)

def idx_list(i, N, r):
    return [((i-1 + d) % N) + 1 for d in range(-r, r+1)]

A = idx_list(ia, Na, r_ab)
B = idx_list(ib, Nb, r_ab)
C = idx_list(ic, Nc, r_c)

# 用第0帧固定mask（成员固定）
ref = read(INFILE, index=0)
frac = ref.get_scaled_positions(wrap=True)

mask = np.zeros(len(ref), dtype=bool)
for a in A:
    a0, a1 = (a-1)/Na, a/Na
    for b in B:
        b0, b1 = (b-1)/Nb, b/Nb
        for c in C:
            c0, c1 = (c-1)/Nc, c/Nc
            mask |= (
                (frac[:,0] >= a0) & (frac[:,0] < a1) &
                (frac[:,1] >= b0) & (frac[:,1] < b1) &
                (frac[:,2] >= c0) & (frac[:,2] < c1)
            )

# 所有帧沿用同一组 index，不再重算mask
frames = read(INFILE, index=":")
out_frames = [atoms[mask] for atoms in frames]

write(OUT, out_frames)
print("Written:", OUT, "frames:", len(out_frames), "atoms/frame:", len(out_frames[0]))
