import numpy as np
from ase.io import read

traj = read("Frankel_md_77K.traj", index=":")

def frame_info(i):
    a = traj[i]
    mom = a.get_momenta()
    vel = a.get_velocities()
    # mom 可能是 None；vel 可能是 None 或者数组
    mom_max = None if mom is None else float(np.abs(mom).max())
    vel_max = None if vel is None else float(np.abs(vel).max())
    return mom_max, vel_max

for i in [0, 1, 2, -3, -2, -1]:
    mom_max, vel_max = frame_info(i)
    print(i, "mom_max=", mom_max, "vel_max=", vel_max, "n_atoms=", len(traj[i]))
