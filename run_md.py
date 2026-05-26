import sys
from ase.io import read, write, Trajectory
from ase.md import Langevin
from ase import units
from mace.calculators import MACECalculator
from ase.md import MDLogger

cif_file = sys.argv[1]
gpu_id = sys.argv[2]

name = cif_file.replace(".cif", "")

atoms = read(cif_file)

calc = MACECalculator(
    model_paths="MACE-MP-0a.model",
    device=f"cuda:{gpu_id}"
)
atoms.calc = calc

timestep = 0.5 * units.fs
temperature = 298
friction = 0.0005 / units.fs

dyn = Langevin(
    atoms,
    timestep=timestep,
    temperature_K=temperature,
    friction=friction,
)

traj = Trajectory(f"{name}.traj", "w", atoms)
dyn.attach(traj.write, interval=10)

dyn.attach(MDLogger(dyn, atoms, f"{name}.log", header=True, stress=False, peratom=False), interval=10)

dyn.run(80000)

write(f"{name}_final.cif", atoms)

traj_frames = read(f"{name}.traj", index=":")
write(f"{name}.xyz", traj_frames)

print(f"MD finished: {cif_file}")
