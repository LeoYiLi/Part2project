from ase.io import read, write, Trajectory
from ase.md import Langevin
from ase import units
from mace.calculators import MACECalculator
from ase.md import MDLogger
from ase.md.npt import NPT

atoms = read("YBa2Cu3O7_553.cif")

calc = MACECalculator(
    model_path="/Users/liyi/Desktop/Part 2 project/MACE-MP-0a.model",
    device="cpu"
)
atoms.set_calculator(calc)

temperatures = [20, 77, 300]
nsteps = 50000  # 50 ps

for T in temperatures:
    dyn = Langevin(
        atoms,
        timestep=1.0 * units.fs,
        temperature_K=T,
        friction=0.01 / units.fs,
    )

    traj = Trajectory(f"yba553_md_{T}K.traj", "w", atoms)
    dyn.attach(traj.write, interval=100)

    dyn.attach(
        MDLogger(
            dyn, atoms, f"md_{T}K.log",
            header=True, stress=False, peratom=False
        ),
        interval=100
    )

    dyn.run(nsteps)

traj = Trajectory("yba553_md_30K.traj", "w", atoms)
dyn.attach(traj.write, interval=10)   

dyn.attach(MDLogger(dyn, atoms, "md.log", header=True, stress=False, peratom=False), interval=10)

dyn.run(5000)


write("YBa2Cu3O7_553_md.cif", atoms)

print("MD simulation finished!")