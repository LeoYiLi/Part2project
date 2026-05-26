from ase.io import read, write, Trajectory
from ase.md import Langevin
from ase import units
from mace.calculators import MACECalculator
from ase.md import MDLogger

atoms = read("YBa2Cu3O7_553_Frenkel4.cif")

calc = MACECalculator(
    model_paths="/home/leo/Part2project/Part2project/MACE-MP-0a.model",
    device="cuda"
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

traj = Trajectory("YBa2Cu3O7_553_Frenkel4.traj", "w", atoms)
dyn.attach(traj.write, interval=1)  
# initially 50 

dyn.attach(MDLogger(dyn, atoms, "YBa2Cu3O7_553_Frenkel4.log", header=True, stress=False, peratom=False), interval=1)

dyn.run(5000)
#initially 100000


write("Frankel_md.cif", atoms)

traj_frames = read("YBa2Cu3O7_553_Frenkel3.traj", index=":")
write("YBa2Cu3O7_553_Frenkel3.xyz", traj_frames)

print("MD simulation finished!")
