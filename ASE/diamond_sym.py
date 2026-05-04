from ase import Atoms
from ase.io import read, write
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from ase.constraints import UnitCellFilter, FixSymmetry
import spglib
import re
import numpy as np

USE_SYMMETRY = True

atoms0 = read("diamond.cif")

calc = MACECalculator(
    model_path="/Users/liyi/Desktop/Part 2 project/MACE-MP-0a.model",
    device="cpu"
)

atoms0.calc = calc

if USE_SYMMETRY:
    lattice = np.array(atoms0.get_cell(), dtype=np.float64, order='C')
    frac = np.mod(atoms0.get_scaled_positions(), 1.0).astype(np.float64, copy=False, order='C')
    types = np.asarray(atoms0.get_atomic_numbers(), dtype=np.intc)
    cell = (lattice, frac, types)
    
    sym = spglib.get_symmetry(cell, symprec=1e-5, angle_tolerance=-1.0)
    rotations = sym["rotations"]
    translations = sym["translations"]
    spacegroup = spglib.spglib.get_spacegroup_type_from_symmetry(rotations=rotations, translations=translations, lattice=lattice)
    
    print(f"Initial space group: {spacegroup.number}, {spacegroup.international_short}")
    print(f"Number of symmetry operations: {len(rotations)}")
    
    fix_symmetry = FixSymmetry(atoms0, symprec=1e-5)
    atoms0.set_constraint([fix_symmetry])
    print("Symmetry constraints applied")
else:
    print("Running optimization without symmetry constraints")

ucf = UnitCellFilter(atoms0, mask=[1,1,1,0,0,0])
dyn = BFGS(ucf, logfile="-", trajectory="vc.traj")
dyn.run(fmax=0.005, steps=1000)

write("diamond_relaxed_sym_MACE.cif", atoms0)

E = atoms0.get_potential_energy()
print("Final energy (eV):", E)
print("Per-atom energy (eV/atom):", E / len(atoms0))

path = "diamond_relaxed_sym_MACE.cif"

atoms = read(path, index=0)
lattice = np.array(atoms.get_cell(), dtype=np.float64, order='C')
frac = np.mod(atoms.get_scaled_positions(), 1.0).astype(np.float64, copy=False, order='C')
types = np.asarray(atoms.get_atomic_numbers(), dtype=np.intc)
cell = (lattice, frac, types)

sym = spglib.get_symmetry(cell, symprec=1e-5, angle_tolerance=-1.0)
rotations = sym["rotations"]
translations = sym["translations"]
spacegroup = spglib.spglib.get_spacegroup_type_from_symmetry(rotations=rotations, translations=translations, lattice=lattice)

n_ops = len(rotations)

filename = "symmetry_operation.txt"
with open(filename, "w") as f:
    f.write(f"# Symmetry operations from spglib\n")
    f.write(f"# Total number: {n_ops}\n\n")
    for i, (rot, trans) in enumerate(zip(rotations, translations), start=1):
        for row in rot:
            f.write(" " + " ".join(f"{x:3d}" for x in row) + "\n")
        f.write("\n")
    f.write(f"Space group: {spacegroup.number}, {spacegroup.international_short}")

new_name = spacegroup.international_short
new_number = spacegroup.number

with open(path, "r") as f:
    text = f.read()
pattern = r"^(_space_group_name_H-M_alt\s+)(['\"])(.*)(['\"])"
new_text = re.sub(pattern, rf"\1'{new_name}'", text, flags=re.MULTILINE)
with open(path, "w") as f:
    f.write(new_text)

with open(path, "r") as f:
    text1 = f.read()
pattern1 = r"^(_space_group_IT_number\s+)(['\"]?)(\d+)(['\"]?)"
new_text1 = re.sub(pattern1, rf"\1'{new_number}'", text1, flags=re.MULTILINE)
with open(path, "w") as f:
    f.write(new_text1)

print(f"Space group: {spacegroup.number}, {spacegroup.international_short}")
