"""
Example 2: Al adatom diffusion on Al(100)  hollow → bridge → hollow

This is the canonical EMT/NEB benchmark used throughout the ASE documentation.
EMT gives a barrier of ~0.4 eV for this process.

Run:
    python examples/al_diffusion_emt.py

Requires: ase (pip install ase)
"""

import numpy as np
from ase import Atoms
from ase.build import fcc100, add_adsorbate
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from nebwalk import NEB, linear_interpolate

# ---------------------------------------------------------------------------
# Build slab and relax endpoints
# ---------------------------------------------------------------------------

def make_slab_with_adatom(hollow_x_offset=0.0):
    """
    Build an Al(100) slab (2×2, 3 layers) with an Al adatom at a hollow site.
    Returns the Atoms object with EMT calculator attached.
    """
    slab = fcc100("Al", size=(2, 2, 3), vacuum=10.0)
    add_adsorbate(slab, "Al", height=1.7, position="hollow")

    # Fix bottom two layers (tag > 1 → lower layers)
    constraint = FixAtoms(mask=[atom.tag > 1 for atom in slab])
    slab.set_constraint(constraint)
    slab.calc = EMT()
    return slab


print("Building and relaxing initial state (hollow site A) ...")
initial = make_slab_with_adatom()
opt = BFGS(initial, logfile=None)
opt.run(fmax=0.01)
e_initial = initial.get_potential_energy()
print(f"  E_initial = {e_initial:.4f} eV")

# Final state: adatom shifted by (a/2, 0) to the adjacent hollow site
final = initial.copy()
final.calc = EMT()
# The adatom is the last atom added
adatom_idx = len(final) - 1
a_lat = 4.05  # Al lattice constant (Å)
final.positions[adatom_idx] += np.array([a_lat / 2, 0.0, 0.0])

print("Relaxing final state (hollow site B) ...")
opt = BFGS(final, logfile=None)
opt.run(fmax=0.01)
e_final = final.get_potential_energy()
print(f"  E_final   = {e_final:.4f} eV")

# ---------------------------------------------------------------------------
# NEB
# ---------------------------------------------------------------------------

N_IMAGES = 5
print(f"\nBuilding NEB path with {N_IMAGES} intermediate images ...")
images = linear_interpolate(initial, final, n_images=N_IMAGES)
for img in images:
    img.set_constraint(FixAtoms(mask=[atom.tag > 1 for atom in img]))
    img.calc = EMT()

neb = NEB(images, k=0.10, climb=True, climb_delay=50)
converged = neb.optimize(fmax=0.05, max_steps=300)

print(f"\nConverged : {converged}")
print(f"Forward barrier : {neb.get_barrier():.4f} eV")

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

neb.plot("al_diffusion_profile.png", title="Al adatom diffusion on Al(100) [EMT]")
neb.save_csv("al_diffusion_profile.csv")
neb.save_trajectory("al_diffusion.traj")
