"""
Cu vacancy migration in FCC Cu — EMT calculator.

Reference barrier: ~0.70 eV (DFT-PBE).
EMT underestimates this (~0.50 eV expected).

Usage
-----
    python examples/cu_vacancy_emt.py
"""

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS

from nebwalk import NEB, linear_interpolate

# ---------------------------------------------------------------------------
# 1. Build supercell and create vacancy
# ---------------------------------------------------------------------------

cu_bulk = bulk("Cu", "fcc", a=3.615, cubic=True)
supercell = cu_bulk.repeat(2)        # 2x2x2 → 32 atoms

vacancy_index = 0
vacancy_pos = supercell.positions[vacancy_index].copy()
del supercell[vacancy_index]         # 31 atoms remain

# ---------------------------------------------------------------------------
# 2. Find nearest neighbour to vacancy site (will migrate into vacancy)
# ---------------------------------------------------------------------------

distances = np.linalg.norm(supercell.positions - vacancy_pos, axis=1)
nn_index = int(np.argmin(distances))
nn_dist = distances[nn_index]
print(f"Migrating atom : index {nn_index}")
print(f"NN distance    : {nn_dist:.4f} Å  (expected ~2.556 Å for FCC Cu)")

# ---------------------------------------------------------------------------
# 3. Build initial and final states — no constraints
# ---------------------------------------------------------------------------

initial = supercell.copy()
initial.calc = EMT()

final = supercell.copy()
final.positions[nn_index] = vacancy_pos   # atom jumps into vacancy
final.calc = EMT()

# Relax both endpoints fully
BFGS(initial, logfile=None).run(fmax=0.02)
BFGS(final,   logfile=None).run(fmax=0.02)

print(f"\nInitial energy : {initial.get_potential_energy():.4f} eV")
print(f"Final energy   : {final.get_potential_energy():.4f} eV")
print(f"ΔE             : "
      f"{final.get_potential_energy() - initial.get_potential_energy():.5f} eV")

# ---------------------------------------------------------------------------
# 4. Build NEB path
# ---------------------------------------------------------------------------

n_images = 7
images = linear_interpolate(initial, final, n_images=n_images)
for img in images:
    img.calc = EMT()

# ---------------------------------------------------------------------------
# 5. Run NEB
# ---------------------------------------------------------------------------

neb = NEB(images, k=0.1, climb=True, climb_delay=30)
converged = neb.optimize(fmax=0.05, max_steps=300)

# ---------------------------------------------------------------------------
# 6. Output
# ---------------------------------------------------------------------------

barrier = neb.get_barrier()
ref = 0.70
print(f"\nConverged      : {converged}")
print(f"Forward barrier: {barrier:.4f} eV")
print(f"Reference (DFT): ~{ref} eV")
print(f"Error          : {abs(barrier - ref) / ref * 100:.1f}%")

neb.plot("cu_vacancy_emt_profile.png",
         title="Cu vacancy migration in FCC Cu [EMT, CI-NEB]")
neb.save_csv("cu_vacancy_emt_profile.csv")
neb.save_trajectory("cu_vacancy_emt.traj")
