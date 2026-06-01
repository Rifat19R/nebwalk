"""
Ag vacancy migration in FCC Ag — EMT calculator.

Reference barrier: ~0.66 eV (DFT-PBE).

Usage
-----
    python examples/ag_vacancy_emt.py
"""

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS

from nebwalk import NEB, linear_interpolate


# ---------------------------------------------------------------------------
# 1. Build supercell and create vacancy
# ---------------------------------------------------------------------------

ag_bulk = bulk("Ag", "fcc", a=4.085, cubic=True)
supercell = ag_bulk.repeat(2)        # 2x2x2 → 32 atoms

vacancy_index = 0
vacancy_pos = supercell.positions[vacancy_index].copy()
del supercell[vacancy_index]         # 31 atoms remain

# ---------------------------------------------------------------------------
# 2. Find nearest neighbour to vacancy site
# ---------------------------------------------------------------------------

distances = np.linalg.norm(supercell.positions - vacancy_pos, axis=1)
nn_index = int(np.argmin(distances))
nn_dist = distances[nn_index]
print(f"Migrating atom : index {nn_index}")
print(f"NN distance    : {nn_dist:.4f} Å  (expected ~2.889 Å for FCC Ag)")

# ---------------------------------------------------------------------------
# 3. Build and relax endpoints
# ---------------------------------------------------------------------------

initial = supercell.copy()
initial.calc = EMT()

final = supercell.copy()
final.positions[nn_index] = vacancy_pos
final.calc = EMT()

BFGS(initial, logfile=None).run(fmax=0.02)
BFGS(final,   logfile=None).run(fmax=0.02)

print(f"\nInitial energy : {initial.get_potential_energy():.4f} eV")
print(f"Final energy   : {final.get_potential_energy():.4f} eV")
print(f"ΔE             : "
      f"{final.get_potential_energy() - initial.get_potential_energy():.5f} eV")

# ---------------------------------------------------------------------------
# 4. Build NEB path
# ---------------------------------------------------------------------------

images = linear_interpolate(initial, final, n_images=7)
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
ref = 0.66
print(f"\nConverged      : {converged}")
print(f"Forward barrier: {barrier:.4f} eV")
print(f"Reference (DFT): ~{ref} eV")
print(f"Error          : {abs(barrier - ref) / ref * 100:.1f}%")

neb.plot("ag_vacancy_emt_profile.png",
         title="Ag vacancy migration in FCC Ag [EMT, CI-NEB]")
neb.save_csv("ag_vacancy_emt_profile.csv")
neb.save_trajectory("ag_vacancy_emt.traj")
