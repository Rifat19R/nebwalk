"""
Example 3: Ethane C–C torsional barrier using Egret-1t

Reaction : staggered C2H6 → [eclipsed TS] → staggered C2H6  (120° CH3 rotation)
Expected barrier : ~0.126 eV  (12.2 kJ/mol, experimental)

Two fixes vs the naive version:
  1. All images get an explicit 20 Å cubic cell + centered.
     MACE/matscipy requires a cell to build the neighbour list.
     Without it, matscipy crashes with "Failed to allocate seed array".
  2. Images are pre-built at exact dihedral angles (60°→180°), NOT by
     linear Cartesian interpolation.  Cartesian interpolation for a torsion
     moves atoms through each other → fmax > 10^6 eV/Å → FIRE diverges.

Run:
    python examples/ethane_egret.py
"""

import numpy as np
import torch
from ase import Atoms
from ase.optimize import BFGS
from mace.calculators import MACECalculator
from nebwalk import NEB

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
MODEL_PATH = "EGRET_1T.model"
DEVICE     = "cuda" if torch.cuda.is_available() else "cpu"
N_IMAGES   = 7          # intermediate images (7 gives good resolution)
K_SPRING   = 0.10       # eV/Å²

# ---------------------------------------------------------------------------
# Geometry parameters (D3d ethane)
# ---------------------------------------------------------------------------
CC   = 0.770                                  # half C-C bond length (Å)
r_CH = 1.090                                  # C-H bond length (Å)
_dz  = abs(r_CH * np.cos(np.radians(111.2))) # H z-offset from C = 0.394 Å
_rl  = r_CH * np.sin(np.radians(111.2))      # H lateral radius  = 1.016 Å

CELL = np.diag([20.0, 20.0, 20.0])           # 20 Å vacuum box — required for MACE


def ethane(phi_C2_deg: float) -> Atoms:
    """
    Return a C2H6 Atoms object with an explicit cell, centered in the box.

    phi_C2_deg : dihedral angle of C2 methyl relative to C1 methyl (degrees).
        60°  → staggered  (start / end)
        120° → eclipsed   (approximate TS)
    """
    phi = np.radians(phi_C2_deg)
    pos = [
        [0.0, 0.0, -CC],  # C1
        [0.0, 0.0, +CC],  # C2
    ]
    for k in range(3):                        # C1 methyls at 0°, 120°, 240°
        a = k * (2 * np.pi / 3)
        pos.append([_rl * np.cos(a), _rl * np.sin(a), -CC - _dz])
    for k in range(3):                        # C2 methyls starting at phi_C2
        a = phi + k * (2 * np.pi / 3)
        pos.append([_rl * np.cos(a), _rl * np.sin(a), +CC + _dz])

    atoms = Atoms("C2H6", positions=pos)
    atoms.set_cell(CELL)
    atoms.pbc = False
    atoms.center()                            # move molecule to cell centre
    return atoms


# ---------------------------------------------------------------------------
# Calculator factory — one instance per image (MACE is not thread-safe)
# ---------------------------------------------------------------------------
def make_calc():
    kw = dict(device=DEVICE, default_dtype="float32")
    try:
        return MACECalculator(model_paths=MODEL_PATH, **kw)
    except TypeError:
        return MACECalculator(model_path=MODEL_PATH, **kw)


# ---------------------------------------------------------------------------
# Relax endpoints
# ---------------------------------------------------------------------------
print(f"Device   : {DEVICE}")
print(f"Model    : {MODEL_PATH}\n")

print("Relaxing initial state (staggered, phi = 60°) ...")
initial = ethane(60.0)
initial.calc = make_calc()
opt = BFGS(initial, logfile=None)
opt.run(fmax=0.01)
E0 = initial.get_potential_energy()
print(f"  E_initial = {E0:.4f} eV")

print("Relaxing final state   (staggered, phi = 180°) ...")
final = ethane(180.0)
final.calc = make_calc()
opt = BFGS(final, logfile=None)
opt.run(fmax=0.01)
E1 = final.get_potential_energy()
print(f"  E_final   = {E1:.4f} eV")
print(f"  ΔE        = {E1 - E0:.5f} eV  (should be ~0 by symmetry)")

# ---------------------------------------------------------------------------
# Pre-build images at evenly-spaced dihedral angles
#
# Linear Cartesian interpolation would move atoms through each other for a
# torsional coordinate.  Instead we construct each image at a known angle:
#   phi = 60°, 74°, 88°, ..., 166°, 180°
# ---------------------------------------------------------------------------
print(f"\nBuilding {N_IMAGES} intermediate images at exact dihedral angles ...")
phi_values = np.linspace(60.0, 180.0, N_IMAGES + 2)   # includes endpoints

images = []
for i, phi in enumerate(phi_values):
    if i == 0:
        img = initial           # already relaxed endpoint
    elif i == N_IMAGES + 1:
        img = final             # already relaxed endpoint
    else:
        img = ethane(phi)       # pre-built at correct dihedral
        img.calc = make_calc()
    images.append(img)

print("  Phi values (°):", [f"{p:.1f}" for p in phi_values])

# ---------------------------------------------------------------------------
# NEB optimisation
# ---------------------------------------------------------------------------
neb = NEB(images, k=K_SPRING, climb=True, climb_delay=60)
converged = neb.optimize(fmax=0.05, max_steps=400)

print(f"\nConverged       : {converged}")
print(f"Forward barrier : {neb.get_barrier():.4f} eV")
print(f"Literature      : ~0.126 eV  (12.2 kJ/mol)")

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
neb.plot("ethane_egret_profile.png",
         title="Ethane C–C torsion [Egret-1t, CI-NEB]")
neb.save_csv("ethane_egret_profile.csv")
neb.save_trajectory("ethane_egret.traj")

print("\nDone.  View trajectory with:  ase gui ethane_egret.traj")