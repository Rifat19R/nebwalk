"""Au vacancy migration in FCC Au with EMT.

Reference context:
    FCC Au vacancy migration barriers are commonly reported near 0.68 eV from
    DFT-PBE literature.  This EMT example is a calculator sanity benchmark for
    nebwalk, not a publication-quality Au diffusion calculation.

Run:
    python examples/au_vacancy_emt.py
"""

from __future__ import annotations

import numpy as np
from ase.build import bulk
from ase.calculators.emt import EMT
from ase.optimize import BFGS

from nebwalk import NEB, linear_interpolate

NAME = "au_vacancy_emt"
REFERENCE_BARRIER_EV = 0.68


def make_calc() -> EMT:
    return EMT()


def make_relaxed_endpoints():
    full = bulk("Au", "fcc", a=4.08, cubic=True).repeat((2, 2, 2))
    vacancy_pos = full.positions[0].copy()

    base = full.copy()
    del base[0]

    distances = np.linalg.norm(base.positions - vacancy_pos, axis=1)
    nn_index = int(np.argmin(distances))
    nn_distance = float(distances[nn_index])

    initial = base.copy()
    initial.calc = make_calc()

    final = base.copy()
    final.positions[nn_index] = vacancy_pos
    final.calc = make_calc()

    BFGS(initial, logfile=None).run(fmax=0.02)
    BFGS(final, logfile=None).run(fmax=0.02)

    print(f"Migrating atom : index {nn_index}")
    print(f"NN distance    : {nn_distance:.4f} Angstrom")
    endpoint_delta = final.get_potential_energy() - initial.get_potential_energy()
    print(f"Endpoint dE    : {endpoint_delta:.6f} eV")
    return initial, final


def main() -> None:
    initial, final = make_relaxed_endpoints()

    images = linear_interpolate(initial, final, n_images=7)
    for image in images:
        image.calc = make_calc()

    neb = NEB(images, k=0.1, climb=True, climb_delay=30)
    converged = neb.optimize(fmax=0.05, max_steps=300)

    barrier = neb.get_barrier()
    error = abs(barrier - REFERENCE_BARRIER_EV) / REFERENCE_BARRIER_EV * 100.0

    print(f"\nConverged      : {converged}")
    print(f"Forward barrier: {barrier:.4f} eV")
    print(f"Reference (DFT): ~{REFERENCE_BARRIER_EV:.2f} eV")
    print(f"Error          : {error:.1f}%")

    neb.plot(
        f"{NAME}_profile.png",
        title="Au vacancy migration in FCC Au [EMT, CI-NEB]",
    )
    neb.save_csv(f"{NAME}_profile.csv")
    neb.save_trajectory(f"{NAME}.traj")


if __name__ == "__main__":
    main()
