"""Tests for NEB user-facing helpers."""

from __future__ import annotations

from ase import Atoms
from ase.calculators.emt import EMT
from ase.io import write

from nebwalk import NEB, NEBRunConfig, run_neb_calculation


def test_from_trajectory_attaches_fresh_calculators(tmp_path):
    images = [
        Atoms("Al", positions=[[0.0, 0.0, 0.0]]),
        Atoms("Al", positions=[[0.1, 0.0, 0.0]]),
        Atoms("Al", positions=[[0.2, 0.0, 0.0]]),
    ]
    path = tmp_path / "restart.traj"
    write(path, images)

    neb = NEB.from_trajectory(str(path), calculator_factory=EMT)

    assert len(neb.images) == 3
    assert all(img.calc is not None for img in neb.images)
    assert len({id(img.calc) for img in neb.images}) == 3


def test_run_neb_calculation_engine_smoke():
    initial = Atoms("Al", positions=[[0.0, 0.0, 0.0]])
    final = Atoms("Al", positions=[[0.1, 0.0, 0.0]])
    result = run_neb_calculation(
        initial,
        final,
        calculator_factory=EMT,
        config=NEBRunConfig(
            n_images=1,
            interpolation="linear",
            max_steps=1,
            fmax=10.0,
            verbose=False,
        ),
    )
    assert result.converged
    assert result.neb.get_barrier() == result.barrier
