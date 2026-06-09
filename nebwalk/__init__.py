"""Public API for nebwalk."""

from __future__ import annotations

from .engine import NEBRunConfig, NEBRunResult, run_neb_calculation
from .forces import compute_neb_forces, variable_spring_constants
from .interpolate import geodesic_interpolate, idpp_interpolate, linear_interpolate
from .neb import NEB
from .output import plot_energy_profile, save_csv, save_trajectory
from .qe import QEParams, make_qe_factory, validate_qe_setup

__version__ = "0.6.0"

__all__ = [
    "NEB",
    "NEBRunConfig",
    "NEBRunResult",
    "run_neb_calculation",
    "linear_interpolate",
    "idpp_interpolate",
    "geodesic_interpolate",
    "compute_neb_forces",
    "variable_spring_constants",
    "plot_energy_profile",
    "save_csv",
    "save_trajectory",
    "QEParams",
    "make_qe_factory",
    "validate_qe_setup",
]
