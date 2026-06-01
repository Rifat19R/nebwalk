"""
nebwalk – minimal Nudged Elastic Band (NEB) implementation for ASE.

Public API
----------
NEB                       Main class. Wraps path images and runs FIRE optimisation.
linear_interpolate        Build an initial path by Cartesian linear interpolation.
idpp_interpolate          Build a path using Image Dependent Pair Potential (IDPP).
compute_neb_forces        Low-level force computation (useful for custom loops).
variable_spring_constants Compute energy-weighted spring constants.
plot_energy_profile       Standalone energy-profile plotter.
save_csv                  Write image energies to CSV.
save_trajectory           Write all images to an ASE .traj file.
"""

from .neb import NEB
from .interpolate import linear_interpolate, idpp_interpolate
from .forces import compute_neb_forces, variable_spring_constants
from .output import plot_energy_profile, save_csv, save_trajectory

__version__ = "0.4.0"
__all__ = [
    "NEB",
    "linear_interpolate",
    "idpp_interpolate",
    "compute_neb_forces",
    "variable_spring_constants",
    "plot_energy_profile",
    "save_csv",
    "save_trajectory",
]
