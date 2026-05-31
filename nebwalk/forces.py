"""
NEB force computation.

References
----------
Improved tangent:
    Henkelman & Jónsson, J. Chem. Phys. 113, 9978 (2000). DOI: 10.1063/1.1323224
Climbing image:
    Henkelman, Uberuaga, Jónsson, J. Chem. Phys. 113, 9901 (2000). DOI: 10.1063/1.1329672
"""

import numpy as np
from ase.geometry import find_mic


def _mic_disp(dr, cell, pbc):
    """
    Apply the Minimum Image Convention (MIC) to inter-image displacement vectors.

    For non-periodic systems (all pbc=False), returns dr unchanged.
    For periodic systems, wraps each atomic displacement to the shortest
    vector consistent with the periodic boundary, preventing spurious
    spring forces when atoms cross cell boundaries.

    Parameters
    ----------
    dr : ndarray, shape (N_atoms, 3)
        Raw displacement vectors (R_j - R_i) in Cartesian coordinates.
    cell : ase.cell.Cell or ndarray, shape (3, 3)
        Unit cell vectors.
    pbc : array-like of bool, shape (3,)
        Periodic boundary conditions for each cell direction.

    Returns
    -------
    dr_mic : ndarray, shape (N_atoms, 3)
        MIC-corrected displacement vectors.
    """
    if not np.any(pbc):
        return dr
    dr_mic, _ = find_mic(dr, cell, pbc)
    return dr_mic


def _improved_tangent(dr_fwd, dr_bwd, d_fwd, d_bwd, E_prev, E_curr, E_next):
    """
    Improved tangent estimate for NEB.

    Selects between forward, backward, or energy-weighted bisection based on
    local energy ordering. This avoids kinks at saddle points.

    Parameters
    ----------
    dr_fwd : ndarray, shape (N, 3)
        R_{i+1} - R_i (MIC-corrected, unnormalized).
    dr_bwd : ndarray, shape (N, 3)
        R_i - R_{i-1} (MIC-corrected, unnormalized).
    d_fwd, d_bwd : float
        Euclidean norms of dr_fwd, dr_bwd.
    E_prev, E_curr, E_next : float
        Potential energies of images i-1, i, i+1.

    Returns
    -------
    tau : ndarray, shape (N, 3)
        Normalised tangent vector.
    """
    if E_next > E_curr > E_prev:
        return dr_fwd / d_fwd

    if E_prev > E_curr > E_next:
        return dr_bwd / d_bwd

    # Local extremum: energy-weighted bisection (Eqs. 10-11, Henkelman 2000)
    dE_fwd = abs(E_next - E_curr)
    dE_bwd = abs(E_prev - E_curr)
    dE_max = max(dE_fwd, dE_bwd)
    dE_min = min(dE_fwd, dE_bwd)

    if E_next >= E_prev:
        tau = dr_fwd * dE_max + dr_bwd * dE_min
    else:
        tau = dr_fwd * dE_min + dr_bwd * dE_max

    norm = np.linalg.norm(tau)
    if norm < 1e-12:
        return dr_fwd / d_fwd
    return tau / norm


def compute_neb_forces(images, k, climb=False, climb_index=None):
    """
    Compute NEB forces for all images.

    Supports both gas-phase (pbc=False) and periodic (pbc=True) systems.
    For periodic images, the Minimum Image Convention (MIC) is applied to
    inter-image displacements so that atoms crossing cell boundaries are
    handled correctly.

    Parameters
    ----------
    images : list of ase.Atoms
        Full path, each with a calculator attached.
    k : float
        Spring constant in eV/Å².
    climb : bool
        If True, apply climbing image force to image at ``climb_index``.
    climb_index : int or None
        Index of the climbing image (highest-energy movable image).

    Returns
    -------
    forces : list of ndarray, each shape (N_atoms, 3)
        NEB forces. Endpoints are zero arrays.
    """
    n = len(images)
    energies = [img.get_potential_energy() for img in images]

    forces = []
    for i in range(n):
        if i == 0 or i == n - 1:
            forces.append(np.zeros_like(images[i].positions))
            continue

        pbc  = images[i].pbc
        cell = images[i].cell

        # MIC-corrected inter-image displacements
        dr_fwd = _mic_disp(images[i + 1].positions - images[i].positions, cell, pbc)
        dr_bwd = _mic_disp(images[i].positions - images[i - 1].positions, cell, pbc)

        d_fwd = np.linalg.norm(dr_fwd)
        d_bwd = np.linalg.norm(dr_bwd)

        if d_fwd < 1e-12 or d_bwd < 1e-12:
            forces.append(np.zeros_like(images[i].positions))
            continue

        tau = _improved_tangent(
            dr_fwd, dr_bwd, d_fwd, d_bwd,
            energies[i - 1], energies[i], energies[i + 1],
        )

        f_pot = images[i].get_forces()  # ASE returns -∇E directly

        if climb and i == climb_index:
            proj    = (f_pot * tau).sum()
            f_total = f_pot - 2.0 * proj * tau
        else:
            f_spring    = k * (d_fwd - d_bwd) * tau
            proj        = (f_pot * tau).sum()
            f_pot_perp  = f_pot - proj * tau
            f_total     = f_spring + f_pot_perp

        forces.append(f_total)

    return forces
