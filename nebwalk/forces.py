"""
NEB force computation.

References
----------
Improved tangent:
    Henkelman & Jónsson, J. Chem. Phys. 113, 9978 (2000). DOI: 10.1063/1.1323224
Climbing image:
    Henkelman, Uberuaga, Jónsson, J. Chem. Phys. 113, 9901 (2000). DOI: 10.1063/1.1329672
Variable spring constants:
    Lindh, Helgaker et al., J. Phys. Chem. 100, 9979 (1996). DOI: 10.1021/jp953468a
"""

import numpy as np
from ase.geometry import find_mic


def _mic_disp(dr, cell, pbc):
    """Apply MIC to inter-image displacement. No-op for non-periodic systems."""
    if not np.any(pbc):
        return dr
    dr_mic, _ = find_mic(dr, cell, pbc)
    return dr_mic


def _improved_tangent(dr_fwd, dr_bwd, d_fwd, d_bwd, E_prev, E_curr, E_next):
    """Improved tangent estimate (Henkelman & Jónsson 2000, Eqs. 8-11)."""
    if E_next > E_curr > E_prev:
        return dr_fwd / d_fwd

    if E_prev > E_curr > E_next:
        return dr_bwd / d_bwd

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


def variable_spring_constants(energies, k_max, k_min):
    """
    Compute energy-weighted spring constants for each inter-image spring.

    Springs near the saddle point (high energy) get k_max; springs far from
    it (low energy) get k_min. This concentrates images near the transition
    state, giving better barrier resolution without adding more images.

    Formula (per spring i connecting images i and i+1):
        E_spring[i] = max(E[i], E[i+1])
        k[i] = k_max - (k_max - k_min) * (E_ref - E_spring[i]) / (E_ref - E_low)
    where E_ref = max(all energies), E_low = min(endpoint energies).
    Values are clipped to [k_min, k_max].

    Parameters
    ----------
    energies : sequence of float, length N
        Potential energies of all images including endpoints.
    k_max : float
        Spring constant at the saddle point (eV/Å²).
    k_min : float
        Spring constant far from the saddle point (eV/Å²). Must be < k_max.

    Returns
    -------
    k : ndarray, shape (N-1,)
        Spring constant for each inter-image spring.
    """
    E = np.array(energies, dtype=float)
    n = len(E)

    E_spring = np.maximum(E[:-1], E[1:])   # (N-1,)
    E_ref = E.max()
    E_low = min(E[0], E[-1])

    dE = E_ref - E_low
    if dE < 1e-10:
        return np.full(n - 1, k_max)

    k = k_max - (k_max - k_min) * (E_ref - E_spring) / dE
    return np.clip(k, k_min, k_max)


def compute_neb_forces(images, k, climb=False, climb_index=None,
                       energies=None, forces=None):
    """
    Compute NEB forces for all images.

    Parameters
    ----------
    images : list of ase.Atoms
        Full path, each with a calculator attached.
    k : float or ndarray of shape (N-1,)
        Spring constant(s) in eV/Å². Scalar → uniform springs. Array →
        one value per inter-image spring (use variable_spring_constants).
    climb : bool
        If True, apply climbing image force to image at ``climb_index``.
    climb_index : int or None
        Index of the climbing image (highest-energy movable image).
    energies : list of float or None
        Pre-computed potential energies for all images. If None, energies
        are fetched from each image's calculator.
    forces : list of ndarray or None, length N
        Pre-computed forces for all images. None entries trigger a
        calculator call for that image. Pass pre-computed values (from
        parallel evaluation) to avoid redundant calculator calls.

    Returns
    -------
    forces_out : list of ndarray, each shape (N_atoms, 3)
        NEB forces. Endpoints are zero arrays.
    """
    n = len(images)

    if energies is None:
        energies = [img.get_potential_energy() for img in images]

    if np.isscalar(k):
        k_arr = np.full(n - 1, float(k))
    else:
        k_arr = np.asarray(k, dtype=float)
        if len(k_arr) != n - 1:
            raise ValueError(
                f"k array length {len(k_arr)} must equal n_images-1 = {n - 1}."
            )

    forces_out = []
    for i in range(n):
        if i == 0 or i == n - 1:
            forces_out.append(np.zeros_like(images[i].positions))
            continue

        pbc  = images[i].pbc
        cell = images[i].cell

        dr_fwd = _mic_disp(images[i + 1].positions - images[i].positions, cell, pbc)
        dr_bwd = _mic_disp(images[i].positions - images[i - 1].positions, cell, pbc)

        d_fwd = np.linalg.norm(dr_fwd)
        d_bwd = np.linalg.norm(dr_bwd)

        if d_fwd < 1e-12 or d_bwd < 1e-12:
            forces_out.append(np.zeros_like(images[i].positions))
            continue

        tau = _improved_tangent(
            dr_fwd, dr_bwd, d_fwd, d_bwd,
            energies[i - 1], energies[i], energies[i + 1],
        )

        # Use pre-computed forces if available, otherwise call calculator
        if forces is not None and forces[i] is not None:
            f_pot = forces[i]
        else:
            f_pot = images[i].get_forces()

        if climb and i == climb_index:
            proj    = (f_pot * tau).sum()
            f_total = f_pot - 2.0 * proj * tau
        else:
            k_fwd      = k_arr[i]
            k_bwd      = k_arr[i - 1]
            f_spring   = (k_fwd * d_fwd - k_bwd * d_bwd) * tau
            proj       = (f_pot * tau).sum()
            f_pot_perp = f_pot - proj * tau
            f_total    = f_spring + f_pot_perp

        forces_out.append(f_total)

    return forces_out
