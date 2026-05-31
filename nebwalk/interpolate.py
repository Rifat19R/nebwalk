"""
Path interpolation utilities for NEB.

Two methods are provided:

linear_interpolate
    Cartesian linear interpolation. Fast but unsuitable for paths involving
    large rotations or torsional motion — atoms may pass through each other.
    Supports PBC via the minimum image convention.

idpp_interpolate
    Image Dependent Pair Potential (IDPP) interpolation. Minimises a
    weighted pairwise distance objective to produce a chemically sensible
    initial path. Recommended for organic molecules and any system where
    linear interpolation produces steric clashes.

    Note: PBC support for IDPP is not yet implemented. Use linear_interpolate
    for periodic systems until MIC support is added in a future release.

Reference for IDPP:
    Smidstrup, Pedersen, Stokbro, Jónsson,
    J. Chem. Phys. 141, 214106 (2014). DOI: 10.1063/1.4878664
"""

import numpy as np
from scipy.optimize import minimize
from ase.geometry import find_mic


# ---------------------------------------------------------------------------
# Linear interpolation (MIC-aware, unchanged)
# ---------------------------------------------------------------------------

def linear_interpolate(start, end, n_images):
    """
    Create linearly interpolated images between start and end.

    The calculator is NOT copied to intermediate images; attach one before
    calling ``NEB.optimize()``.

    Parameters
    ----------
    start : ase.Atoms
        Initial state (endpoint, not modified).
    end : ase.Atoms
        Final state (endpoint, not modified).
    n_images : int
        Number of intermediate images (excluding endpoints).
        Total path length will be n_images + 2.

    Returns
    -------
    images : list of ase.Atoms, length n_images + 2
        [start_copy, img_1, ..., img_n, end_copy].

    Raises
    ------
    ValueError
        If start and end have different numbers of atoms or species.
    """
    if len(start) != len(end):
        raise ValueError(
            f"start and end must have the same number of atoms "
            f"({len(start)} vs {len(end)})."
        )
    if list(start.symbols) != list(end.symbols):
        raise ValueError(
            "start and end must have the same atomic species in the same order."
        )
    if n_images < 1:
        raise ValueError("n_images must be >= 1.")

    dr_total = end.positions - start.positions
    if np.any(start.pbc):
        dr_total, _ = find_mic(dr_total, start.cell, start.pbc)

    images = [start.copy()]

    for m in range(1, n_images + 1):
        alpha = m / (n_images + 1)
        img = start.copy()
        img.set_positions(start.positions + alpha * dr_total)
        images.append(img)

    images.append(end.copy())
    return images


# ---------------------------------------------------------------------------
# IDPP helpers (private)
# ---------------------------------------------------------------------------

def _pairwise(positions):
    """
    Compute pairwise difference vectors and distances.

    Returns
    -------
    diff : (N, N, 3)  diff[i, j] = positions[i] - positions[j]
    dist : (N, N)     Euclidean distances; diagonal is 0.0
    """
    diff = positions[:, None, :] - positions[None, :, :]   # (N, N, 3)
    dist = np.sqrt((diff ** 2).sum(axis=-1))               # (N, N)
    return diff, dist


def _idpp_obj_and_grad(flat_pos, n_atoms, d_target, weights):
    """
    IDPP objective and analytical gradient for one image.

    Objective (Smidstrup et al., eq. 2):
        S = sum_{i<j} w_ij * (d_ij - d_target_ij)^2
          = 0.5 * sum_{i != j} w_ij * (d_ij - d_target_ij)^2

    Gradient (analytical, accounting for both (i,j) and (j,i) contributions):
        dS/dr_k = 2 * sum_{j != k} w_kj * delta_kj / d_kj * (r_k - r_j)

    where delta_ij = d_ij - d_target_ij.

    Parameters
    ----------
    flat_pos : (3N,) current atomic positions
    n_atoms  : int
    d_target : (N, N) target distance matrix (fixed throughout)
    weights  : (N, N) = 1/d_target^4, zeros on diagonal

    Returns
    -------
    obj  : float
    grad : (3N,)
    """
    positions = flat_pos.reshape(n_atoms, 3)
    diff, dist = _pairwise(positions)

    # Avoid division by zero on diagonal (diagonal is excluded via weights=0)
    dist_safe = np.where(dist > 1e-12, dist, 1.0)

    delta = dist - d_target                                    # (N, N)

    # Objective: upper-triangle sum = 0.5 * full symmetric sum
    obj = 0.5 * np.sum(weights * delta ** 2)

    # Gradient: both (k,j) and (j,k) terms contribute equally → factor of 2
    factor = 2.0 * weights * delta / dist_safe                # (N, N)
    np.fill_diagonal(factor, 0.0)
    grad = np.einsum('ij,ijk->ik', factor, diff)              # (N, 3)

    return obj, grad.ravel()


# ---------------------------------------------------------------------------
# IDPP interpolation (public)
# ---------------------------------------------------------------------------

def idpp_interpolate(start, end, n_images, max_iter=500, tol=1e-6):
    """
    IDPP interpolation between start and end.

    Generates a smoother, more chemically sensible initial path than linear
    interpolation by preserving pairwise distance structure across images.
    Recommended whenever the reaction involves significant torsional motion
    or risk of atomic overlap along a linear path.

    Algorithm
    ---------
    1. Generate linearly interpolated starting positions for each image.
    2. For image k, define target pairwise distances by linear interpolation
       between the endpoint distance matrices:
           d_target[k] = d_start + (k / (N-1)) * (d_end - d_start)
    3. Optimise each image independently (L-BFGS-B) to minimise the IDPP
       objective — a 1/d^4-weighted sum of squared distance deviations.
    4. Endpoints are returned as copies of start/end and are never moved.

    Parameters
    ----------
    start : ase.Atoms
        Initial state (endpoint, not modified).
    end : ase.Atoms
        Final state (endpoint, not modified).
    n_images : int
        Number of intermediate images (excluding endpoints).
    max_iter : int
        Maximum L-BFGS-B iterations per image (default 500).
        L-BFGS-B is appropriate here — the IDPP objective IS a conservative
        potential, unlike the NEB force field.
    tol : float
        Convergence tolerance for L-BFGS-B (default 1e-6).

    Returns
    -------
    images : list of ase.Atoms, length n_images + 2
        No calculator is attached to any image.

    Raises
    ------
    ValueError
        If start and end have different numbers of atoms or species,
        or if n_images < 1.
    NotImplementedError
        If start has periodic boundary conditions (PBC). Use
        linear_interpolate for periodic systems.

    Notes
    -----
    PBC support is not yet implemented for IDPP. This will be added
    alongside minimum image convention (MIC) support in a future release.
    """
    if len(start) != len(end):
        raise ValueError(
            f"start and end must have the same number of atoms "
            f"({len(start)} vs {len(end)})."
        )
    if list(start.symbols) != list(end.symbols):
        raise ValueError(
            "start and end must have the same atomic species in the same order."
        )
    if n_images < 1:
        raise ValueError("n_images must be >= 1.")
    if np.any(start.pbc):
        raise NotImplementedError(
            "IDPP interpolation does not yet support periodic boundary conditions. "
            "Use linear_interpolate for periodic systems."
        )

    n_total = n_images + 2
    n_atoms = len(start)

    # Step 1: linear starting guess (re-use existing function)
    images = linear_interpolate(start, end, n_images)

    # Step 2: endpoint pairwise distance matrices
    _, d_start = _pairwise(start.positions)
    _, d_end   = _pairwise(end.positions)

    # Step 3: optimise each intermediate image independently
    for k in range(1, n_total - 1):
        alpha    = k / (n_total - 1)
        d_target = d_start + alpha * (d_end - d_start)

        # Weights: 1/d_target^4; zero on diagonal and for zero-distance pairs
        with np.errstate(divide='ignore', invalid='ignore'):
            weights = np.where(d_target > 1e-10,
                               1.0 / d_target ** 4,
                               0.0)
        np.fill_diagonal(weights, 0.0)

        x0 = images[k].positions.ravel().copy()

        result = minimize(
            _idpp_obj_and_grad,
            x0,
            args=(n_atoms, d_target, weights),
            method='L-BFGS-B',
            jac=True,
            options={'maxiter': max_iter,
                     'ftol': tol,
                     'gtol': tol},
        )

        images[k].set_positions(result.x.reshape(n_atoms, 3))

    return images