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
    Supports PBC via the minimum image convention.

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
    Pairwise difference vectors and distances (non-periodic).

    Returns
    -------
    diff : (N, N, 3)  diff[i, j] = positions[i] - positions[j]
    dist : (N, N)     Euclidean distances; diagonal is 0.0
    """
    diff = positions[:, None, :] - positions[None, :, :]   # (N, N, 3)
    dist = np.sqrt((diff ** 2).sum(axis=-1))               # (N, N)
    return diff, dist


def _pairwise_mic(positions, cell, pbc):
    """
    MIC-aware pairwise difference vectors and distances for periodic systems.

    Applies find_mic to all (N*N) displacement vectors in a single batch
    call — no Python loops over pairs.

    Returns
    -------
    diff : (N, N, 3)  MIC displacement diff[i, j] = mic(r_i - r_j)
    dist : (N, N)     MIC distances; diagonal is 0.0
    """
    n = len(positions)
    diff_raw = positions[:, None, :] - positions[None, :, :]   # (N, N, 3)

    # Batch MIC: reshape to (N², 3), apply find_mic, reshape back
    diff_mic_flat, _ = find_mic(diff_raw.reshape(-1, 3), cell, pbc)
    diff = diff_mic_flat.reshape(n, n, 3)
    dist = np.sqrt((diff ** 2).sum(axis=-1))
    return diff, dist


def _idpp_obj_and_grad(flat_pos, n_atoms, d_target, weights, cell=None, pbc=None):
    """
    IDPP objective and analytical gradient for one image.

    Objective (Smidstrup et al., eq. 2):
        S = sum_{i<j} w_ij * (d_ij - d_target_ij)^2
          = 0.5 * sum_{i != j} w_ij * (d_ij - d_target_ij)^2

    Gradient (analytical, accounting for both (i,j) and (j,i) contributions):
        dS/dr_k = 2 * sum_{j != k} w_kj * delta_kj / d_kj * mic(r_k - r_j)

    For periodic systems, raw displacements are replaced by MIC displacements.
    The MIC image of a pair does not change under the small perturbations used
    by L-BFGS-B, so the gradient formula is identical to the non-periodic case.

    Parameters
    ----------
    flat_pos : (3N,) current atomic positions
    n_atoms  : int
    d_target : (N, N) target distance matrix (fixed throughout)
    weights  : (N, N) = 1/d_target^4, zeros on diagonal
    cell     : ase.cell.Cell or (3, 3) array, optional
    pbc      : (3,) bool array, optional

    Returns
    -------
    obj  : float
    grad : (3N,)
    """
    positions = flat_pos.reshape(n_atoms, 3)

    use_mic = (cell is not None) and (pbc is not None) and np.any(pbc)
    if use_mic:
        diff, dist = _pairwise_mic(positions, cell, pbc)
    else:
        diff, dist = _pairwise(positions)

    dist_safe = np.where(dist > 1e-12, dist, 1.0)
    delta = dist - d_target

    obj    = 0.5 * np.sum(weights * delta ** 2)
    factor = 2.0 * weights * delta / dist_safe
    np.fill_diagonal(factor, 0.0)
    grad   = np.einsum('ij,ijk->ik', factor, diff)

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

    Supports periodic boundary conditions via the minimum image convention.

    Algorithm
    ---------
    1. Generate linearly interpolated starting positions for each image
       (MIC-aware for periodic systems).
    2. For image k, define target pairwise distances by linear interpolation
       between the endpoint MIC distance matrices:
           d_target[k] = d_start + (k / (N-1)) * (d_end - d_start)
    3. Optimise each image independently (L-BFGS-B) to minimise the IDPP
       objective using MIC distances where PBC is active.
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

    n_total  = n_images + 2
    n_atoms  = len(start)
    use_mic  = bool(np.any(start.pbc))
    cell     = start.cell if use_mic else None
    pbc      = start.pbc  if use_mic else None

    # Step 1: linear starting guess (MIC-aware via linear_interpolate)
    images = linear_interpolate(start, end, n_images)

    # Step 2: endpoint pairwise distance matrices
    if use_mic:
        _, d_start = _pairwise_mic(start.positions, cell, pbc)
        _, d_end   = _pairwise_mic(end.positions,   cell, pbc)
    else:
        _, d_start = _pairwise(start.positions)
        _, d_end   = _pairwise(end.positions)

    # Step 3: optimise each intermediate image independently
    for k in range(1, n_total - 1):
        alpha    = k / (n_total - 1)
        d_target = d_start + alpha * (d_end - d_start)

        with np.errstate(divide='ignore', invalid='ignore'):
            weights = np.where(d_target > 1e-10,
                               1.0 / d_target ** 4,
                               0.0)
        np.fill_diagonal(weights, 0.0)

        x0 = images[k].positions.ravel().copy()

        result = minimize(
            _idpp_obj_and_grad,
            x0,
            args=(n_atoms, d_target, weights, cell, pbc),
            method='L-BFGS-B',
            jac=True,
            options={'maxiter': max_iter,
                     'ftol': tol,
                     'gtol': tol},
        )

        images[k].set_positions(result.x.reshape(n_atoms, 3))

    return images