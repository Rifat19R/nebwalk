"""
FIRE optimizer for NEB.

NEB forces are not the gradient of any scalar function — the spring and
perpendicular projections break the conservative field assumption required
by L-BFGS-B. FIRE is robust to this.

Parallel image evaluation uses threads (not processes). PyTorch releases the
GIL during C++ computation, so threads give real parallelism for MACE and
other compiled calculators. Process-based parallelism is avoided because
pickling ASE calculator objects is unreliable across calculator types.

Reference: Bitzek, Koskinen, Gähler, Moseler, Gumbsch,
    PRL 97, 170201 (2006). DOI: 10.1103/PhysRevLett.97.170201
"""

import numpy as np
from concurrent.futures import ThreadPoolExecutor
from ase.constraints import FixAtoms
from .forces import compute_neb_forces, variable_spring_constants


def _find_climb_index(energies):
    """Return index of the highest-energy movable image (1-indexed in full path)."""
    return int(np.argmax(energies[1:-1])) + 1


def _get_free_mask(img):
    """Return boolean mask (n_atoms,): True = atom is free to move."""
    mask = np.ones(len(img), dtype=bool)
    for c in img.constraints:
        if isinstance(c, FixAtoms):
            mask[c.index] = False
    return mask


def _eval_image(img):
    """
    Evaluate energy and forces for a single image.

    This function is the unit of work dispatched to each thread.
    It must be defined at module level so ThreadPoolExecutor can reference it.

    Parameters
    ----------
    img : ase.Atoms with calculator attached

    Returns
    -------
    (energy, forces) : (float, ndarray of shape (N_atoms, 3))
    """
    E = img.get_potential_energy()
    F = img.get_forces()
    return E, F


def _eval_all(images, n_workers):
    """
    Evaluate energy and forces for all images.

    Endpoints (first and last) are always evaluated sequentially — they do
    not move and could be cached, but the cost is negligible compared to
    intermediate images.

    Parameters
    ----------
    images : list of ase.Atoms
    n_workers : int
        Number of threads. 1 = sequential (no thread overhead).
        n_workers > 1 = parallel evaluation of intermediate images.

    Returns
    -------
    energies : list of float, length N
    forces   : list of ndarray or None, length N
        None at endpoints (forces not needed there).
    """
    n = len(images)

    # Endpoint energies (sequential — endpoints never move)
    E_start = images[0].get_potential_energy()
    E_end   = images[-1].get_potential_energy()

    movable = images[1:-1]

    if n_workers == 1 or len(movable) <= 1:
        # Sequential path — no thread overhead
        results = [_eval_image(img) for img in movable]
    else:
        # Parallel path — threads release GIL for MACE/PyTorch
        with ThreadPoolExecutor(max_workers=min(n_workers, len(movable))) as ex:
            results = list(ex.map(_eval_image, movable))

    energies_mid = [r[0] for r in results]
    forces_mid   = [r[1] for r in results]

    energies = [E_start] + energies_mid + [E_end]
    forces   = [None]    + forces_mid   + [None]
    return energies, forces


def fire_optimize(
    images,
    k,
    fmax=0.05,
    max_steps=500,
    climb=False,
    climb_delay=100,
    k_min=None,
    n_workers=1,
    verbose=True,
    dt_start=0.01,
    dt_max=0.10,
    alpha_start=0.10,
    N_min=5,
    f_inc=1.10,
    f_dec=0.50,
    f_alpha=0.99,
):
    """
    Run FIRE optimisation on the NEB path.

    Images are modified in-place via ``set_positions``.
    FixAtoms constraints on intermediate images are respected.

    Parameters
    ----------
    images : list of ase.Atoms
        Full path including fixed endpoints.
    k : float
        Spring constant (uniform) or maximum spring constant (variable), eV/Å².
    fmax : float
        Convergence criterion: max |F| on any free atom < fmax eV/Å.
    max_steps : int
        Hard cap on iterations.
    climb : bool
        Enable CI-NEB after ``climb_delay`` steps.
    climb_delay : int
        Minimum steps before activating CI.
    k_min : float or None
        If given, activates variable spring constants. k is then k_max.
    n_workers : int
        Number of threads for parallel image evaluation (default 1 = sequential).
        Set to the number of intermediate images for maximum speedup.
        Effective for calculators that release the GIL (MACE, MACE-MP-0).
        No benefit for pure-Python calculators (EMT) due to the GIL.
    verbose : bool
        Print convergence info every 10 steps.
    dt_start, dt_max : float
        Initial and maximum time step.
    alpha_start, f_inc, f_dec, f_alpha, N_min : float/int
        Standard FIRE hyper-parameters.

    Returns
    -------
    converged : bool
    steps_taken : int
    history : list of dict
        Keys: ``step``, ``fmax``, ``energies``, ``k_springs``.
    """
    n_movable = len(images) - 2
    n_atoms   = len(images[0])

    free_masks = [_get_free_mask(img) for img in images[1:-1]]
    v = np.zeros((n_movable, n_atoms, 3))

    dt      = dt_start
    alpha   = alpha_start
    n_pos   = 0
    history = []
    climb_active = False
    climb_index  = None
    use_variable_k = k_min is not None

    if verbose and n_workers > 1:
        print(f"  Parallel evaluation: {n_workers} threads "
              f"({len(images) - 2} intermediate images)")

    for step in range(max_steps):

        # Single evaluation pass — energy and forces for all images
        energies, forces_cache = _eval_all(images, n_workers)

        # Variable spring constants (uses energies already computed)
        if use_variable_k:
            k_curr = variable_spring_constants(energies, k_max=k, k_min=k_min)
        else:
            k_curr = k

        # CI-NEB activation
        if climb and not climb_active and step >= climb_delay:
            climb_active = True
            climb_index  = _find_climb_index(energies)
            if verbose:
                print(f"  CI-NEB active at step {step}, image {climb_index} "
                      f"(E = {energies[climb_index]:.4f} eV)")

        # NEB forces — forces_cache avoids all redundant calculator calls
        neb_forces = compute_neb_forces(
            images, k_curr,
            climb=climb_active,
            climb_index=climb_index,
            energies=energies,
            forces=forces_cache,
        )
        f = np.array(neb_forces[1:-1])

        for j in range(n_movable):
            f[j][~free_masks[j]] = 0.0
            v[j][~free_masks[j]] = 0.0

        fmax_curr = float(np.sqrt((f**2).sum(axis=-1)).max())
        k_springs = k_curr if use_variable_k else np.full(len(images) - 1, k)
        history.append({
            'step':     step,
            'fmax':     fmax_curr,
            'energies': energies,
            'k_springs': (k_springs.tolist()
                          if hasattr(k_springs, 'tolist') else k_springs),
        })

        if verbose and step % 10 == 0:
            barrier = max(energies) - energies[0]
            print(f"  Step {step:4d}  fmax = {fmax_curr:.4f} eV/Å  "
                  f"barrier = {barrier:.4f} eV")

        if fmax_curr < fmax:
            if verbose:
                print(f"  Converged in {step} steps  "
                      f"fmax = {fmax_curr:.5f} eV/Å")
            return True, step, history

        P = float((f * v).sum())

        f_norm = float(np.linalg.norm(f))
        v_norm = float(np.linalg.norm(v))
        if f_norm > 1e-10:
            v = (1.0 - alpha) * v + alpha * (v_norm / f_norm) * f

        if P > 0.0:
            n_pos += 1
            if n_pos >= N_min:
                dt    = min(dt * f_inc, dt_max)
                alpha *= f_alpha
        else:
            n_pos  = 0
            v[:]   = 0.0
            dt    *= f_dec
            alpha  = alpha_start

        v += f * dt
        for j in range(n_movable):
            v[j][~free_masks[j]] = 0.0

        for j, img in enumerate(images[1:-1]):
            img.set_positions(img.positions + v[j] * dt)

    if verbose:
        print(f"  Warning: not converged after {max_steps} steps. "
              f"fmax = {fmax_curr:.4f} eV/Å")
    return False, max_steps, history
