"""
FIRE optimizer for NEB.

Why FIRE and not L-BFGS-B via scipy:
    NEB forces are NOT the gradient of any single scalar function (the spring
    and perpendicular projections break the conservative field assumption).
    L-BFGS-B line-search can therefore stall or diverge.  FIRE treats the NEB
    forces as a velocity-dependent damping term and is robust to this.

Reference: Bitzek, Koskinen, Gähler, Moseler, Gumbsch,
    PRL 97, 170201 (2006). DOI: 10.1103/PhysRevLett.97.170201
"""

import numpy as np
from ase.constraints import FixAtoms
from .forces import compute_neb_forces


def _find_climb_index(images):
    """Return index of the highest-energy movable image."""
    energies = [img.get_potential_energy() for img in images[1:-1]]
    return int(np.argmax(energies)) + 1  # offset for endpoint


def _get_free_mask(img):
    """
    Return boolean mask of shape (n_atoms,): True = atom is free to move.
    Reads FixAtoms constraints from the image.
    """
    mask = np.ones(len(img), dtype=bool)
    for c in img.constraints:
        if isinstance(c, FixAtoms):
            mask[c.index] = False
    return mask


def fire_optimize(
    images,
    k,
    fmax=0.05,
    max_steps=500,
    climb=False,
    climb_delay=100,
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
        Full path including fixed endpoints.  All images must have a
        calculator attached before calling this function.
    k : float
        Spring constant, eV/Å².
    fmax : float
        Convergence criterion: max |F| component on any FREE atom in any
        movable image < fmax.
    max_steps : int
        Hard cap on iterations.
    climb : bool
        Enable climbing-image NEB (CI-NEB) after ``climb_delay`` steps.
    climb_delay : int
        Minimum steps before activating CI.
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
        Keys: ``step``, ``fmax``, ``energies``.
    """
    n_movable = len(images) - 2
    n_atoms = len(images[0])

    # Per-image boolean masks: True = atom is free
    free_masks = [_get_free_mask(img) for img in images[1:-1]]  # (n_movable, n_atoms)

    # Velocity array: shape (n_movable, n_atoms, 3); zero on fixed atoms always
    v = np.zeros((n_movable, n_atoms, 3))

    dt = dt_start
    alpha = alpha_start
    n_pos = 0          # consecutive steps with P > 0
    history = []
    climb_active = False
    climb_index = None

    for step in range(max_steps):

        # --- Activate climbing image ---
        if climb and not climb_active and step >= climb_delay:
            climb_active = True
            climb_index = _find_climb_index(images)
            if verbose:
                print(f"  CI-NEB active at step {step}, image {climb_index} "
                      f"(E = {images[climb_index].get_potential_energy():.4f} eV)")

        # --- Compute NEB forces ---
        neb_forces = compute_neb_forces(
            images, k,
            climb=climb_active,
            climb_index=climb_index,
        )
        # (n_movable, n_atoms, 3)
        f = np.array(neb_forces[1:-1])

        # --- Apply FixAtoms constraints: zero forces (and velocities) on fixed atoms ---
        for j in range(n_movable):
            f[j][~free_masks[j]] = 0.0
            v[j][~free_masks[j]] = 0.0

        # --- Convergence check (only on free atoms) ---
        fmax_curr = float(np.max(np.abs(f)))  # fixed atoms already zeroed
        energies = [img.get_potential_energy() for img in images]
        history.append({'step': step, 'fmax': fmax_curr, 'energies': energies})

        if verbose and step % 10 == 0:
            barrier = max(energies) - energies[0]
            print(f"  Step {step:4d}  fmax = {fmax_curr:.4f} eV/Å  "
                  f"barrier = {barrier:.4f} eV")

        if fmax_curr < fmax:
            if verbose:
                print(f"  Converged in {step} steps  "
                      f"fmax = {fmax_curr:.5f} eV/Å")
            return True, step, history

        # --- FIRE step ---
        P = float((f * v).sum())

        f_norm = float(np.linalg.norm(f))
        v_norm = float(np.linalg.norm(v))
        if f_norm > 1e-10:
            v = (1.0 - alpha) * v + alpha * (v_norm / f_norm) * f

        if P > 0.0:
            n_pos += 1
            if n_pos >= N_min:
                dt = min(dt * f_inc, dt_max)
                alpha *= f_alpha
        else:
            n_pos = 0
            v[:] = 0.0
            dt *= f_dec
            alpha = alpha_start

        # Velocity and position update (semi-implicit Euler)
        v += f * dt
        # Re-zero fixed atoms after acceleration
        for j in range(n_movable):
            v[j][~free_masks[j]] = 0.0

        for j, img in enumerate(images[1:-1]):
            img.set_positions(img.positions + v[j] * dt)

    if verbose:
        print(f"  Warning: not converged after {max_steps} steps. "
              f"fmax = {fmax_curr:.4f} eV/Å")
    return False, max_steps, history
