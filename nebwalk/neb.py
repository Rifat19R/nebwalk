"""
NEB class: user-facing API for nebwalk.
"""

from .forces import compute_neb_forces
from .optimize import fire_optimize
from .output import plot_energy_profile, save_csv, save_trajectory


class NEB:
    """
    Nudged Elastic Band calculator.

    Parameters
    ----------
    images : list of ase.Atoms
        Full path ``[start, img_1, ..., img_n, end]``.
        Every image must have an ASE-compatible calculator attached.
    k : float
        Spring constant in eV/Å² (default 0.1).
        With uniform springs: all springs use this value.
        With variable springs (k_min set): this is the maximum value.
    k_min : float or None
        If given, activates energy-weighted variable spring constants.
        Springs near the saddle point use ``k``; springs far from it use
        ``k_min``. Recommended: k_min = k / 3.
        If None (default), uniform springs are used.
    climb : bool
        If True, activate CI-NEB after ``climb_delay`` steps.
    climb_delay : int
        Number of FIRE steps before CI-NEB activates (default 100).
    n_workers : int
        Number of threads for parallel image evaluation (default 1).
        Set to the number of intermediate images for maximum speedup.
        Effective for calculators that release the GIL during computation
        (MACE, MACE-MP-0, Egret-1t). No benefit for pure-Python calculators
        (EMT) due to the Python GIL.

    Raises
    ------
    ValueError
        If fewer than 3 images are provided, any image is missing a
        calculator, k_min >= k, or n_workers < 1.
    """

    def __init__(self, images, k=0.1, k_min=None, climb=False,
                 climb_delay=100, n_workers=1):
        if len(images) < 3:
            raise ValueError(
                f"Need at least 3 images (2 endpoints + ≥1 intermediate), "
                f"got {len(images)}."
            )
        for i, img in enumerate(images):
            if img.calc is None:
                raise ValueError(
                    f"Image {i} has no calculator. "
                    f"Attach one before constructing NEB."
                )
        if k_min is not None and float(k_min) >= float(k):
            raise ValueError(
                f"k_min ({k_min}) must be strictly less than k ({k})."
            )
        if int(n_workers) < 1:
            raise ValueError(f"n_workers must be >= 1, got {n_workers}.")

        self.images      = images
        self.k           = float(k)
        self.k_min       = float(k_min) if k_min is not None else None
        self.climb       = climb
        self.climb_delay = int(climb_delay)
        self.n_workers   = int(n_workers)
        self.history     = []

    def optimize(self, fmax=0.05, max_steps=500, verbose=True):
        """
        Optimise the NEB path using the FIRE algorithm.

        Parameters
        ----------
        fmax : float
            Convergence criterion in eV/Å.
        max_steps : int
            Hard iteration limit.
        verbose : bool
            Print progress every 10 steps.

        Returns
        -------
        converged : bool
        """
        converged, steps, history = fire_optimize(
            self.images,
            self.k,
            fmax=fmax,
            max_steps=max_steps,
            climb=self.climb,
            climb_delay=self.climb_delay,
            k_min=self.k_min,
            n_workers=self.n_workers,
            verbose=verbose,
        )
        self.history = history
        return converged

    def get_energies(self):
        """Return list of potential energies (eV) for all images."""
        return [img.get_potential_energy() for img in self.images]

    def get_barrier(self):
        """Return forward activation barrier in eV (relative to image 0)."""
        e = self.get_energies()
        return max(e) - e[0]

    def get_spring_constants(self):
        """
        Return the current spring constants as an ndarray of shape (N-1,).

        With uniform springs, all values equal self.k.
        With variable springs, values reflect the current energy profile.
        """
        import numpy as np
        from .forces import variable_spring_constants
        energies = self.get_energies()
        if self.k_min is not None:
            return variable_spring_constants(
                energies, k_max=self.k, k_min=self.k_min
            )
        return np.full(len(self.images) - 1, self.k)

    def plot(self, filename="neb_profile.png", show=False, title=None):
        """Plot and save the energy profile."""
        plot_energy_profile(self.images, filename=filename,
                            show=show, title=title)

    def save_csv(self, filename="neb_profile.csv"):
        """Write energies to CSV."""
        save_csv(self.images, filename=filename)

    def save_trajectory(self, filename="neb.traj"):
        """Write all images to an ASE trajectory file."""
        save_trajectory(self.images, filename=filename)
