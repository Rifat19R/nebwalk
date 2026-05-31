"""
NEB class: user-facing API for nebwalk.
"""

from .forces import compute_neb_forces
from .optimize import fire_optimize
from .output import plot_energy_profile, save_csv, save_trajectory


class NEB:
    """
    Nudged Elastic Band calculator.

    Usage
    -----
    >>> from nebwalk import NEB, linear_interpolate
    >>> from ase.calculators.emt import EMT
    >>>
    >>> images = linear_interpolate(start, end, n_images=7)
    >>> for img in images:
    ...     img.calc = EMT()
    >>> neb = NEB(images, k=0.1, climb=True)
    >>> neb.optimize(fmax=0.05)
    >>> neb.plot("profile.png")

    Parameters
    ----------
    images : list of ase.Atoms
        Full path ``[start, img_1, ..., img_n, end]``.
        Every image must have an ASE-compatible calculator attached.
    k : float
        Spring constant in eV/Å² (default 0.1).
        For stiff bonds or wide barriers, values of 0.1–1.0 are typical.
        Do NOT scale by number of images; k has the same physical meaning
        regardless of path length.
    climb : bool
        If True, activate climbing-image NEB (CI-NEB) after ``climb_delay``
        steps.  CI-NEB drives the highest-energy image to the true saddle
        point.  Do not activate on the first optimisation if the path is
        far from converged.
    climb_delay : int
        Number of FIRE steps before CI-NEB activates (default 100).

    Raises
    ------
    ValueError
        If fewer than 3 images are provided, or if any image is missing a
        calculator.
    """

    def __init__(self, images, k=0.1, climb=False, climb_delay=100):
        if len(images) < 3:
            raise ValueError(
                f"Need at least 3 images (2 endpoints + ≥1 intermediate), "
                f"got {len(images)}."
            )
        for i, img in enumerate(images):
            if img.calc is None:
                raise ValueError(
                    f"Image {i} has no calculator.  Attach one before "
                    f"constructing NEB."
                )
        self.images = images
        self.k = float(k)
        self.climb = climb
        self.climb_delay = int(climb_delay)
        self.history = []

    # ------------------------------------------------------------------
    # Optimisation
    # ------------------------------------------------------------------

    def optimize(self, fmax=0.05, max_steps=500, verbose=True):
        """
        Optimise the NEB path using the FIRE algorithm.

        Parameters
        ----------
        fmax : float
            Convergence criterion: maximum force component on any atom in
            any movable image must be below ``fmax`` eV/Å.
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
            verbose=verbose,
        )
        self.history = history
        return converged

    # ------------------------------------------------------------------
    # Inspection
    # ------------------------------------------------------------------

    def get_energies(self):
        """Return list of potential energies (eV) for all images."""
        return [img.get_potential_energy() for img in self.images]

    def get_barrier(self):
        """
        Return forward activation barrier in eV (relative to image 0).

        Returns
        -------
        Ea : float
        """
        e = self.get_energies()
        return max(e) - e[0]

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def plot(self, filename="neb_profile.png", show=False, title=None):
        """Plot and save the energy profile."""
        plot_energy_profile(self.images, filename=filename, show=show, title=title)

    def save_csv(self, filename="neb_profile.csv"):
        """Write energies to CSV."""
        save_csv(self.images, filename=filename)

    def save_trajectory(self, filename="neb.traj"):
        """Write all images to an ASE trajectory file."""
        save_trajectory(self.images, filename=filename)
