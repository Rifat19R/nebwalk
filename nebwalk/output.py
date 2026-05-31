"""
Output utilities: energy profile plot, CSV export, ASE trajectory writer.
"""

import numpy as np
import matplotlib.pyplot as plt
from ase.io import write


def plot_energy_profile(images, filename="neb_profile.png", show=False, title=None):
    """
    Plot the NEB energy profile and save to file.

    Parameters
    ----------
    images : list of ase.Atoms
        All NEB images (including endpoints).
    filename : str
        Output file path (.png, .pdf, .svg are all supported by matplotlib).
    show : bool
        If True, call ``plt.show()`` after saving.
    title : str or None
        Plot title.  Defaults to "NEB Energy Profile".
    """
    energies = np.array([img.get_potential_energy() for img in images])
    energies_rel = energies - energies[0]

    n = len(energies)
    x = np.linspace(0, 1, n)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(x, energies_rel, "o-", color="royalblue", lw=2, ms=7, zorder=3)
    ax.axhline(0, color="gray", lw=0.5, ls="--")
    ax.set_xlabel("Reaction coordinate (normalised)")
    ax.set_ylabel("Energy (eV, relative to image 0)")
    ax.set_title(title or "NEB Energy Profile")
    ax.grid(True, alpha=0.25)

    # Annotate forward and reverse barriers
    i_max = int(np.argmax(energies_rel))
    Ea_fwd = float(energies_rel[i_max])
    Ea_rev = float(Ea_fwd - energies_rel[-1])

    ax.annotate(
        f"$E_a$ = {Ea_fwd:.3f} eV",
        xy=(x[i_max], Ea_fwd),
        xytext=(0.60, 0.75),
        textcoords="axes fraction",
        arrowprops=dict(arrowstyle="->", color="firebrick"),
        color="firebrick",
        fontsize=10,
    )

    plt.tight_layout()
    plt.savefig(filename, dpi=150)
    print(f"Profile saved → {filename}")
    if show:
        plt.show()
    plt.close(fig)


def save_csv(images, filename="neb_profile.csv"):
    """
    Write image index, absolute energy, and relative energy to CSV.

    Parameters
    ----------
    images : list of ase.Atoms
    filename : str
    """
    energies = np.array([img.get_potential_energy() for img in images])
    energies_rel = energies - energies[0]

    with open(filename, "w") as fh:
        fh.write("image,energy_eV,relative_energy_eV\n")
        for i, (e, er) in enumerate(zip(energies, energies_rel)):
            fh.write(f"{i},{e:.8f},{er:.8f}\n")

    print(f"Energies saved → {filename}")


def save_trajectory(images, filename="neb.traj"):
    """
    Write all images to an ASE trajectory file.

    Can be viewed with ``ase gui neb.traj`` or read back with
    ``ase.io.read('neb.traj', index=':')``.

    Parameters
    ----------
    images : list of ase.Atoms
    filename : str
    """
    write(filename, images)
    print(f"Trajectory saved → {filename}")
