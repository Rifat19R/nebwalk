"""FIRE optimizer for NEB."""

from __future__ import annotations

import logging
from collections.abc import Sequence
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np
from ase import Atoms
from ase.constraints import FixAtoms
from numpy.typing import NDArray

from .forces import compute_neb_forces, validate_calculators, variable_spring_constants

logger = logging.getLogger(__name__)
FloatArray = NDArray[np.float64]


def max_force_magnitude(forces: FloatArray) -> float:
    """Return ASE-style maximum per-atom force magnitude."""
    return float(np.sqrt((forces**2).sum(axis=-1)).max())


def check_convergence(forces: FloatArray, fmax: float) -> bool:
    """Return True if max per-atom force magnitude is below fmax."""
    return max_force_magnitude(forces) < fmax


class FIREOptimizer:
    """FIRE optimizer for non-conservative force fields such as NEB."""

    def __init__(
        self,
        dt: float = 0.01,
        dt_max: float = 0.10,
        n_min: int = 5,
        f_inc: float = 1.10,
        f_dec: float = 0.50,
        alpha_start: float = 0.10,
        f_alpha: float = 0.99,
    ) -> None:
        self.dt_start = float(dt)
        self.dt = float(dt)
        self.dt_max = float(dt_max)
        self.n_min = int(n_min)
        self.f_inc = float(f_inc)
        self.f_dec = float(f_dec)
        self.alpha_start = float(alpha_start)
        self.alpha = float(alpha_start)
        self.f_alpha = float(f_alpha)
        self.n_pos = 0
        self.velocity: FloatArray | None = None

    def reset(self) -> None:
        """Reset velocity and FIRE counters."""
        self.dt = self.dt_start
        self.alpha = self.alpha_start
        self.n_pos = 0
        if self.velocity is not None:
            self.velocity[:] = 0.0

    def step(self, positions: FloatArray, forces: FloatArray) -> FloatArray:
        """Perform one FIRE step and return updated positions."""
        if self.velocity is None or self.velocity.shape != positions.shape:
            self.velocity = np.zeros_like(positions, dtype=float)

        velocity = self.velocity
        power = float((forces * velocity).sum())

        force_norm = float(np.linalg.norm(forces))
        velocity_norm = float(np.linalg.norm(velocity))
        if force_norm > 1e-10:
            velocity[:] = (
                (1.0 - self.alpha) * velocity
                + self.alpha * (velocity_norm / force_norm) * forces
            )

        if power > 0.0:
            self.n_pos += 1
            if self.n_pos >= self.n_min:
                self.dt = min(self.dt * self.f_inc, self.dt_max)
                self.alpha *= self.f_alpha
        else:
            self.n_pos = 0
            velocity[:] = 0.0
            self.dt *= self.f_dec
            self.alpha = self.alpha_start

        velocity += forces * self.dt
        return positions + velocity * self.dt


def _find_climb_index(energies: Sequence[float]) -> int:
    """Return index of highest-energy movable image in full path."""
    return int(np.argmax(energies[1:-1])) + 1


def _get_free_mask(img: Atoms) -> NDArray[np.bool_]:
    """Return boolean mask (n_atoms,): True = atom is free to move."""
    mask = np.ones(len(img), dtype=bool)
    for constraint in img.constraints:
        if isinstance(constraint, FixAtoms):
            mask[constraint.index] = False
    return mask


def _eval_image(img: Atoms) -> tuple[float, FloatArray]:
    """Evaluate energy and forces for one image."""
    energy = float(img.get_potential_energy())
    forces = np.asarray(img.get_forces(), dtype=float)
    return energy, forces


def _warn_if_gpu_calculator(images: Sequence[Atoms], n_workers: int) -> None:
    """Warn when thread-parallel evaluation sees a CUDA calculator."""
    if n_workers <= 1:
        return
    for img in images:
        device = getattr(img.calc, "device", None)
        if device is not None and "cuda" in str(device).lower():
            logger.warning(
                "n_workers=%d with a CUDA calculator detected. "
                "Thread-parallel evaluation is NOT safe with GPU calculators "
                "and may silently corrupt results. Use n_workers=1 for GPU runs.",
                n_workers,
            )
            return


def _eval_all(
    images: Sequence[Atoms],
    n_workers: int,
) -> tuple[list[float], list[FloatArray | None]]:
    """Evaluate energy and forces for all images."""
    e_start = float(images[0].get_potential_energy())
    e_end = float(images[-1].get_potential_energy())
    movable = images[1:-1]

    if n_workers == 1 or len(movable) <= 1:
        results = [_eval_image(img) for img in movable]
    else:
        workers = min(int(n_workers), len(movable))
        with ThreadPoolExecutor(max_workers=workers) as executor:
            results = list(executor.map(_eval_image, movable))

    energies_mid = [result[0] for result in results]
    forces_mid = [result[1] for result in results]
    return [e_start] + energies_mid + [e_end], [None] + forces_mid + [None]


def _apply_free_masks(
    values: FloatArray,
    free_masks: Sequence[NDArray[np.bool_]],
) -> None:
    """Zero constrained atoms in-place."""
    for j, mask in enumerate(free_masks):
        values[j][~mask] = 0.0


def fire_optimize(
    images: Sequence[Atoms],
    k: float,
    fmax: float = 0.05,
    max_steps: int = 500,
    climb: bool = False,
    climb_delay: int = 100,
    k_min: float | None = None,
    n_workers: int = 1,
    verbose: bool = True,
    dt_start: float = 0.01,
    dt_max: float = 0.10,
    alpha_start: float = 0.10,
    N_min: int = 5,
    f_inc: float = 1.10,
    f_dec: float = 0.50,
    f_alpha: float = 0.99,
) -> tuple[bool, int, list[dict[str, Any]]]:
    """Run FIRE optimization on the NEB path.

    Thread-based parallelism is intended for CPU calculators. Use serial
    evaluation (n_workers=1) for CUDA-backed calculators.
    """
    validate_calculators(images)

    n_movable = len(images) - 2
    n_atoms = len(images[0])
    free_masks = [_get_free_mask(img) for img in images[1:-1]]
    optimizer = FIREOptimizer(
        dt=dt_start,
        dt_max=dt_max,
        n_min=N_min,
        f_inc=f_inc,
        f_dec=f_dec,
        alpha_start=alpha_start,
        f_alpha=f_alpha,
    )
    optimizer.velocity = np.zeros((n_movable, n_atoms, 3), dtype=float)

    history: list[dict[str, Any]] = []
    climb_active = False
    climb_index: int | None = None
    use_variable_k = k_min is not None
    fmax_curr = float("inf")

    _warn_if_gpu_calculator(images, n_workers)

    if verbose and n_workers > 1:
        logger.info(
            "Parallel evaluation: %d threads (%d intermediate images)",
            n_workers,
            len(images) - 2,
        )

    for step in range(max_steps):
        energies, forces_cache = _eval_all(images, n_workers)

        if use_variable_k:
            k_curr = variable_spring_constants(energies, k_max=k, k_min=k_min)
        else:
            k_curr = k

        if climb and not climb_active and step >= climb_delay:
            climb_active = True
            climb_index = _find_climb_index(energies)
            if verbose:
                logger.info(
                    "CI-NEB active at step %d, image %d (E = %.4f eV)",
                    step,
                    climb_index,
                    energies[climb_index],
                )

        neb_forces = compute_neb_forces(
            images,
            k_curr,
            climb=climb_active,
            climb_index=climb_index,
            energies=energies,
            forces=forces_cache,
        )
        forces = np.array(neb_forces[1:-1], dtype=float)
        _apply_free_masks(forces, free_masks)
        if optimizer.velocity is not None:
            _apply_free_masks(optimizer.velocity, free_masks)

        fmax_curr = max_force_magnitude(forces)
        k_springs = k_curr if use_variable_k else np.full(len(images) - 1, k)
        history.append(
            {
                "step": step,
                "fmax": fmax_curr,
                "energies": energies,
                "k_springs": (
                    k_springs.tolist()
                    if hasattr(k_springs, "tolist")
                    else k_springs
                ),
            }
        )

        if verbose and step % 10 == 0:
            barrier = max(energies) - energies[0]
            logger.info(
                "Step %4d  fmax = %.4f eV/Ang  barrier = %.4f eV",
                step,
                fmax_curr,
                barrier,
            )

        if check_convergence(forces, fmax):
            if verbose:
                logger.info(
                    "Converged in %d steps  fmax = %.5f eV/Ang",
                    step,
                    fmax_curr,
                )
            return True, step, history

        positions = np.array([img.positions for img in images[1:-1]], dtype=float)
        new_positions = optimizer.step(positions, forces)
        if optimizer.velocity is not None:
            _apply_free_masks(optimizer.velocity, free_masks)

        for j, img in enumerate(images[1:-1]):
            img.set_positions(new_positions[j])

    if verbose:
        logger.warning(
            "Not converged after %d steps. fmax = %.4f eV/Ang",
            max_steps,
            fmax_curr,
        )
    return False, max_steps, history
