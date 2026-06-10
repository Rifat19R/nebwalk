"""Quantum ESPRESSO calculator helpers for nebwalk.

The functions here intentionally keep Quantum ESPRESSO optional. Importing
``nebwalk`` does not import ASE's Espresso calculator or require ``pw.x``;
callers opt in by creating a QE calculator factory.
"""

from __future__ import annotations

import shlex
import shutil
from collections.abc import Callable, Mapping
from dataclasses import dataclass, field
from os import X_OK, access
from pathlib import Path
from typing import Any


def _espresso_placeholder(*args: Any, **kwargs: Any) -> Any:
    """Placeholder replaced lazily by ASE's Espresso class."""
    raise RuntimeError("ASE Espresso calculator was not loaded.")


Espresso = _espresso_placeholder
EspressoProfile = None  # set lazily alongside Espresso (ASE >= 3.23)


@dataclass(frozen=True)
class QEParams:
    """Common Quantum ESPRESSO parameters for NEB image calculations."""

    ecutwfc: float = 40.0
    ecutrho: float = 320.0
    kpts: tuple[int, int, int] = (1, 1, 1)
    koffset: tuple[int, int, int] = (0, 0, 0)
    occupations: str = "smearing"
    smearing: str = "marzari-vanderbilt"
    degauss: float = 0.02
    conv_thr: float = 1.0e-8
    mixing_beta: float = 0.3
    nspin: int = 1
    starting_magnetization: Mapping[str | int, float] | None = None
    extra_control: Mapping[str, Any] = field(default_factory=dict)
    extra_system: Mapping[str, Any] = field(default_factory=dict)
    extra_electrons: Mapping[str, Any] = field(default_factory=dict)


def _command_binaries(command: str) -> list[str]:
    """Return executable names worth checking from a QE command string."""
    tokens = shlex.split(command)
    if not tokens:
        return ["pw.x"]

    binaries = [tokens[0]]
    flags_with_values = {
        "-np",
        "-n",
        "--np",
        "--n",
        "-machinefile",
        "-hostfile",
        "--host",
        "--hostfile",
    }
    positionals: list[str] = []
    skip_next = False
    for token in tokens[1:]:
        if skip_next:
            skip_next = False
            continue
        if token in flags_with_values:
            skip_next = True
            continue
        if token.startswith("-") or token.isdigit() or "=" in token:
            continue
        positionals.append(token)

    for token in reversed(positionals):
        name = Path(token).name
        if token.endswith(".x") or "." in name or "/" in token or "\\" in token:
            binaries.append(token)
            break

    unique: list[str] = []
    for binary in binaries:
        if binary not in unique:
            unique.append(binary)
    return unique


def _binary_is_available(binary: str) -> bool:
    """Check PATH binaries and explicit executable paths."""
    path = Path(binary).expanduser()
    has_path_separator = "/" in binary or "\\" in binary
    if path.is_absolute() or has_path_separator:
        return path.is_file() and access(path, X_OK)
    return shutil.which(binary) is not None


def validate_qe_setup(
    pseudo_dir: str | Path,
    pseudopotentials: Mapping[str, str],
    command: str = "pw.x",
) -> None:
    """Validate QE executable access and required pseudopotential files.

    Raises
    ------
    FileNotFoundError
        If the pseudopotential directory, a requested UPF file, or a required
        executable is unavailable.
    ValueError
        If the pseudopotential mapping is empty.
    """
    pseudo_path = Path(pseudo_dir).expanduser()
    if not pseudo_path.is_dir():
        raise FileNotFoundError(
            f"QE pseudopotential directory not found: {pseudo_path}"
        )
    if not pseudopotentials:
        raise ValueError("At least one pseudopotential must be provided.")

    missing_pseudos = [
        pseudo_name
        for pseudo_name in pseudopotentials.values()
        if not (pseudo_path / pseudo_name).is_file()
    ]
    if missing_pseudos:
        missing = ", ".join(missing_pseudos)
        raise FileNotFoundError(f"Missing QE pseudopotential file(s): {missing}")

    missing_bins = []
    for binary in _command_binaries(command):
        if not _binary_is_available(binary):
            missing_bins.append(binary)
    if missing_bins:
        missing = ", ".join(missing_bins)
        raise FileNotFoundError(f"QE executable not found or not executable: {missing}")


def _species_index(
    species: str | int,
    pseudopotentials: Mapping[str, str],
) -> int:
    """Map chemical symbol or explicit index to QE species index."""
    if isinstance(species, int):
        return species
    symbols = list(pseudopotentials)
    if species not in pseudopotentials:
        known = ", ".join(symbols)
        raise ValueError(f"Magnetization species {species!r} not in pseudos: {known}")
    return symbols.index(species) + 1


def _input_data(
    params: QEParams,
    pseudo_dir: Path,
    outdir: Path,
    pseudopotentials: Mapping[str, str],
) -> dict[str, dict[str, Any]]:
    """Build ASE Espresso input_data with NEB-safe defaults."""
    control: dict[str, Any] = {
        "calculation": "scf",
        "restart_mode": "from_scratch",
        "tprnfor": True,
        "tstress": False,
        "disk_io": "low",
        "pseudo_dir": str(pseudo_dir),
        "outdir": str(outdir),
    }
    control.update(dict(params.extra_control))

    system: dict[str, Any] = {
        "ecutwfc": params.ecutwfc,
        "ecutrho": params.ecutrho,
        "occupations": params.occupations,
        "nspin": params.nspin,
    }
    if params.occupations == "smearing":
        system["smearing"] = params.smearing
        system["degauss"] = params.degauss
    if params.nspin == 2 and params.starting_magnetization:
        for species, value in params.starting_magnetization.items():
            idx = _species_index(species, pseudopotentials)
            system[f"starting_magnetization({idx})"] = value
    system.update(dict(params.extra_system))

    electrons: dict[str, Any] = {
        "conv_thr": params.conv_thr,
        "mixing_beta": params.mixing_beta,
    }
    electrons.update(dict(params.extra_electrons))

    return {"control": control, "system": system, "electrons": electrons}


def make_qe_factory(
    params: QEParams,
    pseudo_dir: str | Path,
    pseudopotentials: Mapping[str, str],
    base_dir: str | Path = "neb_qe_workdir",
    command: str = "pw.x",
) -> Callable[[], Any]:
    """Create a zero-argument factory returning independent QE calculators.

    Each factory call creates a fresh image directory, which avoids ASE/QE file
    collisions during serial and parallel NEB force evaluations.
    """
    global Espresso, EspressoProfile

    if not pseudopotentials:
        raise ValueError("At least one pseudopotential must be provided.")
    if Espresso is _espresso_placeholder:
        from ase.calculators.espresso import Espresso as ASEEspresso

        Espresso = ASEEspresso

    pseudo_path = Path(pseudo_dir).expanduser().resolve()
    base_path = Path(base_dir).expanduser().resolve()
    counter = {"image": 0}

    def factory() -> Any:
        image_idx = counter["image"]
        counter["image"] += 1

        image_dir = base_path / f"image_{image_idx:03d}"
        outdir = image_dir / "tmp"
        outdir.mkdir(parents=True, exist_ok=True)

        if EspressoProfile is not None:
            _profile = EspressoProfile(shlex.split(command))
            return Espresso(
                profile=_profile,
                input_data=_input_data(params, pseudo_path, outdir, pseudopotentials),
                pseudopotentials=dict(pseudopotentials),
                kpts=params.kpts,
                koffset=params.koffset,
                directory=str(image_dir),
            )
        return Espresso(
            input_data=_input_data(params, pseudo_path, outdir, pseudopotentials),
            pseudopotentials=dict(pseudopotentials),
            kpts=params.kpts,
            koffset=params.koffset,
            directory=str(image_dir),
            label="espresso",
            command=f"{command} -in PREFIX.pwi > PREFIX.pwo",
        )

    return factory


__all__ = ["QEParams", "make_qe_factory", "validate_qe_setup"]
