"""MLIP-assisted NEB workflow utilities."""

from __future__ import annotations

import json
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any

from ase.io import write

from .engine import NEBRunConfig, run_neb_calculation
from .selection import select_images


@dataclass(frozen=True)
class MLIPActiveNEBConfig:
    """Configuration for the MLIP-assisted NEB workflow layer."""

    selection_strategy: str = "peak_plus_neighbors"
    n_select: int = 3
    include_endpoints: bool = False
    export_selected: bool = True
    output_dir: str | Path = "nebwalk_mlip_round0"
    export_formats: tuple[str, ...] = ("xyz", "traj", "json")
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class SelectedImage:
    """Metadata for an image selected for higher-level refinement."""

    index: int
    energy: float
    relative_energy: float
    reason: str = "peak_plus_neighbors"


@dataclass(frozen=True)
class MLIPActiveNEBResult:
    """Result bundle for an MLIP-assisted NEB run."""

    neb_result: Any
    selected_images: tuple[SelectedImage, ...]
    selected_indices: tuple[int, ...]
    output_dir: Path | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    @property
    def mlip_barrier(self) -> float:
        return float(self.neb_result.barrier)

    @property
    def barrier(self) -> float:
        return self.mlip_barrier


def _nebwalk_version() -> str:
    try:
        return importlib_metadata.version("nebwalk")
    except importlib_metadata.PackageNotFoundError:
        return "0.7.0"


def _validate_export_formats(export_formats: Sequence[str]) -> tuple[str, ...]:
    formats = tuple(export_formats)
    supported = {"xyz", "traj", "json"}
    unsupported = sorted(set(formats) - supported)
    if unsupported:
        raise ValueError(f"unsupported export format(s): {', '.join(unsupported)}")
    return formats


def export_selected_images(
    images: Sequence[Any],
    selected: Sequence[SelectedImage],
    output_dir: str | Path,
    export_formats: Sequence[str] = ("xyz", "traj", "json"),
    metadata: dict[str, Any] | None = None,
) -> Path:
    """Export selected NEB images and selection metadata."""
    formats = _validate_export_formats(export_formats)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    selected_atoms = [images[item.index].copy() for item in selected]

    if "xyz" in formats:
        for item, atoms in zip(selected, selected_atoms):
            write(out / f"selected_image_{item.index:02d}.xyz", atoms)

    if "traj" in formats:
        write(out / "selected_images.traj", selected_atoms)

    if "json" in formats:
        payload: dict[str, Any] = {
            "schema": "nebwalk.selected_images.v1",
            "selection_strategy": selected[0].reason if selected else None,
            "selected_indices": [item.index for item in selected],
            "selected_images": [asdict(item) for item in selected],
            "notes": (
                "Selected images are intended for higher-level DFT/QE refinement."
            ),
        }
        if metadata is not None:
            payload["metadata"] = metadata
        with (out / "selected_images.json").open("w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")

    readme = out / "README.md"
    readme.write_text(
        "# nebwalk selected images\n\n"
        "These images were selected from an MLIP-assisted NEB path using the "
        "`peak_plus_neighbors` strategy.\n\n"
        "They are intended for DFT/QE refinement, single-point validation, or "
        "later active-learning labeling.\n\n"
        "This folder does not contain a complete DFT workflow by itself.\n",
        encoding="utf-8",
    )
    return out


def run_mlip_assisted_neb(
    initial: Any,
    final: Any,
    mlip_calculator_factory: Any,
    neb_config: NEBRunConfig | None = None,
    active_config: MLIPActiveNEBConfig | None = None,
) -> MLIPActiveNEBResult:
    """Run an MLIP-assisted NEB and export barrier-sensitive images."""
    neb_cfg = neb_config or NEBRunConfig()
    active_cfg = active_config or MLIPActiveNEBConfig()

    neb_result = run_neb_calculation(
        initial=initial,
        final=final,
        calculator_factory=mlip_calculator_factory,
        config=neb_cfg,
    )

    images = neb_result.neb.images
    energies = [float(energy) for energy in neb_result.neb.get_energies()]
    selected_indices = select_images(
        energies=energies,
        strategy=active_cfg.selection_strategy,
        n_select=active_cfg.n_select,
        include_endpoints=active_cfg.include_endpoints,
    )
    reference_energy = energies[0]
    selected = tuple(
        SelectedImage(
            index=idx,
            energy=energies[idx],
            relative_energy=energies[idx] - reference_energy,
            reason=active_cfg.selection_strategy,
        )
        for idx in selected_indices
    )

    result_metadata = {
        "nebwalk_version": _nebwalk_version(),
        "stage": "mlip_assisted_neb",
        **active_cfg.metadata,
    }
    output_dir = None
    if active_cfg.export_selected:
        output_dir = export_selected_images(
            images=images,
            selected=selected,
            output_dir=active_cfg.output_dir,
            export_formats=active_cfg.export_formats,
            metadata=result_metadata,
        )

    return MLIPActiveNEBResult(
        neb_result=neb_result,
        selected_images=selected,
        selected_indices=tuple(selected_indices),
        output_dir=output_dir,
        metadata=result_metadata,
    )
