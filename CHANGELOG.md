# Changelog

## v0.5.0 - 2026-06-08

- Fix ASE-style convergence criterion using max per-atom force magnitude.
- Add calculator isolation validation to prevent ASE cache corruption.
- Replace package `print()` output with module logging.
- Refactor FIRE state into `FIREOptimizer`.
- Add restart support with `NEB.from_trajectory()`.
- Add reverse barrier and reaction energy helpers.
- Add shared NEB engine utilities.
- Add regularized geodesic-style interpolation.
- Add GitHub Actions CI and PyPI publish workflow.
- Add Quantum ESPRESSO example template for native Linux.
- Highlight Mg HCP vacancy benchmark: 0.508 eV vs ~0.52 eV DFT-PBE.
