# Changelog

All notable changes to **nebwalk** are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.7.1] — 2026-06-12

### Fixed
- Improved-tangent bisection branch now weights unit vectors by energy
  differences per Henkelman-Jonsson Eq. 10, avoiding displacement-length bias
  at extremum images.
- Variable spring constants now reference the global path minimum, preserving
  spring ordering when intermediates lie below both endpoints.
- CUDA-backed calculators now emit a warning when used with thread-parallel
  image evaluation because CUDA force evaluation from multiple Python threads
  can silently corrupt results.
- Exported selected-images README now reports the actual selection strategy.
- Example benchmark runner now defaults to the full suite, avoids duplicate
  safe-mode runs, documents all QE pseudopotentials, uses `np=4` with Open MPI
  oversubscription for WSL, and cleans stale QE workdirs before QE reruns.
- QE example runner now preserves a user-provided `ESPRESSO_COMMAND` and falls
  back to an absolute `pw.x` path, avoiding PATH-dependent failures in ASE
  subprocesses.

### Docs
- Tangent fix may shift converged barriers slightly for paths with uneven image
  spacing. Benchmark table values will be re-validated in a follow-up run.

---

## [0.7.0] — 2026-06-12

### Added
- `nebwalk.active` module with `run_mlip_assisted_neb()` high-level workflow.
- `MLIPActiveNEBConfig`, `MLIPActiveNEBResult`, and `SelectedImage` dataclasses.
- `nebwalk.selection` module with `peak_plus_neighbors` image selection strategy.
- Selected-image export to `.xyz`, `.traj`, and `.json` formats.
- Examples for MLIP-assisted NEB using EMT and MACE template.
- Tests for selection logic and active workflow exports.

### Changed
- Test suite expanded: 94 → 111 tests.

### Notes
- Introduces an active-learning-ready MLIP-assisted NEB workflow.
- Does not yet implement: uncertainty-guided selection, automatic MLIP
  retraining, adaptive image insertion/removal, or QE failed-image recovery.

---

## [0.6.0] — 2026-06-10

### Added
- **Geodesic interpolation** (`geodesic_interpolate`) — mass-weighted Cartesian
  geodesic path; superior to IDPP for large conformational changes.
- **Quantum ESPRESSO interface** (`nebwalk.qe`): `QEParams`, `make_qe_factory`,
  `validate_qe_setup` — generate and validate QE PWSCF input files for
  DFT-level NEB workflows.
- `.gitattributes` for consistent cross-platform line endings.
- `docs/` placeholder for future Sphinx documentation.
- New validated examples (14 systems total):
  - Cu adatom diffusion on Cu(100) / EMT (4.6% vs DFT-PBE)
  - Ni adatom diffusion on Ni(100) / EMT (12% vs DFT-PBE)
  - Li vacancy migration in Li₂O / MACE-MP-0 (1.4% vs DFT-GGA, ref. 0.28 eV)
  - Mg vacancy migration in MgO / MACE-MP-0 (2.5% vs DFT-PBE)
  - H diffusion on Cu(111) / MACE-MP-0 (~16%, documented surface limitation)

### Changed
- Test suite expanded: 68 → 94 tests (QE interface coverage added).
- Repo root cleaned: model files, logs, outputs added to `.gitignore`.

---

## [0.5.0] — 2026-06-08

### Added
- Improved tangent estimate (Henkelman & Jónsson, J. Chem. Phys. 113, 9978, 2000).
- IDPP interpolation (Smidstrup et al., J. Chem. Phys. 141, 214106, 2014).
- MIC-aware periodic boundary conditions.
- Variable spring constants (Lindh et al., Chem. Phys. Lett. 241, 423, 1995).
- FIRE optimizer (Bitzek et al., PRL 97, 170201, 2006).
- High-level API: `NEBRunConfig`, `run_neb_calculation`, `NEBRunResult`.
- Thread-based parallel image evaluation.

### Fixed
- Convergence criterion changed to per-atom force magnitude (norm-based).
  Confirmed by ethane barrier shift: 0.108 → 0.113 eV.

### Removed
- L-BFGS-B optimizer (inappropriate for non-conservative NEB force field).

---

## [0.4.0] — 2026-06-01

### Added
- Core NEB: spring forces, perpendicular force projection, upwind tangent.
- Climbing Image NEB (CI-NEB).
- Linear interpolation (MIC-aware).
- ASE compatibility, energy profile plot, CSV, `.traj` output.
- GitHub Actions CI, MIT license.
