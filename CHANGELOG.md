# Changelog

All notable changes to **nebwalk** are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.7.0] - MLIP-assisted NEB

### Added
- Added `nebwalk.active` with `run_mlip_assisted_neb()`.
- Added `MLIPActiveNEBConfig`, `MLIPActiveNEBResult`, and `SelectedImage`.
- Added `nebwalk.selection` with `peak_plus_neighbors` image selection.
- Added selected-image export to `.xyz`, `.traj`, and `.json`.
- Added examples for MLIP-assisted NEB using EMT and a MACE template.
- Added tests for selection and active workflow exports.

### Notes
- This release introduces an active-learning-ready MLIP-assisted workflow.
- It does not yet implement uncertainty-guided selection, automatic MLIP
  retraining, adaptive image insertion/removal, or QE failed-image recovery.

---

## [0.6.0] — 2026-06-10

### Added
- **Geodesic interpolation** (`geodesic_interpolate`) — mass-weighted Cartesian
  geodesic path; superior to IDPP for large conformational changes
- **Quantum ESPRESSO interface** (`nebwalk.qe`): `QEParams`, `make_qe_factory`,
  `validate_qe_setup` — generate and validate QE PWSCF input files for
  DFT-level NEB workflows
- `.gitattributes` for consistent cross-platform line endings
- `docs/` placeholder for future Sphinx documentation
- New validated examples (14 systems total):
  - Cu adatom diffusion on Cu(100) / EMT (4.6% vs DFT-PBE)
  - Ni adatom diffusion on Ni(100) / EMT (12% vs DFT-PBE)
  - Li vacancy migration in Li2O / MACE-MP-0 (1.4% vs DFT-GGA, ref. 0.28 eV)
  - Mg vacancy migration in MgO / MACE-MP-0 (2.5% vs DFT-PBE)
  - H diffusion on Cu(111) / MACE-MP-0 (~16%, documented surface limitation)

### Changed
- Test suite expanded: 68 -> 94 tests (QE interface coverage added)
- Repo root cleaned: model files, logs, outputs added to .gitignore

---

## [0.5.0] — FILL-FROM-GIT

### Added
- Improved tangent estimate (Henkelman & Jonsson, J. Chem. Phys. 113, 9978, 2000)
- IDPP interpolation (Smidstrup et al., J. Chem. Phys. 141, 214106, 2014)
- MIC-aware periodic boundary conditions
- Variable spring constants (Lindh et al., Chem. Phys. Lett. 241, 423, 1995)
- FIRE optimizer (Bitzek et al., PRL 97, 170201, 2006)
- High-level API: NEBRunConfig, run_neb_calculation, NEBRunResult
- Thread-based parallel image evaluation

### Fixed
- Convergence criterion: per-atom force magnitude (norm-based)
  Confirmed by ethane barrier shift 0.108 -> 0.113 eV

### Removed
- L-BFGS-B optimizer (inappropriate for non-conservative NEB force field)

---

## [0.4.0] — FILL-FROM-GIT

### Added
- Core NEB: spring forces, perpendicular force projection, upwind tangent
- Climbing Image NEB (CI-NEB)
- Linear interpolation (MIC-aware)
- ASE compatibility, energy profile plot, CSV, .traj output
- GitHub Actions CI, MIT license
