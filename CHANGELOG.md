# Changelog

All notable changes to nebwalk are documented here.
Versions follow [Semantic Versioning](https://semver.org/).

---

## [0.6.0] — 2025

### Added
- `nebwalk.qe` module: Quantum ESPRESSO interface via `QEParams`,
  `make_qe_factory()`, and `validate_qe_setup()`
- Per-image directory isolation for QE — prevents wavefunction file
  conflicts during parallel NEB evaluation
- `examples/template_qe_neb.py`: annotated QE + nebwalk template
- `tests/test_qe.py`: 26 unit tests for QE interface (no pw.x required)

### Notes
- Requires ASE >= 3.22 for `ase.calculators.espresso.Espresso`
- No new required dependencies — QE support is built on ASE's built-in Espresso

---

## [0.5.0] — 2025

### Added
- Improved tangent estimate (Henkelman & Jónsson 2000, full weighted scheme)
- IDPP interpolation (Smidstrup et al. 2014) with PBC/MIC support
- Geodesic interpolation for large structural changes
- Variable spring constants (energy-weighted, Lindh et al. 1996)
- FIRE optimizer — replaces L-BFGS-B (NEB forces are non-conservative)
- Minimum Image Convention (MIC) for periodic bulk calculations
- High-level API: `NEBRunConfig`, `run_neb_calculation`, `NEBRunResult`
- Parallel image evaluation via `concurrent.futures.ThreadPoolExecutor`
- 94 tests across forces, interpolation, MIC, variable springs, QE interface

### Validated systems (14 total)
- Morse H3 collinear (4%), ethane C–C torsion/Egret-1t (10%)
- FCC vacancy series EMT: Cu 7.9%, Ni 5.3%, Pd 7.8%, Ag 3.3%
- Surface adatom EMT: Cu/Cu(100) 4.6%, Ni/Ni(100) 12%
- Bulk oxide MACE-MP-0: Li vacancy Li₂O 1.4%, Mg vacancy MgO 2.5%
- HCP MACE-MP-0: Mg vacancy 2%

### Fixed
- Convergence criterion: per-atom force magnitude (norm), not max component
- Atom ordering consistency for vacancy endpoint construction

---

## [0.4.0] — 2025

### Added
- Initial release: NEB class, linear interpolation, CI-NEB, basic FIRE
- EMT and Egret-1t calculator support
- Energy profile plot, CSV, and trajectory output
