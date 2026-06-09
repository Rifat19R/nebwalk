# Changelog

All notable changes to **nebwalk** are documented here.  
Format: [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).  
Versioning: [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
  - Li vacancy migration in Li₂O / MACE-MP-0 (1.4% vs DFT-GGA, ref. 0.28 eV)
  - Mg vacancy migration in MgO / MACE-MP-0 (2.5% vs DFT-PBE)
  - H diffusion on Cu(111) / MACE-MP-0 (∼16%, documented surface limitation)

### Changed
- Test suite expanded: 68 → 94 tests (QE interface coverage added)
- Repo root cleaned: model files, logs, and computed outputs added to `.gitignore`

---

## [0.5.0] — FILL-FROM-GIT

> `git log --format="%as" -- pyproject.toml | head -1`

### Added
- **Improved tangent estimate** (Henkelman & Jónsson, J. Chem. Phys. **113**, 9978, 2000):
  correct handling of local extrema; replaces simple upwind tangent
- **IDPP interpolation** (Smidstrup et al., J. Chem. Phys. **141**, 214106, 2014):
  preserves short-range bonding geometry; recommended over linear for torsional paths
- **MIC-aware periodic boundary conditions**: minimum image convention for
  tangent and spring-force displacements
- **Variable spring constants** (Lindh et al., Chem. Phys. Lett. **241**, 423, 1995):
  energy-weighted springs concentrate images near the saddle point
- **FIRE optimizer** (Bitzek et al., PRL **97**, 170201, 2006): replaces L-BFGS-B;
  robust to non-conservative NEB force field
- **High-level API**: `NEBRunConfig`, `run_neb_calculation`, `NEBRunResult`
- **Thread-based parallel image evaluation** via `concurrent.futures.ThreadPoolExecutor`

### Fixed
- Convergence criterion: switched from max absolute force component to per-atom
  force magnitude (norm-based). Confirmed by ethane barrier shift 0.108 → 0.113 eV

### Removed
- L-BFGS-B optimizer (inappropriate for non-conservative NEB force field)

---

## [0.4.0] — FILL-FROM-GIT

> `git log --format="%as" -- pyproject.toml | tail -1`

### Added
- Core NEB: spring forces, perpendicular force projection, upwind tangent
- Climbing Image NEB (CI-NEB) with configurable activation delay
- Linear interpolation (MIC-aware)
- ASE compatibility: any `Calculator`-attached `Atoms` list
- Output: energy profile plot (`.png`), CSV, ASE `.traj` trajectory
- GitHub Actions CI (Python 3.9–3.12), MIT license

---

## References

1. H. Jónsson, G. Mills, K.W. Jacobsen, *Classical and Quantum Dynamics in
   Condensed Phase Simulations*, World Scientific (1998).
2. G. Henkelman, H. Jónsson, J. Chem. Phys. **113**, 9978 (2000).
   DOI: [10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
3. G. Henkelman, B.P. Uberuaga, H. Jónsson, J. Chem. Phys. **113**, 9901 (2000).
   DOI: [10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
4. E. Bitzek et al., PRL **97**, 170201 (2006).
   DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)
5. S. Smidstrup et al., J. Chem. Phys. **141**, 214106 (2014).
   DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)
6. R. Lindh et al., Chem. Phys. Lett. **241**, 423 (1995).
