# nebwalk

[![CI](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml/badge.svg)](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/nebwalk.svg)](https://pypi.org/project/nebwalk/)
[![Python](https://img.shields.io/badge/python-3.9%20%7C%203.10%20%7C%203.11%20%7C%203.12-blue)](https://pypi.org/project/nebwalk/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

**nebwalk** is a lightweight, transparent, ASE-compatible Nudged Elastic Band
(NEB/CI-NEB) library for minimum-energy paths, transition-state barriers, and
diffusion mechanisms in atomistic simulations.

It is designed for fast prototyping with classical and machine-learned
calculators, while remaining compatible with DFT backends through ASE.

```bash
pip install nebwalk
```

Current source version: **v0.6.0**.

> Note for maintainers: if PyPI still shows an older release after pushing this
> README, publish the v0.6.0 distribution with the commands in the release
> section below. The source package metadata is already set to `0.6.0`.

---

## Why nebwalk?

`nebwalk` focuses on one task: making NEB workflows simple, inspectable, and
calculator-agnostic.

It is useful when you want to:

- build and optimize minimum-energy paths;
- run standard NEB or climbing-image NEB;
- switch between EMT, MACE-MP-0, Egret-1t, Quantum ESPRESSO, VASP, or any ASE-compatible calculator;
- prototype surface diffusion, vacancy migration, adsorbate hopping, and molecular conformational changes;
- export profiles as plots, CSV files, and ASE trajectories;
- keep the NEB implementation transparent enough to inspect, test, and modify.

This is not a black-box workflow manager. It is a compact research-code layer
around ASE calculators and atomic structures.

---

## Features in v0.6.0

- Standard NEB and **climbing-image NEB (CI-NEB)**.
- Improved tangent estimate following Henkelman and Jónsson.
- Perpendicular potential-force projection and spring-force projection.
- **FIRE optimizer**, suitable for the non-conservative NEB force field.
- Linear, IDPP, and regularized geodesic-style interpolation.
- Minimum-image-convention-aware interpolation and NEB forces for periodic systems.
- Variable spring constants to concentrate images near high-energy regions.
- High-level `run_neb_calculation()` API with calculator factories.
- Thread-based parallel image evaluation.
- Restart from ASE `.traj` files with fresh calculator instances.
- Energy-profile plotting, CSV export, and `.traj` output.
- Quantum ESPRESSO helper layer through ASE calculator construction.

---

## Installation

### Stable installation from PyPI

```bash
pip install nebwalk
```

Optional MACE support:

```bash
pip install "nebwalk[mace]"
```

Check the installed version:

```bash
python -c "import nebwalk; print(nebwalk.__version__ if hasattr(nebwalk, '__version__') else 'installed')"
```

### Latest source from GitHub

```bash
pip install git+https://github.com/Rifat19R/nebwalk.git
```

### Development installation

```bash
git clone https://github.com/Rifat19R/nebwalk.git
cd nebwalk
pip install -e ".[dev]"
pytest tests/ -v
```

---

## Quick start

```python
from ase.calculators.emt import EMT
from nebwalk import NEBRunConfig, run_neb_calculation

initial = ...  # relaxed ase.Atoms endpoint
final = ...    # relaxed ase.Atoms endpoint

config = NEBRunConfig(
    n_images=7,
    interpolation="idpp",
    k=0.1,
    k_min=0.033,
    climb=True,
    climb_delay=50,
    fmax=0.05,
    max_steps=500,
)

result = run_neb_calculation(
    initial=initial,
    final=final,
    calculator_factory=lambda: EMT(),
    config=config,
)

print(f"Converged       : {result.converged}")
print(f"Barrier         : {result.barrier:.3f} eV")
print(f"Reverse barrier : {result.reverse_barrier:.3f} eV")
print(f"Reaction energy : {result.reaction_energy:.3f} eV")

result.neb.plot("profile.png")
result.neb.save_csv("profile.csv")
result.neb.save_trajectory("path.traj")
```

The calculator factory must return a fresh calculator instance. Do not share one
ASE calculator object across all images.

---

## Example with MACE-MP-0

```python
from mace.calculators import mace_mp
from nebwalk import NEBRunConfig, run_neb_calculation

def make_calc():
    return mace_mp(
        model="medium",
        dispersion=False,
        default_dtype="float64",
        device="cpu",  # use "cuda" if available
    )

config = NEBRunConfig(
    n_images=5,
    interpolation="idpp",
    k=0.5,
    climb=True,
    climb_delay=50,
    fmax=0.03,
    max_steps=1000,
)

result = run_neb_calculation(initial, final, make_calc, config)
```

MACE-MP-0 is a PBE-level foundation model. It can be excellent for rapid
prototyping, but absolute barriers should be validated against DFT or experiment
for the specific chemistry.

---

## What materials can nebwalk handle?

`nebwalk` operates on ASE `Atoms` objects, so the practical materials space is
defined by the calculator you attach.

Currently demonstrated or directly supported workflow classes include:

| Class | Examples |
|---|---|
| Molecular paths | H3 Morse benchmark, ethane torsion |
| Surface diffusion | Al/Cu/Ni/Pd/Ru slab adsorbate hopping, H on Cu(111) |
| Bulk vacancy migration | fcc metals, hcp Mg, simple oxides |
| Ionic migration | MgO, Li2O-style vacancy/interstitial paths |
| 2D/slab systems | MXenes, MAX phases, graphene-derived slabs, catalytic surfaces |
| DFT-backed systems | Quantum ESPRESSO via ASE calculator helpers |

The code is calculator-agnostic. The scientific reliability of any result still
depends on the chosen calculator, pseudopotentials, k-points, slab thickness,
coverage, spin state, and convergence settings.

---

## Benchmark spotlight: H diffusion on Cu(111)

A reproducible benchmark is included in:

```text
benchmarks/h_cu111_mace_mp/
```

The measured result uses MACE-MP-0 with a 4×4×4 Cu(111) slab, one H adatom
(1/16 ML), 15 Å vacuum, bottom two Cu layers fixed, IDPP interpolation, and
CI-NEB.

### Elementary fcc → hcp hop

| Quantity | Result |
|---|---:|
| Calculator | MACE-MP-0 small |
| Path | fcc hollow → bridge-like TS → hcp hollow |
| Internal images | 5 |
| Converged | True |
| Wall time | 2.67 min |
| Forward barrier | 125.92 meV |
| Reverse barrier | 136.67 meV |
| Reaction energy | −10.75 meV |
| TS image | 3 |

The profile is a clean single-saddle elementary hop:

```text
image 00:   +0.000 meV   fcc initial
image 01:  +44.614 meV
image 02: +111.234 meV
image 03: +125.917 meV   transition state
image 04: +104.770 meV
image 05:  +33.294 meV
image 06:  -10.749 meV   hcp final
```

### fcc → fcc full hop

A longer fcc-to-fcc hop also converged. It shows the expected two-saddle
sequence:

```text
fcc → bridge TS → hcp hollow → bridge TS → fcc
```

| Quantity | Result |
|---|---:|
| Calculator | MACE-MP-0 small |
| Internal images | 7 |
| Converged | True |
| Wall time | 6.42 min |
| Forward barrier | 125.56 meV |
| Reverse barrier | 125.41 meV |
| Reaction energy | +0.146 meV |

This confirms endpoint equivalence and path symmetry.

### Interpretation

The Cu(111) benchmark validates the NEB workflow: endpoint handling,
minimum-image interpolation, CI-NEB optimization, image-energy reporting, and
profile export. The absolute MACE barrier should not be overinterpreted as an
experimental activation energy. Experiments on H/Cu(111) report macroscopic
diffusion behavior affected by thermal activation, quantum effects, coverage,
and surface morphology. MACE-MP-0 is trained to reproduce PBE-level energetics,
so a lower static flat-terrace barrier is expected.

---

## Other examples

```bash
python examples/morse_h3.py
python examples/al_diffusion_emt.py
python examples/ethane_egret.py
python examples/al_vacancy_macemp.py
python examples/mg_vacancy_macemp.py
python examples/al_diffusion_qe.py
```

Some calculators are intentionally optional. Egret-1t model files and MACE model
weights are not distributed with this repository.

---

## Testing

```bash
pip install -e ".[test]"
pytest tests/ -v
```

The test suite covers interpolation, NEB force projection, tangent construction,
minimum-image convention handling, variable springs, parallel image evaluation,
restart helpers, and calculator-factory workflows.

---

## Roadmap

Short-term priorities:

1. Keep PyPI, GitHub tags, and source metadata synchronized for v0.6.x releases.
2. Add H/Cu(111) benchmark scripts, static output CSVs, and profile plots.
3. Add one DFT-backed Quantum ESPRESSO benchmark with small inputs and clear cost warnings.
4. Add post-NEB transition-state refinement using a dimer method.
5. Build a documentation site with theory, API usage, calculator setup, and benchmarks.

Long-term priorities:

- dimer/RFO transition-state refinement;
- imaginary-mode validation utilities;
- better benchmark provenance files;
- larger measured MACE-medium/large CPU/GPU parallel-scaling benchmarks;
- more surface diffusion and ionic migration examples.

---

## Release checklist for maintainers

Use this when publishing a new release:

```bash
python -m pip install --upgrade build twine
rm -rf dist/ build/ *.egg-info
python -m build
python -m twine check dist/*
python -m twine upload dist/*
git tag -a v0.6.0 -m "nebwalk v0.6.0"
git push origin main --tags
```

After release:

```bash
pip install --upgrade nebwalk
python -c "import nebwalk; print(nebwalk.__version__ if hasattr(nebwalk, '__version__') else 'installed')"
```

---

## References

1. H. Jónsson, G. Mills, K. W. Jacobsen, *Nudged Elastic Band Method for Finding Minimum Energy Paths of Transitions*, World Scientific, 1998.
2. G. Henkelman and H. Jónsson, *Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points*, J. Chem. Phys. 113, 9978 (2000).
3. G. Henkelman, B. P. Uberuaga, and H. Jónsson, *A climbing image nudged elastic band method for finding saddle points and minimum energy paths*, J. Chem. Phys. 113, 9901 (2000).
4. E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, and P. Gumbsch, *Structural Relaxation Made Simple*, Phys. Rev. Lett. 97, 170201 (2006).
5. S. Smidstrup, A. Pedersen, K. Stokbro, and H. Jónsson, *Improved initial guess for minimum energy path calculations*, J. Chem. Phys. 141, 214106 (2014).
6. MACE-MP-0 documentation and foundation model papers for ML-potential context.
7. H/Cu(111) diffusion literature should be treated as macroscopic diffusion reference, not a one-to-one static NEB barrier target.

---

## License

MIT License. See [LICENSE](LICENSE).
