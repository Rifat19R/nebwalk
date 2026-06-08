# nebwalk

[![CI](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml/badge.svg)](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml)
![Python](https://img.shields.io/badge/python-3.9%20%7C%203.11%20%7C%203.12-blue)
![License](https://img.shields.io/badge/license-MIT-green)

Minimal, correct Python implementation of the **Nudged Elastic Band (NEB)** method for finding minimum energy paths (MEPs) and transition states. Works with any ASE-compatible calculator â€” EMT for testing, [Egret-1](https://rowansci.com) or MACE-MP-0 for production.

**GitHub:** [github.com/Rifat19R/nebwalk](https://github.com/Rifat19R/nebwalk)

**Result spotlight:** HCP Mg vacancy migration with MACE-MP-0 gives **0.508 eV**
against **~0.52 eV DFT-PBE** (about **2% error**). This is the cleanest
validated result in the repo and the best lead example for outreach.

---

## Features

- Improved tangent estimate (Henkelman & JÃ³nsson 2000)
- Spring forces + perpendicular potential force projection
- Climbing Image NEB (CI-NEB) for true saddle-point location
- **FIRE optimizer** â€” robust to the non-conservative NEB force field
- **IDPP interpolation** (Smidstrup et al. 2014) â€” chemically sensible initial paths, avoids atomic clashes in torsional reactions
- **Regularized geodesic-style interpolation** - IDPP plus short-range repulsion for severe atom-overlap paths
- **Minimum Image Convention (MIC)** â€” correct tangents and spring forces for periodic systems
- **Variable spring constants** â€” energy-weighted springs concentrate images near the saddle point
- **Parallel image evaluation** â€” thread-based; beneficial for large supercells (>50 atoms) or GPU-accelerated MACE. Not recommended for small molecules on CPU where thread overhead dominates.
- Energy profile plot, CSV export, ASE `.traj` output
- Restart from ASE `.traj` with fresh calculator factories
- Local editable install â€” no compiled extensions. PyPI release is pending.

---

## Installation

```bash
# PyPI release pending: pip install nebwalk is not available yet.
pip install .
# or, for development:
pip install -e ".[test]"
```

For Egret-1 / MACE support:

```bash
pip install -e ".[mace]"
```

---

## Quick start

```python
from ase.calculators.emt import EMT
from nebwalk import NEB, idpp_interpolate

# 1. Define endpoints (relaxed structures)
initial = ...   # ase.Atoms
final   = ...   # ase.Atoms

# 2. Build interpolated path (IDPP recommended over linear for molecules)
images = idpp_interpolate(initial, final, n_images=7)
for img in images:
    img.calc = EMT()

# 3. Run NEB with variable spring constants
neb = NEB(images, k=0.1, k_min=0.033, climb=True)
neb.optimize(fmax=0.05)

# 4. Output
print(f"Barrier: {neb.get_barrier():.3f} eV")
neb.plot("profile.png")
neb.save_csv("profile.csv")
neb.save_trajectory("path.traj")
```

Visualise the trajectory: `ase gui path.traj`

---

## With Egret-1t (MACE) + parallel evaluation

```python
from mace.calculators import MACECalculator
from nebwalk import NEB, idpp_interpolate

# Download EGRET_1T.model from https://rowansci.com â€” do not commit to repo
images = idpp_interpolate(initial, final, n_images=9)
for img in images:
    img.calc = MACECalculator(
        model_paths="EGRET_1T.model",
        device="cpu",
        default_dtype="float32",
    )

# n_workers=9: one thread per intermediate image
# Use for large periodic systems (>50 atoms) or GPU; for small molecules on CPU, n_workers=1 is faster
neb = NEB(images, k=0.1, k_min=0.033, climb=True, n_workers=9)
neb.optimize(fmax=0.05)
```

> **Note:** Egret-1 does not compute stress/virial; only atomic forces are needed
> for NEB, so this is fine.
> The model file is not distributed with this repository. Download it from
> [rowansci.com](https://rowansci.com).

---

## Periodic systems (bulk diffusion, vacancy migration)

```python
from ase.build import bulk
from ase.calculators.emt import EMT
from nebwalk import NEB, linear_interpolate

al = bulk("Al", "fcc", a=4.05, cubic=True).repeat(2)
del al[0]   # create vacancy

start = al.copy(); start.calc = EMT()
end   = al.copy(); end.positions[0] += [2.025, 0, 0]; end.calc = EMT()

# linear_interpolate is MIC-aware for periodic systems
# idpp_interpolate also supports PBC
images = linear_interpolate(start, end, n_images=7)
for img in images:
    img.calc = EMT()

neb = NEB(images, k=0.1, k_min=0.033, climb=True)
neb.optimize(fmax=0.05)
```

---

## Theory

### Interpolation

**Linear interpolation** (`linear_interpolate`) divides the displacement vector
from start to end into equal steps. Fast, but fails for torsional reactions â€”
atoms can pass through each other, producing unphysical initial paths.

**IDPP interpolation** (`idpp_interpolate`) minimises a weighted pairwise
distance objective across images (Smidstrup et al. 2014):

```
S = Î£_{i<j} w_ij (d_ij - d_target_ij)Â²     w_ij = 1/d_target_ij^4
```

This preserves short-range bonding structure throughout the path. Recommended
for any system with significant torsional motion or risk of atomic overlap.
Both methods support PBC via the minimum image convention.

Reference: Smidstrup, Pedersen, Stokbro, JÃ³nsson, J. Chem. Phys. **141**,
214106 (2014). DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)

**Regularized geodesic-style interpolation** (`geodesic_interpolate`) is a self-contained
regularized IDPP approximation. It starts from the MIC-aware Cartesian path,
separates exact atom overlaps deterministically, and optimizes an IDPP objective
with a short-range repulsive penalty. This is not a full redundant-internal
coordinate implementation of Zhu et al.; it is a robust, dependency-free first
line of defense against pathological Cartesian paths.

---

### Tangent estimate (improved, Henkelman & JÃ³nsson 2000)

Let Ï„âº = R_{i+1} âˆ’ Ráµ¢ and Ï„â» = Ráµ¢ âˆ’ R_{i-1}.

| Condition | Tangent |
|-----------|---------|
| E_{i+1} > Eáµ¢ > E_{i-1} | Ï„âº / \|Ï„âº\| |
| E_{i-1} > Eáµ¢ > E_{i+1} | Ï„â» / \|Ï„â»\| |
| Local extremum (E_{i+1} â‰¥ E_{i-1}) | (Ï„âº Î”E_max + Ï„â» Î”E_min) normalised |
| Local extremum (E_{i-1} > E_{i+1}) | (Ï„âº Î”E_min + Ï„â» Î”E_max) normalised |

where Î”E_max = max(|E_{i+1}âˆ’Eáµ¢|, |E_{i-1}âˆ’Eáµ¢|) and Î”E_min = min(â€¦).
For periodic systems, Ï„âº and Ï„â» use MIC displacements.

---

### NEB force

```
F_spring = k * (|R_{i+1} - R_i| - |R_i - R_{i-1}|) * Ï„
F_perp   = F_pot - (F_pot Â· Ï„) Ï„          # F_pot = ASE get_forces() = -âˆ‡E
F_NEB    = F_spring + F_perp
```

With variable spring constants, forward and backward spring constants can
differ per image:

```
F_spring = (k_fwd * |R_{i+1} - R_i| - k_bwd * |R_i - R_{i-1}|) * Ï„
```

---

### Variable spring constants

Energy-weighted springs (Lindh et al. 1996) assign k_max to springs adjacent
to the saddle point and k_min to springs far from it:

```
k[i] = k_max - (k_max - k_min) * (E_ref - E_spring[i]) / (E_ref - E_low)
```

where E_spring[i] = max(E[i], E[i+1]), E_ref = max(all energies),
E_low = min(endpoint energies). Values are clipped to [k_min, k_max].

This concentrates images near the transition state, improving barrier
resolution without adding more images. Recommended ratio: k_min = k / 3.

---

### Climbing image

For the highest-energy movable image (after `climb_delay` steps):

```
F_CI = F_pot - 2 (F_pot Â· Ï„) Ï„
```

The spring is removed and the tangent component is inverted, driving the
image to the true saddle point.

---

### Convergence criterion

```
max per-atom |F_NEB| < fmax   (eV/Ã…, default 0.05)
```

---

### Optimizer: FIRE

FIRE (Fast Inertial Relaxation Engine) is used instead of L-BFGS-B because
NEB forces are **not the gradient of any single energy function** â€” the spring
projection and perpendicular projection break the conservative-field assumption
required for L-BFGS-B line search. FIRE's velocity-damping approach is
robust to this.

Reference: Bitzek et al., PRL **97**, 170201 (2006).
DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

---

### Parallel image evaluation

At each FIRE step, energy and forces for all intermediate images are
independent â€” they do not communicate until spring forces are assembled.
`nebwalk` exploits this with `concurrent.futures.ThreadPoolExecutor`.

Threads are used (not processes) because pickling ASE calculator objects is
unreliable across calculator types.

Measured on a 4-core CPU (WSL Ubuntu, Egret-1t, C2H6, 7 images):

| Calculator | Sequential | Parallel (n_workers=7) | Speedup |
|------------|------------|------------------------|---------|
| EMT        | 17.1 ms    | 6.5 ms                 | 2.6Ã—    |
| Egret-1t   | 3.8 ms     | 8.3 ms                 | 0.45Ã—   |

For small molecules on CPU, per-image compute (~0.5 ms) is faster than thread
overhead â€” use `n_workers=1`. Parallelism is beneficial when per-image compute
is large: heavy MACE models, large supercells (>50 atoms), or GPU evaluation.

---

## API

### `NEB(images, k=0.1, k_min=None, climb=False, climb_delay=100, n_workers=1)`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `images` | â€” | List of `ase.Atoms`, all with calculators attached |
| `k` | 0.1 | Spring constant (eV/Ã…Â²). With variable springs: maximum value. |
| `k_min` | None | Minimum spring constant for variable springs. `None` = uniform springs. Recommended: k / 3. |
| `climb` | False | Enable CI-NEB |
| `climb_delay` | 100 | FIRE steps before CI activates |
| `n_workers` | 1 | Threads for parallel image evaluation. Beneficial for large supercells (>50 atoms) or GPU-accelerated MACE. Use default (1) for small molecules on CPU. |

**Methods:**

| Method | Returns | Description |
|--------|---------|-------------|
| `optimize(fmax=0.05, max_steps=500, verbose=True)` | `bool` | Run FIRE. Returns True if converged. |
| `run(fmax=0.05, max_steps=500, verbose=True)` | `bool` | Alias for `optimize`, convenient for restarted runs. |
| `get_energies()` | `list[float]` | Potential energies for all images (eV) |
| `get_barrier()` | `float` | Forward barrier relative to image 0 (eV) |
| `get_reverse_barrier()` | `float` | Reverse barrier relative to final image (eV). |
| `get_reaction_energy()` | `float` | Final minus initial energy (eV). |
| `get_spring_constants()` | `ndarray` | Current spring constants, shape (N-1,) |
| `plot(filename, show, title)` | â€” | Save energy profile plot |
| `save_csv(filename)` | â€” | Write energies to CSV |
| `save_trajectory(filename)` | â€” | Write all images to ASE `.traj` |

---

### `linear_interpolate(start, end, n_images)`

Cartesian linear interpolation. MIC-aware for periodic systems.
Returns `[start_copy, img_1, ..., img_n, end_copy]` â€” no calculators
attached to intermediate images.

### `idpp_interpolate(start, end, n_images, max_iter=500, tol=1e-6)`

IDPP interpolation. MIC-aware for periodic systems.
Returns `[start_copy, img_1, ..., img_n, end_copy]` â€” no calculators
attached to intermediate images. Recommended over `linear_interpolate`
for torsional reactions and any path with risk of atomic overlap.

### `geodesic_interpolate(start, end, n_images, min_distance=0.75)`

Regularized geodesic-style interpolation. Use when Cartesian interpolation
creates severe atom overlaps or exact atom crossings. This method is
dependency-free and MIC-aware.

### `NEB.from_trajectory(path, calculator_factory, **neb_kwargs)`

Restart an interrupted run from an ASE `.traj` file. Calculators are not stored
in trajectory files, so `calculator_factory` must return one fresh calculator
instance per image.

```python
from ase.calculators.emt import EMT
from nebwalk import NEB

neb = NEB.from_trajectory("neb.traj", calculator_factory=EMT, climb=True)
neb.run(fmax=0.05)
```

### `variable_spring_constants(energies, k_max, k_min)`

Returns a `(N-1,)` array of spring constants given the current energy
profile. Called automatically inside `NEB` when `k_min` is set. Exposed
publicly for use in custom optimisation loops.

### `compute_neb_forces(images, k, climb, climb_index, energies, forces)`

Low-level force function. Returns a list of `(N_atoms, 3)` arrays.
`energies` and `forces` can be passed as pre-computed values to avoid
redundant calculator calls (used internally by the parallel evaluator).

---

## Running tests

```bash
pip install -e ".[test]"
pytest tests/ -v
```

Current suite: **76 tests** across forces, interpolation, MIC, variable
springs, parallel evaluation, restart helpers, and the shared engine.

---

## Examples

```bash
python examples/morse_h3.py           # collinear H+H2, Morse potential (self-contained)
python examples/al_diffusion_emt.py   # Al adatom diffusion on Al(100), EMT
python examples/ethane_egret.py       # ethane C-C torsion barrier, Egret-1t
python examples/al_vacancy_macemp.py  # Al vacancy migration in bulk Al, MACE-MP-0
python examples/ni_vacancy_emt.py     # Ni vacancy migration, EMT
python examples/pd_vacancy_emt.py     # Pd vacancy migration, EMT
python examples/ag_vacancy_emt.py     # Ag vacancy migration, EMT
python examples/pt_vacancy_emt.py     # Pt vacancy migration, EMT (see footnote)
python examples/mg_vacancy_macemp.py  # Mg vacancy in HCP Mg, MACE-MP-0
python examples/verify_egret.py       # sanity check: confirm Egret-1t loads correctly
```

### Validated results

| Example | Calculator | Barrier | Reference | Error |
|---------|------------|---------|-----------|-------|
| **Mg vacancy (HCP)** | **MACE-MP-0** | **0.508 eV** | **~0.52 eV (DFT-PBE)** | **2%** |
| Morse H3 | Morse (analytical) | 0.200 eV | 0.193 eV (exact) | 4% |
| Al adatom diffusion | EMT | 0.237 eV | ~0.40 eV (DFT-PBE) | finite-size slab |
| Ethane Câ€“C torsion | Egret-1t | 0.113 eV | 0.126 eV (exp.) | 10% |
| Al vacancy migration | MACE-MP-0 | 0.508 eV | 0.61 eV (DFT-PBE) | 17%â€  |
| Cu vacancy | EMT | 0.755 eV | ~0.70 eV (DFT-PBE) | 7.9% |
| Ni vacancy | EMT | 1.095 eV | ~1.04 eV (DFT-PBE) | 5.3% |
| Pd vacancy | EMT | 0.839 eV | ~0.91 eV (DFT-PBE) | 7.8% |
| Ag vacancy | EMT | 0.682 eV | ~0.66 eV (DFT-PBE) | 3.3% |
| Pt vacancy | EMT | 0.971 eV | ~1.49 eV (DFT-PBE) | 34.8%â€¡ |

â€  MACE-MP-0 systematically underestimates vacancy migration barriers by
10â€“20%. This is a known model limitation, not a nebwalk bug.

â€¡ EMT does not capture relativistic effects significant in Pt. The NEB
itself converged cleanly in 60 steps â€” the error is from the calculator.

---

## Release status

GitHub source install is ready. PyPI package name `nebwalk` is not published
yet, so `pip install nebwalk` will fail until the first PyPI release is made.
The next release target is `v0.5.0`.

Quantum ESPRESSO example support is documented in `examples/al_diffusion_qe.py`.
It is intentionally guarded because reliable QE production runs require native
Linux, a configured QE executable, and validated pseudopotentials.

---

## References

1. H. JÃ³nsson, G. Mills, K.W. Jacobsen, *Nudged Elastic Band Method for
   Finding Minimum Energy Paths of Transitions*, in *Classical and Quantum
   Dynamics in Condensed Phase Simulations*, World Scientific, 1998.
2. G. Henkelman, H. JÃ³nsson, *Improved tangent estimate in the nudged
   elastic band method for finding minimum energy paths and saddle points*,
   J. Chem. Phys. **113**, 9978 (2000).
   DOI: [10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
3. G. Henkelman, B.P. Uberuaga, H. JÃ³nsson, *A climbing image nudged
   elastic band method for finding saddle points and minimum energy paths*,
   J. Chem. Phys. **113**, 9901 (2000).
   DOI: [10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
4. E. Bitzek, P. Koskinen, F. GÃ¤hler, M. Moseler, P. Gumbsch, *Structural
   Relaxation Made Simple*, PRL **97**, 170201 (2006).
   DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)
5. S. Smidstrup, A. Pedersen, K. Stokbro, H. JÃ³nsson, *Improved initial
   guess for minimum energy path calculations*, J. Chem. Phys. **141**,
   214106 (2014).
   DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)
6. R. Lindh, A. Bernhardsson, G. KarlstrÃ¶m, P.-Ã…. Malmqvist, *On the use of
   a Hessian model function in molecular geometry optimizations*,
   Chem. Phys. Lett. **241**, 423 (1995). (Variable spring constant scheme.)

---

## License

MIT Â© Md. Rifat Khandaker
