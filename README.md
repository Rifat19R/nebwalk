# nebwalk

Minimal, correct Python implementation of the **Nudged Elastic Band (NEB)** method for finding minimum energy paths (MEPs) and transition states. Works with any ASE-compatible calculator — EMT for testing, [Egret-1](https://rowansci.com) or MACE-MP-0 for production.

**GitHub:** [github.com/Rifat19R/nebwalk](https://github.com/Rifat19R/nebwalk)

---

## Features

- Improved tangent estimate (Henkelman & Jónsson 2000)
- Spring forces + perpendicular potential force projection
- Climbing Image NEB (CI-NEB) for true saddle-point location
- **FIRE optimizer** — robust to the non-conservative NEB force field
- **IDPP interpolation** (Smidstrup et al. 2014) — chemically sensible initial paths, avoids atomic clashes in torsional reactions
- **Minimum Image Convention (MIC)** — correct tangents and spring forces for periodic systems
- **Variable spring constants** — energy-weighted springs concentrate images near the saddle point
- **Parallel image evaluation** — thread-based, reduces wall time for expensive calculators (MACE, Egret-1t)
- Energy profile plot, CSV export, ASE `.traj` output
- Single `pip install .` — no compiled extensions

---

## Installation

```bash
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

# Download EGRET_1T.model from https://rowansci.com — do not commit to repo
images = idpp_interpolate(initial, final, n_images=9)
for img in images:
    img.calc = MACECalculator(
        model_paths="EGRET_1T.model",
        device="cpu",
        default_dtype="float32",
    )

# n_workers=9: one thread per intermediate image
# Egret-1t (PyTorch) releases the GIL → genuine parallel speedup
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
from start to end into equal steps. Fast, but fails for torsional reactions —
atoms can pass through each other, producing unphysical initial paths.

**IDPP interpolation** (`idpp_interpolate`) minimises a weighted pairwise
distance objective across images (Smidstrup et al. 2014):

```
S = Σ_{i<j} w_ij (d_ij - d_target_ij)²     w_ij = 1/d_target_ij^4
```

This preserves short-range bonding structure throughout the path. Recommended
for any system with significant torsional motion or risk of atomic overlap.
Both methods support PBC via the minimum image convention.

Reference: Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. **141**,
214106 (2014). DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)

---

### Tangent estimate (improved, Henkelman & Jónsson 2000)

Let τ⁺ = R_{i+1} − Rᵢ and τ⁻ = Rᵢ − R_{i-1}.

| Condition | Tangent |
|-----------|---------|
| E_{i+1} > Eᵢ > E_{i-1} | τ⁺ / \|τ⁺\| |
| E_{i-1} > Eᵢ > E_{i+1} | τ⁻ / \|τ⁻\| |
| Local extremum (E_{i+1} ≥ E_{i-1}) | (τ⁺ ΔE_max + τ⁻ ΔE_min) normalised |
| Local extremum (E_{i-1} > E_{i+1}) | (τ⁺ ΔE_min + τ⁻ ΔE_max) normalised |

where ΔE_max = max(|E_{i+1}−Eᵢ|, |E_{i-1}−Eᵢ|) and ΔE_min = min(…).
For periodic systems, τ⁺ and τ⁻ use MIC displacements.

---

### NEB force

```
F_spring = k * (|R_{i+1} - R_i| - |R_i - R_{i-1}|) * τ
F_perp   = F_pot - (F_pot · τ) τ          # F_pot = ASE get_forces() = -∇E
F_NEB    = F_spring + F_perp
```

With variable spring constants, forward and backward spring constants can
differ per image:

```
F_spring = (k_fwd * |R_{i+1} - R_i| - k_bwd * |R_i - R_{i-1}|) * τ
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
F_CI = F_pot - 2 (F_pot · τ) τ
```

The spring is removed and the tangent component is inverted, driving the
image to the true saddle point.

---

### Convergence criterion

```
max |F_NEB|_component < fmax   (eV/Å, default 0.05)
```

---

### Optimizer: FIRE

FIRE (Fast Inertial Relaxation Engine) is used instead of L-BFGS-B because
NEB forces are **not the gradient of any single energy function** — the spring
projection and perpendicular projection break the conservative-field assumption
required for L-BFGS-B line search. FIRE's velocity-damping approach is
robust to this.

Reference: Bitzek et al., PRL **97**, 170201 (2006).
DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

---

### Parallel image evaluation

At each FIRE step, energy and forces for all intermediate images are
independent — they do not communicate until spring forces are assembled.
`nebwalk` exploits this with `concurrent.futures.ThreadPoolExecutor`.

Threads are used (not processes) because pickling ASE calculator objects is
unreliable across calculator types. PyTorch-based calculators (MACE, Egret-1t)
release the GIL during C++ computation, so threads give genuine parallelism.
Pure-Python calculators (EMT) see no speedup from threading but are not
harmed by it.

---

## API

### `NEB(images, k=0.1, k_min=None, climb=False, climb_delay=100, n_workers=1)`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `images` | — | List of `ase.Atoms`, all with calculators attached |
| `k` | 0.1 | Spring constant (eV/Å²). With variable springs: maximum value. |
| `k_min` | None | Minimum spring constant for variable springs. `None` = uniform springs. Recommended: k / 3. |
| `climb` | False | Enable CI-NEB |
| `climb_delay` | 100 | FIRE steps before CI activates |
| `n_workers` | 1 | Threads for parallel image evaluation. Set to n_images for maximum speedup with MACE/Egret-1t. |

**Methods:**

| Method | Returns | Description |
|--------|---------|-------------|
| `optimize(fmax=0.05, max_steps=500, verbose=True)` | `bool` | Run FIRE. Returns True if converged. |
| `get_energies()` | `list[float]` | Potential energies for all images (eV) |
| `get_barrier()` | `float` | Forward barrier relative to image 0 (eV) |
| `get_spring_constants()` | `ndarray` | Current spring constants, shape (N-1,) |
| `plot(filename, show, title)` | — | Save energy profile plot |
| `save_csv(filename)` | — | Write energies to CSV |
| `save_trajectory(filename)` | — | Write all images to ASE `.traj` |

---

### `linear_interpolate(start, end, n_images)`

Cartesian linear interpolation. MIC-aware for periodic systems.
Returns `[start_copy, img_1, ..., img_n, end_copy]` — no calculators
attached to intermediate images.

### `idpp_interpolate(start, end, n_images, max_iter=500, tol=1e-6)`

IDPP interpolation. MIC-aware for periodic systems.
Returns `[start_copy, img_1, ..., img_n, end_copy]` — no calculators
attached to intermediate images. Recommended over `linear_interpolate`
for torsional reactions and any path with risk of atomic overlap.

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

Current suite: **68 tests** across forces, interpolation, MIC, variable
springs, and parallel evaluation.

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
| Morse H3 | Morse (analytical) | 0.200 eV | 0.193 eV (exact) | 4% |
| Al adatom diffusion | EMT | 0.237 eV | ~0.40 eV (DFT-PBE) | finite-size slab |
| Ethane C–C torsion | Egret-1t | 0.113 eV | 0.126 eV (exp.) | 10% |
| Al vacancy migration | MACE-MP-0 | 0.508 eV | 0.61 eV (DFT-PBE) | 17%† |
| Mg vacancy (HCP) | MACE-MP-0 | 0.508 eV | ~0.52 eV (DFT-PBE) | 2% |
| Cu vacancy | EMT | 0.755 eV | ~0.70 eV (DFT-PBE) | 7.9% |
| Ni vacancy | EMT | 1.095 eV | ~1.04 eV (DFT-PBE) | 5.3% |
| Pd vacancy | EMT | 0.839 eV | ~0.91 eV (DFT-PBE) | 7.8% |
| Ag vacancy | EMT | 0.682 eV | ~0.66 eV (DFT-PBE) | 3.3% |
| Pt vacancy | EMT | 0.971 eV | ~1.49 eV (DFT-PBE) | 34.8%‡ |

† MACE-MP-0 systematically underestimates vacancy migration barriers by
10–20%. This is a known model limitation, not a nebwalk bug.

‡ EMT does not capture relativistic effects significant in Pt. The NEB
itself converged cleanly in 60 steps — the error is from the calculator.

---

## References

1. H. Jónsson, G. Mills, K.W. Jacobsen, *Nudged Elastic Band Method for
   Finding Minimum Energy Paths of Transitions*, in *Classical and Quantum
   Dynamics in Condensed Phase Simulations*, World Scientific, 1998.
2. G. Henkelman, H. Jónsson, *Improved tangent estimate in the nudged
   elastic band method for finding minimum energy paths and saddle points*,
   J. Chem. Phys. **113**, 9978 (2000).
   DOI: [10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
3. G. Henkelman, B.P. Uberuaga, H. Jónsson, *A climbing image nudged
   elastic band method for finding saddle points and minimum energy paths*,
   J. Chem. Phys. **113**, 9901 (2000).
   DOI: [10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
4. E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, P. Gumbsch, *Structural
   Relaxation Made Simple*, PRL **97**, 170201 (2006).
   DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)
5. S. Smidstrup, A. Pedersen, K. Stokbro, H. Jónsson, *Improved initial
   guess for minimum energy path calculations*, J. Chem. Phys. **141**,
   214106 (2014).
   DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)
6. R. Lindh, A. Bernhardsson, G. Karlström, P.-Å. Malmqvist, *On the use of
   a Hessian model function in molecular geometry optimizations*,
   Chem. Phys. Lett. **241**, 423 (1995). (Variable spring constant scheme.)

---

## License

MIT © Md. Rifat Khandaker
