# nebwalk

[!\[CI](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml/badge.svg)](https://github.com/Rifat19R/nebwalk/actions/workflows/ci.yml/badge.svg)
[!\[Python](https://img.shields.io/badge/python-3.9%20%7C%203.11%20%7C%203.12-blue)](https://img.shields.io/badge/python-3.9%20%7C%203.11%20%7C%203.12-blue)
[!\[License](https://img.shields.io/badge/license-MIT-green)](https://img.shields.io/badge/license-MIT-green)

Minimal, correct Python implementation of the **Nudged Elastic Band (NEB)** method for finding minimum energy paths (MEPs) and transition states. Works with any ASE-compatible calculator — EMT for testing, [Egret-1](https://rowansci.com) or MACE-MP-0 for production.

**GitHub:** [github.com/Rifat19R/nebwalk](https://github.com/Rifat19R/nebwalk)

\---

## Features

* Improved tangent estimate (Henkelman \& Jónsson 2000)
* Spring forces + perpendicular potential force projection
* Climbing Image NEB (CI-NEB) for true saddle-point location
* **FIRE optimizer** — robust to the non-conservative NEB force field
* **IDPP interpolation (Smidstrup et al. 2014) — chemically sensible initial paths, avoids atomic clashes in torsional reactions**
* **Geodesic interpolation — Riemannian-geometry path for large structural changes where IDPP may generate unphysical intermediates**
* **Minimum Image Convention (MIC)** — correct tangents and spring forces for periodic systems
* **Variable spring constants** — energy-weighted springs concentrate images near the saddle point
* **High-level API** (`run\_neb\_calculation`) — single function call with calculator factory pattern; handles interpolation, calculator attachment, and NEB in one step
* **Parallel image evaluation** — thread-based; beneficial for large supercells (>50 atoms) or GPU-accelerated MACE. Not recommended for small molecules on CPU where thread overhead dominates.
* Energy profile plot, CSV export, ASE `.traj` output
* Local editable install — no compiled extensions. PyPI release is pending.

\---

## Installation

```
# PyPI release pending: pip install nebwalk is not available yet.
pip install .
# or, for development:
pip install -e ".\[test]"
```

For Egret-1 / MACE support:

```
pip install -e ".\[mace]"
```

\---

## Quick start

```python
from ase.calculators.emt import EMT
from nebwalk import NEB, idpp\_interpolate

# 1. Define endpoints (relaxed structures)
initial = ...   # ase.Atoms
final   = ...   # ase.Atoms

# 2. Build interpolated path (IDPP recommended over linear for molecules)
images = idpp\_interpolate(initial, final, n\_images=7)
for img in images:
    img.calc = EMT()

# 3. Run NEB with variable spring constants
neb = NEB(images, k=0.1, k\_min=0.033, climb=True)
neb.optimize(fmax=0.05)

# 4. Output
print(f"Barrier: {neb.get\_barrier():.3f} eV")
neb.plot("profile.png")
neb.save\_csv("profile.csv")
neb.save\_trajectory("path.traj")
```

Visualise the trajectory: `ase gui path.traj`

\---

## High-level API (run\_neb\_calculation)

`run\_neb\_calculation` wraps interpolation, calculator attachment, and NEB
optimisation into a single function call. Particularly useful for periodic
systems where each image needs an independent calculator instance.

```python
from ase.calculators.emt import EMT
from nebwalk import run\_neb\_calculation, NEBRunConfig

initial = ...   # ase.Atoms, relaxed
final   = ...   # ase.Atoms, relaxed

config = NEBRunConfig(
    n\_images  = 7,
    climb     = True,
    k         = 0.1,
    k\_min     = 0.033,
    fmax      = 0.05,
)

result = run\_neb\_calculation(
    initial            = initial,
    final              = final,
    calculator\_factory = lambda: EMT(),   # called once per image
    config             = config,
)

print(f"Converged : {result.converged}")
print(f"Barrier   : {result.barrier:.3f} eV")
result.neb.plot("profile.png")
result.neb.save\_trajectory("path.traj")
```

The `calculator\_factory` must be a zero-argument callable returning a fresh,
independent calculator. Sharing a single calculator across images causes ASE
cache corruption.

\---

## With Egret-1t (MACE) + parallel evaluation

```python
from mace.calculators import MACECalculator
from nebwalk import NEB, idpp\_interpolate

# Download EGRET\_1T.model from https://rowansci.com — do not commit to repo
images = idpp\_interpolate(initial, final, n\_images=9)
for img in images:
    img.calc = MACECalculator(
        model\_paths="EGRET\_1T.model",
        device="cpu",
        default\_dtype="float32",
    )

# n\_workers=9: one thread per intermediate image
# Use for large periodic systems (>50 atoms) or GPU; for small molecules on CPU, n\_workers=1 is faster
neb = NEB(images, k=0.1, k\_min=0.033, climb=True, n\_workers=9)
neb.optimize(fmax=0.05)
```

> \*\*Note:\*\* Egret-1 does not compute stress/virial; only atomic forces are needed
> for NEB, so this is fine.
> The model file is not distributed with this repository. Download it from \[rowansci.com](https://rowansci.com).

\---

## Periodic systems (bulk diffusion, vacancy migration)

```python
from ase.build import bulk
from ase.calculators.emt import EMT
from nebwalk import NEB, linear\_interpolate

al = bulk("Al", "fcc", a=4.05, cubic=True).repeat(2)
del al\[0]   # create vacancy

start = al.copy(); start.calc = EMT()
end   = al.copy(); end.positions\[0] += \[2.025, 0, 0]; end.calc = EMT()

# linear\_interpolate is MIC-aware for periodic systems
# idpp\_interpolate also supports PBC
images = linear\_interpolate(start, end, n\_images=7)
for img in images:
    img.calc = EMT()

neb = NEB(images, k=0.1, k\_min=0.033, climb=True)
neb.optimize(fmax=0.05)
```

\---

## Theory

### Interpolation

**Linear interpolation** (`linear\_interpolate`) divides the displacement vector
from start to end into equal steps. Fast, but fails for torsional reactions —
atoms can pass through each other, producing unphysical initial paths.

**IDPP interpolation** (`idpp\_interpolate`) minimises a weighted pairwise
distance objective across images (Smidstrup et al. 2014):

```
S = Σ\_{i<j} w\_ij (d\_ij - d\_target\_ij)²     w\_ij = 1/d\_target\_ij^4
```

This preserves short-range bonding structure throughout the path. Recommended
for any system with significant torsional motion or risk of atomic overlap.
Both methods support PBC via the minimum image convention.

Reference: Smidstrup, Pedersen, Stokbro, Jónsson, J. Chem. Phys. **141**,
214106 (2014). DOI: [10.1063/1.4878664](https://doi.org/10.1063/1.4878664)

\---

### Tangent estimate (improved, Henkelman \& Jónsson 2000)

Let τ⁺ = R\_{i+1} − Rᵢ and τ⁻ = Rᵢ − R\_{i-1}.

|Condition|Tangent|
|-|-|
|E\_{i+1} > Eᵢ > E\_{i-1}|τ⁺ / \|τ⁺\||
|E\_{i-1} > Eᵢ > E\_{i+1}|τ⁻ / \|τ⁻\||
|Local extremum (E\_{i+1} ≥ E\_{i-1})|(τ⁺ ΔE\_max + τ⁻ ΔE\_min) normalised|
|Local extremum (E\_{i-1} > E\_{i+1})|(τ⁺ ΔE\_min + τ⁻ ΔE\_max) normalised|

where ΔE\_max = max(|E\_{i+1}−Eᵢ|, |E\_{i-1}−Eᵢ|) and ΔE\_min = min(…).
For periodic systems, τ⁺ and τ⁻ use MIC displacements.

\---

### NEB force

```
F\_spring = k \* (|R\_{i+1} - R\_i| - |R\_i - R\_{i-1}|) \* τ
F\_perp   = F\_pot - (F\_pot · τ) τ          # F\_pot = ASE get\_forces() = -∇E
F\_NEB    = F\_spring + F\_perp
```

With variable spring constants, forward and backward spring constants can
differ per image:

```
F\_spring = (k\_fwd \* |R\_{i+1} - R\_i| - k\_bwd \* |R\_i - R\_{i-1}|) \* τ
```

\---

### Variable spring constants

Energy-weighted springs (Lindh et al. 1996) assign k\_max to springs adjacent
to the saddle point and k\_min to springs far from it:

```
k\[i] = k\_max - (k\_max - k\_min) \* (E\_ref - E\_spring\[i]) / (E\_ref - E\_low)
```

where E\_spring\[i] = max(E\[i], E\[i+1]), E\_ref = max(all energies),
E\_low = min(endpoint energies). Values are clipped to \[k\_min, k\_max].

This concentrates images near the transition state, improving barrier
resolution without adding more images. Recommended ratio: k\_min = k / 3.

\---

### Climbing image

For the highest-energy movable image (after `climb\_delay` steps):

```
F\_CI = F\_pot - 2 (F\_pot · τ) τ
```

The spring is removed and the tangent component is inverted, driving the
image to the true saddle point.

\---

### Convergence criterion

```
max per-atom |F\_NEB| < fmax   (eV/Å, default 0.05)
```

\---

### Optimizer: FIRE

FIRE (Fast Inertial Relaxation Engine) is used instead of L-BFGS-B because
NEB forces are **not the gradient of any single energy function** — the spring
projection and perpendicular projection break the conservative-field assumption
required for L-BFGS-B line search. FIRE's velocity-damping approach is
robust to this.

Reference: Bitzek et al., PRL **97**, 170201 (2006).
DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

\---

### Parallel image evaluation

At each FIRE step, energy and forces for all intermediate images are
independent — they do not communicate until spring forces are assembled.
`nebwalk` exploits this with `concurrent.futures.ThreadPoolExecutor`.

Threads are used (not processes) because pickling ASE calculator objects is
unreliable across calculator types.

Measured on a 4-core CPU (WSL Ubuntu, Egret-1t, C2H6, 7 images):

|Calculator|Sequential|Parallel (n\_workers=7)|Speedup|
|-|-|-|-|
|EMT|17.1 ms|6.5 ms|2.6×|
|Egret-1t|3.8 ms|8.3 ms|0.45×|

For small molecules on CPU, per-image compute (\~0.5 ms) is faster than thread
overhead — use `n\_workers=1`. Parallelism is beneficial when per-image compute
is large: heavy MACE models, large supercells (>50 atoms), or GPU evaluation.

\---

## API

### `NEB(images, k=0.1, k\_min=None, climb=False, climb\_delay=100, n\_workers=1)`

|Parameter|Default|Description|
|-|-|-|
|`images`|—|List of `ase.Atoms`, all with calculators attached|
|`k`|0.1|Spring constant (eV/Å²). With variable springs: maximum value.|
|`k\_min`|None|Minimum spring constant for variable springs. `None` = uniform springs. Recommended: k / 3.|
|`climb`|False|Enable CI-NEB|
|`climb\_delay`|100|FIRE steps before CI activates|
|`n\_workers`|1|Threads for parallel image evaluation. Beneficial for large supercells (>50 atoms) or GPU-accelerated MACE. Use default (1) for small molecules on CPU.|

**Methods:**

|Method|Returns|Description|
|-|-|-|
|`optimize(fmax=0.05, max\_steps=500, verbose=True)`|`bool`|Run FIRE. Returns True if converged.|
|`get\_energies()`|`list\[float]`|Potential energies for all images (eV)|
|`get\_barrier()`|`float`|Forward barrier relative to image 0 (eV)|
|`get\_spring\_constants()`|`ndarray`|Current spring constants, shape (N-1,)|
|`plot(filename, show, title)`|—|Save energy profile plot|
|`save\_csv(filename)`|—|Write energies to CSV|
|`save\_trajectory(filename)`|—|Write all images to ASE `.traj`|

\---

### `NEBRunConfig`

Configuration dataclass for `run\_neb\_calculation`. All parameters have
defaults; construct with only the values you want to override.

```python
from nebwalk import NEBRunConfig

config = NEBRunConfig(
    n\_images      = 7,
    interpolation = "idpp",
    k             = 0.1,
    k\_min         = None,
    climb         = False,
    climb\_delay   = 100,
    n\_workers     = 1,
    fmax          = 0.05,
    max\_steps     = 500,
    verbose       = True,
)
```

|Parameter|Type|Default|Description|
|-|-|-|-|
|`n\_images`|`int`|`7`|Number of intermediate images (endpoints excluded)|
|`interpolation`|`str`|`"idpp"`|`"linear"`, `"idpp"`, or `"geodesic"`. IDPP recommended for torsional reactions; geodesic for large structural changes|
|`k`|`float`|`0.1`|Spring constant (eV/Å²). With variable springs: maximum value.|
|`k\_min`|`float \| None`|`None`|Minimum spring constant for variable springs. `None` = uniform springs. Recommended: `k / 3`.|
|`climb`|`bool`|`False`|Enable CI-NEB for true saddle-point location|
|`climb\_delay`|`int`|`100`|FIRE steps before the climbing image activates|
|`n\_workers`|`int`|`1`|Threads for parallel image evaluation. Use `1` for small molecules on CPU.|
|`fmax`|`float`|`0.05`|Convergence criterion: max per-atom force magnitude (eV/Å)|
|`max\_steps`|`int`|`500`|Maximum FIRE steps before giving up|
|`verbose`|`bool`|`True`|Print per-step progress|

\---

### `run\_neb\_calculation(initial, final, calculator\_factory, config=None, prepare\_images=None)`

High-level NEB runner. Handles interpolation, calculator attachment, and
optimisation. Returns a `NEBRunResult`.

|Parameter|Type|Description|
|-|-|-|
|`initial`|`Atoms`|Relaxed initial state. No calculator required.|
|`final`|`Atoms`|Relaxed final state. No calculator required.|
|`calculator\_factory`|`Callable\[\[], Any]`|Zero-argument callable that returns a fresh calculator. Called once per image. **Must return independent instances** — sharing causes ASE cache corruption.|
|`config`|`NEBRunConfig \| None`|Configuration object. `None` = default `NEBRunConfig()`.|
|`prepare\_images`|`Callable\[\[list\[Atoms]], None] \| None`|Optional hook called on all images after interpolation, before NEB. Use to apply constraints or custom per-image setup.|

**Example with MACE-MP-0:**

```python
from mace.calculators import mace\_mp
from nebwalk import run\_neb\_calculation, NEBRunConfig

def make\_calc():
    return mace\_mp(model="small", dispersion=False, default\_dtype="float64")

config = NEBRunConfig(n\_images=7, climb=True, k=0.1, k\_min=0.033)

result = run\_neb\_calculation(
    initial            = initial,
    final              = final,
    calculator\_factory = make\_calc,
    config             = config,
)
```

\---

### `NEBRunResult`

Returned by `run\_neb\_calculation`.

|Attribute|Type|Description|
|-|-|-|
|`neb`|`NEB`|The NEB object after optimisation. Access images via `result.neb.images`, plot via `result.neb.plot(...)`.|
|`converged`|`bool`|`True` if the `fmax` criterion was met|
|`barrier`|`float`|Forward barrier relative to image 0 (eV)|
|`reverse\_barrier`|`float`|Reverse barrier relative to the last image (eV)|
|`reaction\_energy`|`float`|ΔE = E\_final − E\_initial (eV). Should be \~0 for equivalent sites.|

\---

### `linear\_interpolate(start, end, n\_images)`

Cartesian linear interpolation. MIC-aware for periodic systems.
Returns `\[start\_copy, img\_1, ..., img\_n, end\_copy]` — no calculators
attached to intermediate images.

### `idpp\_interpolate(start, end, n\_images, max\_iter=500, tol=1e-6)`

IDPP interpolation. MIC-aware for periodic systems.
Returns `\[start\_copy, img\_1, ..., img\_n, end\_copy]` — no calculators
attached to intermediate images. Recommended over `linear\_interpolate` for
torsional reactions and any path with risk of atomic overlap.

### `geodesic\_interpolate(start, end, n\_images)`

Geodesic interpolation on the Riemannian manifold of atomic configurations.

Handles large structural changes where IDPP may generate unphysical intermediates.

MIC-aware for periodic systems.

Returns `\[start\_copy, img\_1, ..., img\_n, end\_copy]` — no calculators

attached to intermediate images.`compute\_neb\_forces(images, k, climb, climb\_index, energies, forces)`

Low-level force function. Returns a list of `(N\_atoms, 3)` arrays. `energies`
and `forces` can be passed as pre-computed values to avoid redundant calculator
calls (used internally by the parallel evaluator).

\---

## Running tests

```
pip install -e ".\[test]"
pytest tests/ -v
```

Current suite: **94 tests** across forces, interpolation, MIC, variable
springs, parallel evaluation, restart helpers, the shared engine, and QE
calculator helpers.

\---

## Examples

```
python examples/morse\_h3.py                # collinear H+H2, Morse potential (self-contained)
python examples/al\_diffusion\_emt.py        # Al adatom diffusion on Al(100), EMT
python examples/cu\_adatom\_cu100\_emt.py     # Cu adatom diffusion on Cu(100), EMT
python examples/ni\_adatom\_ni100\_emt.py     # Ni adatom diffusion on Ni(100), EMT
python examples/ethane\_egret.py            # ethane C-C torsion barrier, Egret-1t
python examples/al\_vacancy\_macemp.py       # Al vacancy migration in bulk Al, MACE-MP-0
python examples/mg\_vacancy\_macemp.py       # Mg vacancy in HCP Mg, MACE-MP-0
python examples/li2o\_vacancy\_macemp.py     # Li vacancy migration in bulk Li2O, MACE-MP-0
python examples/mg\_vacancy\_mgo\_macemp.py   # Mg vacancy migration in bulk MgO, MACE-MP-0
python examples/cu\_vacancy\_emt.py          # Cu vacancy migration, EMT
python examples/ni\_vacancy\_emt.py          # Ni vacancy migration, EMT
python examples/pd\_vacancy\_emt.py          # Pd vacancy migration, EMT
python examples/ag\_vacancy\_emt.py          # Ag vacancy migration, EMT
python examples/pt\_vacancy\_emt.py          # Pt vacancy migration, EMT (see footnote)
python examples/verify\_egret.py            # sanity check: confirm Egret-1t loads correctly
```

### Validated results

|Example|Calculator|Barrier|Reference|Error|
|-|-|-|-|-|
|Morse H3|Morse (analytical)|0.200 eV|0.193 eV (exact)|4%|
|Al adatom diffusion|EMT|0.237 eV|\~0.40 eV (DFT-PBE)|finite-size slab|
|Cu adatom Cu(100)|EMT|0.418 eV|\~0.40 eV (DFT-LDA/exp.)|4.6%|
|Ni adatom Ni(100)|EMT|0.555 eV|\~0.63 eV (DFT-GGA)|12.0%§|
|Ethane C–C torsion|Egret-1t|0.113 eV|0.126 eV (exp.)|10%|
|Al vacancy migration|MACE-MP-0|0.508 eV|0.61 eV (DFT-PBE)|17%†|
|Mg vacancy (HCP Mg)|MACE-MP-0|0.508 eV|\~0.52 eV (DFT-PBE)|2%|
|Li vacancy (Li₂O)|MACE-MP-0|0.284 eV|\~0.28 eV (DFT-GGA)|1.4%|
|Mg vacancy (MgO)|MACE-MP-0|2.254 eV|\~2.20 eV (DFT-PBE)|2.5%|
|Cu vacancy|EMT|0.755 eV|\~0.70 eV (DFT-PBE)|7.9%|
|Ni vacancy|EMT|1.095 eV|\~1.04 eV (DFT-PBE)|5.3%|
|Pd vacancy|EMT|0.839 eV|\~0.91 eV (DFT-PBE)|7.8%|
|Ag vacancy|EMT|0.682 eV|\~0.66 eV (DFT-PBE)|3.3%|
|Pt vacancy|EMT|0.971 eV|\~1.49 eV (DFT-PBE)|34.8%‡|

† MACE-MP-0 systematically underestimates vacancy migration barriers by
10–20%. This is a known model limitation, not a nebwalk bug.

‡ EMT does not capture relativistic effects significant in Pt. The NEB
itself converged cleanly in 60 steps — the error is from the calculator.

§ EMT underestimates Ni adatom barriers due to the strong d-band character
at the Ni(100) saddle-point geometry, which effective-medium theory does not
resolve. Profile shape and convergence are correct.

\---

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

\---

## License

MIT © Md. Rifat Khandaker

