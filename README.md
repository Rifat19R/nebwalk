# nebwalk

Minimal, correct Python implementation of the **Nudged Elastic Band (NEB)** method for finding minimum energy paths (MEPs) and transition states. Works with any ASE-compatible calculator — EMT for testing, [Egret-1](https://rowansci.com) for production.

**GitHub:** [github.com/Rifat19R/nebwalk](https://github.com/Rifat19R/nebwalk)

---

## Features

- Improved tangent estimate (Henkelman & Jónsson 2000)
- Spring forces + perpendicular potential force projection
- Climbing Image NEB (CI-NEB) for true saddle-point location
- **FIRE optimizer** — robust to the non-conservative NEB force field (unlike L-BFGS-B)
- Minimum Image Convention (MIC) for periodic systems
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
from nebwalk import NEB, linear_interpolate

# 1. Define endpoints (relaxed)
initial = ...   # ase.Atoms
final   = ...   # ase.Atoms

# 2. Build interpolated path
images = linear_interpolate(initial, final, n_images=7)
for img in images:
    img.calc = EMT()

# 3. Run NEB
neb = NEB(images, k=0.1, climb=True)
neb.optimize(fmax=0.05)

# 4. Output
print(f"Barrier: {neb.get_barrier():.3f} eV")
neb.plot("profile.png")
neb.save_csv("profile.csv")
neb.save_trajectory("path.traj")
```

Visualise the trajectory: `ase gui path.traj`

---

## With Egret-1 (MACE)

```python
from mace.calculators import MACECalculator

# Download EGRET_1T.model from https://rowansci.com — do not commit to repo
images = linear_interpolate(initial, final, n_images=9)
for img in images:
    img.calc = MACECalculator(
        model_paths="EGRET_1T.model",
        device="cpu",
        default_dtype="float32",
    )

neb = NEB(images, k=0.1, climb=True)
neb.optimize(fmax=0.05)
```

> **Note:** Egret-1 does not compute stress/virial; only atomic forces are needed for NEB, so this is fine.  
> The model file is not distributed with this repository. Download it from [rowansci.com](https://rowansci.com).

---

## Theory

### Tangent estimate (improved, Henkelman & Jónsson 2000)

Let τ⁺ = R_{i+1} − Rᵢ and τ⁻ = Rᵢ − R_{i-1}.

| Condition | Tangent |
|-----------|---------|
| E_{i+1} > Eᵢ > E_{i-1} | τ⁺ / \|τ⁺\| |
| E_{i-1} > Eᵢ > E_{i+1} | τ⁻ / \|τ⁻\| |
| Local extremum (E_{i+1} ≥ E_{i-1}) | (τ⁺ ΔE_max + τ⁻ ΔE_min) normalised |
| Local extremum (E_{i-1} > E_{i+1}) | (τ⁺ ΔE_min + τ⁻ ΔE_max) normalised |

where ΔE_max = max(\|E_{i+1}−Eᵢ\|, \|E_{i-1}−Eᵢ\|) and ΔE_min = min(…).

### NEB force

```
F_spring = k * (|R_{i+1} - R_i| - |R_i - R_{i-1}|) * τ
F_perp   = F_pot - (F_pot · τ) τ          # F_pot = ASE get_forces() = -∇E
F_NEB    = F_spring + F_perp
```

### Climbing image

For the highest-energy movable image (after `climb_delay` steps):

```
F_CI = F_pot - 2 (F_pot · τ) τ
```

The spring is removed and the tangent component is inverted, driving the image to the true saddle point.

### Convergence criterion

```
max |F_NEB|_component < fmax   (eV/Å, default 0.05)
```

### Optimizer: FIRE

FIRE (Fast Inertial Relaxation Engine) is used instead of scipy L-BFGS-B because NEB forces are **not the gradient of any single energy function** — the spring projection and perpendicular projection break the conservative-field assumption required for L-BFGS-B line search. FIRE's velocity-damping approach is robust to this.

Reference: Bitzek et al., PRL 97, 170201 (2006).

---

## API

### `NEB(images, k=0.1, climb=False, climb_delay=100)`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `images` | — | List of `ase.Atoms`, all with calculators |
| `k` | 0.1 | Spring constant (eV/Å²). Do **not** scale by N_images. |
| `climb` | False | Enable CI-NEB |
| `climb_delay` | 100 | Steps before CI activates |

**Methods:**
- `optimize(fmax=0.05, max_steps=500, verbose=True)` → `bool`
- `get_energies()` → `list[float]`
- `get_barrier()` → `float` (eV, relative to image 0)
- `plot(filename, show, title)`
- `save_csv(filename)`
- `save_trajectory(filename)`

### `linear_interpolate(start, end, n_images)`

Returns `[start_copy, img_1, ..., img_n, end_copy]`.  
No calculators are attached to intermediate images — attach one before calling `NEB`.

### `compute_neb_forces(images, k, climb, climb_index)`

Low-level force function. Returns a list of `(N_atoms, 3)` arrays.  
Useful for custom optimisation loops.

---

## Running tests

```bash
pip install -e ".[test]"
pytest tests/ -v
```

---

## Examples

```bash
python examples/morse_h3.py          # collinear H+H2, Morse potential (self-contained)
python examples/al_diffusion_emt.py  # Al adatom diffusion on Al(100), EMT
python examples/ethane_egret.py      # ethane torsion barrier, Egret-1t
python examples/al_vacancy_macemp.py # Al vacancy migration in bulk Al, MACE-MP-0
python examples/verify_egret.py      # sanity check: confirm Egret-1t loads correctly
```

### Validated results

| Example | Calculator | Barrier | Reference | Error |
|---------|-----------|---------|-----------|-------|
| Morse H3 | Morse (analytical) | 0.200 eV | 0.193 eV (analytical) | 4% |
| Al adatom diffusion | EMT | 0.237 eV | ~0.40 eV (DFT) | finite-size slab |
| Ethane C–C torsion | Egret-1t | 0.113 eV | 0.126 eV (exp.) | 10% |
| Al vacancy migration | MACE-MP-0 | 0.508 eV | 0.61 eV (DFT-PBE) | 17%* |

*MACE-MP-0 small systematically underestimates vacancy migration barriers by 10–20%. This is a known model limitation, not a nebwalk bug.

---

## References

1. H. Jónsson, G. Mills, K.W. Jacobsen, *Nudged Elastic Band Method for Finding Minimum Energy Paths of Transitions*, in *Classical and Quantum Dynamics in Condensed Phase Simulations*, World Scientific, 1998.
2. G. Henkelman, H. Jónsson, *Improved tangent estimate in the nudged elastic band method for finding minimum energy paths and saddle points*, J. Chem. Phys. **113**, 9978 (2000). DOI: [10.1063/1.1323224](https://doi.org/10.1063/1.1323224)
3. G. Henkelman, B.P. Uberuaga, H. Jónsson, *A climbing image nudged elastic band method for finding saddle points and minimum energy paths*, J. Chem. Phys. **113**, 9901 (2000). DOI: [10.1063/1.1329672](https://doi.org/10.1063/1.1329672)
4. E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, P. Gumbsch, *Structural Relaxation Made Simple*, PRL **97**, 170201 (2006). DOI: [10.1103/PhysRevLett.97.170201](https://doi.org/10.1103/PhysRevLett.97.170201)

---

## License

MIT © Md. Rifat Khandaker
