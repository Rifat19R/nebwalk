# H diffusion on Cu(111) with MACE-MP-0

This benchmark demonstrates a surface-diffusion workflow using `nebwalk`
CI-NEB and MACE-MP-0.

The benchmark contains two measured paths:

1. **Elementary fcc → hcp hop**
   - one dominant bridge-like transition state;
   - best README-facing benchmark because it is clean and easy to interpret.

2. **Full fcc → fcc hop**
   - two equivalent bridge-like transition states;
   - shallow hcp intermediate;
   - useful symmetry check.

## System

| Parameter | Value |
|---|---:|
| Surface | Cu(111) |
| Slab | 4×4×4 |
| Cu lattice parameter | 3.6147 Å |
| Vacuum | 15 Å |
| H coverage | 1/16 ML |
| Fixed atoms | bottom two Cu layers |
| Endpoint relaxation fmax | 0.015 eV/Å |
| NEB fmax | 0.03 eV/Å |
| Interpolation | IDPP |
| NEB mode | CI-NEB |
| Calculator | MACE-MP-0 small |
| Precision | float64 |

## Results

### Elementary fcc → hcp

| Quantity | Result |
|---|---:|
| Converged | True |
| Wall time | 2.67 min |
| Forward barrier | 125.92 meV |
| Reverse barrier | 136.67 meV |
| Reaction energy | −10.75 meV |
| TS image | 3 |

Profile:

```text
image 00:   +0.000 meV
image 01:  +44.614 meV
image 02: +111.234 meV
image 03: +125.917 meV  <-- TS
image 04: +104.770 meV
image 05:  +33.294 meV
image 06:  -10.749 meV
```

### Full fcc → fcc

| Quantity | Result |
|---|---:|
| Converged | True |
| Wall time | 6.42 min |
| Forward barrier | 125.56 meV |
| Reverse barrier | 125.41 meV |
| Reaction energy | +0.146 meV |

Profile:

```text
image 00:   +0.000 meV
image 01:  +75.999 meV
image 02: +125.559 meV  <-- TS
image 03:  +75.126 meV
image 04:   -9.635 meV
image 05:  +86.956 meV
image 06: +125.405 meV
image 07:  +69.358 meV
image 08:   +0.146 meV
```

The fcc → fcc path is interpreted as:

```text
fcc → bridge TS → hcp hollow → bridge TS → fcc
```

## Scientific interpretation

This benchmark validates the `nebwalk` workflow rather than claiming exact
experimental reproduction. The path is converged, smooth, and internally
consistent. Endpoint symmetry in the fcc→fcc case is excellent: the reaction
energy is only 0.146 meV and the forward/reverse barriers differ by only
0.15 meV.

The measured barrier is lower than macroscopic experimental H/Cu(111) diffusion
activation energies because this is a static 0 K flat-terrace MACE/PBE-like
model. Experiments include thermal activation, quantum effects, coverage,
defects, and surface morphology. MACE-MP-0 is a foundation model trained on
DFT-like data and should be treated as a fast prototyping calculator unless
system-specific validation is performed.

## Files

```text
results/cu111_fcc_to_hcp_profile.csv
results/cu111_fcc_to_hcp_metadata.json
results/cu111_fcc_to_hcp_profile.png
results/cu111_fcc_to_fcc_profile.csv
results/cu111_fcc_to_fcc_metadata.json
results/cu111_fcc_to_fcc_profile.png
scripts/run_h_cu111_fcc_to_hcp.sh
```

## Re-running

Install MACE support first:

```bash
pip install "nebwalk[mace]"
```

Then run:

```bash
bash benchmarks/h_cu111_mace_mp/scripts/run_h_cu111_fcc_to_hcp.sh
```

The run script is intentionally written in terminal-heredoc style so each step
can also be copied and executed manually.
