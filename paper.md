---
title: 'nebwalk: a calculator-agnostic nudged elastic band library with MLIP-assisted barrier pre-screening'
tags:
  - Python
  - computational chemistry
  - nudged elastic band
  - transition state theory
  - machine learning interatomic potentials
  - MACE
  - Quantum ESPRESSO
  - materials simulation
  - active learning
authors:
  - name: Md. Rifat Khandaker
    orcid: 0000-0003-0520-4654
    affiliation: 1
affiliations:
  - name: Department of Chemical Engineering, Dhaka University of Engineering and Technology, Gazipur 1707, Bangladesh
    index: 1
date: 13 June 2026
bibliography: paper.bib
---

# Summary

nebwalk is a Python library for nudged elastic band (NEB) calculations. It
takes two endpoint structures — a reactant and a product — builds an
interpolated path between them, and optimizes the path to locate the
transition state and report the activation barrier.

The library implements the improved tangent NEB and climbing-image NEB
(CI-NEB) of @henkelman2000a and @henkelman2000b, the FIRE optimizer of
@bitzek2006, and IDPP path interpolation of @smidstrup2014. Spring constants
can be uniform or varied along the path to concentrate images near the
transition state. Periodic systems are handled through ASE's minimum-image
convention [@larsen2017].

nebwalk supports four calculator backends without any code change: the EMT
classical calculator for fast local testing, machine-learned interatomic
potentials (MACE-MP-0 [@batatia2024], Egret-1t) for affordable pre-screening,
and Quantum ESPRESSO [@giannozzi2009] for DFT-level results. The same
configuration object and the same function call work with all of them.

A key feature is the MLIP-assisted two-stage workflow. This runs a fast
MACE or Egret NEB first, identifies the images closest to the transition
state using an energy-based selection strategy, and exports those images
for DFT single-point validation or active-learning labeling. This reduces
the number of expensive DFT calculations needed to characterize a barrier
from a full NEB run to three to five targeted image evaluations.

A reproducibility module records a self-contained bundle for every completed
calculation: input structures, the full configuration, results, trajectory,
software environment, SHA-256 checksums of all files, and a rerun script
that reproduces the original calculation from scratch.

# Statement of need

Finding transition states is a routine task in computational materials
research: vacancy migration in metals and oxides, surface diffusion in
heterogeneous catalysis, conformational changes in molecular systems, and
elementary reaction steps in electrocatalysis. The NEB method [@jonsson1998]
is the standard approach, but most implementations are tied to one specific
code. QE's neb.x requires pw.x. VASP's NEB implementation requires VASP.
Neither can accept a machine-learned potential without significant wrapping.

ASE [@larsen2017] provides a built-in NEB module, but its default optimizer
(L-BFGS) is technically inappropriate for the NEB force field, because NEB
forces are not the gradient of any single potential energy function
[@bitzek2006]. Beyond the optimizer, ASE's NEB is not structured for
switching calculators mid-workflow.

The rise of universal machine-learned potentials such as MACE-MP-0
[@batatia2024] has made MLIP-level NEB practical. A single MACE NEB run
on a 50-atom supercell takes minutes on a CPU, compared to hours for a
DFT NEB with the same system. However, MLIP barriers carry systematic
errors relative to DFT that depend on the system and calculator. The
standard response — run the full NEB at the DFT level — is expensive and
often unnecessary. What researchers need is a way to use the MLIP result
to identify which images to evaluate with DFT, rather than evaluating all
of them.

nebwalk addresses this directly. Its MLIP-assisted workflow runs the cheap
NEB first, selects the three to five images that carry the most information
about the barrier shape, and exports them in a format ready for DFT
single-point calculations or inclusion in an active-learning loop. This
is useful for any researcher who wants to combine the speed of MLIPs with
the accuracy of DFT without writing their own orchestration code.

The library is also small enough to inspect and modify. The NEB force
equations in the source code map directly to the published expressions in
@henkelman2000a, which makes it useful as a reference implementation for
teaching and methods development.

# Implementation

nebwalk is implemented in Python and depends on ASE [@larsen2017], NumPy,
SciPy, and Matplotlib. MACE and Egret are optional at runtime. Quantum
ESPRESSO support goes through ASE's Espresso interface and requires a
working pw.x binary and UPF pseudopotential files.

The NEB force routine applies the improved tangent estimate of @henkelman2000a,
computes the perpendicular potential force, and adds the projected spring
force. The climbing-image modification [@henkelman2000b] activates at the
highest-energy image after a configurable delay. The FIRE optimizer
[@bitzek2006] manages a shared velocity array across all intermediate images,
which is the correct approach for the non-conservative NEB force field. For
periodic systems, all inter-image displacements are computed using ASE's
`find_mic` to handle atoms that cross cell boundaries.

The MLIP-assisted workflow is implemented in `nebwalk.active`. The function
`run_mlip_assisted_neb()` accepts any ASE-compatible calculator factory,
runs the NEB, and then calls the image selection module. The default
strategy, `peak_plus_neighbors`, selects the highest-energy intermediate
image and its immediate path neighbors. Selected images are exported as
`.xyz`, `.traj`, and `.json` files with full metadata. The selection
strategy is a named, swappable component so that users can implement and
register alternative strategies without modifying the core library.

The reproducibility module in `nebwalk.reproduce` provides `save_bundle()`,
which takes the result of a completed NEB run and writes a directory
containing: the initial and final structures in extended XYZ format, the
full NEB configuration as a JSON file, the results and per-step convergence
history, the complete path trajectory, a capture of the Python environment
via `pip freeze`, SHA-256 checksums of all files, and a rerun script that
reads the stored configuration and re-executes the calculation. The
calculator itself is not serialized — the rerun script contains a marked
placeholder where the user fills in their calculator factory. This is an
intentional design choice: arbitrary Python callables cannot be reliably
serialized, and storing a parameter dictionary alongside a clear placeholder
is more honest and more portable than storing source code.

The package includes 132 tests covering interpolation, force projection,
tangent construction, minimum-image handling, variable springs, the
MLIP-assisted workflow, reproducibility bundles, and the QE interface.
It is distributed on PyPI and tested continuously via GitHub Actions.

# Validation

Nineteen validated system-calculator combinations across four backends have been checked against
DFT-PBE or experimental reference values. A representative selection is
shown below. Where errors are large — Pt vacancy with EMT (34.8%), Al
vacancy with MACE-MP-0 (17%) — these reflect documented limitations of the
calculator, not of the NEB optimization. In every case the NEB converged
cleanly and the endpoint energy difference for symmetric hops was below
1 meV. Full results, per-system notes, and calculator-specific caveats
are in the repository README.

| System | Calculator | Barrier | Reference | Error |
|--------|-----------|---------|-----------|-------|
| Al vacancy / FCC Al | QE/PBE | 0.564 eV | ~0.61 eV | 7.6% |
| W vacancy / BCC W | QE/PBE | 1.561 eV | ~1.66 eV | 5.9% |
| Mo vacancy / BCC Mo | QE/PBE | 1.281 eV | ~1.35 eV | 5.1% |
| Li vacancy / Li₂O | MACE-MP-0 | 0.284 eV | ~0.28 eV | 1.4% |
| Mg vacancy / MgO | MACE-MP-0 | 2.254 eV | ~2.20 eV | 2.5% |
| Ethane C–C torsion | Egret-1t | 0.113 eV | 0.126 eV | 10% |

# Acknowledgements

The author thanks Prof. Dr. Mohammad Asaduzzaman Chowdhury (DUET) for
supervision and research support. Development was assisted by Claude
(Anthropic) and OpenAI Codex for code review, debugging, and documentation
guidance. All scientific decisions, algorithm choices, validation strategy,
and implementation responsibility belong to the author.

# References
