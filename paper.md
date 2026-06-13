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

nebwalk is a Python library for nudged elastic band (NEB) calculations. The
user provides two relaxed structures, a reactant and a product. nebwalk builds
an interpolated path between them, optimizes that path, and reports the
minimum-energy path along with the forward and reverse activation barriers.

The library implements the improved tangent estimate [@henkelman2000a] and the
climbing-image variant [@henkelman2000b], the FIRE optimizer [@bitzek2006], and
IDPP path interpolation [@smidstrup2014]. Spring constants can be held constant
or varied along the path so that images cluster near the saddle point. For
periodic systems, all displacements between images use the minimum-image
convention through ASE [@larsen2017].

The same code runs with four kinds of calculator and no changes beyond the
calculator factory: EMT for fast local tests, the MACE-MP-0 [@batatia2024] and
Egret-1t machine-learned potentials for cheap screening, and Quantum ESPRESSO
[@giannozzi2009] for DFT. Each QE image runs in its own working directory so
that wavefunction files from different images do not overwrite each other.

nebwalk also includes a two-stage workflow for combining machine-learned
potentials with DFT. It runs the NEB once with a fast potential, picks out the
few images that matter most for the barrier, and writes them to disk ready for
DFT single-point checks or for an active-learning dataset. A separate module
saves a full record of any run, including structures, settings, results,
trajectory, file checksums, and a script to repeat the calculation.

# Statement of need

Transition states show up in a lot of materials work: vacancies moving through
a crystal, atoms hopping across a surface, molecules changing shape, and
elementary steps in catalysis. The NEB method [@jonsson1998] is the usual way
to find them, but most NEB code is locked to one engine. The neb.x tool in
Quantum ESPRESSO only drives pw.x. The NEB inside VASP only works with VASP.
Getting a machine-learned potential into either of these takes real effort.

ASE [@larsen2017] ships its own NEB, which is more flexible, but its default
optimizer is L-BFGS. That choice is not ideal here, because the NEB force is
not the gradient of a single energy function and quasi-Newton methods assume
that it is [@bitzek2006]. nebwalk uses FIRE everywhere instead, which does not
make that assumption.

The other gap is more recent. Universal machine-learned potentials like
MACE-MP-0 [@batatia2024] are now fast enough that running a full NEB on a
50-atom cell takes minutes on a laptop CPU, where the same calculation in DFT
takes hours. The catch is that the barrier from a universal potential carries a
systematic error against DFT, and how large that error is depends on the
chemistry. The common fix, rerunning the whole NEB in DFT, is expensive and
usually overkill. What is actually needed is a way to use the cheap result to
decide which images deserve a DFT calculation.

nebwalk does this. It runs the cheap NEB, selects the three to five images that
carry the most information about the barrier, and exports them in a form ready
for DFT or for active learning. Anyone who wants the speed of a machine-learned
potential together with the accuracy of DFT can use this without writing their
own glue code.

The library is small on purpose. The force routine maps line for line onto the
equations in @henkelman2000a, so it is also useful for teaching the method or
for testing changes to it.

# Implementation

nebwalk is pure Python and depends on ASE [@larsen2017], NumPy, SciPy, and
Matplotlib. MACE and Egret are optional and are only imported if the user's own
calculator factory imports them. Quantum ESPRESSO support goes through ASE's
Espresso interface and needs a working pw.x and UPF pseudopotential files.

The force routine computes the improved tangent [@henkelman2000a], projects out
the component of the potential force along the tangent, and adds the spring
force along it. The climbing image [@henkelman2000b] turns on at the
highest-energy image after a set number of steps. FIRE [@bitzek2006] carries a
single velocity array across all the moving images. For periodic cells, every
displacement between neighboring images is wrapped with ASE's `find_mic` so
that an atom crossing a cell boundary does not create a false jump in the path.

The two-stage workflow lives in `nebwalk.active`. The function
`run_mlip_assisted_neb()` takes any ASE calculator factory, runs the NEB, and
hands the energies to the selection module. The default rule,
`peak_plus_neighbors`, keeps the highest-energy image and the images on either
side of it. The selected images are written out as `.xyz`, `.traj`, and `.json`
with their metadata. The selection rule is a named component, so a different
rule can be added without touching the core code.

The reproducibility module is `nebwalk.reproduce`, with one main function,
`save_bundle()`. It writes the two endpoints in extended XYZ, the full
configuration as JSON, the results and the per-step convergence history, the
whole path trajectory, the output of `pip freeze`, SHA-256 checksums for every
file, and a rerun script. The calculator is deliberately not saved as code. A
Python calculator factory can be a lambda or a closure, and those cannot be
serialized in a reliable way. Instead the user passes a plain dictionary
describing the calculator, which is stored as JSON, and the rerun script leaves
a clearly marked spot for the user to put the factory back. This is more honest
than trying to capture arbitrary source code and more likely to actually work
on another machine.

The test suite has 132 tests. They cover interpolation, the force projection,
the tangent construction, minimum-image handling, variable springs, the
two-stage workflow, the reproducibility bundles, and the QE interface. The
package is on PyPI and is tested on every push through GitHub Actions.

# Validation

Nineteen system-calculator combinations across the four backends were checked
against DFT-PBE or experimental references. A few are shown below; the full set
with per-system notes is in the repository README. Where the error is large,
such as the Pt vacancy with EMT (34.8%) or the Al vacancy with MACE-MP-0 (17%),
the cause is the calculator, not the NEB. EMT has no relativistic terms for Pt,
and MACE-MP-0 is known to underbind vacancy migration barriers. In every case
the path converged and the two endpoints came out equal to within 1 meV for
symmetric hops, which is the internal check that the optimization itself is
sound.

| System | Calculator | Barrier | Reference | Error |
|--------|-----------|---------|-----------|-------|
| Al vacancy / FCC Al | QE/PBE | 0.564 eV | ~0.61 eV | 7.6% |
| W vacancy / BCC W | QE/PBE | 1.561 eV | ~1.66 eV | 5.9% |
| Mo vacancy / BCC Mo | QE/PBE | 1.281 eV | ~1.35 eV | 5.1% |
| Li vacancy / Li₂O | MACE-MP-0 | 0.284 eV | ~0.28 eV | 1.4% |
| Mg vacancy / MgO | MACE-MP-0 | 2.254 eV | ~2.20 eV | 2.5% |
| Ethane C–C torsion | Egret-1t | 0.113 eV | 0.126 eV | 10% |

The three QE results cover two crystal structures, FCC aluminium and BCC
tungsten and molybdenum, and stay within about 8% of the DFT references. The
remaining error is mostly finite-size: the supercells are small, and a larger
cell would bring the barriers closer to the reference values.

# Acknowledgements

I thank Prof. Dr. Mohammad Asaduzzaman Chowdhury (DUET) for supervision and
support. I used Claude (Anthropic) and OpenAI Codex during development for code
review, debugging, and help with the documentation. The scientific choices, the
algorithm and validation work, and responsibility for the code are mine.

# References
