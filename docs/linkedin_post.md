# LinkedIn Post Draft

I built `nebwalk`, a minimal Python Nudged Elastic Band implementation for
ASE-compatible calculators.

Best validation result so far: **HCP Mg vacancy migration with MACE-MP-0 gives
0.508 eV**, compared with **~0.52 eV from DFT-PBE** -- about **2% error**.
That is the strongest headline result because it checks a real HCP materials
diffusion pathway, not just a toy molecule.

Other validated cases include:

- Morse H3 exchange: 0.200 eV vs 0.193 eV analytical reference
- Ethane C-C torsion with Egret-1t: 0.113 eV vs 0.126 eV experiment
- Al vacancy with MACE-MP-0: 0.508 eV vs 0.61 eV DFT-PBE
- EMT vacancy migration across Cu, Ni, Pd, Ag, and Pt

The package now includes CI-NEB, improved tangents, FIRE optimization, MIC/PBC
support, IDPP and regularized geodesic-style interpolation, variable springs,
parallel image evaluation, restart from `.traj`, and a QE/ASE template for
DFT-level workflows.

GitHub: https://github.com/Rifat19R/nebwalk

#computationalchemistry #materialsscience #opensource #NEB #DFT
