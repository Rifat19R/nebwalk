"""Ag vacancy migration in FCC Ag with MACE-MP-0.

Reference barrier: ~0.66 eV DFT-PBE. Calculator: MACE-MP-0.
"""

from vacancy_benchmark_suite import main

if __name__ == "__main__":
    main("ag", "mace")
