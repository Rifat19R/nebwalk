"""Si vacancy migration in diamond Si with MACE-MP-0.

Reference barrier: ~0.45 eV DFT-PBE order-of-magnitude. Calculator: MACE-MP-0.
"""

from vacancy_benchmark_suite import main

if __name__ == "__main__":
    main("si", "mace")
