"""Al vacancy migration in FCC Al with EMT.

Reference barrier: ~0.61 eV DFT-PBE. Calculator: EMT.
"""

from vacancy_benchmark_suite import main

if __name__ == "__main__":
    main("al", "emt")
