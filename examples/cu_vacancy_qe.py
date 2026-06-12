"""Cu vacancy migration QE/PBE setup or run.

Set NEBWALK_RUN_QE=1, ESPRESSO_COMMAND="mpirun --oversubscribe -np 4 pw.x",
ESPRESSO_PSEUDO, and CU_PSEUDO to run QE.
"""

from vacancy_benchmark_suite import main

if __name__ == "__main__":
    main("cu", "qe")
