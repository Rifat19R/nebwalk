#!/usr/bin/env bash
set -euo pipefail

cd /mnt/d/Rifat_kh/nebwalk_universal

if [[ -d venv ]]; then
  # shellcheck disable=SC1091
  source venv/bin/activate
fi

MODE="${1:-all}"

run() {
  echo
  echo ">>> python $1"
  python "$1"
}

run_safe() {
  run examples/active_selection_logic.py
  run examples/active_al_vacancy_emt.py
  run examples/au_vacancy_emt.py
  run examples/au_vacancy_macemp.py
}

run_emt() {
  run examples/al_vacancy_emt.py
  run examples/cu_vacancy_emt.py
  run examples/ag_vacancy_emt.py
  run examples/ni_vacancy_emt.py
  run examples/pd_vacancy_emt.py
  run examples/au_vacancy_emt.py
}

run_mace() {
  run examples/al_vacancy_mace.py
  run examples/cu_vacancy_mace.py
  run examples/ag_vacancy_mace.py
  run examples/ni_vacancy_mace.py
  run examples/pd_vacancy_mace.py
  run examples/au_vacancy_macemp.py
}

run_dry() {
  for item in \
    al:emt cu:emt ag:emt ni:emt pd:emt au:emt \
    al:mace cu:mace ag:mace ni:mace pd:mace au:mace \
    al:qe cu:qe ag:qe w:qe mo:qe si:qe
  do
    material="${item%%:*}"
    backend="${item##*:}"
    echo
    echo ">>> python examples/vacancy_benchmark_suite.py ${material} ${backend} --dry-run"
    python examples/vacancy_benchmark_suite.py "${material}" "${backend}" --dry-run
  done
}

run_qe() {
  export NEBWALK_RUN_QE=1
  export ESPRESSO_PSEUDO=/mnt/d/Rifat_kh/SSSP_1.3.0_PBE_efficiency
  export ESPRESSO_COMMAND="mpirun --oversubscribe -np 4 pw.x"
  export NEBWALK_QE_CLEAN=1

  export AL_PSEUDO=Al.pbe-n-kjpaw_psl.1.0.0.UPF
  export CU_PSEUDO=Cu.paw.z_11.ld1.psl.v1.0.0-low.upf
  export AG_PSEUDO=Ag_ONCV_PBE-1.0.oncvpsp.upf
  export W_PSEUDO=W_pbe_v1.2.uspp.F.UPF
  export MO_PSEUDO=Mo_ONCV_PBE-1.0.oncvpsp.upf
  export SI_PSEUDO=Si.pbe-n-rrkjus_psl.1.0.0.UPF

  run examples/al_vacancy_qe_suite.py
  run examples/cu_vacancy_qe.py
  run examples/ag_vacancy_qe.py
  run examples/w_vacancy_qe.py
  run examples/mo_vacancy_qe.py
  run examples/si_vacancy_qe.py
}

run_qe_targets() {
  export NEBWALK_RUN_QE=1
  export ESPRESSO_PSEUDO=/mnt/d/Rifat_kh/SSSP_1.3.0_PBE_efficiency
  export ESPRESSO_COMMAND="mpirun --oversubscribe -np 4 pw.x"
  export NEBWALK_QE_CLEAN=1

  export W_PSEUDO=W_pbe_v1.2.uspp.F.UPF
  export MO_PSEUDO=Mo_ONCV_PBE-1.0.oncvpsp.upf
  export SI_PSEUDO=Si.pbe-n-rrkjus_psl.1.0.0.UPF

  run examples/w_vacancy_qe.py
  run examples/mo_vacancy_qe.py
  run examples/si_vacancy_qe.py
}

run_all_local() {
  run_safe
  run examples/al_vacancy_emt.py
  run examples/cu_vacancy_emt.py
  run examples/ag_vacancy_emt.py
  run examples/ni_vacancy_emt.py
  run examples/pd_vacancy_emt.py
  run examples/al_vacancy_mace.py
  run examples/cu_vacancy_mace.py
  run examples/ag_vacancy_mace.py
  run examples/ni_vacancy_mace.py
  run examples/pd_vacancy_mace.py
}

case "${MODE}" in
  safe)
    run_safe
    ;;
  dry)
    run_dry
    ;;
  emt)
    run_emt
    ;;
  mace)
    run_mace
    ;;
  qe)
    run_qe
    ;;
  qe-targets)
    run_qe_targets
    ;;
  all)
    run_all_local
    ;;
  full)
    run_all_local
    run_qe
    ;;
  *)
    echo "Usage: bash run.sh [all|safe|dry|emt|mace|qe|qe-targets|full]" >&2
    echo "  all : local validation only (safe + EMT + MACE), no QE" >&2
    echo "  qe-targets: real QE for W/Mo/Si vacancy tests" >&2
    echo "  full: local validation plus real QE calculations" >&2
    exit 2
    ;;
esac
