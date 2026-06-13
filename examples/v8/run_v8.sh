#!/usr/bin/env bash
set -uo pipefail

cd /mnt/d/Rifat_kh/nebwalk_universal

if [[ -d venv ]]; then
  # shellcheck disable=SC1091
  source venv/bin/activate
fi

scripts=(
  examples/v8_001_au100_adatom_emt.py
  examples/v8_002_al100_adatom_emt.py
  examples/v8_003_pd100_adatom_emt.py
  examples/v8_004_au_fcc_vacancy_emt.py
  examples/v8_005_fe_bcc_vacancy_emt.py
  examples/v8_006_cu_fcc_interstitial_emt.py
  examples/v8_007_ag111_adatom_emt.py
  examples/v8_008_al_fcc_divacancy_emt.py
  examples/v8_009_cao_ca_vacancy_mace.py
  examples/v8_010_sro_sr_vacancy_mace.py
  examples/v8_011_bao_ba_vacancy_mace.py
  examples/v8_012_na2o_na_vacancy_mace.py
  examples/v8_013_zno_zn_vacancy_mace.py
  examples/v8_014_tio2_o_vacancy_mace.py
  examples/v8_015_ceo2_ce_vacancy_mace.py
  examples/v8_016_lif_li_vacancy_mace.py
  examples/v8_017_tic_ti_vacancy_mace.py
  examples/v8_018_tin_n_vacancy_mace.py
  examples/v8_019_tic_c_vacancy_mace.py
  examples/v8_020_sic_si_vacancy_mace.py
  examples/v8_021_aln_al_vacancy_mace.py
  examples/v8_022_mo2c_mo_vacancy_mace.py
  examples/v8_023_licoo2_li_vacancy_mace.py
  examples/v8_024_li3po4_li_diffusion_mace.py
  examples/v8_025_nacoo2_na_diffusion_mace.py
  examples/v8_026_mg2si_mg_vacancy_mace.py
  examples/v8_027_h_pd_interstitial_mace.py
  examples/v8_028_h_ni_interstitial_mace.py
  examples/v8_029_n_fe100_surface_mace.py
  examples/v8_030_n_ru0001_surface_mace.py
  examples/v8_031_h_pt111_surface_mace.py
  examples/v8_032_o_cu111_surface_mace.py
  examples/v8_033_co_cu100_surface_mace.py
  examples/v8_034_n_mo110_surface_mace.py
  examples/v8_035_ethane_torsion_egret.py
  examples/v8_036_propane_torsion_egret.py
  examples/v8_037_dme_torsion_egret.py
  examples/v8_038_acetaldehyde_methyl_torsion_egret.py
)

failed=0
skipped=0
passed=0

for script in "${scripts[@]}"; do
  echo
  echo ">>> python ${script}"
  python "${script}"
  status=$?
  if [[ ${status} -eq 0 ]]; then
    passed=$((passed + 1))
  elif [[ ${status} -eq 77 ]]; then
    skipped=$((skipped + 1))
    echo ">>> skipped: ${script}"
  else
    failed=$((failed + 1))
    echo ">>> failed (${status}): ${script}" >&2
    break
  fi
done

echo
echo "v8 run summary: passed=${passed} skipped=${skipped} failed=${failed}"

if [[ ${failed} -ne 0 ]]; then
  exit 1
fi
