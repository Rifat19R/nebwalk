"""v0.8.0 reproducibility benchmark inputs.

These examples are intentionally isolated from the older benchmark scripts.
Each case writes to examples/v8_outputs/<case>/ and uses a unique
reproducibility bundle directory.

Reference barriers are literature/database guide values for benchmark
comparison, not fitted targets. Calculator error can be large for EMT and for
MACE-MP-0 outside its training domain.
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

import numpy as np
from ase import Atom, Atoms
from ase.build import bcc100, bcc110, bulk, fcc100, fcc111, hcp0001, molecule
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import BFGS
from ase.spacegroup import crystal

from nebwalk import NEBRunConfig, run_neb_calculation

CalcKind = Literal["emt", "mace", "egret"]
BuildKind = Literal[
    "surface_hop",
    "bulk_vacancy",
    "bulk_interstitial",
    "divacancy",
    "ionic_vacancy",
    "h_interstitial",
    "molecule_torsion",
]

ROOT = Path(__file__).resolve().parent
OUTPUT_ROOT = ROOT / "v8_outputs"
SKIP_EXIT = 77


@dataclass(frozen=True)
class BenchmarkSpec:
    key: str
    title: str
    calc: CalcKind
    build: BuildKind
    reference_barrier: float
    reference_source: str
    structure_source: str
    priority: str
    symbol: str = ""
    host: str = ""
    adsorbate: str = ""
    crystal_type: str = ""
    lattice_a: float = 0.0
    lattice_c: float | None = None
    repeat: tuple[int, int, int] = (2, 2, 2)
    surface: str = ""
    migrating_species: str = ""
    n_images: int = 5
    fmax: float = 0.05
    max_steps: int = 250
    climb_delay: int = 40
    k: float = 0.1
    k_min: float | None = None
    relax_fmax: float = 0.05
    relax_steps: int = 80
    notes: str = ""


SYSTEMS: dict[str, BenchmarkSpec] = {
    "v8_001_au100_adatom_emt": BenchmarkSpec(
        key="v8_001_au100_adatom_emt",
        title="Au adatom hop on Au(100)",
        calc="emt",
        build="surface_hop",
        symbol="Au",
        host="Au",
        adsorbate="Au",
        surface="fcc100",
        lattice_a=4.078,
        reference_barrier=0.08,
        reference_source="DFT surface-diffusion literature guide value",
        structure_source="Experimental FCC Au lattice, CRC/NIST scale",
        priority="high",
        n_images=5,
        k=0.3,
    ),
    "v8_002_al100_adatom_emt": BenchmarkSpec(
        key="v8_002_al100_adatom_emt",
        title="Al adatom hop on Al(100)",
        calc="emt",
        build="surface_hop",
        symbol="Al",
        host="Al",
        adsorbate="Al",
        surface="fcc100",
        lattice_a=4.05,
        reference_barrier=0.04,
        reference_source="DFT surface-diffusion literature guide value",
        structure_source="Experimental FCC Al lattice, CRC/NIST scale",
        priority="high",
        n_images=5,
        k=0.3,
    ),
    "v8_003_pd100_adatom_emt": BenchmarkSpec(
        key="v8_003_pd100_adatom_emt",
        title="Pd adatom hop on Pd(100)",
        calc="emt",
        build="surface_hop",
        symbol="Pd",
        host="Pd",
        adsorbate="Pd",
        surface="fcc100",
        lattice_a=3.89,
        reference_barrier=0.18,
        reference_source="DFT surface-diffusion literature guide value",
        structure_source="Experimental FCC Pd lattice, CRC/NIST scale",
        priority="high",
        n_images=5,
        k=0.3,
    ),
    "v8_004_au_fcc_vacancy_emt": BenchmarkSpec(
        key="v8_004_au_fcc_vacancy_emt",
        title="FCC Au vacancy migration",
        calc="emt",
        build="bulk_vacancy",
        symbol="Au",
        crystal_type="fcc",
        lattice_a=4.078,
        reference_barrier=0.68,
        reference_source="DFT-PBE vacancy-migration literature guide value",
        structure_source="Experimental FCC Au lattice, CRC/NIST scale",
        priority="medium",
        n_images=5,
    ),
    "v8_005_fe_bcc_vacancy_emt": BenchmarkSpec(
        key="v8_005_fe_bcc_vacancy_emt",
        title="BCC Fe vacancy migration",
        calc="emt",
        build="bulk_vacancy",
        symbol="Fe",
        crystal_type="bcc",
        lattice_a=2.866,
        reference_barrier=0.65,
        reference_source="DFT-GGA ferromagnetic BCC Fe vacancy literature",
        structure_source="Experimental BCC Fe lattice, CRC/NIST scale",
        priority="medium",
        notes="ASE EMT has no Fe parameterization; this case skips cleanly.",
    ),
    "v8_006_cu_fcc_interstitial_emt": BenchmarkSpec(
        key="v8_006_cu_fcc_interstitial_emt",
        title="FCC Cu self-interstitial hop",
        calc="emt",
        build="bulk_interstitial",
        symbol="Cu",
        crystal_type="fcc",
        lattice_a=3.615,
        reference_barrier=0.09,
        reference_source="DFT Cu self-interstitial guide value",
        structure_source="Experimental FCC Cu lattice, CRC/NIST scale",
        priority="medium",
        n_images=5,
    ),
    "v8_007_ag111_adatom_emt": BenchmarkSpec(
        key="v8_007_ag111_adatom_emt",
        title="Ag adatom hollow-to-hollow hop on Ag(111)",
        calc="emt",
        build="surface_hop",
        symbol="Ag",
        host="Ag",
        adsorbate="Ag",
        surface="fcc111",
        lattice_a=4.085,
        reference_barrier=0.06,
        reference_source="DFT Ag(111) adatom-diffusion literature guide value",
        structure_source="Experimental FCC Ag lattice, CRC/NIST scale",
        priority="low",
        n_images=5,
        k=0.3,
    ),
    "v8_008_al_fcc_divacancy_emt": BenchmarkSpec(
        key="v8_008_al_fcc_divacancy_emt",
        title="FCC Al divacancy second-neighbour migration",
        calc="emt",
        build="divacancy",
        symbol="Al",
        crystal_type="fcc",
        lattice_a=4.05,
        reference_barrier=0.54,
        reference_source="DFT Al divacancy-migration literature guide value",
        structure_source="Experimental FCC Al lattice, CRC/NIST scale",
        priority="low",
        n_images=5,
    ),
    "v8_009_cao_ca_vacancy_mace": BenchmarkSpec(
        key="v8_009_cao_ca_vacancy_mace",
        title="CaO Ca-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="CaO",
        migrating_species="Ca",
        crystal_type="rocksalt",
        lattice_a=4.81,
        reference_barrier=2.2,
        reference_source="DFT-GGA alkaline-earth oxide vacancy literature",
        structure_source="Materials Project / experimental rocksalt CaO",
        priority="high",
    ),
    "v8_010_sro_sr_vacancy_mace": BenchmarkSpec(
        key="v8_010_sro_sr_vacancy_mace",
        title="SrO Sr-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="SrO",
        migrating_species="Sr",
        crystal_type="rocksalt",
        lattice_a=5.16,
        reference_barrier=2.0,
        reference_source="DFT-GGA alkaline-earth oxide vacancy literature",
        structure_source="Materials Project / experimental rocksalt SrO",
        priority="medium",
    ),
    "v8_011_bao_ba_vacancy_mace": BenchmarkSpec(
        key="v8_011_bao_ba_vacancy_mace",
        title="BaO Ba-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="BaO",
        migrating_species="Ba",
        crystal_type="rocksalt",
        lattice_a=5.52,
        reference_barrier=1.8,
        reference_source="DFT-GGA alkaline-earth oxide vacancy literature",
        structure_source="Materials Project / experimental rocksalt BaO",
        priority="medium",
    ),
    "v8_012_na2o_na_vacancy_mace": BenchmarkSpec(
        key="v8_012_na2o_na_vacancy_mace",
        title="Na2O Na-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="Na2O",
        migrating_species="Na",
        crystal_type="antifluorite",
        lattice_a=5.55,
        reference_barrier=0.30,
        reference_source="DFT-GGA alkali oxide diffusion literature",
        structure_source="Materials Project / ICSD antifluorite Na2O",
        priority="high",
    ),
    "v8_013_zno_zn_vacancy_mace": BenchmarkSpec(
        key="v8_013_zno_zn_vacancy_mace",
        title="Wurtzite ZnO Zn-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="ZnO",
        migrating_species="Zn",
        crystal_type="wurtzite",
        lattice_a=3.25,
        lattice_c=5.207,
        reference_barrier=1.2,
        reference_source="DFT-PBE+U ZnO vacancy-migration literature",
        structure_source="Materials Project / experimental wurtzite ZnO",
        priority="medium",
    ),
    "v8_014_tio2_o_vacancy_mace": BenchmarkSpec(
        key="v8_014_tio2_o_vacancy_mace",
        title="Rutile TiO2 O-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="TiO2",
        migrating_species="O",
        crystal_type="rutile",
        lattice_a=4.594,
        lattice_c=2.959,
        reference_barrier=1.1,
        reference_source="DFT-PBE rutile TiO2 oxygen-vacancy literature",
        structure_source="Materials Project / ICSD rutile TiO2",
        priority="medium",
        notes="Ti3+ physics can challenge universal MLIPs.",
    ),
    "v8_015_ceo2_ce_vacancy_mace": BenchmarkSpec(
        key="v8_015_ceo2_ce_vacancy_mace",
        title="Fluorite CeO2 Ce-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="CeO2",
        migrating_species="Ce",
        crystal_type="fluorite",
        lattice_a=5.411,
        reference_barrier=0.9,
        reference_source="DFT-PBE+U ceria defect-migration literature",
        structure_source="Materials Project / ICSD fluorite CeO2",
        priority="medium",
        notes="Ce3+ localization can challenge universal MLIPs.",
    ),
    "v8_016_lif_li_vacancy_mace": BenchmarkSpec(
        key="v8_016_lif_li_vacancy_mace",
        title="LiF Li-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="LiF",
        migrating_species="Li",
        crystal_type="rocksalt",
        lattice_a=4.03,
        reference_barrier=0.45,
        reference_source="DFT-LDA LiF vacancy-diffusion literature",
        structure_source="Materials Project / experimental rocksalt LiF",
        priority="medium",
    ),
    "v8_017_tic_ti_vacancy_mace": BenchmarkSpec(
        key="v8_017_tic_ti_vacancy_mace",
        title="TiC Ti-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="TiC",
        migrating_species="Ti",
        crystal_type="rocksalt",
        lattice_a=4.33,
        reference_barrier=2.8,
        reference_source="DFT carbide vacancy-migration literature",
        structure_source="Materials Project / experimental rocksalt TiC",
        priority="medium",
    ),
    "v8_018_tin_n_vacancy_mace": BenchmarkSpec(
        key="v8_018_tin_n_vacancy_mace",
        title="TiN N-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="TiN",
        migrating_species="N",
        crystal_type="rocksalt",
        lattice_a=4.24,
        reference_barrier=1.9,
        reference_source="DFT TiN nitrogen-vacancy literature",
        structure_source="Materials Project / experimental rocksalt TiN",
        priority="high",
    ),
    "v8_019_tic_c_vacancy_mace": BenchmarkSpec(
        key="v8_019_tic_c_vacancy_mace",
        title="TiC C-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="TiC",
        migrating_species="C",
        crystal_type="rocksalt",
        lattice_a=4.33,
        reference_barrier=2.1,
        reference_source="DFT carbide vacancy-migration literature",
        structure_source="Materials Project / experimental rocksalt TiC",
        priority="medium",
    ),
    "v8_020_sic_si_vacancy_mace": BenchmarkSpec(
        key="v8_020_sic_si_vacancy_mace",
        title="3C-SiC Si-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="SiC",
        migrating_species="Si",
        crystal_type="zincblende",
        lattice_a=4.36,
        reference_barrier=2.4,
        reference_source="DFT SiC vacancy-migration literature",
        structure_source="Materials Project / experimental 3C-SiC",
        priority="medium",
    ),
    "v8_021_aln_al_vacancy_mace": BenchmarkSpec(
        key="v8_021_aln_al_vacancy_mace",
        title="Wurtzite AlN Al-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="AlN",
        migrating_species="Al",
        crystal_type="wurtzite",
        lattice_a=3.112,
        lattice_c=4.982,
        reference_barrier=2.6,
        reference_source="DFT nitride vacancy-migration literature",
        structure_source="Materials Project / experimental wurtzite AlN",
        priority="high",
    ),
    "v8_022_mo2c_mo_vacancy_mace": BenchmarkSpec(
        key="v8_022_mo2c_mo_vacancy_mace",
        title="Hexagonal Mo2C Mo-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="Mo2C",
        migrating_species="Mo",
        crystal_type="mo2c",
        lattice_a=3.01,
        lattice_c=4.73,
        reference_barrier=2.0,
        reference_source="DFT Mo2C vacancy-migration literature",
        structure_source="Materials Project / ICSD hexagonal beta-Mo2C",
        priority="medium",
    ),
    "v8_023_licoo2_li_vacancy_mace": BenchmarkSpec(
        key="v8_023_licoo2_li_vacancy_mace",
        title="Layered LiCoO2 Li-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="LiCoO2",
        migrating_species="Li",
        crystal_type="layered_r3m",
        lattice_a=2.815,
        lattice_c=14.05,
        reference_barrier=0.35,
        reference_source="Van der Ven and Ceder, PRB/ESL first-principles Li diffusion",
        structure_source="Materials Project / ICSD layered R-3m LiCoO2",
        priority="high",
    ),
    "v8_024_li3po4_li_diffusion_mace": BenchmarkSpec(
        key="v8_024_li3po4_li_diffusion_mace",
        title="Li3PO4 Li migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="Li3PO4",
        migrating_species="Li",
        crystal_type="li3po4",
        lattice_a=6.12,
        lattice_c=10.49,
        reference_barrier=0.55,
        reference_source="DFT-GGA Li3PO4 ion-conduction literature",
        structure_source="Materials Project / ICSD gamma-Li3PO4 proxy",
        priority="medium",
    ),
    "v8_025_nacoo2_na_diffusion_mace": BenchmarkSpec(
        key="v8_025_nacoo2_na_diffusion_mace",
        title="Layered NaCoO2 Na migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="NaCoO2",
        migrating_species="Na",
        crystal_type="layered_r3m",
        lattice_a=2.84,
        lattice_c=15.59,
        reference_barrier=0.45,
        reference_source="DFT-GGA layered NaCoO2 diffusion literature",
        structure_source="Materials Project / ICSD layered NaCoO2",
        priority="medium",
    ),
    "v8_026_mg2si_mg_vacancy_mace": BenchmarkSpec(
        key="v8_026_mg2si_mg_vacancy_mace",
        title="Mg2Si Mg-vacancy migration",
        calc="mace",
        build="ionic_vacancy",
        symbol="Mg2Si",
        migrating_species="Mg",
        crystal_type="antifluorite_mg2si",
        lattice_a=6.35,
        reference_barrier=0.80,
        reference_source="DFT-PBE Mg2Si vacancy-diffusion literature",
        structure_source="Materials Project / ICSD antifluorite Mg2Si",
        priority="medium",
    ),
    "v8_027_h_pd_interstitial_mace": BenchmarkSpec(
        key="v8_027_h_pd_interstitial_mace",
        title="H interstitial diffusion in FCC Pd",
        calc="mace",
        build="h_interstitial",
        symbol="Pd",
        host="Pd",
        crystal_type="fcc",
        lattice_a=3.89,
        reference_barrier=0.23,
        reference_source="DFT plus experimental H/Pd diffusion literature",
        structure_source="Experimental FCC Pd lattice, CRC/NIST scale",
        priority="high",
    ),
    "v8_028_h_ni_interstitial_mace": BenchmarkSpec(
        key="v8_028_h_ni_interstitial_mace",
        title="H interstitial diffusion in FCC Ni",
        calc="mace",
        build="h_interstitial",
        symbol="Ni",
        host="Ni",
        crystal_type="fcc",
        lattice_a=3.524,
        reference_barrier=0.41,
        reference_source="DFT plus experimental H/Ni diffusion literature",
        structure_source="Experimental FCC Ni lattice, CRC/NIST scale",
        priority="high",
    ),
    "v8_029_n_fe100_surface_mace": BenchmarkSpec(
        key="v8_029_n_fe100_surface_mace",
        title="N adatom hop on Fe(100)",
        calc="mace",
        build="surface_hop",
        symbol="Fe",
        host="Fe",
        adsorbate="N",
        surface="bcc100",
        lattice_a=2.866,
        reference_barrier=0.80,
        reference_source="Norskov et al., N/Fe(100) STM plus DFT diffusion study",
        structure_source="Experimental BCC Fe lattice, CRC/NIST scale",
        priority="high",
        k=0.3,
    ),
    "v8_030_n_ru0001_surface_mace": BenchmarkSpec(
        key="v8_030_n_ru0001_surface_mace",
        title="N adatom hop on Ru(0001)",
        calc="mace",
        build="surface_hop",
        symbol="Ru",
        host="Ru",
        adsorbate="N",
        surface="hcp0001",
        lattice_a=2.706,
        lattice_c=4.282,
        reference_barrier=1.3,
        reference_source="DFT Ru(0001) nitrogen-adatom diffusion literature",
        structure_source="Experimental HCP Ru lattice, CRC/NIST scale",
        priority="high",
        k=0.3,
    ),
    "v8_031_h_pt111_surface_mace": BenchmarkSpec(
        key="v8_031_h_pt111_surface_mace",
        title="H fcc-to-hcp hop on Pt(111)",
        calc="mace",
        build="surface_hop",
        symbol="Pt",
        host="Pt",
        adsorbate="H",
        surface="fcc111",
        lattice_a=3.924,
        reference_barrier=0.12,
        reference_source="DFT/experiment HER H/Pt(111) diffusion literature",
        structure_source="Experimental FCC Pt lattice, CRC/NIST scale",
        priority="medium",
        k=0.3,
    ),
    "v8_032_o_cu111_surface_mace": BenchmarkSpec(
        key="v8_032_o_cu111_surface_mace",
        title="O adatom hop on Cu(111)",
        calc="mace",
        build="surface_hop",
        symbol="Cu",
        host="Cu",
        adsorbate="O",
        surface="fcc111",
        lattice_a=3.615,
        reference_barrier=0.25,
        reference_source="DFT O/Cu(111) diffusion literature",
        structure_source="Experimental FCC Cu lattice, CRC/NIST scale",
        priority="medium",
        k=0.3,
    ),
    "v8_033_co_cu100_surface_mace": BenchmarkSpec(
        key="v8_033_co_cu100_surface_mace",
        title="CO diffusion on Cu(100)",
        calc="mace",
        build="surface_hop",
        symbol="Cu",
        host="Cu",
        adsorbate="CO",
        surface="fcc100",
        lattice_a=3.615,
        reference_barrier=0.06,
        reference_source="DFT CO/Cu(100) surface-diffusion literature",
        structure_source="Experimental FCC Cu lattice, CRC/NIST scale",
        priority="medium",
        k=0.3,
    ),
    "v8_034_n_mo110_surface_mace": BenchmarkSpec(
        key="v8_034_n_mo110_surface_mace",
        title="N adatom hop on Mo(110)",
        calc="mace",
        build="surface_hop",
        symbol="Mo",
        host="Mo",
        adsorbate="N",
        surface="bcc110",
        lattice_a=3.147,
        reference_barrier=0.55,
        reference_source="DFT N/Mo(110) diffusion literature",
        structure_source="Experimental BCC Mo lattice, CRC/NIST scale",
        priority="medium",
        k=0.3,
    ),
    "v8_035_ethane_torsion_egret": BenchmarkSpec(
        key="v8_035_ethane_torsion_egret",
        title="Ethane C-C torsion",
        calc="egret",
        build="molecule_torsion",
        symbol="C2H6",
        reference_barrier=0.113,
        reference_source="nebwalk v0.8.0 validated Egret-1t ethane torsion",
        structure_source="Standard gas-phase ethane geometry",
        priority="high",
        n_images=7,
        k=0.10,
        k_min=0.033,
        max_steps=400,
        climb_delay=60,
    ),
    "v8_036_propane_torsion_egret": BenchmarkSpec(
        key="v8_036_propane_torsion_egret",
        title="Propane C-C internal rotation",
        calc="egret",
        build="molecule_torsion",
        symbol="C3H8",
        reference_barrier=0.14,
        reference_source="MP2/internal-rotation literature guide value",
        structure_source="ASE molecule plus torsion endpoint construction",
        priority="medium",
        n_images=7,
        k=0.10,
        k_min=0.033,
    ),
    "v8_037_dme_torsion_egret": BenchmarkSpec(
        key="v8_037_dme_torsion_egret",
        title="Dimethyl ether C-O internal rotation",
        calc="egret",
        build="molecule_torsion",
        symbol="CH3OCH3",
        reference_barrier=0.11,
        reference_source="MP2/internal-rotation literature guide value",
        structure_source="Constructed DME gas-phase geometry",
        priority="medium",
        n_images=7,
        k=0.10,
        k_min=0.033,
    ),
    "v8_038_acetaldehyde_methyl_torsion_egret": BenchmarkSpec(
        key="v8_038_acetaldehyde_methyl_torsion_egret",
        title="Acetaldehyde methyl internal rotation",
        calc="egret",
        build="molecule_torsion",
        symbol="CH3CHO",
        reference_barrier=0.05,
        reference_source="MP2/internal-rotation literature guide value",
        structure_source="ASE molecule plus torsion endpoint construction",
        priority="medium",
        n_images=7,
        k=0.10,
        k_min=0.033,
    ),
}


def main(key: str | None = None) -> None:
    if key is None:
        if len(sys.argv) != 2 or sys.argv[1] not in SYSTEMS:
            names = "\n".join(f"  {name}" for name in SYSTEMS)
            raise SystemExit(
                f"Usage: python {Path(__file__).name} CASE\n\nCases:\n{names}"
            )
        key = sys.argv[1]
    try:
        run_case(SYSTEMS[key])
    except SkipCase as exc:
        print(f"SKIP: {exc}")
        raise SystemExit(SKIP_EXIT) from exc


def run_case(spec: BenchmarkSpec) -> None:
    out_dir = OUTPUT_ROOT / spec.key
    repro_dir = out_dir / "repro_bundle"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 72)
    print(spec.title)
    print(f"Case              : {spec.key}")
    print(f"Calculator        : {spec.calc.upper()}")
    print(f"Reference barrier : ~{spec.reference_barrier:.3f} eV")
    print(f"Reference source  : {spec.reference_source}")
    print(f"Structure source  : {spec.structure_source}")
    if spec.notes:
        print(f"Note              : {spec.notes}")
    print("=" * 72)

    initial, final = build_endpoints(spec)
    write_starter_files(spec, out_dir, initial, final)

    if os.environ.get("NEBWALK_V8_DRY_RUN") == "1":
        print(f"Dry run only. Starter structures written to {out_dir}")
        return

    calc_factory = make_calculator_factory(spec)
    relax_endpoint(initial, calc_factory, spec)
    relax_endpoint(final, calc_factory, spec)

    config = NEBRunConfig(
        n_images=spec.n_images,
        interpolation="idpp",
        k=spec.k,
        k_min=spec.k_min,
        climb=True,
        climb_delay=spec.climb_delay,
        n_workers=1,
        fmax=spec.fmax,
        max_steps=spec.max_steps,
        verbose=True,
    )
    result = run_neb_calculation(
        initial,
        final,
        calc_factory,
        config=config,
        reproduce_dir=repro_dir,
        calc_params=calc_params(spec),
    )
    write_outputs(spec, out_dir, result)
    print_summary(spec, result, repro_dir)


def make_calculator_factory(spec: BenchmarkSpec) -> Callable[[], object]:
    if spec.calc == "emt":
        check_emt_supported(spec)
        return EMT
    if spec.calc == "mace":
        return make_mace_factory()
    return make_egret_factory()


def check_emt_supported(spec: BenchmarkSpec) -> None:
    symbols: set[str] = set()
    for label in (spec.symbol, spec.host, spec.adsorbate):
        if label == "CO":
            symbols.update({"C", "O"})
        elif label:
            symbols.add(label)
    unsupported = sorted(sym for sym in symbols if sym and sym not in EMT_SYMBOLS)
    if unsupported:
        raise SkipCase(
            f"ASE EMT has no parameterization for: {', '.join(unsupported)}"
        )


def make_mace_factory() -> Callable[[], object]:
    try:
        from mace.calculators import mace_mp
    except ImportError as exc:
        raise SkipCase("MACE is not installed in this environment.") from exc

    model = os.environ.get("NEBWALK_V8_MACE_MODEL", "small")
    device = os.environ.get("NEBWALK_MACE_DEVICE", "cpu")
    dtype = os.environ.get("NEBWALK_MACE_DTYPE", "float32")

    def factory():
        return mace_mp(
            model=model,
            dispersion=False,
            default_dtype=dtype,
            device=device,
        )

    return factory


def make_egret_factory() -> Callable[[], object]:
    try:
        from mace.calculators import MACECalculator
    except ImportError as exc:
        raise SkipCase("MACE/Egret calculator support is not installed.") from exc

    model_path = Path(os.environ.get("EGRET_MODEL_PATH", "EGRET_1T.model"))
    if not model_path.exists():
        raise SkipCase(f"Egret model file not found: {model_path}")
    device = os.environ.get(
        "EGRET_DEVICE",
        os.environ.get("NEBWALK_MACE_DEVICE", "cpu"),
    )
    dtype = os.environ.get("EGRET_DTYPE", "float32")

    def factory():
        kwargs = {"device": device, "default_dtype": dtype}
        try:
            return MACECalculator(model_paths=str(model_path), **kwargs)
        except TypeError:
            return MACECalculator(model_path=str(model_path), **kwargs)

    return factory


def build_endpoints(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    if spec.build == "surface_hop":
        return build_surface_hop(spec)
    if spec.build == "bulk_vacancy":
        return build_bulk_vacancy(spec)
    if spec.build == "bulk_interstitial":
        return build_bulk_interstitial(spec)
    if spec.build == "divacancy":
        return build_divacancy(spec)
    if spec.build == "ionic_vacancy":
        return build_ionic_vacancy(spec)
    if spec.build == "h_interstitial":
        return build_h_interstitial(spec)
    if spec.build == "molecule_torsion":
        return build_molecule_torsion(spec)
    raise ValueError(f"Unknown build type: {spec.build}")


def build_surface_hop(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    slab = make_slab(spec)
    fix_bottom_layers(slab)
    z_top = float(slab.positions[:, 2].max())
    x_step = float(np.linalg.norm(slab.cell[0]) / 3.0)
    y_step = float(np.linalg.norm(slab.cell[1]) / 3.0)
    pos_a = np.array([0.5 * x_step, 0.5 * y_step, z_top + adsorbate_height(spec)])
    pos_b = np.array([1.5 * x_step, 0.5 * y_step, z_top + adsorbate_height(spec)])
    return add_adsorbate(slab, spec.adsorbate, pos_a), add_adsorbate(
        slab,
        spec.adsorbate,
        pos_b,
    )


def make_slab(spec: BenchmarkSpec) -> Atoms:
    if spec.surface == "fcc100":
        return fcc100(spec.host, size=(3, 3, 4), a=spec.lattice_a, vacuum=10.0)
    if spec.surface == "fcc111":
        return fcc111(spec.host, size=(3, 3, 4), a=spec.lattice_a, vacuum=10.0)
    if spec.surface == "bcc100":
        return bcc100(spec.host, size=(3, 3, 4), a=spec.lattice_a, vacuum=10.0)
    if spec.surface == "bcc110":
        return bcc110(spec.host, size=(3, 3, 4), a=spec.lattice_a, vacuum=10.0)
    if spec.surface == "hcp0001":
        return hcp0001(
            spec.host,
            size=(3, 3, 4),
            a=spec.lattice_a,
            c=spec.lattice_c,
            vacuum=10.0,
        )
    raise ValueError(f"Unsupported surface: {spec.surface}")


def add_adsorbate(slab: Atoms, adsorbate: str, pos: np.ndarray) -> Atoms:
    atoms = slab.copy()
    if adsorbate == "CO":
        atoms.append(Atom("C", position=pos))
        atoms.append(Atom("O", position=pos + np.array([0.0, 0.0, 1.15])))
    else:
        atoms.append(Atom(adsorbate, position=pos))
    atoms.info.clear()
    return atoms


def adsorbate_height(spec: BenchmarkSpec) -> float:
    return {"H": 1.0, "N": 1.2, "O": 1.2, "CO": 1.85}.get(spec.adsorbate, 1.8)


def fix_bottom_layers(slab: Atoms, n_layers: int = 2) -> None:
    z_values = np.unique(np.round(slab.positions[:, 2], 3))
    z_cut = np.sort(z_values)[max(0, n_layers - 1)]
    slab.set_constraint(FixAtoms(mask=slab.positions[:, 2] <= z_cut + 1e-3))


def build_bulk_vacancy(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    atoms = make_bulk(spec)
    return vacancy_hop(atoms, spec.symbol)


def build_bulk_interstitial(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    atoms = make_bulk(spec)
    a = spec.lattice_a
    pos_a = np.array([0.5 * a, 0.5 * a, 0.5 * a])
    pos_b = np.array([1.0 * a, 0.5 * a, 0.5 * a])
    initial = atoms.copy()
    final = atoms.copy()
    initial.append(Atom(spec.symbol, position=pos_a))
    final.append(Atom(spec.symbol, position=pos_b))
    return initial, final


def build_divacancy(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    full = make_bulk(spec)
    vacancy_pos = full.positions[0].copy()
    distances = np.linalg.norm(full.positions - vacancy_pos, axis=1)
    remove = sorted([0, int(np.argsort(distances)[2])], reverse=True)
    base = full.copy()
    for idx in remove:
        del base[idx]
    distances = np.linalg.norm(base.positions - vacancy_pos, axis=1)
    mover = int(np.argmin(distances))
    initial = base.copy()
    final = base.copy()
    final.positions[mover] = vacancy_pos
    return initial, final


def build_ionic_vacancy(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    atoms = make_bulk(spec)
    return vacancy_hop(atoms, spec.migrating_species)


def build_h_interstitial(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    atoms = make_bulk(spec)
    a = spec.lattice_a
    pos_a = np.array([1.5 * a, 1.5 * a, 0.5 * a])
    pos_b = np.array([1.0 * a, 1.0 * a, 0.5 * a])
    initial = atoms.copy()
    final = atoms.copy()
    initial.append(Atom("H", position=pos_a))
    final.append(Atom("H", position=pos_b))
    return initial, final


def vacancy_hop(atoms: Atoms, migrating_species: str) -> tuple[Atoms, Atoms]:
    indices = [
        idx for idx, atom in enumerate(atoms) if atom.symbol == migrating_species
    ]
    if len(indices) < 2:
        raise ValueError(
            f"Need at least two {migrating_species} atoms for vacancy hop."
        )
    vacancy_index = indices[0]
    vacancy_pos = atoms.positions[vacancy_index].copy()
    base = atoms.copy()
    del base[vacancy_index]
    candidates = [
        idx for idx, atom in enumerate(base) if atom.symbol == migrating_species
    ]
    distances = np.linalg.norm(base.positions[candidates] - vacancy_pos, axis=1)
    mover = candidates[int(np.argmin(distances))]
    initial = base.copy()
    final = base.copy()
    final.positions[mover] = vacancy_pos
    return initial, final


def make_bulk(spec: BenchmarkSpec) -> Atoms:
    a = spec.lattice_a
    c = spec.lattice_c
    if spec.crystal_type in {"fcc", "bcc"}:
        atoms = bulk(spec.symbol, spec.crystal_type, a=a, cubic=True)
    elif spec.crystal_type == "rocksalt":
        atoms = bulk(spec.symbol, "rocksalt", a=a, cubic=True)
    elif spec.crystal_type == "zincblende":
        atoms = bulk(spec.symbol, "zincblende", a=a, cubic=True)
    elif spec.crystal_type == "wurtzite":
        atoms = bulk(spec.symbol, "wurtzite", a=a, c=c)
    elif spec.crystal_type == "rutile":
        atoms = crystal(
            ("Ti", "O"),
            basis=[(0.0, 0.0, 0.0), (0.305, 0.305, 0.0)],
            spacegroup=136,
            cellpar=[a, a, c, 90, 90, 90],
        )
    elif spec.crystal_type == "fluorite":
        atoms = crystal(
            ("Ce", "O"),
            basis=[(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)],
            spacegroup=225,
            cellpar=[a, a, a, 90, 90, 90],
        )
    elif spec.crystal_type == "antifluorite":
        atoms = crystal(
            ("O", "Na"),
            basis=[(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)],
            spacegroup=225,
            cellpar=[a, a, a, 90, 90, 90],
        )
    elif spec.crystal_type == "antifluorite_mg2si":
        atoms = crystal(
            ("Si", "Mg"),
            basis=[(0.0, 0.0, 0.0), (0.25, 0.25, 0.25)],
            spacegroup=225,
            cellpar=[a, a, a, 90, 90, 90],
        )
    elif spec.crystal_type == "mo2c":
        atoms = crystal(
            ("Mo", "C"),
            basis=[(0.333, 0.667, 0.25), (0.0, 0.0, 0.0)],
            spacegroup=194,
            cellpar=[a, a, c, 90, 90, 120],
        )
    elif spec.crystal_type == "layered_r3m":
        alkali = "Li" if spec.symbol.startswith("Li") else "Na"
        atoms = crystal(
            (alkali, "Co", "O"),
            basis=[(0.0, 0.0, 0.0), (0.0, 0.0, 0.5), (0.0, 0.0, 0.24)],
            spacegroup=166,
            cellpar=[a, a, c, 90, 90, 120],
        )
    elif spec.crystal_type == "li3po4":
        atoms = make_li3po4_proxy(a, c or 10.49)
    else:
        raise ValueError(f"Unsupported crystal type: {spec.crystal_type}")
    atoms = atoms.repeat(spec.repeat)
    atoms.pbc = True
    return atoms


def make_li3po4_proxy(a: float, c: float) -> Atoms:
    cell = [a, a, c]
    positions = [
        (0.0, 0.0, 0.0),
        (0.25, 0.25, 0.15),
        (0.75, 0.25, 0.15),
        (0.25, 0.75, 0.15),
        (0.5, 0.5, 0.5),
        (0.62, 0.5, 0.5),
        (0.38, 0.5, 0.5),
        (0.5, 0.62, 0.5),
        (0.5, 0.38, 0.5),
    ]
    symbols = ["P", "Li", "Li", "Li", "O", "O", "O", "O", "O"]
    atoms = Atoms(symbols=symbols, scaled_positions=positions, cell=cell, pbc=True)
    return atoms


def build_molecule_torsion(spec: BenchmarkSpec) -> tuple[Atoms, Atoms]:
    if spec.symbol == "C2H6":
        initial = ethane(60.0)
        final = ethane(180.0)
    elif spec.symbol == "C3H8":
        initial = molecule("C3H8")
        final = rotate_subset(initial, [0, 3, 4, 5], 120.0)
    elif spec.symbol == "CH3OCH3":
        initial = dimethyl_ether()
        final = rotate_subset(initial, [2, 6, 7, 8], 120.0)
    elif spec.symbol == "CH3CHO":
        initial = molecule("CH3CHO")
        final = rotate_subset(initial, [0, 4, 5, 6], 120.0)
    else:
        raise ValueError(f"Unsupported molecule: {spec.symbol}")
    for atoms in (initial, final):
        atoms.set_cell(np.diag([20.0, 20.0, 20.0]))
        atoms.pbc = False
        atoms.center()
    return initial, final


def ethane(phi_deg: float) -> Atoms:
    cc = 0.770
    r_ch = 1.090
    dz = abs(r_ch * np.cos(np.radians(111.2)))
    rl = r_ch * np.sin(np.radians(111.2))
    phi = np.radians(phi_deg)
    pos = [[0.0, 0.0, -cc], [0.0, 0.0, cc]]
    for idx in range(3):
        angle = idx * 2.0 * np.pi / 3.0
        pos.append([rl * np.cos(angle), rl * np.sin(angle), -cc - dz])
    for idx in range(3):
        angle = phi + idx * 2.0 * np.pi / 3.0
        pos.append([rl * np.cos(angle), rl * np.sin(angle), cc + dz])
    return Atoms("C2H6", positions=pos)


def dimethyl_ether() -> Atoms:
    atoms = Atoms(
        "COC" + "H6",
        positions=[
            [-1.42, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.42, 0.0, 0.0],
            [-1.78, 1.02, 0.0],
            [-1.78, -0.51, 0.88],
            [-1.78, -0.51, -0.88],
            [1.78, 1.02, 0.0],
            [1.78, -0.51, 0.88],
            [1.78, -0.51, -0.88],
        ],
    )
    return atoms


def rotate_subset(atoms: Atoms, indices: list[int], angle_deg: float) -> Atoms:
    out = atoms.copy()
    origin = out.positions[indices[0]].copy()
    axis = np.array([1.0, 0.0, 0.0])
    theta = np.radians(angle_deg)
    rot = rotation_matrix(axis, theta)
    for idx in indices[1:]:
        out.positions[idx] = origin + rot @ (out.positions[idx] - origin)
    return out


def rotation_matrix(axis: np.ndarray, theta: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    ux, uy, uz = axis
    c = np.cos(theta)
    s = np.sin(theta)
    return np.array(
        [
            [
                c + ux * ux * (1 - c),
                ux * uy * (1 - c) - uz * s,
                ux * uz * (1 - c) + uy * s,
            ],
            [
                uy * ux * (1 - c) + uz * s,
                c + uy * uy * (1 - c),
                uy * uz * (1 - c) - ux * s,
            ],
            [
                uz * ux * (1 - c) - uy * s,
                uz * uy * (1 - c) + ux * s,
                c + uz * uz * (1 - c),
            ],
        ]
    )


def relax_endpoint(
    atoms: Atoms,
    calc_factory: Callable[[], object],
    spec: BenchmarkSpec,
) -> None:
    atoms.calc = calc_factory()
    BFGS(atoms, logfile=None).run(fmax=spec.relax_fmax, steps=spec.relax_steps)


def write_starter_files(
    spec: BenchmarkSpec,
    out_dir: Path,
    initial: Atoms,
    final: Atoms,
) -> None:
    from ase.io import write

    start = initial.copy()
    end = final.copy()
    start.info.clear()
    end.info.clear()
    write(out_dir / f"{spec.key}_initial.xyz", start, format="extxyz")
    write(out_dir / f"{spec.key}_final.xyz", end, format="extxyz")


def write_outputs(spec: BenchmarkSpec, out_dir: Path, result) -> None:
    result.neb.plot(str(out_dir / f"{spec.key}_profile.png"), title=spec.title)
    result.neb.save_csv(str(out_dir / f"{spec.key}_profile.csv"))
    result.neb.save_trajectory(str(out_dir / f"{spec.key}_path.traj"))


def print_summary(spec: BenchmarkSpec, result, repro_dir: Path) -> None:
    barrier = float(result.barrier)
    err = abs(barrier - spec.reference_barrier) / spec.reference_barrier * 100.0
    print()
    print(f"Converged       : {result.converged}")
    print(f"Forward barrier : {barrier:.4f} eV")
    print(f"Reference       : ~{spec.reference_barrier:.4f} eV")
    print(f"Error           : {err:.1f}%")
    print(f"Repro bundle    : {repro_dir}")


def calc_params(spec: BenchmarkSpec) -> dict[str, object]:
    payload: dict[str, object] = {
        "type": spec.calc.upper(),
        "case": spec.key,
        "reference_barrier_eV": spec.reference_barrier,
        "reference_source": spec.reference_source,
        "structure_source": spec.structure_source,
        "priority": spec.priority,
    }
    if spec.calc == "mace":
        payload.update(
            {
                "model": os.environ.get("NEBWALK_V8_MACE_MODEL", "small"),
                "device": os.environ.get("NEBWALK_MACE_DEVICE", "cpu"),
            }
        )
    if spec.calc == "egret":
        payload.update(
            {
                "model_path": os.environ.get("EGRET_MODEL_PATH", "EGRET_1T.model"),
                "device": os.environ.get("EGRET_DEVICE", "cpu"),
            }
        )
    return payload


class SkipCase(RuntimeError):
    """A scientifically unsupported or unavailable optional benchmark."""


EMT_SYMBOLS = {"H", "C", "N", "O", "Al", "Ni", "Cu", "Pd", "Ag", "Pt", "Au"}


if __name__ == "__main__":
    main()
