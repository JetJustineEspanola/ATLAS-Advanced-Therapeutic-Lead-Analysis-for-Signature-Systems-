#!/usr/bin/env python3
"""
ATLAS — Stage 04S Target-Supported Molecular Docking

Strict docking implementation for the current 04R shortlist.

Default experimental structure plan
-----------------------------------
CYP3A4  P08684  -> PDB 3NXU, ritonavir-bound, 2.00 Å
DPP4    P27487  -> PDB 1X70, sitagliptin-bound
NPSR1   Q6W5P4  -> no experimental structure is used by default; candidate is skipped

Important guardrails
--------------------
- Docking score is not proof of binding, efficacy, or trastuzumab-resistance reversal.
- The network-selected target is used only after 04R target validation.
- CYP3A4 heme (HEM) is retained in the receptor.
- Binding boxes are derived from co-crystallized ligand coordinates.
- Experimental structures are preferred. NPSR1 is not silently replaced by a predicted model.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from vina import Vina

try:
    from openbabel import openbabel as ob
    OPENBABEL_PYTHON_AVAILABLE = True
except Exception:
    ob = None
    OPENBABEL_PYTHON_AVAILABLE = False


PROJECT_ROOT = Path(__file__).resolve().parents[1]

INPUT_SHORTLIST = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "final_prioritization"
    / "ATLAS_docking_shortlist.csv"
)

OUT_DIR = (
    PROJECT_ROOT
    / "results"
    / "cmap"
    / "docking"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_RESULTS = OUT_DIR / "ATLAS_docking_results.csv"
OUT_SUMMARY = OUT_DIR / "ATLAS_docking_summary.csv"
OUT_META = OUT_DIR / "ATLAS_docking_metadata.json"


STRUCTURE_PLAN = {
    "CYP3A4": {
        "uniprot": "P08684",
        "pdb_id": "3NXU",
        "reference_ligand": "RIT",
        "keep_hetero": {"HEM"},
        "box_padding": 6.0,
        "minimum_box": 22.0,
        "note": "Human CYP3A4; heme retained; ritonavir defines active-site box.",
    },
    "DPP4": {
        "uniprot": "P27487",
        "pdb_id": "1X70",
        "reference_ligand": "715",
        "keep_hetero": set(),
        "box_padding": 6.0,
        "minimum_box": 22.0,
        "note": "Human DPP4; co-crystallized sitagliptin defines active-site box.",
    },
}


def header(text: str) -> None:
    print("\n" + "=" * 78, flush=True)
    print(text, flush=True)
    print("=" * 78, flush=True)


def clean_text(x: Any) -> str:
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def slugify(x: str) -> str:
    x = re.sub(r"[^A-Za-z0-9._-]+", "_", x.strip())
    return x.strip("_") or "item"


def atomic_csv(df: pd.DataFrame, path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(payload: dict[str, Any], path: Path) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def command_exists(name: str) -> bool:
    return shutil.which(name) is not None


def run_cmd(cmd: list[str], cwd: Path | None = None) -> tuple[int, str]:
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    return p.returncode, p.stdout


def download(url: str, path: Path, timeout: int = 30) -> None:
    if path.exists() and path.stat().st_size > 100:
        return
    r = requests.get(
        url,
        timeout=timeout,
        headers={"User-Agent": "ATLAS-04S/1.0"},
    )
    r.raise_for_status()
    path.write_bytes(r.content)


def pdb_atom_fields(line: str) -> dict[str, Any]:
    return {
        "record": line[0:6].strip(),
        "atom_name": line[12:16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21:22].strip(),
        "resseq": line[22:26].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "element": line[76:78].strip() if len(line) >= 78 else "",
    }


def choose_reference_instance(
    pdb_lines: list[str],
    ligand_resname: str,
) -> tuple[str, str, list[str]]:
    """
    Pick the co-crystal ligand instance with the largest atom count.
    Returns chain, residue number, ligand PDB lines.
    """
    groups: dict[tuple[str, str], list[str]] = {}

    for line in pdb_lines:
        if not line.startswith("HETATM"):
            continue
        f = pdb_atom_fields(line)
        if f["resname"].upper() != ligand_resname.upper():
            continue
        key = (f["chain"], f["resseq"])
        groups.setdefault(key, []).append(line)

    if not groups:
        raise RuntimeError(
            f"Could not find co-crystal ligand {ligand_resname} in PDB."
        )

    key = max(groups, key=lambda k: len(groups[k]))
    return key[0], key[1], groups[key]


def write_receptor_and_reference(
    raw_pdb: Path,
    receptor_pdb: Path,
    ref_pdb: Path,
    reference_ligand: str,
    keep_hetero: set[str],
) -> dict[str, Any]:
    lines = raw_pdb.read_text(errors="replace").splitlines()

    chain, ref_resseq, ligand_lines = choose_reference_instance(
        lines,
        reference_ligand,
    )

    receptor_lines: list[str] = []

    for line in lines:
        if line.startswith("ATOM"):
            f = pdb_atom_fields(line)
            if f["chain"] == chain:
                receptor_lines.append(line)
        elif line.startswith("HETATM"):
            f = pdb_atom_fields(line)
            if f["chain"] == chain and f["resname"] in keep_hetero:
                receptor_lines.append(line)

    receptor_lines.append("END")
    receptor_pdb.write_text("\n".join(receptor_lines) + "\n")

    ref_pdb.write_text(
        "\n".join(ligand_lines + ["END"]) + "\n"
    )

    coords = np.array([
        [pdb_atom_fields(x)["x"],
         pdb_atom_fields(x)["y"],
         pdb_atom_fields(x)["z"]]
        for x in ligand_lines
    ])

    center = coords.mean(axis=0)
    span = coords.max(axis=0) - coords.min(axis=0)

    return {
        "chain": chain,
        "reference_resseq": ref_resseq,
        "reference_atom_n": int(len(ligand_lines)),
        "center": center,
        "span": span,
    }


def ensure_ligand_3d(input_sdf: Path, output_sdf: Path) -> None:
    supplier = Chem.SDMolSupplier(
        str(input_sdf),
        removeHs=False,
        sanitize=True,
    )
    mol = next((m for m in supplier if m is not None), None)

    if mol is None:
        raise RuntimeError(f"RDKit could not parse {input_sdf}")

    mol = Chem.AddHs(mol, addCoords=True)

    needs_embed = mol.GetNumConformers() == 0
    if not needs_embed:
        conf = mol.GetConformer()
        coords = np.array([
            list(conf.GetAtomPosition(i))
            for i in range(mol.GetNumAtoms())
        ])
        if np.nanstd(coords) < 1e-3:
            needs_embed = True

    if needs_embed:
        params = AllChem.ETKDGv3()
        params.randomSeed = 20260831
        status = AllChem.EmbedMolecule(mol, params)
        if status != 0:
            raise RuntimeError("RDKit ETKDG embedding failed.")

    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        pass

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(mol)
    writer.close()


def fetch_candidate_sdf(
    row: pd.Series,
    path: Path,
) -> str:
    cid = clean_text(row.get("pubchem_cid"))
    name = clean_text(row.get("pert_iname"))

    urls = []

    if cid and cid.lower() != "nan":
        cid = re.sub(r"\.0$", "", cid)
        urls.append(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
            f"cid/{cid}/record/SDF?record_type=3d"
        )

    if name:
        urls.append(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
            f"name/{requests.utils.quote(name)}/record/SDF?record_type=3d"
        )

    last_error = ""

    for url in urls:
        try:
            r = requests.get(
                url,
                timeout=30,
                headers={"User-Agent": "ATLAS-04S/1.0"},
            )
            if r.status_code == 200 and len(r.content) > 100:
                path.write_bytes(r.content)
                return url
            last_error = f"HTTP {r.status_code}"
        except Exception as e:
            last_error = str(e)

    raise RuntimeError(
        f"Could not retrieve PubChem SDF for {name}: {last_error}"
    )


def prepare_ligand(
    raw_sdf: Path,
    prepared_sdf: Path,
    pdbqt: Path,
) -> str:
    ensure_ligand_3d(raw_sdf, prepared_sdf)

    code, out = run_cmd([
        "mk_prepare_ligand.py",
        "-i", str(prepared_sdf),
        "-o", str(pdbqt),
    ])

    if code != 0 or not pdbqt.exists():
        raise RuntimeError(
            "Meeko ligand preparation failed:\n" + out[-3000:]
        )

    return out



def openbabel_convert(
    input_path: Path,
    output_path: Path,
    in_format: str,
    out_format: str,
    add_hydrogens: bool = False,
) -> None:
    """
    Convert molecular files using OpenBabel Python bindings.
    This avoids requiring the system `obabel` executable.
    """
    if not OPENBABEL_PYTHON_AVAILABLE:
        raise RuntimeError(
            "OpenBabel Python bindings are not available. "
            "Install with: python -m pip install openbabel-wheel"
        )

    conv = ob.OBConversion()
    if not conv.SetInAndOutFormats(in_format, out_format):
        raise RuntimeError(
            f"OpenBabel could not configure {in_format}->{out_format}"
        )

    mol = ob.OBMol()

    if not conv.ReadFile(mol, str(input_path)):
        raise RuntimeError(
            f"OpenBabel failed to read {input_path}"
        )

    if add_hydrogens:
        mol.AddHydrogens()

    if not conv.WriteFile(mol, str(output_path)):
        raise RuntimeError(
            f"OpenBabel failed to write {output_path}"
        )



def rigidify_pdbqt_text(text: str) -> str:
    """
    Convert an OpenBabel ligand-style PDBQT into rigid-receptor PDBQT.

    OpenBabel may emit ROOT/BRANCH/TORSDOF records because it represents the
    molecule as a torsion tree. AutoDock Vina rejects those records in a rigid
    receptor. For a rigid receptor we keep atom/heteroatom records and
    non-torsion remarks, and discard torsion-tree control records.
    """
    forbidden_prefixes = (
        "ROOT",
        "ENDROOT",
        "BRANCH",
        "ENDBRANCH",
        "TORSDOF",
    )

    out = []

    for line in text.splitlines():
        stripped = line.strip()

        if not stripped:
            continue

        if stripped.startswith(forbidden_prefixes):
            continue

        # Keep PDBQT atom records plus harmless metadata.
        if (
            line.startswith("ATOM")
            or line.startswith("HETATM")
            or line.startswith("REMARK")
            or line.startswith("TER")
        ):
            out.append(line)

    if not any(
        x.startswith("ATOM") or x.startswith("HETATM")
        for x in out
    ):
        raise RuntimeError(
            "Rigidified receptor PDBQT contains no atoms."
        )

    return "\n".join(out) + "\n"


def validate_rigid_receptor_pdbqt(path: Path) -> None:
    """
    Fail early if a receptor PDBQT still contains ligand torsion-tree tags.
    """
    text = path.read_text(errors="replace")

    bad = [
        tag for tag in [
            "ROOT",
            "ENDROOT",
            "BRANCH",
            "ENDBRANCH",
            "TORSDOF",
        ]
        if re.search(rf"(?m)^\s*{tag}\b", text)
    ]

    if bad:
        raise RuntimeError(
            "Rigid receptor PDBQT still contains torsion-tree tags: "
            + ", ".join(bad)
        )

    atom_n = sum(
        1 for line in text.splitlines()
        if line.startswith("ATOM") or line.startswith("HETATM")
    )

    if atom_n == 0:
        raise RuntimeError(
            "Rigid receptor PDBQT contains zero atoms."
        )


def prepare_receptor(
    receptor_pdb: Path,
    receptor_pdbqt: Path,
    contains_critical_cofactor: bool,
) -> tuple[str, str]:
    """
    Prepare a rigid receptor PDBQT accepted by AutoDock Vina.

    Strategy
    --------
    1. For ordinary protein receptors, try Meeko with ProDy and allow bad/
       incomplete residues to be skipped rather than aborting.
    2. For heme-containing CYP3A4, or if Meeko fails, use OpenBabel Python.
    3. OpenBabel can emit ligand torsion-tree tags (ROOT/BRANCH/TORSDOF).
       Those tags are stripped because the receptor is rigid.
    """
    logs = []

    if not contains_critical_cofactor:
        meeko_attempts = [
            [
                "mk_prepare_receptor.py",
                "-i", str(receptor_pdb),
                "--allow_bad_res",
                "--write_pdbqt", str(receptor_pdbqt),
            ],
            [
                "mk_prepare_receptor.py",
                "--read_pdb", str(receptor_pdb),
                "--allow_bad_res",
                "--write_pdbqt", str(receptor_pdbqt),
            ],
        ]

        for attempt_i, cmd in enumerate(meeko_attempts, start=1):
            # Remove stale output so success is unambiguous.
            if receptor_pdbqt.exists():
                receptor_pdbqt.unlink()

            code, out = run_cmd(cmd)
            logs.append(
                f"MEEKO ATTEMPT {attempt_i}:\n{out}"
            )

            if code == 0 and receptor_pdbqt.exists():
                try:
                    validate_rigid_receptor_pdbqt(
                        receptor_pdbqt
                    )
                    return (
                        "MEEKO",
                        "\n".join(logs),
                    )
                except Exception as e:
                    logs.append(
                        "MEEKO OUTPUT VALIDATION FAILED: "
                        + str(e)
                    )

    if not OPENBABEL_PYTHON_AVAILABLE:
        raise RuntimeError(
            "OpenBabel Python bindings are required for receptor fallback. "
            "Run: python -m pip install openbabel-wheel"
        )

    temp_pdbqt = receptor_pdbqt.with_suffix(
        ".openbabel_raw.pdbqt"
    )

    try:
        if temp_pdbqt.exists():
            temp_pdbqt.unlink()

        openbabel_convert(
            receptor_pdb,
            temp_pdbqt,
            "pdb",
            "pdbqt",
            add_hydrogens=True,
        )

        raw_text = temp_pdbqt.read_text(
            errors="replace"
        )

        rigid_text = rigidify_pdbqt_text(
            raw_text
        )

        receptor_pdbqt.write_text(
            rigid_text,
            encoding="utf-8",
        )

        validate_rigid_receptor_pdbqt(
            receptor_pdbqt
        )

        logs.append(
            "OPENBABEL_PYTHON: converted receptor to PDBQT; "
            "removed ROOT/BRANCH/TORSDOF torsion-tree records "
            "for rigid-receptor compatibility."
        )

    except Exception as e:
        raise RuntimeError(
            "OpenBabel rigid receptor preparation failed: "
            + str(e)
        )

    return (
        "OPENBABEL_PYTHON_RIGIDIFIED",
        "\n".join(logs),
    )


def run_vina(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    out_pdbqt: Path,
    center: np.ndarray,
    size: np.ndarray,
    exhaustiveness: int,
    num_modes: int,
) -> tuple[list[float], str]:
    """
    Run AutoDock Vina through its installed Python API.

    The PyPI `vina` package exposes the Vina engine as a Python module but
    normally does not install a `vina` shell executable.
    """
    try:
        validate_rigid_receptor_pdbqt(receptor_pdbqt)

        v = Vina(sf_name="vina")
        v.set_receptor(str(receptor_pdbqt))
        v.set_ligand_from_file(str(ligand_pdbqt))

        v.compute_vina_maps(
            center=[float(x) for x in center],
            box_size=[float(x) for x in size],
        )

        v.dock(
            exhaustiveness=int(exhaustiveness),
            n_poses=int(num_modes),
        )

        v.write_poses(
            str(out_pdbqt),
            n_poses=int(num_modes),
            overwrite=True,
        )

        energies = np.asarray(
            v.energies(n_poses=int(num_modes)),
            dtype=float,
        )

        if energies.ndim == 1:
            energies = energies.reshape(1, -1)

        scores = (
            energies[:, 0].astype(float).tolist()
            if energies.size
            else []
        )

        log = (
            "AutoDock Vina Python API\n"
            f"receptor={receptor_pdbqt}\n"
            f"ligand={ligand_pdbqt}\n"
            f"center={[float(x) for x in center]}\n"
            f"box_size={[float(x) for x in size]}\n"
            f"exhaustiveness={exhaustiveness}\n"
            f"num_modes={num_modes}\n"
            f"energies={energies.tolist() if energies.size else []}\n"
        )

        return scores, log

    except Exception as e:
        raise RuntimeError(
            "AutoDock Vina Python API failed: " + str(e)
        )


def try_reference_redocking(
    ref_pdb: Path,
    receptor_pdbqt: Path,
    workdir: Path,
    center: np.ndarray,
    size: np.ndarray,
    exhaustiveness: int,
) -> dict[str, Any]:
    """
    Attempt reference-ligand redocking as a protocol sanity check.
    RMSD is attempted only if Open Babel can reconstruct both crystal and pose.
    """
    result = {
        "reference_redocking_attempted": False,
        "reference_best_affinity_kcal_mol": np.nan,
        "reference_redocking_rmsd_A": np.nan,
        "reference_redocking_status": "NOT_ATTEMPTED",
    }

    if not OPENBABEL_PYTHON_AVAILABLE:
        result["reference_redocking_status"] = (
            "SKIPPED_NO_OPENBABEL_PYTHON"
        )
        return result

    crystal_sdf = workdir / "reference_crystal.sdf"

    try:
        openbabel_convert(
            ref_pdb,
            crystal_sdf,
            "pdb",
            "sdf",
            add_hydrogens=True,
        )
    except Exception:
        result["reference_redocking_status"] = (
            "CRYSTAL_CONVERSION_FAILED"
        )
        return result

    prepared_sdf = workdir / "reference_prepared.sdf"
    ref_pdbqt = workdir / "reference_ligand.pdbqt"

    try:
        # Preserve crystal geometry; only add Hs.
        supplier = Chem.SDMolSupplier(
            str(crystal_sdf),
            removeHs=False,
            sanitize=True,
        )
        mol = next((m for m in supplier if m is not None), None)
        if mol is None:
            raise RuntimeError("RDKit could not parse reference ligand.")

        mol = Chem.AddHs(mol, addCoords=True)
        w = Chem.SDWriter(str(prepared_sdf))
        w.write(mol)
        w.close()

        code, out = run_cmd([
            "mk_prepare_ligand.py",
            "-i", str(prepared_sdf),
            "-o", str(ref_pdbqt),
        ])

        if code != 0:
            raise RuntimeError(out)

        pose = workdir / "reference_redocked.pdbqt"
        scores, vina_log = run_vina(
            receptor_pdbqt,
            ref_pdbqt,
            pose,
            center,
            size,
            exhaustiveness,
            1,
        )

        result["reference_redocking_attempted"] = True
        if scores:
            result["reference_best_affinity_kcal_mol"] = scores[0]

        docked_sdf = workdir / "reference_redocked.sdf"

        conversion_ok = True
        try:
            openbabel_convert(
                pose,
                docked_sdf,
                "pdbqt",
                "sdf",
                add_hydrogens=False,
            )
        except Exception:
            conversion_ok = False

        if conversion_ok and docked_sdf.exists():
            ref_sup = Chem.SDMolSupplier(
                str(prepared_sdf),
                removeHs=True,
            )
            dock_sup = Chem.SDMolSupplier(
                str(docked_sdf),
                removeHs=True,
            )

            ref_mol = next((m for m in ref_sup if m is not None), None)
            dock_mol = next((m for m in dock_sup if m is not None), None)

            if (
                ref_mol is not None
                and dock_mol is not None
                and ref_mol.GetNumAtoms() == dock_mol.GetNumAtoms()
            ):
                try:
                    rms = AllChem.GetBestRMS(
                        ref_mol,
                        dock_mol,
                    )
                    result["reference_redocking_rmsd_A"] = float(rms)
                except Exception:
                    pass

        rms = result["reference_redocking_rmsd_A"]

        if pd.notna(rms):
            result["reference_redocking_status"] = (
                "PASS_RMSD_LE_2A"
                if rms <= 2.0
                else "CAUTION_RMSD_GT_2A"
            )
        else:
            result["reference_redocking_status"] = (
                "COMPLETED_RMSD_UNAVAILABLE"
            )

    except Exception as e:
        result["reference_redocking_status"] = (
            "FAILED: " + str(e)[:250]
        )

    return result


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()

    p.add_argument(
        "--exhaustiveness",
        type=int,
        default=16,
    )
    p.add_argument(
        "--num-modes",
        type=int,
        default=9,
    )
    p.add_argument(
        "--max-candidates",
        type=int,
        default=5,
    )
    p.add_argument(
        "--skip-redocking",
        action="store_true",
    )

    return p.parse_args()


def main() -> int:
    args = parse_args()

    header("ATLAS — Stage 04S Target-Supported Molecular Docking")

    if not INPUT_SHORTLIST.exists():
        print(
            f"ERROR: 04R docking shortlist not found:\n{INPUT_SHORTLIST}",
            flush=True,
        )
        return 1

    missing_tools = [
        x for x in [
            "mk_prepare_ligand.py",
            "mk_prepare_receptor.py",
        ]
        if not command_exists(x)
    ]

    if missing_tools:
        print("\nMissing required Meeko command-line tools:", flush=True)
        for x in missing_tools:
            print(f"  - {x}", flush=True)

        print(
            "\nInstall inside the ATLAS venv with:\n"
            "  python -m pip install -U meeko gemmi\n",
            flush=True,
        )
        return 2

    if not OPENBABEL_PYTHON_AVAILABLE:
        print(
            "\nOpenBabel Python bindings are missing.\n"
            "Install them without pacman using:\n"
            "  python -m pip install -U openbabel-wheel\n",
            flush=True,
        )
        return 2

    print("Vina backend: Python API", flush=True)
    print("OpenBabel backend: Python bindings", flush=True)

    shortlist = pd.read_csv(INPUT_SHORTLIST).head(
        args.max_candidates
    )

    print(f"\n04R docking candidates: {len(shortlist)}", flush=True)
    print(f"Vina exhaustiveness: {args.exhaustiveness}", flush=True)
    print(f"Vina modes: {args.num_modes}", flush=True)

    records: list[dict[str, Any]] = []

    for i, row in shortlist.iterrows():
        compound = clean_text(row.get("pert_iname"))
        target = clean_text(row.get("validated_target_symbol"))
        uniprot = clean_text(row.get("validated_uniprot_accession"))

        header(f"DOCKING {compound} -> {target}")

        base_record = {
            "docking_rank": row.get("docking_rank", i + 1),
            "pert_iname": compound,
            "target_symbol": target,
            "target_uniprot": uniprot,
            "04r_integrated_score": row.get(
                "integrated_prioritization_score",
                np.nan,
            ),
        }

        plan = STRUCTURE_PLAN.get(target)

        if plan is None or plan["uniprot"] != uniprot:
            print(
                "SKIPPED: no vetted experimental docking structure "
                "configured for this target.",
                flush=True,
            )
            records.append({
                **base_record,
                "docking_status": "SKIPPED_NO_VETTED_EXPERIMENTAL_STRUCTURE",
                "pdb_id": "",
                "best_affinity_kcal_mol": np.nan,
                "protocol_validation": "NOT_APPLICABLE",
                "interpretation": (
                    "Candidate retained as exploratory; no predicted "
                    "structure was silently substituted."
                ),
            })
            continue

        workdir = OUT_DIR / (
            f"{int(row.get('docking_rank', i + 1)):02d}_"
            f"{slugify(compound)}_{target}_{plan['pdb_id']}"
        )
        workdir.mkdir(parents=True, exist_ok=True)

        try:
            # ---------------------------------------------------------
            # PDB download and binding-site definition
            # ---------------------------------------------------------
            pdb_id = plan["pdb_id"]
            raw_pdb = workdir / f"{pdb_id}.pdb"

            download(
                f"https://files.rcsb.org/download/{pdb_id}.pdb",
                raw_pdb,
            )

            receptor_pdb = workdir / "receptor_clean.pdb"
            ref_pdb = workdir / (
                f"reference_{plan['reference_ligand']}.pdb"
            )

            site = write_receptor_and_reference(
                raw_pdb,
                receptor_pdb,
                ref_pdb,
                plan["reference_ligand"],
                plan["keep_hetero"],
            )

            center = np.asarray(site["center"], dtype=float)
            span = np.asarray(site["span"], dtype=float)

            size = np.maximum(
                span + 2.0 * float(plan["box_padding"]),
                float(plan["minimum_box"]),
            )
            size = np.minimum(size, 30.0)

            print(
                f"PDB: {pdb_id}; protein chain: {site['chain']}; "
                f"reference ligand: {plan['reference_ligand']}",
                flush=True,
            )
            print(
                "Box center: "
                f"{center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f}",
                flush=True,
            )
            print(
                "Box size: "
                f"{size[0]:.1f} x {size[1]:.1f} x {size[2]:.1f} Å",
                flush=True,
            )

            # ---------------------------------------------------------
            # Receptor
            # ---------------------------------------------------------
            receptor_pdbqt = workdir / "receptor.pdbqt"
            receptor_method, receptor_log = prepare_receptor(
                receptor_pdb,
                receptor_pdbqt,
                contains_critical_cofactor=(
                    "HEM" in plan["keep_hetero"]
                ),
            )
            (workdir / "receptor_prep.log").write_text(
                receptor_log,
                encoding="utf-8",
            )

            # ---------------------------------------------------------
            # Candidate ligand
            # ---------------------------------------------------------
            ligand_raw = workdir / "candidate_pubchem.sdf"
            ligand_source = fetch_candidate_sdf(row, ligand_raw)

            ligand_prepared = workdir / "candidate_prepared.sdf"
            ligand_pdbqt = workdir / "candidate.pdbqt"

            ligand_log = prepare_ligand(
                ligand_raw,
                ligand_prepared,
                ligand_pdbqt,
            )
            (workdir / "ligand_prep.log").write_text(
                ligand_log,
                encoding="utf-8",
            )

            # ---------------------------------------------------------
            # Reference redocking
            # ---------------------------------------------------------
            if args.skip_redocking:
                redock = {
                    "reference_redocking_attempted": False,
                    "reference_best_affinity_kcal_mol": np.nan,
                    "reference_redocking_rmsd_A": np.nan,
                    "reference_redocking_status": "SKIPPED_BY_USER",
                }
            else:
                redock = try_reference_redocking(
                    ref_pdb,
                    receptor_pdbqt,
                    workdir,
                    center,
                    size,
                    args.exhaustiveness,
                )

            print(
                "Protocol validation: "
                f"{redock['reference_redocking_status']}",
                flush=True,
            )

            # ---------------------------------------------------------
            # Candidate docking
            # ---------------------------------------------------------
            pose_out = workdir / "candidate_docked.pdbqt"

            scores, vina_log = run_vina(
                receptor_pdbqt,
                ligand_pdbqt,
                pose_out,
                center,
                size,
                args.exhaustiveness,
                args.num_modes,
            )

            (workdir / "vina.log").write_text(
                vina_log,
                encoding="utf-8",
            )

            best_score = scores[0] if scores else np.nan

            print(
                f"Best Vina affinity: {best_score:.3f} kcal/mol"
                if pd.notna(best_score)
                else "No Vina score parsed.",
                flush=True,
            )

            if scores:
                score_df = pd.DataFrame({
                    "mode": range(1, len(scores) + 1),
                    "affinity_kcal_mol": scores,
                })
                score_df.to_csv(
                    workdir / "vina_modes.csv",
                    index=False,
                )

            records.append({
                **base_record,
                "docking_status": "COMPLETED",
                "pdb_id": pdb_id,
                "protein_chain": site["chain"],
                "reference_ligand": plan["reference_ligand"],
                "receptor_preparation_method": receptor_method,
                "critical_cofactor_retained": (
                    "HEM" in plan["keep_hetero"]
                ),
                "box_center_x": center[0],
                "box_center_y": center[1],
                "box_center_z": center[2],
                "box_size_x": size[0],
                "box_size_y": size[1],
                "box_size_z": size[2],
                "best_affinity_kcal_mol": best_score,
                "vina_mode_n": len(scores),
                "reference_best_affinity_kcal_mol": redock[
                    "reference_best_affinity_kcal_mol"
                ],
                "reference_redocking_rmsd_A": redock[
                    "reference_redocking_rmsd_A"
                ],
                "protocol_validation": redock[
                    "reference_redocking_status"
                ],
                "ligand_source": ligand_source,
                "workdir": str(workdir),
                "interpretation": (
                    "Computational structural-compatibility evidence only; "
                    "not proof of binding or trastuzumab-resistance reversal."
                ),
            })

        except Exception as e:
            print(f"FAILED: {e}", flush=True)
            records.append({
                **base_record,
                "docking_status": "FAILED",
                "pdb_id": plan["pdb_id"],
                "best_affinity_kcal_mol": np.nan,
                "protocol_validation": "FAILED",
                "interpretation": str(e)[:500],
            })

    # -----------------------------------------------------------------
    # Outputs
    # -----------------------------------------------------------------
    results = pd.DataFrame(records)

    if not results.empty:
        results["_status_order"] = results["docking_status"].map({
            "COMPLETED": 0,
            "SKIPPED_NO_VETTED_EXPERIMENTAL_STRUCTURE": 1,
            "FAILED": 2,
        }).fillna(9)

        results = results.sort_values(
            ["_status_order", "best_affinity_kcal_mol"],
            ascending=[True, True],
            na_position="last",
        ).drop(columns="_status_order")

    atomic_csv(results, OUT_RESULTS)

    if results.empty:
        summary = pd.DataFrame()
    else:
        summary = (
            results.groupby("docking_status", dropna=False)
            .agg(
                candidate_count=("pert_iname", "size"),
                median_best_affinity_kcal_mol=(
                    "best_affinity_kcal_mol",
                    "median",
                ),
            )
            .reset_index()
        )

    atomic_csv(summary, OUT_SUMMARY)

    metadata = {
        "stage": "04S",
        "implementation": "v3_rigid_receptor_pdbqt_fix",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input": str(INPUT_SHORTLIST),
        "output": str(OUT_RESULTS),
        "exhaustiveness": args.exhaustiveness,
        "num_modes": args.num_modes,
        "structure_plan": {
            k: {
                **v,
                "keep_hetero": sorted(v["keep_hetero"]),
            }
            for k, v in STRUCTURE_PLAN.items()
        },
        "guardrails": [
            "Docking score is not proof of binding.",
            "Docking score is not proof of efficacy.",
            "Docking score is not proof of reversal of trastuzumab resistance.",
            "Experimental structures are preferred to predicted structures.",
            "CYP3A4 heme is retained.",
            "Binding-site boxes are derived from co-crystal ligand coordinates.",
        ],
        "next_stage": "04T_ADMET_structural_assessment",
    }

    atomic_json(metadata, OUT_META)

    header("STAGE 04S SUMMARY")

    if results.empty:
        print("No docking results were produced.", flush=True)
    else:
        show = [
            c for c in [
                "docking_rank",
                "pert_iname",
                "target_symbol",
                "pdb_id",
                "docking_status",
                "best_affinity_kcal_mol",
                "reference_redocking_rmsd_A",
                "protocol_validation",
                "critical_cofactor_retained",
            ]
            if c in results.columns
        ]
        print(results[show].to_string(index=False), flush=True)

    header("STAGE 04S COMPLETE")
    print("\nOutputs:", flush=True)
    print(f"  {OUT_RESULTS}", flush=True)
    print(f"  {OUT_SUMMARY}", flush=True)
    print(f"  {OUT_META}", flush=True)
    print(
        "\nPer-pair receptor, ligand, pose, and log files are under:",
        flush=True,
    )
    print(f"  {OUT_DIR}", flush=True)
    print("\nNext: 04T — ADMET / structural assessment", flush=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
