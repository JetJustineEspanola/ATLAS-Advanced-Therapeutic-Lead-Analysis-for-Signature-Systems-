#!/usr/bin/env python3

from __future__ import annotations

import io
import json
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[1]

STRICT_CORE = (
    ROOT
    / "results/external_validation/three_dataset_strict_core_genes.csv"
)

DRUG_TARGETS = (
    ROOT
    / "results/cmap/drug_targets/ATLAS_CMap_drug_target_pairs.csv"
)

TARGET_RESISTANCE = (
    ROOT
    / "results/cmap/network_integration/ATLAS_target_resistance_gene_links.csv"
)

DOCKING_RESULTS = (
    ROOT
    / "results/cmap/docking/ATLAS_docking_results.csv"
)

OUT = ROOT / "results/publication_remaining_data"

GO_OUT = OUT / "go_enrichment"
STRING_OUT = OUT / "string_ppi"
NETWORK_OUT = OUT / "drug_target_resistance_network"
DOCKING_OUT = OUT / "docking_interactions"

for p in [OUT, GO_OUT, STRING_OUT, NETWORK_OUT, DOCKING_OUT]:
    p.mkdir(parents=True, exist_ok=True)


# ------------------------------------------------------------
# Utilities
# ------------------------------------------------------------

def require(path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required file missing: {path}")


def save_json(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, default=str))


for p in [STRICT_CORE, DRUG_TARGETS, TARGET_RESISTANCE, DOCKING_RESULTS]:
    require(p)


core = pd.read_csv(STRICT_CORE)

if "gene_symbol" not in core.columns:
    raise RuntimeError("gene_symbol column not found in strict-core file.")

genes = (
    core["gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
)

genes = sorted(set(g for g in genes if g))

up_genes = sorted(
    set(
        core.loc[
            core["replication_direction"].astype(str).eq("UP_IN_RESISTANT"),
            "gene_symbol",
        ]
        .dropna()
        .astype(str)
    )
)

down_genes = sorted(
    set(
        core.loc[
            core["replication_direction"].astype(str).eq("DOWN_IN_RESISTANT"),
            "gene_symbol",
        ]
        .dropna()
        .astype(str)
    )
)

print("=" * 72)
print("ATLAS PUBLICATION DATA COMPLETION")
print("=" * 72)
print(f"Strict-core genes : {len(genes)}")
print(f"UP in resistance  : {len(up_genes)}")
print(f"DOWN in resistance: {len(down_genes)}")


# ------------------------------------------------------------
# 1. Save gene sets
# ------------------------------------------------------------

(OUT / "strict_core_all_genes.txt").write_text("\n".join(genes) + "\n")
(OUT / "strict_core_up_genes.txt").write_text("\n".join(up_genes) + "\n")
(OUT / "strict_core_down_genes.txt").write_text("\n".join(down_genes) + "\n")

core.to_csv(OUT / "strict_core_242_genes.csv", index=False)


# ------------------------------------------------------------
# 2. g:Profiler GO enrichment
# ------------------------------------------------------------

GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"


def run_gprofiler(gene_list, label):
    print(f"\n[g:Profiler] {label}: {len(gene_list)} genes")

    payload = {
        "organism": "hsapiens",
        "query": gene_list,
        "sources": ["GO:BP", "GO:CC", "GO:MF"],
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "no_evidences": False,
    }

    r = requests.post(
        GPROFILER_URL,
        json=payload,
        timeout=180,
    )
    r.raise_for_status()

    data = r.json()

    rows = []

    for x in data.get("result", []):
        intersections = x.get(
            "intersections",
            x.get("intersection", [])
        )

        if isinstance(intersections, list):
            intersection_string = ";".join(
                str(v) for v in intersections
            )
        else:
            intersection_string = str(intersections)

        rows.append(
            {
                "gene_set": label,
                "source": x.get("source"),
                "term_id": x.get("native"),
                "term_name": x.get("name"),
                "adjusted_p_value": x.get("p_value"),
                "significant": x.get("significant"),
                "term_size": x.get("term_size"),
                "query_size": x.get("query_size"),
                "intersection_size": x.get("intersection_size"),
                "effective_domain_size": x.get(
                    "effective_domain_size"
                ),
                "precision": x.get("precision"),
                "recall": x.get("recall"),
                "intersection_genes": intersection_string,
            }
        )

    return pd.DataFrame(rows)


go_frames = []

for label, glist in [
    ("STRICT_CORE_ALL", genes),
    ("UP_IN_RESISTANT", up_genes),
    ("DOWN_IN_RESISTANT", down_genes),
]:
    try:
        go_frames.append(run_gprofiler(glist, label))
    except Exception as exc:
        print(f"WARNING: g:Profiler failed for {label}: {exc}")


if go_frames:
    go = pd.concat(go_frames, ignore_index=True)

    go = go.sort_values(
        ["gene_set", "source", "adjusted_p_value"],
        ascending=[True, True, True],
    )

    go.to_csv(
        GO_OUT / "ATLAS_GO_enrichment_all.csv",
        index=False,
    )

    top_go = (
        go.groupby(["gene_set", "source"], group_keys=False)
        .head(15)
        .copy()
    )

    top_go.to_csv(
        GO_OUT / "ATLAS_GO_enrichment_top15.csv",
        index=False,
    )

    for source in ["GO:BP", "GO:CC", "GO:MF"]:
        subset = go[go["source"] == source]
        subset.to_csv(
            GO_OUT
            / f"ATLAS_{source.replace(':', '_')}_enrichment.csv",
            index=False,
        )

    print(f"[OK] GO enrichment: {len(go)} terms")
else:
    print("[WARNING] No GO enrichment data generated.")


# ------------------------------------------------------------
# 3. STRING PPI network
# ------------------------------------------------------------

STRING_URL = "https://string-db.org/api/tsv/network"

print("\n[STRING] Requesting high-confidence human PPI network...")

string_payload = {
    "identifiers": "\r".join(genes),
    "species": 9606,
    "required_score": 700,
    "network_type": "functional",
    "caller_identity": "ATLAS_publication_pipeline",
}

try:
    r = requests.post(
        STRING_URL,
        data=string_payload,
        timeout=300,
    )
    r.raise_for_status()

    ppi = pd.read_csv(
        io.StringIO(r.text),
        sep="\t",
    )

    ppi.to_csv(
        STRING_OUT / "ATLAS_STRING_PPI_edges.csv",
        index=False,
    )

    # Build node-degree statistics
    if {
        "preferredName_A",
        "preferredName_B",
        "score",
    }.issubset(ppi.columns):

        a = ppi[["preferredName_A", "score"]].rename(
            columns={"preferredName_A": "gene"}
        )

        b = ppi[["preferredName_B", "score"]].rename(
            columns={"preferredName_B": "gene"}
        )

        both = pd.concat([a, b], ignore_index=True)

        degree = (
            both.groupby("gene")
            .agg(
                degree=("gene", "size"),
                mean_STRING_score=("score", "mean"),
                sum_STRING_score=("score", "sum"),
            )
            .reset_index()
        )

        all_nodes = pd.DataFrame({"gene": genes})

        nodes = all_nodes.merge(
            degree,
            how="left",
            on="gene",
        )

        nodes["degree"] = nodes["degree"].fillna(0).astype(int)
        nodes["mean_STRING_score"] = (
            nodes["mean_STRING_score"].fillna(0)
        )
        nodes["sum_STRING_score"] = (
            nodes["sum_STRING_score"].fillna(0)
        )

        direction_map = (
            core.drop_duplicates("gene_symbol")
            .set_index("gene_symbol")["replication_direction"]
            .to_dict()
        )

        nodes["direction"] = nodes["gene"].map(direction_map)

        nodes = nodes.sort_values(
            ["degree", "sum_STRING_score"],
            ascending=False,
        )

        nodes.to_csv(
            STRING_OUT / "ATLAS_STRING_PPI_nodes.csv",
            index=False,
        )

        nodes.head(25).to_csv(
            STRING_OUT
            / "ATLAS_STRING_top25_network_central_genes.csv",
            index=False,
        )

        print(
            f"[OK] STRING: {len(ppi)} edges; "
            f"{(nodes.degree > 0).sum()} connected genes"
        )

    else:
        print(
            "WARNING: STRING response did not contain "
            "expected columns."
        )

except Exception as exc:
    print(f"WARNING: STRING request failed: {exc}")


# ------------------------------------------------------------
# 4. Publication-ready compound → target → resistance edges
# ------------------------------------------------------------

targets = pd.read_csv(DRUG_TARGETS)
links = pd.read_csv(TARGET_RESISTANCE)

compound_target = (
    targets[
        [
            "pert_iname",
            "target_pref_name",
            "target_accessions",
        ]
    ]
    .dropna(subset=["pert_iname"])
    .drop_duplicates()
)

compound_target["source"] = compound_target["pert_iname"]
compound_target["target"] = compound_target["target_accessions"]
compound_target["edge_type"] = "COMPOUND_TARGET"

ct_edges = compound_target[
    [
        "source",
        "target",
        "edge_type",
        "pert_iname",
        "target_pref_name",
        "target_accessions",
    ]
].copy()


# Existing target-resistance links
tr_edges = links.copy()

tr_edges["source"] = tr_edges["target_symbol"]
tr_edges["target"] = tr_edges["resistance_gene"]
tr_edges["edge_type"] = "TARGET_RESISTANCE_NETWORK"

tr_edges.to_csv(
    NETWORK_OUT
    / "ATLAS_target_resistance_edges_publication.csv",
    index=False,
)

ct_edges.to_csv(
    NETWORK_OUT
    / "ATLAS_compound_target_edges_publication.csv",
    index=False,
)

print(
    f"\n[OK] Existing drug-network data exported: "
    f"{len(ct_edges)} compound-target edges and "
    f"{len(tr_edges)} target-resistance edges."
)


# ------------------------------------------------------------
# 5. Docking complex creation
# ------------------------------------------------------------

docking = pd.read_csv(DOCKING_RESULTS)


def find_workdir(compound, target, pdb_id):
    row = docking[
        docking["pert_iname"].astype(str).str.lower().eq(
            compound.lower()
        )
        & docking["target_symbol"].astype(str).eq(target)
        & docking["pdb_id"].astype(str).eq(pdb_id)
    ]

    if not row.empty:
        raw = row.iloc[0].get("workdir")

        if pd.notna(raw):
            p = Path(str(raw))
            if p.exists():
                return p

    patterns = [
        f"*_{compound}_{target}_{pdb_id}",
        f"*{compound}*{target}*{pdb_id}*",
    ]

    docking_root = ROOT / "results/cmap/docking"

    for pattern in patterns:
        hits = sorted(docking_root.glob(pattern))
        if hits:
            return hits[0]

    raise FileNotFoundError(
        f"Could not locate docking workdir for "
        f"{compound}/{target}/{pdb_id}"
    )


def extract_first_vina_pose(src: Path, dst: Path):
    text = src.read_text().splitlines()

    if not any(x.startswith("MODEL") for x in text):
        dst.write_text("\n".join(text) + "\n")
        return

    output = []
    inside = False

    for line in text:
        if line.startswith("MODEL"):
            if inside:
                break
            inside = True
            output.append(line)
            continue

        if inside:
            output.append(line)

            if line.startswith("ENDMDL"):
                break

    dst.write_text("\n".join(output) + "\n")


def convert_ligand_pdbqt_to_pdb(src: Path, dst: Path):
    obabel = shutil.which("obabel")

    if not obabel:
        raise RuntimeError(
            "Open Babel command 'obabel' not found. "
            "Install it with: sudo pacman -S openbabel"
        )

    subprocess.run(
        [
            obabel,
            "-ipdbqt",
            str(src),
            "-opdb",
            "-O",
            str(dst),
        ],
        check=True,
    )


def make_complex(receptor: Path, ligand: Path, complex_out: Path):
    receptor_lines = receptor.read_text().splitlines()
    ligand_lines = ligand.read_text().splitlines()

    max_serial = 0
    clean_receptor = []

    for line in receptor_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                max_serial = max(
                    max_serial,
                    int(line[6:11].strip())
                )
            except Exception:
                pass

        if line.startswith("END"):
            continue

        clean_receptor.append(line)

    ligand_atoms = []
    serial = max_serial + 1

    for line in ligand_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue

        line = line.ljust(80)

        atom_name = line[12:16]
        altloc = line[16:17]
        xyz_rest = line[26:]

        # Force ligand residue identity to LIG chain Z residue 1.
        newline = (
            f"HETATM{serial:5d} "
            f"{atom_name}{altloc}"
            f"LIG Z"
            f"{1:4d}"
            f"{xyz_rest}"
        )

        ligand_atoms.append(newline[:80])
        serial += 1

    complex_out.write_text(
        "\n".join(
            clean_receptor
            + ["TER"]
            + ligand_atoms
            + ["END"]
        )
        + "\n"
    )


DOCKING_CASES = [
    ("sitagliptin", "DPP4", "1X70"),
    ("clofibrate", "CYP3A4", "3NXU"),
]


prepared_complexes = []

for compound, target, pdb_id in DOCKING_CASES:
    print(
        f"\n[Docking interactions] "
        f"{compound} → {target} ({pdb_id})"
    )

    try:
        workdir = find_workdir(
            compound,
            target,
            pdb_id,
        )

        receptor = workdir / "receptor_clean.pdb"
        pose = workdir / "candidate_docked.pdbqt"

        if not receptor.exists():
            raise FileNotFoundError(receptor)

        if not pose.exists():
            raise FileNotFoundError(pose)

        case_out = (
            DOCKING_OUT
            / f"{compound}_{target}_{pdb_id}"
        )

        case_out.mkdir(parents=True, exist_ok=True)

        first_pose = case_out / "best_pose.pdbqt"
        ligand_pdb = case_out / "best_pose_ligand.pdb"
        complex_pdb = case_out / "complex_best_pose.pdb"

        extract_first_vina_pose(pose, first_pose)

        convert_ligand_pdbqt_to_pdb(
            first_pose,
            ligand_pdb,
        )

        make_complex(
            receptor,
            ligand_pdb,
            complex_pdb,
        )

        prepared_complexes.append(
            {
                "compound": compound,
                "target": target,
                "pdb_id": pdb_id,
                "complex": str(complex_pdb),
                "output_dir": str(case_out),
            }
        )

        print(f"[OK] Complex: {complex_pdb}")

    except Exception as exc:
        print(
            f"WARNING: Could not prepare "
            f"{compound}/{target}: {exc}"
        )


# ------------------------------------------------------------
# 6. PLIP residue interaction analysis
# ------------------------------------------------------------

def plip_command():
    cli = shutil.which("plip")

    if cli:
        return [cli]

    try:
        import plip  # noqa: F401
        return [
            sys.executable,
            "-m",
            "plip.plipcmd",
        ]
    except Exception:
        return None


def flatten_plip_xml(xml_file, compound, target, pdb_id):
    root = ET.parse(xml_file).getroot()

    records = []

    for site in root.findall(".//bindingsite"):
        identifiers = site.find("identifiers")

        if identifiers is None:
            continue

        hetid = identifiers.findtext("hetid")
        chain = identifiers.findtext("chain")
        position = identifiers.findtext("position")

        # Only analyze our docked ligand.
        if hetid != "LIG":
            continue

        interactions = site.find("interactions")

        if interactions is None:
            continue

        for interaction_group in list(interactions):
            interaction_type = interaction_group.tag

            for interaction in list(interaction_group):
                record = {
                    "compound": compound,
                    "target": target,
                    "pdb_id": pdb_id,
                    "ligand_hetid": hetid,
                    "ligand_chain": chain,
                    "ligand_position": position,
                    "interaction_type": interaction_type,
                }

                for child in list(interaction):
                    if len(list(child)) == 0:
                        record[child.tag] = child.text

                records.append(record)

    return records


plip_base = plip_command()
all_interactions = []

if plip_base:
    print("\n[PLIP] PLIP detected. Running residue analysis...")

    for case in prepared_complexes:
        complex_pdb = Path(case["complex"])

        plip_out = (
            Path(case["output_dir"])
            / "plip"
        )

        plip_out.mkdir(parents=True, exist_ok=True)

        command = (
            plip_base
            + [
                "-f",
                str(complex_pdb),
                "-x",
                "-t",
                "-o",
                str(plip_out),
            ]
        )

        print("Running:", " ".join(command))

        try:
            subprocess.run(
                command,
                check=True,
            )

            xml_candidates = list(
                plip_out.rglob("report.xml")
            )

            if not xml_candidates:
                xml_candidates = list(
                    plip_out.rglob("*.xml")
                )

            if xml_candidates:
                records = flatten_plip_xml(
                    xml_candidates[0],
                    case["compound"],
                    case["target"],
                    case["pdb_id"],
                )

                all_interactions.extend(records)

                if records:
                    pd.DataFrame(records).to_csv(
                        Path(case["output_dir"])
                        / "PLIP_interactions.csv",
                        index=False,
                    )

                print(
                    f"[OK] PLIP: {len(records)} "
                    f"interaction records"
                )
            else:
                print(
                    "WARNING: PLIP completed but "
                    "no XML report was found."
                )

        except Exception as exc:
            print(
                f"WARNING: PLIP failed for "
                f"{case['compound']}: {exc}"
            )

else:
    status = """
PLIP was not detected.

The docking complexes HAVE been prepared successfully.

Install PLIP and rerun this script to generate:
- hydrogen bonds
- hydrophobic contacts
- salt bridges
- pi interactions
- interacting residue numbers
- interaction distances

Try:

    python -m pip install plip

If PLIP installation fails because of Open Babel,
the generated complex_best_pose.pdb files can also
be analyzed using the PLIP Docker image.
""".strip()

    (DOCKING_OUT / "PLIP_NOT_INSTALLED.txt").write_text(
        status + "\n"
    )

    print("\n[WARNING] PLIP is not installed.")
    print(
        "Complex PDB files were still prepared. "
        "Install PLIP and rerun."
    )


if all_interactions:
    pd.DataFrame(all_interactions).to_csv(
        DOCKING_OUT
        / "ATLAS_PLIP_all_docking_interactions.csv",
        index=False,
    )


# ------------------------------------------------------------
# 7. Final manifest
# ------------------------------------------------------------

manifest = {
    "strict_core_gene_count": len(genes),
    "up_in_resistant": len(up_genes),
    "down_in_resistant": len(down_genes),
    "string_required_score": 700,
    "string_species": 9606,
    "go_sources": [
        "GO:BP",
        "GO:CC",
        "GO:MF",
    ],
    "docking_cases": [
        {
            "compound": c,
            "target": t,
            "pdb_id": p,
        }
        for c, t, p in DOCKING_CASES
    ],
    "prepared_complexes": prepared_complexes,
    "plip_interaction_records": len(all_interactions),
}

save_json(
    manifest,
    OUT / "publication_remaining_data_manifest.json",
)


print("\n" + "=" * 72)
print("DONE")
print("=" * 72)
print(f"Output directory:\n{OUT}")
print("\nExpected outputs:")
print("  go_enrichment/")
print("  string_ppi/")
print("  drug_target_resistance_network/")
print("  docking_interactions/")
print("  publication_remaining_data_manifest.json")
