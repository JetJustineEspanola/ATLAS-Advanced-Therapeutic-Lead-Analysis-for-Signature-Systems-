#!/usr/bin/env python3
"""
ATLAS — Stage 00D1 v2
Phenotype Metadata Audit with token-aware hint detection.

Fixes false positives from substring matching such as:
- "res" inside unrelated words
- "tr" inside longer tokens
- "hr" inside longer tokens
- "wt" inside longer tokens

Also writes a compact distinct-label file to make study-specific phenotype
mapping easier.

Outputs
-------
data/enriched/phenotype_audit_summary_v2.csv
data/enriched/phenotype_audit_samples_v2.csv
data/enriched/phenotype_distinct_labels.csv
data/enriched/phenotype_manual_mapping_template.csv
"""

from __future__ import annotations

from pathlib import Path
import re
import sys
from collections import Counter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import pandas as pd

INFILE = PROJECT_ROOT / "data" / "enriched" / "sample_metadata.csv"
OUTDIR = PROJECT_ROOT / "data" / "enriched"
OUTDIR.mkdir(parents=True, exist_ok=True)

STOPWORDS = {
    "sample","samples","cell","cells","rna","seq","rna-seq","rnaseq","replicate",
    "rep","control","treated","treatment","experiment","library","human","breast",
    "cancer","line","with","from","and","the","for","of","in","to","a","an"
}

# Long phrases use normal substring matching.
RESISTANCE_PHRASES = [
    "trastuzumab-resistant",
    "trastuzumab resistant",
    "herceptin-resistant",
    "herceptin resistant",
    "acquired resistance",
    "drug resistant",
    "resistant clone",
]
CONTROL_PHRASES = [
    "trastuzumab-sensitive",
    "trastuzumab sensitive",
    "parental",
    "wild-type",
    "wild type",
    "untreated",
    "vehicle",
    "naive",
    "naïve",
]

# Short labels must match token boundaries exactly.
RESISTANCE_TOKENS = {
    "resistant", "resistance", "refractory", "adapted",
    "bt474r", "skbr3r", "res",
}
CONTROL_TOKENS = {
    "sensitive", "parental", "control", "vehicle",
    "wt", "naive", "naïve",
}

# These are ambiguous and should only be surfaced, never treated as phenotype
# evidence by themselves.
AMBIGUOUS_SHORT_TOKENS = {"tr", "hr"}


def clean(x):
    if x is None:
        return ""
    try:
        if pd.isna(x):
            return ""
    except Exception:
        pass
    return str(x).strip()


def combined_text(r):
    fields = [
        "title",
        "source_name",
        "characteristics",
        "treatment",
        "description",
        "cell_line",
        "biological_group",
    ]
    return " | ".join(clean(r.get(c)) for c in fields if clean(r.get(c)))


def lexical_tokens(text: str) -> list[str]:
    # Preserve useful labels like BT474-R, T15, SKBR3R.
    return re.findall(r"[A-Za-z0-9]+(?:[-_/][A-Za-z0-9]+)*", text.lower())


def top_tokens(text):
    toks = lexical_tokens(text)
    return [t for t in toks if t not in STOPWORDS and not t.isdigit()]


def exact_token_hits(text: str, wanted: set[str]) -> list[str]:
    toks = set(lexical_tokens(text))
    return sorted(toks & wanted)


def phrase_hits(text: str, phrases: list[str]) -> list[str]:
    low = text.lower()
    return sorted({p for p in phrases if p in low})


def main():
    if not INFILE.exists():
        print(f"ERROR: missing {INFILE}")
        return 1

    df = pd.read_csv(INFILE)

    summary_rows = []
    preview_rows = []
    distinct_rows = []
    manual_rows = []

    for dataset_id, g in df.groupby("dataset_id", dropna=False):
        texts = [combined_text(r) for _, r in g.iterrows()]

        counter = Counter()
        res_hits = Counter()
        ctl_hits = Counter()
        ambiguous_hits = Counter()

        for t in texts:
            counter.update(top_tokens(t))

            for h in phrase_hits(t, RESISTANCE_PHRASES):
                res_hits[h] += 1
            for h in exact_token_hits(t, RESISTANCE_TOKENS):
                res_hits[h] += 1

            for h in phrase_hits(t, CONTROL_PHRASES):
                ctl_hits[h] += 1
            for h in exact_token_hits(t, CONTROL_TOKENS):
                ctl_hits[h] += 1

            for h in exact_token_hits(t, AMBIGUOUS_SHORT_TOKENS):
                ambiguous_hits[h] += 1

        resolved = g.get(
            "resistance_status",
            pd.Series(["UNRESOLVED"] * len(g))
        ).value_counts(dropna=False)

        summary_rows.append({
            "dataset_id": dataset_id,
            "sample_n": len(g),
            "currently_resistant_n": int(resolved.get("RESISTANT", 0)),
            "currently_sensitive_parental_n": int(
                resolved.get("SENSITIVE_OR_PARENTAL", 0)
            ),
            "currently_ambiguous_n": int(resolved.get("AMBIGUOUS", 0)),
            "currently_unresolved_n": int(resolved.get("UNRESOLVED", 0)),
            "resistance_hint_counts": " | ".join(
                f"{k}:{v}" for k, v in res_hits.most_common()
            ),
            "control_hint_counts": " | ".join(
                f"{k}:{v}" for k, v in ctl_hits.most_common()
            ),
            "ambiguous_short_token_counts": " | ".join(
                f"{k}:{v}" for k, v in ambiguous_hits.most_common()
            ),
            "top_metadata_tokens": " | ".join(
                f"{tok}:{n}" for tok, n in counter.most_common(25)
            ),
        })

        # all samples for inspectability
        for _, r in g.iterrows():
            row_text = combined_text(r)
            preview_rows.append({
                "dataset_id": dataset_id,
                "sample_id": clean(r.get("sample_id")),
                "title": clean(r.get("title")),
                "source_name": clean(r.get("source_name")),
                "characteristics": clean(r.get("characteristics")),
                "treatment": clean(r.get("treatment")),
                "description": clean(r.get("description")),
                "current_resistance_status": clean(r.get("resistance_status")),
                "cell_line": clean(r.get("cell_line")),
                "replicate_type": clean(r.get("replicate_type")),
                "resistance_hints": " | ".join(
                    phrase_hits(row_text, RESISTANCE_PHRASES)
                    + exact_token_hits(row_text, RESISTANCE_TOKENS)
                ),
                "control_hints": " | ".join(
                    phrase_hits(row_text, CONTROL_PHRASES)
                    + exact_token_hits(row_text, CONTROL_TOKENS)
                ),
                "ambiguous_short_tokens": " | ".join(
                    exact_token_hits(row_text, AMBIGUOUS_SHORT_TOKENS)
                ),
            })

        # Distinct labels are more useful than scanning 565 raw rows.
        label_cols = ["title", "source_name", "characteristics", "treatment"]
        seen = Counter()
        examples = {}
        for _, r in g.iterrows():
            label = " | ".join(
                f"{c}={clean(r.get(c))}" for c in label_cols if clean(r.get(c))
            )
            if label:
                seen[label] += 1
                examples[label] = clean(r.get("sample_id"))

        for label, n in seen.most_common(100):
            distinct_rows.append({
                "dataset_id": dataset_id,
                "sample_count_with_label": n,
                "example_sample_id": examples[label],
                "metadata_label": label,
            })

        manual_rows.append({
            "dataset_id": dataset_id,
            "include_for_validation": "",
            "resistant_regex": "",
            "sensitive_parental_regex": "",
            "exclude_regex": "",
            "cell_line_regex": "",
            "replicate_group_regex": "",
            "notes": "",
        })

    summary = pd.DataFrame(summary_rows).sort_values(
        ["currently_unresolved_n", "sample_n"],
        ascending=[False, False],
    )
    previews = pd.DataFrame(preview_rows)
    distinct = pd.DataFrame(distinct_rows)
    manual = pd.DataFrame(manual_rows)

    p1 = OUTDIR / "phenotype_audit_summary_v2.csv"
    p2 = OUTDIR / "phenotype_audit_samples_v2.csv"
    p3 = OUTDIR / "phenotype_distinct_labels.csv"
    p4 = OUTDIR / "phenotype_manual_mapping_template.csv"

    summary.to_csv(p1, index=False)
    previews.to_csv(p2, index=False)
    distinct.to_csv(p3, index=False)
    manual.to_csv(p4, index=False)

    print("=" * 78)
    print("ATLAS — 00D1 v2 PHENOTYPE METADATA AUDIT")
    print("=" * 78)
    print(f"Datasets audited: {len(summary)}")
    print(f"Samples audited: {len(df)}")

    print("\nAudit summary:")
    cols = [
        "dataset_id",
        "sample_n",
        "currently_resistant_n",
        "currently_sensitive_parental_n",
        "currently_unresolved_n",
        "resistance_hint_counts",
        "control_hint_counts",
        "ambiguous_short_token_counts",
    ]
    print(summary[cols].head(20).to_string(index=False))

    print("\nOutputs:")
    print(p1)
    print(p2)
    print(p3)
    print(p4)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
