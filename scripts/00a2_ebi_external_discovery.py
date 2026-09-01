#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, re, time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
DISC = ROOT / "data" / "discovery"
DISC.mkdir(parents=True, exist_ok=True)
EXISTING = DISC / "dataset_candidates.csv"
OUT = DISC / "ebi_external_candidates.csv"
XREF = DISC / "ebi_cross_reference_map.csv"

BASE = "https://www.ebi.ac.uk/biostudies/api/v1"
UA = "ATLAS-EBI-discovery/1.0"
QUERIES = [
    '"HER2" AND trastuzumab AND resistance',
    'ERBB2 AND trastuzumab AND resistance',
    '"trastuzumab resistant" AND breast',
    'BT474 AND trastuzumab',
    'SKBR3 AND trastuzumab',
    '"anti-HER2" AND resistance AND breast',
    '"HER2-positive" AND breast AND resistance',
    'TGF-beta AND HER2 AND resistance',
]
PATS = {
    "GEO": re.compile(r"\bGSE\d+\b", re.I),
    "ENA_STUDY": re.compile(r"\b(?:ERP|SRP|DRP)\d+\b", re.I),
    "ENA_PROJECT": re.compile(r"\bPRJ(?:EB|NA|DB)\d+\b", re.I),
    "BIOSAMPLE": re.compile(r"\b(?:SAMEA|SAMN|SAMD)\d+\b", re.I),
    "GXA": re.compile(r"\bE-(?:MTAB|GEOD|TABM|MEXP)-\d+\b", re.I),
}

def getj(url):
    req = Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urlopen(req, timeout=45) as r:
        return json.load(r)

def search(q):
    params = urlencode({"query": q, "page": 1, "pageSize": 100})
    for u in (f"{BASE}/arrayexpress/search?{params}", f"{BASE}/search?{params}"):
        try:
            d = getj(u)
            hits = d.get("hits") or d.get("studies") or d.get("content") or d.get("results") or []
            if "/search?" in u:
                hits = [h for h in hits if str(h.get("accession") or h.get("accno") or h.get("id") or "").upper().startswith("E-")]
            return hits if isinstance(hits, list) else []
        except Exception:
            pass
    return []

def acc(h):
    return str(h.get("accession") or h.get("accno") or h.get("id") or "").strip()

def title(h):
    return str(h.get("title") or h.get("name") or h.get("description") or "").strip()

def flatten(x):
    parts=[]
    def w(v):
        if isinstance(v, dict):
            for z in v.values(): w(z)
        elif isinstance(v, list):
            for z in v: w(z)
        elif isinstance(v, (str,int,float)):
            parts.append(str(v))
    w(x)
    return "\n".join(parts)

def score(text):
    t=text.lower(); s=0
    for pts,needle in [(25,"trastuzumab"),(15,"her2"),(15,"erbb2"),(15,"resistan"),(10,"breast"),(8,"bt474"),(8,"skbr3"),(5,"rna-seq"),(5,"transcript")]:
        if needle in t: s += pts
    return min(s,100)

def discover():
    found={}
    for q in QUERIES:
        print(f"[BioStudies/ArrayExpress] {q}", flush=True)
        hs=search(q)
        print(f"  hits: {len(hs)}", flush=True)
        for h in hs:
            a=acc(h)
            if not a: continue
            found.setdefault(a, {"hit":h,"queries":set()})["queries"].add(q)
        time.sleep(.2)

    rows=[]; xrows=[]
    for i,(a,obj) in enumerate(sorted(found.items()),1):
        print(f"[detail {i}/{len(found)}] {a}", flush=True)
        try: detail=getj(f"{BASE}/studies/{a}")
        except Exception: detail={}
        try: info=getj(f"{BASE}/studies/{a}/info")
        except Exception: info={}
        txt=flatten({"hit":obj["hit"],"detail":detail,"info":info})
        refs={k:sorted({m.upper() for m in p.findall(txt)}) for k,p in PATS.items()}
        r={
            "source":"BIOSTUDIES_ARRAYEXPRESS",
            "source_accession":a,
            "title":title(obj["hit"]) or str(detail.get("title") or ""),
            "discovery_score":score(txt),
            "query_hits":" || ".join(sorted(obj["queries"])),
            "geo_accessions":"|".join(refs["GEO"]),
            "ena_study_accessions":"|".join(refs["ENA_STUDY"]),
            "ena_project_accessions":"|".join(refs["ENA_PROJECT"]),
            "biosample_accessions":"|".join(refs["BIOSAMPLE"]),
            "expression_atlas_accessions":"|".join(refs["GXA"]),
            "biostudies_url":f"https://www.ebi.ac.uk/biostudies/studies/{a}",
            "ftp_link":str(info.get("ftpLink") or ""),
            "relationship_role":"EXTERNAL_DISCOVERY_LEAD",
            "independence_note":"Cross-database mirrors must not be counted as independent validation datasets."
        }
        rows.append(r)
        for k,vals in refs.items():
            for v in vals:
                xrows.append({"biostudies_accession":a,"xref_type":k,"xref_accession":v})
        time.sleep(.1)
    return sorted(rows,key=lambda r:(-r["discovery_score"],r["source_accession"])),xrows

def write(rows,xrows):
    fields=list(rows[0].keys()) if rows else ["source","source_accession","title"]
    with OUT.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=fields); w.writeheader(); w.writerows(rows)
    with XREF.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=["biostudies_accession","xref_type","xref_accession"]); w.writeheader(); w.writerows(xrows)

def merge(rows):
    if not EXISTING.exists():
        print(f"MERGE SKIPPED: {EXISTING} missing"); return
    with EXISTING.open(newline="",encoding="utf-8") as f:
        rd=csv.DictReader(f); old=list(rd); cols=rd.fieldnames or []
    keys={(str(r.get("source","")).upper(),str(r.get("source_accession","")).upper()) for r in old}
    n=0
    for r in rows:
        k=(r["source"].upper(),r["source_accession"].upper())
        if k in keys: continue
        nr={c:"" for c in cols}
        mp={
            "dataset_id":"BIOSTUDIES:" + r["source_accession"],
            "source":r["source"],"source_accession":r["source_accession"],"title":r["title"],
            "url":r["biostudies_url"],"dataset_url":r["biostudies_url"],"landing_url":r["biostudies_url"],
            "query_group":"ebi_external","query":r["query_hits"],"eligibility_score":r["discovery_score"],
            "relationship_role":r["relationship_role"],"relationship_reason":r["independence_note"]
        }
        for kk,v in mp.items():
            if kk in nr: nr[kk]=v
        old.append(nr); keys.add(k); n+=1
    tmp=EXISTING.with_suffix(".csv.tmp")
    with tmp.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(old)
    tmp.replace(EXISTING)
    print(f"Merged {n} new BioStudies/ArrayExpress leads into {EXISTING}")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--merge-existing",action="store_true")
    args=ap.parse_args()
    print("="*88)
    print("ATLAS — 00A2 EMBL-EBI / ARRAYEXPRESS EXTERNAL DISCOVERY")
    print("="*88)
    rows,xrows=discover()
    write(rows,xrows)
    print(f"Unique external leads: {len(rows)}")
    print(f"Cross-reference rows: {len(xrows)}")
    print(OUT); print(XREF)
    if args.merge_existing: merge(rows)
    print("SUCCESS [00a2] External discovery complete")

if __name__=="__main__":
    main()
