"""Fetch measured binding affinities from BindingDB and report the net-new yield per target.

BindingDB is an independent measured-affinity database (curated from papers, patents and PubChem).
It overlaps ChEMBL heavily for well-studied targets, so the only honest way to state the gain is to
standardise both to InChIKey and count how many BindingDB compounds are *not already* in our
ChEMBL-derived endpoint sets. This script fetches, standardises, and writes that comparison. It does
not modify the training data; rebuilding the endpoints is a separate, deliberate step taken only if
the net-new yield is worthwhile.

Affinities (nM) are converted to the same -log10(molar) potency scale as pChEMBL and labelled with
the identical rule (>=6 active, <5 inactive, 5-6 grey zone dropped).

Output: results/tables/bindingdb_yield.csv
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import certifi
import pandas as pd
import requests
import urllib3

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_compound_library import standardise  # noqa: E402

urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[3]
CACHE = ROOT / "data" / "_bindingdb_cache"
CACHE.mkdir(parents=True, exist_ok=True)

# Human UniProt accessions for the eleven ChEMBL targets.
UNIPROT = {
    "AChE": "P22303", "BChE": "P06276", "BACE1": "P56817", "GSK3B": "P49841",
    "MAO_B": "P27338", "MAO_A": "P21397", "D2": "P14416", "A2A": "P29274",
    "HT2A": "P28223", "SERT": "P31645", "hERG": "Q12809",
}
CUTOFF_NM = 1_000_000  # 1 mM ceiling: captures weak binders and inactives (the server rejects 10 mM)


_HEADERS = {"User-Agent": "Mozilla/5.0 (research; BrainSafe QSAR data collection)"}


def _get(url: str) -> requests.Response:
    """GET with the machine's CA bundle, falling back to unverified for hosts not in the store."""
    try:
        return requests.get(url, timeout=120, verify=certifi.where(), headers=_HEADERS)
    except requests.exceptions.SSLError:
        return requests.get(url, timeout=120, verify=False, headers=_HEADERS)


def fetch_target(name: str, uniprot: str) -> list[dict]:
    cache = CACHE / f"{name}.json"
    if cache.exists():
        rows = json.loads(cache.read_text())
        print(f"  [{name}] cache: {len(rows)} affinities")
        return rows
    url = f"https://bindingdb.org/rest/getLigandsByUniprot?uniprot={uniprot}&cutoff={CUTOFF_NM}"
    rows, hit = [], None
    for attempt in range(5):
        r = _get(url)
        if r.status_code == 200 and r.content:
            try:
                resp = r.json().get("getLindsByUniprotResponse", {})
                rows = resp.get("bdb.affinities", []) or []
                hit = resp.get("bdb.hit")
                break
            except ValueError:
                pass
        wait = 5 * (2 ** attempt)  # exponential backoff for rate-limiting: 5,10,20,40,80s
        print(f"  [{name}] attempt {attempt + 1}: empty/non-JSON, waiting {wait}s")
        time.sleep(wait)
    if not rows:
        print(f"  [{name}] no data after retries (not cached, will refetch next run)")
        return rows
    cache.write_text(json.dumps(rows))
    print(f"  [{name}] fetched: {len(rows)} affinities (hit={hit})")
    return rows


def parse_affinity(raw: str) -> float | None:
    """Convert a BindingDB affinity string in nM to a -log10(molar) potency value."""
    s = str(raw).strip().lstrip("><~=").strip()
    try:
        v = float(s)
    except ValueError:
        return None
    if v <= 0:
        return None
    return 9.0 - math.log10(v)  # nM -> pX


def label_from(pvalue: float) -> int:
    return 1 if pvalue >= 6.0 else (0 if pvalue < 5.0 else -1)


def build_bindingdb_compounds(rows: list[dict]) -> pd.DataFrame:
    """Standardise and median-aggregate BindingDB rows into labelled compounds keyed by InChIKey."""
    recs = []
    for a in rows:
        p = parse_affinity(a.get("bdb.affinity"))
        if p is None:
            continue
        csmi, ik = standardise(a.get("bdb.smile"))
        if ik is None:
            continue
        recs.append({"inchikey": ik, "smiles": csmi, "pvalue": p})
    if not recs:
        return pd.DataFrame(columns=["inchikey", "smiles", "pchembl", "label"])
    df = pd.DataFrame(recs)
    g = df.groupby("inchikey").agg(smiles=("smiles", "first"),
                                   pchembl=("pvalue", "median")).reset_index()
    g["label"] = g["pchembl"].apply(label_from)
    return g[g["label"] >= 0].reset_index(drop=True)


def current_inchikeys(endpoint: str) -> set[str]:
    path = ROOT / "data" / "endpoints" / f"{endpoint}.csv"
    keys = set()
    if path.exists():
        for smi in pd.read_csv(path)["smiles"]:
            _, ik = standardise(smi)
            if ik:
                keys.add(ik)
    return keys


def main() -> None:
    summary = []
    for i, (name, uniprot) in enumerate(UNIPROT.items()):
        if i:
            time.sleep(6)  # be gentle with the BindingDB server between targets
        rows = fetch_target(name, uniprot)
        bdb = build_bindingdb_compounds(rows)
        have = current_inchikeys(name)
        new = bdb[~bdb.inchikey.isin(have)]
        summary.append({
            "endpoint": name,
            "current_chembl_compounds": len(have),
            "bindingdb_labelled": len(bdb),
            "net_new_compounds": len(new),
            "net_new_active": int((new.label == 1).sum()),
            "net_new_inactive": int((new.label == 0).sum()),
            "combined_total": len(have) + len(new),
        })
        # persist the standardised BindingDB compounds for the later rebuild step
        bdb.to_csv(CACHE / f"{name}_labelled.csv", index=False)
        print(f"  [{name}] current {len(have)} + net-new {len(new)} "
              f"-> combined {len(have)+len(new)}")
    out = pd.DataFrame(summary)
    (ROOT / "results" / "tables").mkdir(parents=True, exist_ok=True)
    out.to_csv(ROOT / "results" / "tables" / "bindingdb_yield.csv", index=False)
    print("\n=== BindingDB net-new yield ===")
    print(out.to_string(index=False))
    print(f"\nTOTAL net-new compounds: {out.net_new_compounds.sum():,}")
    print("wrote results/tables/bindingdb_yield.csv")


if __name__ == "__main__":
    main()
