"""Audit the MeSH-to-panel mapping used by H6.

The mapping is the only subjective choice in the clinical-indication test, so it is worth showing
exactly what it did rather than asserting that it is reasonable. This lists every distinct MeSH
heading that matched each panel condition, and every discarded heading that mentions the brain or
nervous system, which is where a mapping error would hide.

Read-only. Writes inversion/results/H6_mesh_mapping_audit.csv
"""
from __future__ import annotations

import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
import requests
import urllib3

urllib3.disable_warnings()
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from inv_clinical_indication import MESH_MAP, MESH_EXCLUDE, MIN_PHASE, BASE  # noqa: E402

OUT = ROOT / "inversion" / "results"

# words that mark a heading as plausibly neurological or psychiatric, used only to decide what to
# show as discarded; they play no part in the mapping itself
# "coma" is matched as a whole word only, otherwise every sarcoma is listed as neurological
NEURO_WORDS = ["brain", "cerebr", "neuro", "nerv", "psych", "cognit", "memory", "mental",
               "spinal", "sclerosis", "tremor", "dystonia", "ataxia", "encephal", "myelin",
               "stroke", "delirium", "autis", "tourette", "restless legs"]
NEURO_EXACT = ["coma", "seizure", "convulsion"]


def main():
    matched = defaultdict(Counter)
    discarded = Counter()
    url = f"{BASE}/drug_indication.json?limit=1000"
    page = 0
    while url:
        r = requests.get(url, timeout=90, verify=False)
        if r.status_code != 200:
            break
        j = r.json()
        for d in j.get("drug_indications", []):
            try:
                if float(d.get("max_phase_for_ind") or 0) < MIN_PHASE:
                    continue
            except (TypeError, ValueError):
                continue
            head = d.get("mesh_heading") or ""
            low = head.lower()
            if any(x in low for x in MESH_EXCLUDE):
                discarded[head] += 1
                continue
            hit = next((v for k, v in MESH_MAP.items() if k in low), None)
            if hit:
                matched[hit][head] += 1
            elif any(w in low for w in NEURO_WORDS) or \
                    any(w in low.replace(",", " ").split() for w in NEURO_EXACT):
                discarded[head] += 1
        page += 1
        nxt = j["page_meta"]["next"]
        url = ("https://www.ebi.ac.uk" + nxt) if nxt else None
        time.sleep(0.05)
    print(f"read {page} pages", flush=True)

    rows = []
    for cond in sorted(matched):
        for head, n in matched[cond].most_common():
            rows.append({"panel_condition": cond, "mesh_heading": head, "records": n,
                         "status": "mapped"})
    for head, n in discarded.most_common(60):
        rows.append({"panel_condition": "", "mesh_heading": head, "records": n,
                     "status": "discarded (neurological wording, no clear panel match)"})
    pd.DataFrame(rows).to_csv(OUT / "H6_mesh_mapping_audit.csv", index=False)

    for cond in sorted(matched):
        heads = ", ".join(h for h, _ in matched[cond].most_common(8))
        print(f"\n{cond} ({sum(matched[cond].values())} records, "
              f"{len(matched[cond])} distinct headings)\n    {heads}")
    print(f"\ndiscarded headings with neurological wording, top 25 of {len(discarded)}:")
    for h, n in discarded.most_common(25):
        print(f"  {n:5}  {h}")
    print("\nwrote", OUT / "H6_mesh_mapping_audit.csv")


if __name__ == "__main__":
    main()
