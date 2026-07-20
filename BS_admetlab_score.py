"""
BS_admetlab_score.py
--------------------
Same-test-set comparison (reviewer point 5): BrainSafe vs ADMETlab 3.0 on IDENTICAL
molecules for BBB and hERG.

Input:
  admetlab_testset.csv          -- our 240 test compounds (endpoint, smiles, y, brainsafe_p)  [already created]
  admetlab_results.csv (arg)    -- the CSV downloaded from ADMETlab 3.0 Screening for the same SMILES
                                   (Services -> ADMET Screening -> paste SMILES -> Submit -> "Download as CSV")

It auto-detects ADMETlab's SMILES / BBB / hERG columns, canonicalises SMILES on both
sides, matches, and reports AUROC for BrainSafe and ADMETlab on the shared compounds.
Higher ADMETlab probability = more likely BBB-penetrant / hERG-blocker (positive class),
matching our label convention. Nothing is hand-entered.

Usage:  python BS_admetlab_score.py [admetlab_results.csv]
Output: BS_admetlab_comparison.json ; supplementary/STable16_admetlab_headtohead.csv
"""
import os, sys, re
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
from rdkit import Chem
from sklearn.metrics import roc_auc_score

RESULTS = sys.argv[1] if len(sys.argv) > 1 else "admetlab_results.csv"
if not os.path.exists(RESULTS):
    sys.exit(f"ADMETlab results file not found: {RESULTS}\n"
             "Download it from ADMETlab 3.0 Screening ('Download as CSV') and pass its path.")

def canon(s):
    try:
        m = Chem.MolFromSmiles(str(s)); return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None

test = pd.read_csv("admetlab_testset.csv"); test["c"] = test["smiles"].map(canon)
adm = pd.read_csv(RESULTS)

def find_col(cols, pat, avoid=None):
    for c in cols:
        if re.search(pat, str(c), re.I) and (avoid is None or not re.search(avoid, str(c), re.I)):
            return c
    return None
smi_col = find_col(adm.columns, r"^smiles$") or find_col(adm.columns, r"smiles")
bbb_col = find_col(adm.columns, r"\bbbb\b") or find_col(adm.columns, r"bbb")
herg_col = find_col(adm.columns, r"\bherg\b") or find_col(adm.columns, r"herg")
print(f"ADMETlab columns detected -> SMILES: {smi_col!r} | BBB: {bbb_col!r} | hERG: {herg_col!r}")
if not (smi_col and bbb_col and herg_col):
    print("Available columns:", list(adm.columns)); sys.exit("Could not detect required columns; edit the regex above.")

adm["c"] = adm[smi_col].map(canon)
for col in (bbb_col, herg_col):
    adm[col] = pd.to_numeric(adm[col], errors="coerce")
m = test.merge(adm[["c", bbb_col, herg_col]], on="c", how="inner")
print(f"matched {len(m)} of {len(test)} test compounds to ADMETlab output")

rows = []
for ep, acol in [("BBB", bbb_col), ("hERG", herg_col)]:
    sub = m[(m.endpoint == ep) & m[acol].notna()]
    if sub["y"].nunique() < 2 or len(sub) < 10:
        print(f"{ep}: too few matched/both-class compounds ({len(sub)}) to score"); continue
    bs = roc_auc_score(sub["y"], sub["brainsafe_p"])
    al = roc_auc_score(sub["y"], sub[acol])
    rows.append({"endpoint": ep, "n_matched": int(len(sub)),
                 "BrainSafe_AUROC": round(bs, 3), "ADMETlab_AUROC": round(al, 3),
                 "delta_BrainSafe_minus_ADMETlab": round(bs - al, 3)})
    print(f"{ep:5} n={len(sub):3}  BrainSafe {bs:.3f}  ADMETlab {al:.3f}  Δ {bs-al:+.3f}")

if rows:
    pd.DataFrame(rows).to_csv("supplementary/STable16_admetlab_headtohead.csv", index=False)
    import json
    json.dump({"rows": rows,
               "note": ("Same molecules, same labels; ADMETlab may have trained on some public BBB/hERG "
                        "compounds, which would favour it, whereas BrainSafe scores are strict out-of-fold. "
                        "Interpret accordingly.")},
              open("BS_admetlab_comparison.json", "w"), indent=2)
    print("\nWrote BS_admetlab_comparison.json + STable16_admetlab_headtohead.csv")
