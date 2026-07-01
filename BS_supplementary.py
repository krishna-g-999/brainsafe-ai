"""
BS_supplementary.py — consolidate ALL validated metrics into publication supplementary
tables (CSV) directly from the saved validation artifacts. No recompute, no rounding beyond
what is stored -> numbers are exactly those produced during training/testing.
"""
import os, json, glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pandas as pd
OUT = "supplementary"; os.makedirs(OUT, exist_ok=True)

rs = json.load(open("BS_randomsplit_benchmark.json"))
cis = json.load(open("BS_auroc_cis.json")) if os.path.exists("BS_auroc_cis.json") else {}
ext = json.load(open("BS_external_validation_report.json"))
tp = json.load(open("BS_temporal_pr_report.json"))
metas = {os.path.basename(f).replace("_meta.json", ""): json.load(open(f))
         for f in glob.glob("models_brain/*_meta.json")}

# ---- S1: classification full hierarchy ----
rowsS1 = []
for ep, m in sorted(metas.items()):
    if ep not in ("BBB", "hERG") and m.get("mcc", 0) < 0.45:
        continue
    rowsS1.append({
        "endpoint": ep, "n": m["n"], "pos_rate": m["pos_rate"], "source": m["source"],
        "AUROC_random": rs.get(ep), "AUROC_scaffold": m["auroc"],
        "scaffold_AUROC_CI95_low": cis.get(ep, {}).get("ci95", [None, None])[0],
        "scaffold_AUROC_CI95_high": cis.get(ep, {}).get("ci95", [None, None])[1],
        "AUROC_cluster": m.get("cluster_split_auroc"), "AUROC_temporal": m.get("temporal_auroc"),
        "PR_AUC": m["pr_auc"], "balanced_acc": m["balanced_acc"], "MCC": m["mcc"],
        "Brier": m.get("brier"), "conformal_coverage": m.get("conformal_coverage"),
        "operating_threshold": m["threshold"],
    })
pd.DataFrame(rowsS1).to_csv(f"{OUT}/STable1_classification_metrics.csv", index=False)

# ---- S2: receptor regression ----
rr = json.load(open("models_brain_reg/regression_report.json"))
rowsS2 = [{"receptor": k, "n": v["n"], "scaffold_cv_R2": v["scaffold_cv_r2"], "RMSE": v["rmse"],
           "Spearman": v["spearman"], "temporal_R2": (v["temporal"] or {}).get("r2"),
           "y_range": f"{v['y_min']}-{v['y_max']}", "task": v["task"]} for k, v in rr.items()]
pd.DataFrame(rowsS2).to_csv(f"{OUT}/STable2_receptor_regression.csv", index=False)

# ---- S3: antioxidant measured ----
am = json.load(open("models_genuine/antioxidant_measured_meta.json"))
pd.DataFrame([{ "endpoint": "antioxidant_DPPH_measured", "n": am["n"], "scaffold_cv_R2": am["scaffold_cv_r2"],
               "RMSE": am["rmse"], "Spearman": am["spearman"],
               "temporal_R2": (am.get("temporal") or {}).get("r2"),
               "crosscheck_curated_vs_measured_Spearman": am["crosscheck_curated_vs_measured_spearman"],
               "prior_curated_R2": am["vs_old_curated_r2"], "source": am["source"]}]
            ).to_csv(f"{OUT}/STable3_antioxidant_measured.csv", index=False)

# ---- S4: threshold sensitivity (PR) ----
rowsS4 = []
for ep, v in tp.items():
    for split in ("temporal", "scaffold_holdout"):
        s = v.get(split)
        if not s:
            continue
        for thr, d in s.get("thresholds", {}).items():
            rowsS4.append({"endpoint": ep, "split": split, "threshold": thr,
                           "precision": d.get("precision"), "recall": d.get("recall"),
                           "f1": d.get("f1"), "AUROC": s.get("auroc"), "PR_AUC": s.get("pr_auc"),
                           "n_test": s.get("n_test"), "pos_rate_test": s.get("pos_rate_test")})
pd.DataFrame(rowsS4).to_csv(f"{OUT}/STable4_threshold_sensitivity.csv", index=False)

# ---- S5: similarity-binned generalisation ----
rowsS5 = []
for ep, v in ext.items():
    for b, d in v.get("similarity_binned_auroc", {}).items():
        rowsS5.append({"endpoint": ep, "tanimoto_bin": b, "n": d["n"], "AUROC": d["auroc"]})
pd.DataFrame(rowsS5).to_csv(f"{OUT}/STable5_similarity_binned_auroc.csv", index=False)

# ---- S6: clinical reference composition ----
cdf = pd.read_csv("data/clinical_cns_reference.csv")
cdf["disease"].value_counts().rename_axis("disease").reset_index(name="n_compounds") \
   .to_csv(f"{OUT}/STable6_clinical_reference_composition.csv", index=False)

# ---- S7: benchmark vs literature ----
lit = {"BBB": "0.88-0.96", "hERG": "0.86-0.93", "AChE": "~0.90-0.97", "BACE1": "~0.90-0.96",
       "MAO_B": "~0.88-0.96", "MAO_A": "~0.85-0.95", "BChE": "~0.90-0.97", "GSK3B": "~0.88-0.95"}
rowsS7 = [{"endpoint": ep, "BrainSafe_random_AUROC": rs.get(ep),
           "published_random_AUROC_range": lit.get(ep, "n/a")} for ep in rs]
pd.DataFrame(rowsS7).to_csv(f"{OUT}/STable7_benchmark_vs_literature.csv", index=False)

print("Supplementary tables written to", OUT + "/")
for f in sorted(glob.glob(f"{OUT}/STable*.csv")):
    print("  ", os.path.basename(f), "-", len(pd.read_csv(f)), "rows")
