"""Build the complete, manuscript-ready endpoint tables (markdown + CSV).

Table 1  Core target panel: 13 endpoints, data size, class balance, random and scaffold 10-fold.
Table 2  Decoy-aware binder classifiers: 18 receptor/transporter/kinase endpoints.
Table 3  ADME / exposure layer: 9 endpoints.
Table 4  Error-bar (between-fold variance) decomposition.

Written to manuscript/tables_generated.md and results/tables/*.csv
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
TAB = ROOT / "results" / "tables"
OUT_MD = ROOT / "manuscript" / "tables_generated.md"

TARGET_BIOLOGY = {
    "BBB": ("Blood-brain barrier penetration", "Gate: no CNS effect without brain entry"),
    "AChE": ("Acetylcholinesterase", "Symptomatic Alzheimer therapy (donepezil class)"),
    "BChE": ("Butyrylcholinesterase", "Rises as Alzheimer progresses; selective-inhibitor target"),
    "BACE1": ("Beta-secretase 1", "Rate-limiting step of amyloid-beta generation"),
    "GSK3B": ("Glycogen synthase kinase-3 beta", "Tau hyperphosphorylation; neuroprotection"),
    "MAO_A": ("Monoamine oxidase A", "Serotonin/noradrenaline catabolism; depression"),
    "MAO_B": ("Monoamine oxidase B", "Dopamine catabolism; Parkinson therapy (selegiline)"),
    "hERG": ("hERG potassium channel", "Cardiotoxicity liability; principal safety filter"),
    "D2": ("Dopamine D2 receptor", "Antipsychotic efficacy; motor control"),
    "A2A": ("Adenosine A2A receptor", "Non-dopaminergic Parkinson target (istradefylline)"),
    "HT2A": ("Serotonin 5-HT2A receptor", "Atypical antipsychotics; psychedelics"),
    "SERT": ("Serotonin transporter", "SSRI antidepressant target"),
    "antioxidant_DPPH": ("Radical-scavenging capacity", "Oxidative stress in neurodegeneration"),
}
BINDER_BIOLOGY = {
    "HT1A": ("5-HT1A receptor", "Anxiety, depression (buspirone)"),
    "HT6": ("5-HT6 receptor", "Cognition enhancement in Alzheimer"),
    "HT7": ("5-HT7 receptor", "Mood, circadian rhythm, sleep"),
    "H3": ("Histamine H3 receptor", "Wakefulness (pitolisant), cognition"),
    "DAT": ("Dopamine transporter", "ADHD, addiction, stimulant liability"),
    "NET": ("Noradrenaline transporter", "Depression, ADHD (atomoxetine)"),
    "Sigma1": ("Sigma-1 receptor", "Neuroprotection, ER-stress chaperone"),
    "CB1": ("Cannabinoid CB1 receptor", "Pain, appetite, mood"),
    "OPRK1": ("Kappa-opioid receptor", "Analgesia, dysphoria, mood"),
    "OPRM1": ("Mu-opioid receptor", "Analgesia, addiction liability"),
    "D3": ("Dopamine D3 receptor", "Addiction, Parkinson motor complications"),
    "A1": ("Adenosine A1 receptor", "Neuroprotection, epilepsy, sedation"),
    "a7nAChR": ("Alpha-7 nicotinic receptor", "Cognition, neuroinflammation in Alzheimer"),
    "LRRK2": ("LRRK2 kinase", "Most common genetic cause of Parkinson disease"),
    "D2": ("Dopamine D2 receptor", "Antipsychotic efficacy"),
    "A2A": ("Adenosine A2A receptor", "Non-dopaminergic Parkinson target"),
    "HT2A": ("Serotonin 5-HT2A receptor", "Atypical antipsychotic profile"),
    "SERT": ("Serotonin transporter", "SSRI antidepressant target"),
}
ADME_BIOLOGY = {
    "kpuu": ("Unbound brain/plasma ratio (Kp,uu)", "Free drug available to CNS targets"),
    "logbb": ("Total brain/plasma ratio (logBB)", "Bulk brain distribution"),
    "caco2_permeability": ("Caco-2 permeability", "Passive membrane permeability"),
    "pgp_substrate": ("P-glycoprotein substrate", "Active efflux out of the brain"),
    "pgp_inhibition": ("P-glycoprotein inhibition", "Efflux-mediated drug interactions"),
    "solubility": ("Aqueous solubility (logS)", "Formulation and absorption"),
    "lipophilicity": ("Lipophilicity (logD)", "Permeability/promiscuity balance"),
    "plasma_protein_binding": ("Plasma protein binding", "Determines free fraction"),
    "clearance_hepatocyte": ("Hepatocyte clearance", "Metabolic stability, exposure duration"),
}
ADME_ORDER = ["kpuu", "logbb", "caco2_permeability", "pgp_substrate", "pgp_inhibition",
              "solubility", "lipophilicity", "plasma_protein_binding", "clearance_hepatocyte"]


def md_table(headers, rows):
    out = "| " + " | ".join(headers) + " |\n"
    out += "|" + "|".join(["---"] * len(headers)) + "|\n"
    for r in rows:
        out += "| " + " | ".join(str(c) for c in r) + " |\n"
    return out


def main():
    T1 = pd.read_csv(TAB / "manuscript_T1_endpoints.csv")
    T4 = pd.read_csv(TAB / "manuscript_T4_variance_decomposition.csv")
    adme = pd.read_csv(TAB / "adme_cv_summary.csv")

    md = ["# Generated tables (all numbers computed from saved cross-validation predictions)\n"]

    # ---------------- Table 1: core target panel ----------------
    rows = []
    for ep in ["BBB", "AChE", "BChE", "BACE1", "GSK3B", "MAO_A", "MAO_B", "hERG",
               "D2", "A2A", "HT2A", "SERT", "antioxidant_DPPH"]:
        a = T1[(T1.endpoint == ep) & (T1.split == "random")]
        b = T1[(T1.endpoint == ep) & (T1.split == "scaffold")]
        if a.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        nm, why = TARGET_BIOLOGY[ep]
        if a.task == "classification":
            metric = "AUROC"
            rnd = f"{a.roc_auc_mean:.3f} ± {a.roc_auc_sd:.3f}"
            scf = f"{b.roc_auc_mean:.3f} ± {b.roc_auc_sd:.3f}"
            bal = f"{a.positive_rate:.2f}"
        else:
            metric = "R2"
            rnd = f"{a.r2_mean:.3f} ± {a.r2_sd:.3f}"
            scf = f"{b.r2_mean:.3f} ± {b.r2_sd:.3f}"
            bal = "n/a"
        rows.append([ep, nm, "classification" if a.task == "classification" else "regression",
                     f"{int(a.n_compounds):,}", f"{int(a.n_scaffolds):,}", bal,
                     f"{int(a.median_train_per_fold):,}", f"{int(a.median_test_per_fold):,}",
                     metric, rnd, scf, why])
    md.append("\n## Table 1. Core target panel (13 endpoints), 10-fold cross-validation\n")
    md.append(md_table(["Endpoint", "Target", "Task", "Compounds", "Scaffolds", "Active fraction",
                        "Train/fold", "Test/fold", "Metric", "Random 10-fold", "Scaffold 10-fold",
                        "Why this endpoint"], rows))
    pd.DataFrame(rows, columns=["endpoint", "target", "task", "compounds", "scaffolds", "active_fraction",
                                "train_per_fold", "test_per_fold", "metric", "random_10fold",
                                "scaffold_10fold", "rationale"]).to_csv(TAB / "manuscript_Table1_core.csv", index=False)

    # ---------------- Table 2: decoy-aware binder classifiers ----------------
    rows = []
    for ep in ["D2", "A2A", "HT2A", "SERT", "HT1A", "HT6", "HT7", "H3", "DAT", "NET",
               "Sigma1", "CB1", "OPRK1", "OPRM1", "D3", "A1", "a7nAChR", "LRRK2"]:
        f = ROOT / "models_rf" / f"{ep}_binder_meta.json"
        if not f.exists():
            continue
        m = json.loads(f.read_text())
        nm, why = BINDER_BIOLOGY[ep]
        n_pos = m.get("n_positive", m.get("n_active"))
        n_neg = m.get("n_decoy") or m.get("n_measured_inactive_train") or 0
        auroc = m.get("auroc_vs_heldout_measured_inactives", m.get("auroc_hard_decoys"))
        thr = m.get("decision_threshold", m.get("threshold"))
        sens = m.get("sensitivity_at_threshold")
        rows.append([ep, nm, f"{n_pos:,}" if n_pos else "n/a", f"{n_neg:,}",
                     f"{auroc:.3f}" if auroc else "n/a",
                     f"{thr:.3f}" if thr else "n/a",
                     f"{sens:.3f}" if sens else "n/a", why])
    md.append("\n## Table 2. Binder classifiers validated against held-out measured inactives\n")
    md.append(md_table(["Endpoint", "Target", "Measured binders", "Training negatives",
                        "AUROC (held-out measured inactives)", "Threshold", "Sensitivity",
                        "Why this endpoint"], rows))
    pd.DataFrame(rows, columns=["endpoint", "target", "n_binders", "n_training_negatives",
                                "auroc_heldout_measured_inactives", "threshold", "sensitivity",
                                "rationale"]).to_csv(TAB / "manuscript_Table2_binders.csv", index=False)

    # ---------------- Table 3: ADME ----------------
    rows = []
    for ep in ADME_ORDER:
        a = adme[(adme.endpoint == ep) & (adme.split == "random")]
        b = adme[(adme.endpoint == ep) & (adme.split == "scaffold")]
        if a.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        nm, why = ADME_BIOLOGY[ep]
        if a.task == "classification":
            metric, rnd, scf = "AUROC", f"{a.roc_auc_mean:.3f} ± {a.roc_auc_sd:.3f}", f"{b.roc_auc_mean:.3f} ± {b.roc_auc_sd:.3f}"
        else:
            metric, rnd, scf = "R2", f"{a.r2_mean:.3f} ± {a.r2_sd:.3f}", f"{b.r2_mean:.3f} ± {b.r2_sd:.3f}"
        rows.append([ep, nm, "classification" if a.task == "classification" else "regression",
                     f"{int(a.n):,}", metric, rnd, scf, why])
    md.append("\n## Table 3. ADME / exposure layer (9 endpoints)\n")
    md.append(md_table(["Endpoint", "Property", "Task", "Compounds", "Metric",
                        "Random 10-fold", "Scaffold 10-fold", "Why this endpoint"], rows))
    pd.DataFrame(rows, columns=["endpoint", "property", "task", "compounds", "metric",
                                "random_10fold", "scaffold_10fold",
                                "rationale"]).to_csv(TAB / "manuscript_Table3_adme.csv", index=False)

    # ---------------- Table 4: variance decomposition ----------------
    rows = []
    for ep in T4.endpoint.unique():
        a = T4[(T4.endpoint == ep) & (T4.split == "random")].iloc[0]
        b = T4[(T4.endpoint == ep) & (T4.split == "scaffold")].iloc[0]
        rows.append([ep, f"{a.sd_observed:.4f}", f"{a.sd_sampling:.4f}", f"{a.pct_variance_heterogeneity:.0f}%",
                     f"{b.sd_observed:.4f}", f"{b.sd_sampling:.4f}", f"{b.pct_variance_heterogeneity:.0f}%"])
    md.append("\n## Table 4. Between-fold error-bar decomposition\n")
    md.append("Observed between-fold SD separated into sampling noise (finite test set) and genuine "
              "fold-to-fold heterogeneity, by within-fold bootstrap.\n\n")
    md.append(md_table(["Endpoint", "Random SD", "of which sampling", "heterogeneity share",
                        "Scaffold SD", "of which sampling", "heterogeneity share"], rows))

    # ---------------- Table 5: temporal (future-compound) validation ----------------
    tp = TAB / "rf_temporal.csv"
    if tp.exists():
        t = pd.read_csv(tp)
        rows = []
        for _, x in t.iterrows():
            ep = x.endpoint
            nm = TARGET_BIOLOGY.get(ep, (ep, ""))[0]
            rows.append([ep, nm, int(x.cutoff_year), f"{int(x.n_train):,}", f"{int(x.n_test):,}",
                         "AUROC" if x.metric == "auroc" else "R2", f"{x.score:.3f}"])
        md.append("\n## Table 5. Temporal (future-compound) validation\n")
        md.append("Models are trained only on compounds published before the cutoff year and tested on "
                  "compounds published after it. This is the most demanding regime and the closest analogue "
                  "of prospective use.\n\n")
        md.append(md_table(["Endpoint", "Target", "Cutoff year", "Train", "Test", "Metric", "Score"], rows))
        c = t[t.task == "classification"]; r_ = t[t.task == "regression"]
        md.append(f"\nClassifier endpoints: mean AUROC {c.score.mean():.3f} "
                  f"(range {c.score.min():.3f} to {c.score.max():.3f}). "
                  f"Regression endpoints: mean R2 {r_.score.mean():.3f} "
                  f"(range {r_.score.min():.3f} to {r_.score.max():.3f}).\n")
        t.to_csv(TAB / "manuscript_Table5_temporal.csv", index=False)


    # ---- Table 6: scaffold hold-out prospective validation ----
    sh = TAB / "scaffold_holdout_results.csv"
    if sh.exists():
        d6 = pd.read_csv(sh).sort_values("holdout_recall", ascending=False)
        rows = []
        for _, x in d6.iterrows():
            rows.append([x.target, f"{int(x.train_actives):,}", f"{int(x.holdout_actives):,}",
                         f"{int(x.holdout_scaffolds):,}", f"{x.threshold:.3f}",
                         f"{x.holdout_recall:.3f}",
                         f"[{x.recall_ci95_low:.3f}, {x.recall_ci95_high:.3f}]",
                         "excluded" if x.threshold_collapsed else ""])
        good6 = d6[~d6.threshold_collapsed]
        n6 = int(good6.holdout_actives.sum())
        k6 = int((good6.holdout_recall * good6.holdout_actives).round().sum())
        md.append("\n## Table 6. Prospective validation under a scaffold hold-out\n")
        md.append("Twenty per cent of Bemis-Murcko scaffolds were withheld per target and every model "
                  "retrained on the remainder, so no held-out compound shares a scaffold with anything "
                  "its model saw. Thresholds were recalibrated on held-out negatives and an independent "
                  "background sample. Targets marked excluded produced a threshold at the permitted "
                  "floor, meaning no separation from background chemistry, and do not contribute to "
                  "the pooled estimate.\n\n")
        md.append(md_table(["Target", "Train actives", "Held-out actives", "Held-out scaffolds",
                            "Threshold", "Recall", "95% CI", "Note"], rows))
        md.append(f"\nPooled recall {k6:,}/{n6:,} = {k6/n6:.3f}; median per-target "
                  f"{good6.holdout_recall.median():.3f}; {(good6.holdout_recall>=0.80).sum()} of "
                  f"{len(good6)} targets at or above 0.80.\n")

    OUT_MD.write_text("\n".join(md), encoding="utf-8")
    print("wrote", OUT_MD)
    print("core:", len(pd.read_csv(TAB / 'manuscript_Table1_core.csv')),
          "binders:", len(pd.read_csv(TAB / 'manuscript_Table2_binders.csv')),
          "adme:", len(pd.read_csv(TAB / 'manuscript_Table3_adme.csv')))


if __name__ == "__main__":
    main()
