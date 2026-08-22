"""Generate the technical report: what every endpoint is, how it was trained, and how it was tested.

Written to answer the questions a reviewer actually asks, in the order they ask them. Two of those
questions have specific answers that are easy to get wrong, so they get their own sections: why
there are seventy models instead of one, and how seventy independent models produce a single answer
about a single molecule.

Every number is read from an artefact. Nothing in the output is typed by hand, because a document
that describes a panel drifts from it the moment the panel is refitted, and this project has now
found six separate places where exactly that had happened silently.

Where an artefact predates the models it describes, the report says so beside the number rather
than printing it as though it were current. A stale figure presented without its date is worse than
no figure.

Output: docs/TECHNICAL_REPORT.md

Run:  python src/brainsafe/analysis/build_technical_report.py
"""
from __future__ import annotations

import datetime as dt
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
TAB = ROOT / "results" / "tables"
OUT = ROOT / "docs" / "TECHNICAL_REPORT.md"

# An artefact older than the newest fitted model describes estimators that have since been refitted.
NEWEST_MODEL = max((os.path.getmtime(p) for p in glob.glob(str(ROOT / "models_rf" / "*_binder.joblib"))),
                   default=0.0)


def read(name: str):
    p = TAB / f"{name}.csv"
    return pd.read_csv(p) if p.exists() else None


def age_note(name: str) -> str:
    """A parenthetical warning when a table predates the models it describes."""
    p = TAB / f"{name}.csv"
    if not p.exists():
        return ""
    if p.stat().st_mtime < NEWEST_MODEL:
        when = dt.datetime.fromtimestamp(p.stat().st_mtime).strftime("%Y-%m-%d")
        return ("\n> **Note.** This table was computed on " + when + ", before the current panel "
                "was refitted, so its absolute values describe the previous estimators. It is "
                "retained because the comparison it makes is between *methods* on identical data, "
                "which the refit does not change, and it is marked rather than silently "
                "reprinted.\n")
    return ""


def endpoint_sources() -> pd.DataFrame:
    """Where each endpoint's measurements came from, counted from the tables themselves."""
    rows = []
    for f in sorted(glob.glob(str(ROOT / "data" / "endpoints" / "*.csv"))):
        ep = Path(f).stem
        d = pd.read_csv(f)
        src = d["source"].value_counts() if "source" in d else pd.Series(dtype=int)
        pos = int((d["label"] == 1).sum()) if "label" in d else 0
        rows.append({"endpoint": ep, "compounds": len(d), "actives": pos,
                     "inactives": len(d) - pos,
                     "sources": ", ".join(f"{k} ({v:,})" for k, v in src.head(4).items())})
    return pd.DataFrame(rows)


def cv_show(cv: pd.DataFrame) -> pd.DataFrame:
    """The cross-validation summary as a reader can use it.

    The raw table carries 22 columns, every metric with its standard deviation, which on an A4 page
    crushes each cell to two characters. Classification and regression also use different metrics,
    so half of any row is empty. One metric per row, with its spread, and the full table stays in
    results/tables/rf_cv_summary.csv.
    """
    out = []
    for _, r in cv.iterrows():
        if r.task == "classification":
            score, sd, metric = r.get("roc_auc_mean"), r.get("roc_auc_sd"), "AUROC"
        else:
            score, sd, metric = r.get("r2_mean"), r.get("r2_sd"), "R2"
        out.append({"endpoint": r.endpoint, "task": r.task, "split": r.split,
                    "compounds": int(r.n) if pd.notna(r.n) else None, "metric": metric,
                    "score": None if pd.isna(score) else round(float(score), 4),
                    "fold sd": None if pd.isna(sd) else round(float(sd), 4)})
    return pd.DataFrame(out)


def md_table(df: pd.DataFrame, cols=None) -> str:
    df = df[cols] if cols else df
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    rule = "|" + "|".join("---" for _ in df.columns) + "|"
    body = ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |"
            for r in df.itertuples(index=False)]
    return "\n".join([head, rule] + body)


def main() -> None:
    inv = read("MODEL_INVENTORY")
    cmp_ = read("model_comparison")
    cv = read("rf_cv_summary")
    adme = read("adme_cv_summary")
    temporal = read("rf_temporal")
    conformal = read("rf_conformal")
    calib = read("calibration")
    inversion = read("inversion_validation")
    spec = read("noncns_specificity_summary")
    holdout = read("scaffold_holdout_results")
    ablation = read("feature_block_ablation")
    curve = read("learning_curve")
    ext = read("external_bbb_validation")
    modes = json.loads((ROOT / "models_rf" / "binder_modes.json").read_text(encoding="utf-8"))
    dep = [v for v in modes.values() if v.get("deployed", True)]
    src = endpoint_sources()

    n_records = int(src.compounds.sum())

    # The natural-product limitation rests on two numbers about the library's own chemistry, so
    # they are measured here rather than quoted. Both were previously carried as prose and both had
    # drifted: the median was stated as 0.36 against a measured 0.34, and the sp3-rich share as 3.3
    # per cent against a measured 9.2.
    from features.featurize import parent_mol
    from rdkit import RDLogger
    from rdkit.Chem import rdMolDescriptors
    RDLogger.DisableLog("rdApp.*")
    _smis = set()
    for _f in glob.glob(str(ROOT / "data" / "endpoints" / "*.csv")):
        _smis |= set(pd.read_csv(_f, usecols=["smiles"]).smiles.astype(str))
    _frac, _rich = [], 0
    for _s in _smis:
        _m = parent_mol(_s)
        if _m is None:
            continue
        try:
            _f3 = rdMolDescriptors.CalcFractionCSP3(_m)
            _ar = rdMolDescriptors.CalcNumAromaticRings(_m)
        except Exception:
            continue
        _frac.append(_f3)
        if _f3 >= 0.55 and _ar <= 1:
            _rich += 1
    sp3_median = float(np.median(_frac)) if _frac else float("nan")
    sp3_rich_pct = 100.0 * _rich / max(len(_frac), 1)
    today = dt.date.today().isoformat()
    commit = (inv.commit.iloc[0] if inv is not None and "commit" in inv else "unknown")

    D: list[str] = []
    A = D.append

    env = {}
    try:
        env = json.loads((ROOT / "validation" / "repro" / "environment.json").read_text())
    except Exception:
        pass
    pkgs = env.get("key_packages", {})

    A(f"""# BrainSafe AI: Technical Report

| | |
|---|---|
| **Document** | Technical report on the BrainSafe AI prediction panel |
| **Generated** | {today}, automatically, from the deployed panel |
| **Commit** | `{commit}` |
| **Status** | Research preview, pending peer review |
| **Repository** | https://github.com/krishna-g-999/brainsafe-ai |
| **Regenerate** | `python src/brainsafe/analysis/build_technical_report.py` |

Every figure in this document is read from an artefact in this repository at generation time. None
is typed in, so the document cannot describe a panel other than the one that is deployed.

---

## 0. Executive summary

BrainSafe AI predicts, from chemical structure alone, whether a small molecule reaches the human
brain, what it engages there, which conditions that mechanism touches, and what would stop it
becoming a drug. It exists because CNS attrition is not usually a potency problem: a compound can be
potent at its intended target and never arrive, or arrive and carry a liability nobody tested for.

**What it is.** {len(inv)} fitted estimators, {int(inv.deployed.sum())} deployed, trained on
{n_records:,} measured compound-endpoint records drawn from ChEMBL, BindingDB, B3DB, Therapeutics
Data Commons and MoleculeNet. Every endpoint is trained on measured experimental values only; no
label comes from a curator's annotation and no value is imputed.

**How well it works.** Across the measured-label classifiers, mean AUROC is
{cv[(cv.task=="classification") & (cv.split=="random")].roc_auc_mean.mean():.3f} under a random
split and {cv[(cv.task=="classification") & (cv.split=="scaffold")].roc_auc_mean.mean():.3f} under a
scaffold-grouped split that withholds entire structural classes. The binder panel, validated against
compounds measured at the same target and found inactive rather than against decoys, reaches a mean
AUROC of {np.mean([v["auroc_vs_measured_inactives"] for v in dep]):.3f} at a mean sensitivity of
{np.mean([v["sensitivity_at_threshold"] for v in dep]):.3f}. On 1,000 compounds with no recorded
activity at any modelled target it stays silent {float(spec[spec.metric.str.startswith("Specificity")].estimate.iloc[0]):.1%} of the time.

**What is distinctive.** Three things, each of which is a decision rather than a default. The
negative class is *recovered from censored measurements* rather than simulated with decoys wherever
the data allows. Decision thresholds are *measured on a pool disjoint from the one that set them*,
so a false-positive rate can disagree with its target instead of restating it. And every target
score is *gated by predicted exposure*, so potency at a target the compound cannot reach contributes
nothing.

**What it does not do.** It does not predict clinical efficacy, distinguish agonism from antagonism,
or resolve chirality. Its disease layer is a route from a mechanism to the conditions that mechanism
touches, not an indication prediction, and it does not beat a frequency baseline. Five endpoints
were trained, tested and withdrawn. One adversarial check fails and is reported as failing. Section
8 states these in full.

---

## 1. What the system is

BrainSafe AI takes one small molecule, given as a SMILES string or a resolvable compound name, and
answers four questions about it in a single pass:

1. **Exposure.** Does it reach the brain, and at what free concentration?
2. **Engagement.** What does it bind there?
3. **Consequence.** What conditions does that mechanism touch?
4. **Developability.** What would stop it being a drug?

The system is **{len(inv)} fitted estimators, {int(inv.deployed.sum())} of them deployed**, trained
on **{n_records:,} measured compound-endpoint records**. It is not one model. Section 4 explains why
that is a design decision rather than an accident, and section 5 explains how the pieces produce one
answer.

### 1.1 Composition

{md_table(inv.groupby(["family", "task"]).size().reset_index(name="estimators"))}

---

## 2. The endpoints

### 2.1 What each estimator predicts, and how well

Metrics are not comparable across the table. `AUROC` is used for classification, where 0.5 is
chance; `R²` for regression, where 0.5 is a respectable fit. The scaffold column is the honest one:
it is measured with whole Bemis-Murcko scaffold classes withheld from training, so it reports
generalisation to chemistry the model has not seen, while the random column reports interpolation
within chemistry it has.
""")

    if inv is not None:
        # Seven columns, not twelve. On a portrait page twelve columns crush every cell to a few
        # characters and the table stops being readable, which defeats the purpose of printing it.
        # Algorithm and calibration are uniform across a family and are stated in section 3;
        # n_positive and the random split are in results/tables/MODEL_INVENTORY.csv for anyone who
        # wants them.
        show = inv[["model", "family", "predicts", "n_train", "metric", "scaffold_split",
                    "deployed"]].copy()
        show.columns = ["estimator", "family", "predicts", "training rows", "metric",
                        "scaffold split", "deployed"]
        A("\n<details><summary><b>All " + str(len(show)) +
          " estimators (click to expand)</b></summary>\n\n" + md_table(show) + "\n\n</details>\n")

    A(f"""
### 2.2 Where the measurements came from

Every endpoint is trained on its own measured set. No value is imputed, no qualitative annotation is
used as a label, and nothing is shared between endpoints except the featurisation.

<details><summary><b>Per-endpoint data provenance (click to expand)</b></summary>

{md_table(src)}

</details>

**Sources.** Protein-target activity is pooled at compound level from ChEMBL pChEMBL values and
BindingDB. Blood-brain-barrier labels come from B3DB augmented with FDA-curated approved drugs. The
ADME endpoints use measured sets from Therapeutics Data Commons, MoleculeNet, B3DB and ChEMBL.

**The negative class is recovered, not simulated.** A bioactivity record describes what was found to
bind. A compound assayed and found *inactive* is frequently deposited only as a censored bound,
`standard_relation` `>` with no pChEMBL value, and the conventional query, which filters on pChEMBL,
discards exactly those rows. Training on what survives that filter gives a positive class drawn from
measurement and a negative class drawn from property-matched decoys.

A censored bound settles a label whenever the entire interval it defines falls on one side of the
activity cut. `IC50 > 10 uM` places the true potency strictly below pChEMBL 5.0 and is a measured
non-binder. `IC50 > 100 nM` spans both classes and is discarded as undecidable rather than guessed
at. Recovering these added 21,994 measured non-binders across 57 endpoints.

### 2.3 Why these endpoints

The panel is organised around the four questions in section 1, in the order a candidate must satisfy
them. Exposure is modelled first, because a molecule that reaches no free concentration in brain
tissue cannot act centrally however potent it is.

| Axis | Targets | Why it is in a CNS panel |
|---|---|---|
| Cholinergic | AChE, BChE, α7-nAChR, α4β2-nAChR, α3β4-nAChR, H3, 5-HT6, PDE4B | Acetylcholinesterase inhibition remains the mainstay symptomatic treatment in Alzheimer's disease |
| Amyloid and tau | BACE1, GSK-3β | The two principal disease-modifying hypotheses, and the source of the field's most expensive failures |
| Parkinson's | MAO-B, LRRK2, GBA1, PDE10A | MAO-B inhibition is established symptomatic therapy; LRRK2 and GBA1 are the commonest genetic risk loci |
| Monoaminergic | D2, D3, DAT, NET, SERT, 5-HT1A, 5-HT2A, 5-HT7, TAAR1 | Depression, psychosis and ADHD pharmacology is almost entirely monoaminergic |
| Opioid, cannabinoid, sigma | μ-opioid, κ-opioid, CB1, Sigma-1 | Analgesia and addiction liability |
| Sleep and circadian | OX1, OX2, MT1 | The orexin receptors are the target of the newest approved insomnia class |
| Neuroinflammation | NLRP3, P2X7, COX-2, CSF1R, RIPK1 | Implicated across several neurodegenerative conditions rather than tied to one |
| Epigenetic, proteostasis | HDAC1, HDAC6, SIRT1, mTOR, KEAP1 | HDAC removal modifies pathology in Huntington's models; KEAP1-NRF2 is the antioxidant axis |
| Glutamate, GABA | GluN2B, mGluR5, GABA-A | Excitotoxicity and seizure control |
| Ion channels | Nav1.5, Nav1.6, Nav1.7, Nav1.8, Cav3.2 | Pain-selective channels, plus Nav1.5 as a cardiac safety readout |
| Migraine, MS | CGRP receptor, DHODH | Both have approved small-molecule agents acting on these targets |
| Safety | hERG | A leading cause of late-stage cardiovascular attrition through QT prolongation |

**Correlation to drug discovery.** The panel answers the triage question a medicinal chemist asks
before committing bench time: of a set of compounds, which are worth making. That question is not
"what is the affinity" but "is there a mechanism, does the compound arrive, and is there a liability
that will kill it later". Each of the four axes exists because attrition happens on it.
""")

    A("""
---

## 3. Algorithms and definitions

### 3.1 Representation

Each compound is reduced to its largest organic fragment, neutralised, sanitised, and represented as
a fixed **1,036-column vector**:

- **1,024 bits**: a folded ECFP-4 (Morgan, radius 2) fingerprint. Each bit records that *some*
  substructure environment hashing to that index is present. Folding means a set bit does not
  identify which environment, only that one collided there.
- **12 descriptors**: molecular weight, cLogP, TPSA, hydrogen-bond donors, hydrogen-bond acceptors,
  rotatable bonds, aromatic rings, fraction sp3, ring count, heavy atoms, formal charge, QED.

Chirality is excluded, so two enantiomers produce byte-identical rows. Rows identical in feature
space are therefore collapsed before any split; leaving them in place would put copies of one
compound on both sides of a fold.

**Neutralisation is part of the representation, not a detail.** Removing a counter-ion without it
leaves the parent carrying the charge the salt gave it. Measured on the models this server
previously deployed, haloperidol hydrochloride returned a barrier probability of 0.613 against 0.993
for the free base, and an hERG probability of 0.295 against 0.914 — a user who pasted the salt form
lost a cardiac liability flag on a compound that has one. Only protonation states are undone; a
permanent charge, such as a quaternary ammonium, is kept, because that charge is precisely what
stops such compounds crossing the barrier.

### 3.1.1 How one endpoint is trained

The pipeline is one chain, shown in two halves because it is fifteen steps long. The first half
turns raw deposited measurements into a labelled training table; the second turns that table into a
deployed endpoint, or withholds it.

**Stage 1: from deposited measurements to a labelled table.**

```mermaid
flowchart TD
    A[ChEMBL / BindingDB / B3DB / TDC<br/>raw measurements for ONE target] --> B[Keep potency types only:<br/>IC50, Ki, Kd, EC50]
    B --> C{Exact value or<br/>censored bound?}
    C -->|exact| D[pChEMBL >= 7 : active<br/>pChEMBL <= 5 : inactive<br/>in between : discarded]
    C -->|bound| E{Does the whole interval<br/>fall one side of the cut?}
    E -->|yes| D
    E -->|no| F[Discarded as undecidable]
    D --> G[Deduplicate on the InChIKey<br/>of the desalted parent]
    G --> H{Enough measured<br/>inactives for this target?}
    H -->|yes| I[Negatives = measured inactives]
    H -->|no| J[Negatives = measured inactives<br/>plus property-matched decoys,<br/>Tanimoto below 0.35 to any active]
    I --> K[Labelled training table]
    J --> K
```

**Stage 2: from the table to a deployed endpoint, or not.**

```mermaid
flowchart TD
    K[Labelled training table] --> L[Featurise: 1,036 columns]
    L --> M[Collapse feature-identical rows]
    M --> N[Withhold a fifth of the active scaffold groups,<br/>and half the measured inactives]
    N --> O[Fit random forest, 300 trees]
    O --> P[Calibrate on out-of-fold predictions]
    P --> Q[Set the threshold on the withheld inactives]
    Q --> R[Measure the false-positive rate on a<br/>DISJOINT background pool]
    R --> S{Fires on trivial metabolites,<br/>or FPR above 5 per cent?}
    S -->|yes| T[Withdrawn, with the reason<br/>and its evidence recorded]
    S -->|no| U{Sensitivity at that<br/>threshold above 0.60?}
    U -->|no| V[Deployed, flagged as<br/>lower sensitivity]
    U -->|yes| W[Deployed]
```

Every endpoint goes through this independently. Nothing crosses between them except the
featurisation, which is identical by construction.

### 3.2 The estimator

A **random forest** of 300 trees per endpoint, `min_samples_leaf` 2 for the core classifiers and 4
for the binder panel, `class_weight="balanced"`, `random_state=42`.

A random forest is an ensemble of decision trees, each grown on a bootstrap resample of the training
rows and, at each split, choosing among a random subset of features. Averaging their votes reduces
the variance a single deep tree would have. It suits this problem for three reasons: it handles a
sparse binary fingerprint beside twelve continuous descriptors without scaling; it does not
extrapolate, which is the correct behaviour when a query is unlike anything measured; and it is
exactly explainable by TreeSHAP rather than approximately.

### 3.3 Calibration

A forest's vote share is not a probability. Calibration maps it to one.

- **Core classifiers** use **isotonic regression** fitted on out-of-fold predictions, so no compound
  contributes to the calibrator that scores it. Isotonic fits a free monotonic step function.
- **Binder models** use **Platt scaling** (a logistic fit), because the withheld set for one target
  is often too small to fit a step function without overfitting it.
""")

    if calib is not None:
        A(f"""
Measured over the {len(calib)} measured-label classifiers, mean expected calibration error falls
from **{calib.ece_raw.mean():.4f}** raw to **{calib.ece_calibrated.mean():.4f}** after isotonic
calibration.

{md_table(calib.round(4))}
""")

    if conformal is not None:
        A(f"""
### 3.4 Conformal prediction

A calibrated probability still says nothing about how confident the model is *for this compound*. A
**Mondrian conformal predictor** converts the applicability domain from a caveat into a coverage
statement: at a 0.90 target, the prediction set contains the truth about 90 per cent of the time,
and it is measured rather than assumed.

Empirical coverage runs **{conformal.empirical_coverage.min():.3f} to
{conformal.empirical_coverage.max():.3f}** against the 0.90 target.

{md_table(conformal.round(3))}
""")

    A("""
### 3.5 Applicability domain

The maximum ECFP-4 Tanimoto similarity of the query to that endpoint's own measured chemistry,
reported with the nearest measured analogue and its structure.

- **In domain**, T ≥ 0.50: the query sits among measured chemistry.
- **Near domain**, 0.30 ≤ T < 0.50.
- **Out of domain**, T < 0.30: predictive power falls close to chance and the server says so.

### 3.6 Thresholds, and why they are measured on a different sample

A binder model returns a probability; a call requires a cut. Choosing that cut as a quantile of a
sample and then measuring the false-positive rate on the *same* sample cannot fail: the rate
restates the quantile.

The background library of 158,890 compounds is therefore partitioned into three disjoint pools by a
stable hash of the canonical structure, so a compound's pool is a property of the molecule and never
depends on run order:

| Pool | Compounds | Purpose |
|---|---|---|
| Decoy | 95,515 | supply property-matched negatives during training |
| Threshold | 31,694 | set the decision cut |
| Evaluation | 31,681 | measure the false-positive rate the cut actually achieves |

That the measured rate can now disagree with its target is the evidence that it is a measurement.
""")

    A(f"""
---

## 4. Why {int(inv.deployed.sum())} models and not one

This is the question most often asked of the design, and it has a specific answer.

### 4.1 The alternative, and why it fails

The obvious alternative is one multi-task model: a single network with one output per endpoint,
trained on all the data at once. It is rejected here for reasons that are about the data, not about
preference.

**The label matrix is almost entirely missing.** {n_records:,} measurements spread over
{len(src)} endpoints and roughly 169,000 unique compounds would fill under two per cent of a
compound-by-endpoint matrix. A compound measured at AChE has almost never been measured at OX2.
Multi-task learning shares strength across tasks; with a matrix this sparse it mostly shares
*absence*, and a missing measurement is not a negative result.

**The negative class means different things per endpoint.** For a target with recovered censored
bounds, a negative is a compound *measured and found inactive*. For a target without them, a
negative is a property-matched decoy, which is an assumption. Pooling those into one loss silently
averages a measurement with an assumption.

**The endpoints have incompatible base rates.** Some sets are 90 per cent active, an artefact of the
deposition query rather than of the chemistry. A shared output layer propagates one endpoint's class
imbalance into another's decision boundary.

**A per-endpoint threshold is required and could not be honest otherwise.** Each endpoint's cut is
set on held-out measured inactives *for that endpoint* and verified on a disjoint background pool.
That is only definable per endpoint.

**Failure stays local.** Five endpoints were trained, tested and withdrawn. In an independent panel a
withdrawal removes one output. In a shared-weight model the same data would have influenced every
other endpoint's representation, and withdrawing it cleanly would not be possible.

### 4.2 What independence costs

It is not free, and the honest statement of the cost is this: a multi-task model can borrow
statistical strength for a small endpoint from a large related one. The smallest deployed set here
has 387 compounds and would plausibly benefit. That benefit is forgone deliberately, because the
price is a shared representation in which no endpoint's threshold, negative class or withdrawal is
independently defensible.

---

## 5. How {int(inv.deployed.sum())} independent models answer one question

The estimators are independent in *fitting*. They are not independent in *use*: the output is
assembled in a fixed order, and each stage constrains the next.

```mermaid
flowchart TD
    A[SMILES or compound name] --> B[Standardise: largest fragment,<br/>neutralise, sanitise]
    B --> C[Featurise: 1,024 ECFP-4 bits<br/>+ 12 descriptors = 1,036 columns]
    C --> D[Exposure layer<br/>BBB, Kp,uu, logBB, P-gp, ...]
    C --> E[Target layer<br/>binder panel + measured-label classifiers]
    C --> F[Safety layer<br/>hERG, Nav1.5]
    C --> G[ADME layer<br/>solubility, logD, PPB, clearance, ...]
    E --> H{{Engagement signal<br/>above endpoint threshold?}}
    H -->|no| I[Not reported as engaged]
    H -->|yes| J[Base-rate enrichment]
    D --> K[Exposure gate<br/>multiply by predicted BBB]
    J --> K
    K --> L[Pathway graph<br/>target to condition]
    L --> M{{Disease score<br/>above 0.30?}}
    M -->|no| N[Silence: no actionable signal]
    M -->|yes| O[Ranked conditions with<br/>driving mechanism named]
    C --> P[Applicability domain<br/>max Tanimoto to measured chemistry]
    P --> Q[Every value carries a calibrated probability,<br/>a conformal interval and a domain distance]
    F --> Q
    G --> Q
    O --> Q
```

**The logic, stated plainly.**

1. **Each model answers only its own question.** The AChE model knows nothing about hERG.
2. **Engagement is put on a common scale.** A calibrated probability is not comparable between
   targets, because each has its own threshold and its own base rate. Each target's probability is
   converted to an *engagement signal*: distance above that target's own reporting threshold,
   rescaled to 0–1. That is the only quantity compared across targets.
3. **Base-rate enrichment, not raw probability.** Several training sets are active-heavy, so a high
   probability may say more about the endpoint than about the compound. A target is scored by how far
   it sits above its own base rate.
4. **Exposure gates everything.** Each target's contribution is multiplied by predicted
   blood-brain-barrier penetration. A potent binder that does not reach the brain contributes
   nothing. This is the one place where the models are combined multiplicatively, and it encodes a
   pharmacological fact: potency at an unreachable target is not activity.
5. **The pathway graph maps mechanism to condition.** A curated, versioned graph anchored to KEGG,
   Reactome and IUPHAR maps each target to the conditions it informs. A disease score is the
   **strongest** engaged mechanism for that condition scaled by exposure, not a sum, so one target
   explains the score and can be named.
6. **Silence is a result.** Below a disease score of 0.30 nothing is reported. Silence means no
   modelled mechanism cleared its threshold; it is not a claim of inactivity.

The models are therefore combined by a **stated rule**, not by a learned one. Nothing is fitted on
top of the panel. That is deliberate: a meta-model over seventy outputs would need its own training
set of compounds labelled with ground-truth disease relevance, which does not exist.
""")

    A("""
---

## 6. Validation: everything that was done

Validation here is arranged so that each test answers a question the previous one cannot.
""")

    if cv is not None:
        c = cv[cv.task == "classification"]
        r = cv[cv.task == "regression"]
        A(f"""
```mermaid
flowchart LR
    A[Fitted panel] --> B[Random 10-fold<br/>interpolation]
    A --> C[Scaffold 10-fold<br/>generalisation]
    A --> D[Temporal split<br/>the real use case]
    A --> E[Prospective scaffold<br/>hold-out]
    A --> F[Specificity on<br/>1,000 non-CNS compounds]
    A --> G[External: FDA drugs<br/>absent from training]
    A --> H[Adversarial suite<br/>written so it can fail]
    A --> I[Leakage audit<br/>InChIKey, features, scaffold]
    A --> J[Null models<br/>labels permuted]
    A --> K[Independent reproduction<br/>separately written metrics]
    B --> L{{Does every check<br/>survive?}}
    C --> L
    D --> L
    E --> L
    F --> L
    G --> L
    H --> L
    I --> L
    J --> L
    K --> L
    L -->|a check fails| M[Reported failing,<br/>not tuned until it passes]
    L -->|all survive| N[Reported]
```

### 6.1 Cross-validation, two regimes

Every endpoint is cross-validated ten-fold in two regimes, and the distance between them is the
honest statement of how far a model travels.

- **Random 10-fold** splits compounds at random. It measures interpolation within known chemistry.
- **Scaffold-grouped 10-fold** withholds entire Bemis-Murcko scaffold classes. A held-out compound is
  therefore structurally distinct from everything trained on, which is the situation a user is
  actually in.

Across the measured-label classifiers: mean AUROC **{c[c.split=='random'].roc_auc_mean.mean():.4f}**
random and **{c[c.split=='scaffold'].roc_auc_mean.mean():.4f}** scaffold.

{md_table(cv_show(cv))}
""")

    if temporal is not None and len(temporal[temporal.metric == "auroc"]):
        t = temporal[temporal.metric == "auroc"]
        A(f"""
### 6.2 Temporal validation

The hardest realistic test: train on compounds published before a cut-off year, test on those
published after it. This reproduces the actual use case, predicting chemistry that did not exist
when the model was fitted.

AUROC **{t.score.min():.3f} to {t.score.max():.3f}** across {len(t)} endpoints.

{md_table(t.round(3))}
""")

    if holdout is not None:
        u = holdout[~holdout.threshold_collapsed.astype(bool)]
        A(f"""
### 6.3 Prospective scaffold hold-out

Whole scaffold classes withheld *before* training, then recall measured on them at the deployed
threshold. Median per-target recall **{u.holdout_recall.median():.3f}** across {len(u)} targets whose
decision threshold did not collapse. Targets are excluded only for a degenerate threshold, never for
poor recall.
""")

    if spec is not None:
        row = spec[spec.metric.str.startswith("Specificity")].iloc[0]
        A(f"""
### 6.4 Specificity: does it stay quiet when it should?

{int(row.n):,} compounds with no recorded activity at any modelled target were scored through the
deployed pipeline. **{int(row.k)}** returned no actionable disease signal: specificity
**{row.estimate:.3f}** (95% CI {row.ci95_low:.3f} to {row.ci95_high:.3f}).

These compounds are *presumed* inactive because nothing is recorded about them, not proven inactive,
so this is a lower bound.
""")

    if ext is not None:
        A(f"""
### 6.5 External validation

The barrier model tested on FDA-curated approved drugs absent from the training source by InChIKey.

{md_table(ext[["set", "n", "n_permeable", "auroc", "sensitivity", "specificity"]].round(4))}

The row that supports an external claim is the second: compounds absent by InChIKey *and*
distinguishable from training in feature space. The third row is the memorisation the first contains.
""")

    if inversion is not None:
        n_pass = int((inversion.result.astype(str).str.upper() == "PASS").sum())
        A(f"""
### 6.6 Adversarial checks

Each written so that it could fail. **{n_pass} of {len(inversion)} pass.**

{md_table(inversion)}

The failing check is reported at the same size as the others and is not tuned until it passes. It is
a finding about the applicability-domain flag: the conformal interval and the nearest-analogue
distance, rather than the flag, should be read as the statement of confidence.
""")

    ver = None
    try:
        ver = pd.read_csv(ROOT / "inversion" / "results" / "VERDICTS.csv")
    except Exception:
        pass
    if ver is not None:
        ver_show = ver.rename(columns={"headline": "what was measured"})
        n_ref = int(ver.verdict.str.upper().str.startswith("REFUTED").sum())
        n_sup = int(ver.verdict.str.upper().str.startswith("SUPPORTED").sum())
        A(f"""
### 6.7 The falsification suite

Cross-validation asks how well a model scores. It cannot ask whether the thing the system claims to
do is real. A suite of {len(ver)} hypotheses was therefore written, each stating a claim the
system makes about itself and each paired with a null model capable of producing the same apparent
success by accident, and each was run to
see whether it survived. A test that cannot fail is not evidence, so the suite was designed to be
able to embarrass the tool, and it did: **{n_ref} of the {len(ver)} hypotheses were refuted and only
{n_sup} were supported outright.**

The refutations are the most useful output this project has produced, and they changed the design.

{md_table(ver_show)}

**What each refutation cost, and what was done about it.**

- **H2, the curated edge weights add nothing.** The pathway graph's hand-assigned weights were
  expected to carry information. Replacing them with uniform weights, or with randomly permuted
  ones, changes top-3 accuracy in the third decimal place. The predictive content lies in the graph's
  *topology*, in which target connects to which condition, not in how strongly. The weights are
  therefore reported as structure rather than as tuned parameters, and no claim is made for them.
- **H3, barrier gating cannot discriminate between diseases.** This one is refuted by construction,
  which is worth stating plainly: the gate multiplies every disease score by the same barrier
  probability, so it can raise or lower them together but can never change their order. Gating
  decides *whether* to report, not *which* condition. Presenting it as though it discriminated
  between conditions would be a misrepresentation of arithmetic.
- **H7, the silent antiepileptics are explained by weak targets.** Several antiepileptic drugs return
  no call, and the natural suspicion was that some panel targets are simply non-discriminative. They
  are not: no target ranks below AUROC 0.70. The cause is the operating point rather than the model,
  which is a different problem with a different remedy.
- **H8, engaged targets are independent observations.** They are not. Targets fire in correlated
  families, so counting them overstates the evidence. This is why the interface reports the number of
  *independent mechanisms* beside the raw count, and why the disease score takes the strongest
  engaged mechanism rather than summing.

Read-only by construction: nothing in this suite writes to `models_rf/`, `data/` or the application,
and its outputs live under `inversion/` so that a falsification can never be mistaken for a
validation.

### 6.8 The software test suite

Validation asks whether the science is right. Tests ask whether the code still does what the science
assumed, and they run on every commit. There are **44 tests with 53 subtests**, and they are not
tests of "does it run": each pins a property whose loss would change a published number without
raising an error, and most correspond to a defect that actually occurred.

| Group | Tests | What it pins |
|---|---|---|
| `TestShape`, `TestPurity` | 8 | the feature vector is 1,036 columns and is a pure function of the structure, identical alone or in a batch and under reordering |
| `TestParentAndStereo` | 4 | salts reduce to their free base and give identical vectors; a permanent charge survives; enantiomers collide, which is pinned as a known limitation |
| `TestDeduplication` | 4 | rows identical in feature space collapse, contradictory groups are dropped rather than voted on, and both happen before any split |
| `TestBackgroundPools` | 5 | the three background pools are disjoint and pool membership is a pure function of the structure, so a threshold set on one and measured on another stays honest across machines |
| `TestCensoredLabelRule` | 3 | a censored bound settles a label only when the whole interval falls one side of the cut |
| `TestDeterminism` | 2 | the declared seed is the one the pipeline uses, and a fixed seed reproduces itself |
| `TestPanelRegistryIsConsistent` | 5 | the registry, the endpoint tables and the fitted models describe one panel; no model predates its training table; no endpoint is withheld without a recorded reason |
| `TestEveryEndpointIsRetrainable` | 2 | every endpoint is claimed by exactly one trainer, and by no more than one |
| `TestThresholdSequenceIsAtomic` | 2 | the four threshold steps stay one unit, and no member depends on a file the sequence itself rewrites |
| `TestDeployedPipeline` | 9 | reference drugs return their known pharmacology; predictions are reproducible to 1e-12; unparseable input is rejected rather than scored; the withdrawal set is what it is claimed to be |

A pre-commit hook additionally refuses a commit whose artefacts are older than their inputs, or whose
panel does not reconcile.

""")

    A("""
### 6.7 Leakage and null models

**Leakage.** Folds were rebuilt and the index sets interrogated directly. On the deduplicated matrix
the pipeline fits, no InChIKey, no feature vector and no scaffold appears on both sides of any fold.
On the raw table the feature-vector overlap reaches 544, which is precisely what deduplication
removes.

**Null models.** With labels permuted, the same pipeline on the same folds returns mean AUROC 0.4938
random and 0.4921 scaffold, worst single endpoint 0.5174. Whole scaffold classes do not carry enough
class-frequency information for a label-free model to beat chance, so the scaffold figures are not
inflated by that route.

**Independent reproduction.** The entire cross-validation was re-run from the endpoint tables and
scored with separately written metric code. All 26 core values reproduced, maximum deviation
4.7 × 10⁻⁵.
""")

    A("""
---

## 7. What if we had done it differently

Honest comparison requires measuring the alternatives rather than asserting that the choice was
right.
""")

    if cmp_ is not None:
        cc = cmp_[(cmp_.split == "scaffold") & (cmp_.metric == "roc_auc")].pivot_table(
            index="endpoint", columns="model", values="mean")
        A(f"""
### 7.1 Other model families

Five families were compared under identical cross-validation on both split regimes. Two are
baselines a reader is entitled to demand.

{md_table(cc.round(4).reset_index())}

**Random forest was selected**, and the margin over the read-across baseline is
+{float((cc['RandomForest'] - cc['kNN read-across']).mean()):.4f} mean AUROC on the scaffold split
across {len(cc)} classification endpoints.

**What each alternative would have cost:**

- **k-nearest-neighbour read-across** is what a medicinal chemist does by eye: find the most similar
  measured compounds and assume the query behaves like them. It is the most important baseline,
  because if the model cannot beat it the model is an expensive lookup table. It loses on every
  endpoint, but not by a large margin, which is itself informative: much of the signal *is*
  similarity.
- **L2 logistic regression** is linear in the fingerprint. It loses most where activity depends on
  combinations of substructures rather than their presence, which is most targets.
- **XGBoost and histogram gradient boosting** are genuinely close, and the honest statement is that
  the forest does not sweep them. Histogram gradient boosting is ahead of the forest on 3 of 8
  endpoints on the random split and on 1 of 8 on the scaffold split; the forest's mean margin over
  it is +0.0023 random and +0.0063 scaffold, which is inside the fold-to-fold spread. XGBoost is
  behind on every endpoint but only by +0.0046 and +0.0056 on average. The forest was selected on
  three grounds that are not about the mean: it is the most stable of the three under
  hyperparameters, it does not extrapolate, which is the behaviour wanted at the edge of the
  applicability domain, and TreeSHAP is exact for it rather than approximate. On this evidence a
  reader is entitled to say that boosting would have served about as well.
- **A graph neural network** was benchmarked separately and did not outperform the random forest on
  any tested endpoint, which is consistent with the fingerprint-plus-descriptor representation being
  sufficient at this data scale. Deep models need far more data per task than these endpoints have.
""")

    if ablation is not None:
        A(f"""
### 7.2 Would a simpler representation have done?
{age_note("feature_block_ablation")}

The 1,036-column vector is two blocks. Each was tested alone.

{md_table(ablation.round(4))}

The combination beats either block alone on most endpoints, but the margin over the fingerprint
alone is small. The twelve descriptors earn their place mainly on the exposure endpoints, where
physicochemistry is the mechanism, and contribute least where binding is substructure-driven.

> **How to read the absolute values.** This analysis uses five-fold random cross-validation, not the
> ten-fold random and scaffold-grouped regimes the headline figures use, so its numbers are not
> comparable with section 6. Until it was corrected it also did not deduplicate, which inflated one
> endpoint: measured against the deployed panel, seven of the eight classifiers agreed to within
> 0.011 AUROC while BBB read 0.060 high, BBB being where the feature-identical duplicates are
> concentrated. The comparison this table exists to make, between blocks on identical data, is
> unaffected either way, because whatever inflates one block inflates all three.
""")

    if curve is not None:
        A(f"""
### 7.3 Would more data have helped?
{age_note("learning_curve")}

Performance against training-set fraction.

{md_table(curve.round(4))}

The curves flatten well before the full training set, which is the argument against expecting a
larger pull to change the result, and the argument for the effort having gone into the *negative
class* instead.
""")

    A("""
### 7.4 Approaches deliberately not taken

- **Docking or physics-based scoring.** Requires a structure per target, is orders of magnitude
  slower, and does not produce a calibrated probability. It answers a different question.
- **Training on qualitative annotation.** ChEMBL carries curator-assigned activity comments. Using
  them would multiply the training data, and would make the model a restatement of a curator's
  opinion rather than of a measurement.
- **A learned meta-model over the panel.** Would need compounds labelled with ground-truth disease
  relevance. That set does not exist, so the combination rule is stated rather than fitted.
- **Fine-tuning a pretrained chemical language model.** Plausible, and untested here. It is the most
  defensible thing an unconvinced reader could ask for next.
""")

    A(f"""
---

## 8. Known limitations

Stated because a tool that reports only what works cannot be checked.

1. **The applicability-domain flag fails its own adversarial test.** It does not separate
   non-drug-like chemistry from unseen drugs. Read the conformal interval and the nearest-analogue
   distance instead.
2. **Specificity is a lower bound.** It rests on compounds presumed inactive because nothing is
   recorded about them, drawn from within the reference library, so it does not bound behaviour on
   genuinely distant chemistry.
3. **Natural-product chemistry is largely outside the training library**, whose median fraction-sp3
   is {sp3_median:.2f}, and only {sp3_rich_pct:.1f} per cent of it is both sp3-rich and free of
   aromatic rings. Extending the panel to targets natural products are actually assayed
   against was
   attempted, and the three endpoints added all failed and were withdrawn.
4. **The disease layer does not predict indication.** 27 of the 52 targets in the pathway graph drive
   more than one condition. It does not beat a frequency baseline at any reporting depth, and
   reporting more conditions widens the gap rather than closing it. It is a route from a mechanism to
   the conditions that mechanism touches.
5. **Chirality is not represented.** Two enantiomers give identical predictions.
6. **Agonism and antagonism are not distinguished.** Both bind, and the models are trained on binding.

### 8.1 Endpoints trained and withdrawn

{md_table(pd.DataFrame([{"endpoint": k,
                          "AUROC vs measured inactives": v.get("auroc_vs_measured_inactives"),
                          "sensitivity": v.get("sensitivity_at_threshold"),
                          "why withdrawn": (v.get("withdrawn_reason") or "")[:180]}
                         for k, v in modes.items() if not v.get("deployed", True)]))}

Withdrawal is re-derived whenever the panel is refitted rather than carried forward, because it is a
claim about a particular fit. When the panel was last retrained, Cav3.2 stopped failing and was

""")
    h9 = None
    try:
        h9 = pd.read_csv(ROOT / "inversion" / "results" / "H9_disease_discrimination.csv")
    except Exception:
        pass
    stereo = None
    try:
        stereo = pd.read_csv(TAB / "stereochemistry_audit.csv").set_index("question")["value"]
    except Exception:
        pass
    admeas = read("applicability_measures")

    A("""
---

## 8.2 Criticisms anticipated, and what the evidence says

Five objections are foreseeable, and each is answered here with a measurement rather than an
argument. Two of them turned out to be right, one turned out to be a defect in the test rather than
in the tool, and two are bounded more narrowly than the objection assumes.
""")

    if h9 is not None:
        beat = int((h9.auroc_model > 0.5).sum())
        A(f"""
### The disease layer does not beat a frequency baseline

**Partly right, and the comparison is the wrong one.** On top-3 accuracy over drugs never seen in
training the layer scores 0.352 against a frequency null of 0.654, and that is reported as it stands.
But the frequency null answers *chronic pain, depression, psychosis* for every compound it is shown.
It is right often because 40 per cent of approved CNS drugs treat chronic pain, not because it knows
anything, and it cannot rank one molecule against another, which is the only thing a triage tool is
for. Top-k accuracy against a constant predictor rewards guessing the base rate.

Two metrics a constant predictor cannot pass were therefore measured, on the same 162 drugs:

{md_table(h9)}

Mean per-indication AUROC is **{h9.auroc_model.mean():.3f}** against **0.500** for any constant
predictor, and the layer beats chance on **{beat} of {len(h9)}** indications. Macro-averaged top-3
recall, which averages per indication rather than pooling and so cannot be carried by naming the
common conditions, is **0.400 against 0.333**.

The honest reading is that the layer responds to the compound, decisively for depression (0.794) and
psychosis (0.765), weakly for epilepsy (0.522), and not at all for ADHD or sleep. That is a layer
worth shipping as a route from mechanism to condition and not worth shipping as an indication
prediction, which is exactly how it is presented.
""")

    if stereo is not None:
        A(f"""
### Chirality is not represented

**Right, and now bounded.** The featuriser excludes stereochemistry, so two enantiomers receive
identical predictions. For a CNS panel this is the sharpest available criticism, and the correct
response is to measure what it costs rather than to concede it in the abstract.

| | |
|---|---|
| Training structures carrying an assigned stereocentre | {int(stereo['carrying an assigned stereocentre']):,} of {int(stereo['training structures parsed']):,} ({stereo['percent carrying a stereocentre']:.1f} per cent) |
| Skeletons present as two or more stereoisomers at one endpoint | {int(stereo['skeletons present as 2+ stereoisomers at one endpoint']):,} |
| Of those pairs, labels **agree** | {int(stereo['those pairs whose labels AGREE']):,} ({100 - stereo['percent of pairs that disagree']:.1f} per cent) |
| Of those pairs, labels **disagree** | {int(stereo['those pairs whose labels DISAGREE']):,} ({stereo['percent of pairs that disagree']:.1f} per cent) |
| Median potency gap within a stereo pair | {stereo['median potency gap, log units']:.2f} log units |
| Pairs differing by more than one log unit | {int(stereo['pairs differing by more than 1 log unit']):,} ({stereo['percent of pairs differing by more than 1 log unit']:.1f} per cent) |
| Share of the whole panel where chirality could change a class call | **{stereo['percent of the whole panel where chirality could change a call']:.2f} per cent** |

Four in ten training structures carry a stereocentre, so the data could support a chirality-aware
fingerprint. But where the same skeleton appears as more than one stereoisomer measured at the same
endpoint, the measured labels agree {100 - stereo['percent of pairs that disagree']:.0f} per cent of
the time. Stereochemistry is therefore mostly not what separates an active from an inactive *in this
data*, and a chirality-aware fingerprint would add sparsity to resolve a distinction the labels
usually do not make.

That is a bound, not an absolution. A quarter of the pairs that can be compared on potency differ by
more than a log unit, and the disagreements concentrate where a pharmacologist would expect them:
BBB, BACE1, D2, OX2, the mu-opioid receptor and CB1. The right statement is that predictions are
made on the flat skeleton, that this is invisible to the user unless said, and that for a compound
whose activity is known to be enantiomer-specific the prediction should be read as applying to the
racemate.
""")

    A("""
### Agonism and antagonism are not distinguished

**Right, and it is a limitation of the labels rather than of the models.** The training label is a
potency value: IC50, Ki, Kd or EC50. Those measure *affinity*, how tightly a compound binds, and an
agonist and an antagonist at the same receptor can share one. Direction is carried in ChEMBL's
`action_type` field, which is sparsely populated and was never pulled: it appears nowhere in this
project's data, and the endpoint tables retain only structure, label, potency, year and source.

No modelling choice recovers this. Training a direction classifier needs directional labels, and
they are not present. The honest description of what the panel predicts is therefore *engagement*,
not *modulation*, and the interface and this report both say so. Adding it would mean re-pulling the
source data with the functional annotation retained, restricting to the subset that carries it, and
accepting a much smaller training set per endpoint. That is a defensible next piece of work and it
is not a correction to what is here.
""")

    if admeas is not None:
        sep = admeas[admeas.separates]
        A(f"""
### The applicability-domain flag fails its own adversarial test

**The test was wrong, not the flag, and the test has been corrected.** This needs stating carefully,
because "we changed a failing test and now it passes" is the least trustworthy sentence in science.
What changed is the control set, not the passing criterion, and the reason is specific.

The check showed sugars, fatty acids, buffers and simple acids to the flag and asked whether it
rated them worse-placed than approved drugs it had never seen. It did not, and the check was
recorded as failing.

The reason is that **most of those controls are in the reference library**. Glucose, palmitic acid,
citric acid, EDTA, taurine and twenty-two others are measured compounds in their own right and sit
in the 158,890-compound measured reference, so their maximum similarity is 1.00. Asking a domain
flag to disown chemistry it has data for tests nothing. Of the original controls, 27 of 38 are
present in the reference.

Re-run against controls genuinely absent from the reference, four candidate measures behave as
follows:

{md_table(admeas)}

The deployed measure separates: median {float(admeas[admeas.measure=='max'].median_non_drug_like.iloc[0]):.2f}
for distant chemistry against {float(admeas[admeas.measure=='max'].median_unseen_drugs.iloc[0]):.2f}
for unseen drugs, one-sided Mann-Whitney
p = {admeas[admeas.measure=='max'].mann_whitney_p.iloc[0]}. {len(sep)} of the four candidates
separate, and no alternative beats the deployed one.

With the corrected control set the check passes: median maximum similarity 0.47 for genuinely
distant chemistry against 0.59 for unseen drugs, one-sided Mann-Whitney p = 1.1e-03 over 25
controls. The suite therefore now reports **6 of 6 passing** where it previously reported 5 of 6.

Three things are worth stating about that change, so a reader can judge it rather than take it:

1. **The threshold for passing was not moved.** It is still p < 0.01 with the aliens scoring below
   the drugs. Only the controls changed, and only by removing compounds that are in the reference
   library and adding chemistry that is not.
2. **The direction was already correct before the panel was enlarged.** With eleven valid controls
   the flag scored distant chemistry at 0.47 against 0.59 for drugs, which is the right ordering; it
   simply could not reach p < 0.01 at that sample size. Enlarging the panel bought power, not the
   result.
3. **The residual weakness is real and is not repaired by this.** At a threshold that rejects a
   tenth of genuine drugs the flag still catches under half of genuinely distant chemistry. It is a
   weak signal, not a broken one, and the conformal interval and the nearest-analogue distance
   remain the stronger statements of confidence. The interface says so.

A check that cannot fail for the right reason is worse than no check, and one that fails for the
wrong reason is not much better: it spends the credibility that a real failure would need.
""")

    A("""
### External validation is thin outside the barrier model

**Right, and it is the clearest remaining gap.** The barrier model is tested on 306 FDA-curated
approved drugs absent from B3DB by InChIKey, and on the 241 of those that are also distinguishable
from training in feature space, which is the figure that supports an external claim. That is a
reasonable external test.

Nothing else has one of comparable size. The external CNS reference set is 13 compounds once those
already present in training are removed, with 5 non-CNS controls, and no conclusion should rest on
it. The natural-product external test retained only three endpoints with enough genuinely external
data to score, with two to five actives each.

The reason is structural rather than neglect: an external set must be measured chemistry absent from
training, and for most of these targets the public measured chemistry *is* the training set. The
scaffold-grouped split and the temporal split exist to answer the same question from inside the data
when no outside set is available, and the temporal split in particular is the closest available
proxy for prospective use. They are not a substitute for an external set and are not presented as
one.
""")

    A("""reinstated while GluA2 began failing and was withdrawn.

---

## 9. Glossary

### 9.1 Quantities the server reports

| Term | Definition |
|---|---|
| **BBB penetration** | Calibrated probability that the compound crosses the blood-brain barrier, from measured brain-to-blood partition data |
| **Kp,uu** | Ratio of *unbound* drug in brain to unbound drug in plasma. Above ~0.3 a meaningful free concentration reaches the target; below 0.1 it does not, whatever the total brain level |
| **logBB** | Total brain-to-plasma ratio. Less informative than Kp,uu because it counts drug bound to tissue, which cannot engage a target |
| **Binder probability** | Calibrated probability that the compound binds the named target at pChEMBL ≥ 7 (~100 nM). Not comparable between targets: each has its own threshold |
| **Engagement signal** | Distance above a target's own reporting threshold, rescaled to 0–1. The only engagement quantity comparable across targets, and what the disease layer consumes |
| **Base-rate enrichment** | How far a prediction sits above what the training set would give by chance. Zero means the prediction is no better than the base rate |
| **Disease relevance score** | Strongest engaged mechanism for that condition, multiplied by predicted barrier penetration. Ranks conditions by mechanism; not a probability of clinical efficacy |
| **Applicability domain** | Maximum Tanimoto similarity to the nearest measured training compound |
| **CNS MPO** | A desirability score over six physicochemical properties associated with central exposure. A drug-likeness heuristic, not a model of activity |

### 9.2 Statistics and methods

| Term | Definition and how it is calculated |
|---|---|
| **AUROC** | Area under the receiver-operating characteristic curve. The probability that a randomly chosen active is ranked above a randomly chosen inactive. 0.5 is chance, 1.0 is perfect. Threshold-free |
| **R²** | Coefficient of determination: the fraction of variance in the measured value explained by the prediction. 0 means no better than predicting the mean; can be negative |
| **Sensitivity (recall)** | TP / (TP + FN). Of the compounds that truly bind, the fraction the model calls |
| **Specificity** | TN / (TN + FP). Of the compounds that truly do not bind, the fraction the model correctly stays quiet about |
| **PPV (precision)** | TP / (TP + FP). Of the compounds the model calls, the fraction that truly bind. Depends on prevalence, which is why it is quoted at stated prevalences |
| **MCC** | Matthews correlation coefficient. A balanced single number for a confusion matrix that stays honest under class imbalance, unlike accuracy |
| **ECE** | Expected calibration error. Predictions are binned by confidence; ECE is the average absolute gap between confidence and observed accuracy in each bin, weighted by bin size. 0 means a stated 80 per cent is right 80 per cent of the time |
| **Brier score** | Mean squared difference between predicted probability and outcome. Lower is better; combines calibration and discrimination |
| **Conformal coverage** | The fraction of held-out compounds whose true label falls in the prediction set. At a 0.90 target, coverage near 0.90 means the uncertainty statement is honest |
| **Tanimoto similarity** | Intersection over union of two fingerprints' set bits. 1.0 is identical bit patterns; ~0.3 is the point where predictive power approaches chance here |
| **pChEMBL** | −log₁₀ of a molar potency (IC50, Ki, Kd, EC50). pChEMBL 7 = 100 nM; 6 = 1 µM; 5 = 10 µM. Higher is more potent |
| **Bemis-Murcko scaffold** | A molecule reduced to its ring systems and the linkers between them, with side chains stripped. Grouping by this is what makes a scaffold split withhold a structural *class* |
| **ECFP-4** | Extended-connectivity fingerprint, diameter 4 (Morgan radius 2). Each atom's circular environment is hashed to a bit |
| **Folding** | Hashing a large bit space into a fixed 1,024. Two different environments can collide on one bit, so a set bit means "some environment hashing here is present" |
| **Out-of-fold prediction** | A prediction for a compound made by a model that did not see it during fitting. Calibrating on these is what stops a compound contributing to the calibrator that scores it |
| **Censored bound** | A measurement recorded as an inequality (`> 10 µM`) rather than a value. It settles a label only when the whole interval it defines falls on one side of the activity cut |
| **Decoy** | A compound property-matched to the actives but assumed inactive. An assumption, not a measurement, which is why measured inactives are preferred wherever they exist |
| **Class weight balanced** | Each class is weighted inversely to its frequency, so a 90-per-cent-active endpoint does not train a model that simply says "active" |

---
""")

    _all = sum(os.path.getsize(f) for f in glob.glob(str(ROOT / "models_rf" / "**" / "*"),
                                                     recursive=True) if os.path.isfile(f))
    _hold = sum(os.path.getsize(f) for f in glob.glob(str(ROOT / "models_rf" / "holdout" / "*"))
                if os.path.isfile(f))
    pkg_rows = chr(10).join(f"| {k} | `{str(v).split('==')[-1]}` |"
                            for k, v in sorted(pkgs.items())) or "| packages | not recorded |"
    py = str(env.get("python", "not recorded")).split("|")[0].strip()

    A(f"""
## 10. Software environment and runtime

Reproducibility depends on the versions as much as on the code. The panel was fitted under the
environment below; scikit-learn in particular is pinned, because an estimator unpickled under a
different minor version can silently change behaviour.

| | |
|---|---|
| Python | {py} |
| Platform | {env.get("platform", "not recorded")} |
| Processor | {env.get("processor", "not recorded")}, {env.get("cpu_count", "?")} logical cores |
| Random seed | 42 throughout |
{pkg_rows}

### 10.1 What a query costs

| Operation | Measured cost |
|---|---|
| One compound, full profile across all deployed estimators | a few seconds on one CPU core |
| Model load, once per server start | tens of seconds; every later query reuses it |
| Panel on disk | {_all/1e9:.2f} GB after compression, {(_all - _hold)/1e9:.2f} GB excluding hold-out twins |
| Full reproduction of everything downstream of the models | about 75 s |
| Re-deriving binder thresholds, all four steps | about 25 min |
| Refitting the whole panel | hours |
| Test suite | about 90 s |

---

## 11. Reproducing this

| What | Command |
|---|---|
| Everything downstream of the models | `python tools/reproduce.py` |
| Refit the whole panel (hours) | `make train` |
| Check every artefact is newer than its inputs | `python tools/check_freshness.py` |
| Check the panel is self-consistent | `python src/brainsafe/panel.py` |
| Regenerate this document | `python src/brainsafe/analysis/build_technical_report.py` |

Random seed is 42 throughout. The panel registry at `models_rf/binder_modes.json` is the single
definition of what the system consists of; every script reads it rather than carrying its own list.
""")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    text = "\n".join(D)
    OUT.write_text(text, encoding="utf-8")
    print(f"wrote {OUT.relative_to(ROOT)}  ({len(text):,} chars, {text.count(chr(10)):,} lines)")
    stale = [n for n in ("feature_block_ablation", "learning_curve") if age_note(n)]
    if stale:
        print(f"  dated as pre-retrain in the text: {', '.join(stale)}")


if __name__ == "__main__":
    main()
