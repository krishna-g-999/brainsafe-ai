"""BrainSafe AI, interactive predictor (current RF + ADME/exposure models).

A single-page Streamlit app over the deployed random-forest models: target engagement (calibrated
probabilities), receptor potency, a safety flag, the ADME / exposure panel, and the directly-modelled
unbound brain exposure (K_p,uu) with a free-brain-exposure verdict. Every number is a model output on
measured public data; this is a research triage tool, not a clinical or diagnostic device.

Run:  streamlit run app.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import streamlit as st
from rdkit import Chem
from rdkit.Chem import Draw

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src" / "brainsafe"))
from features.featurize import featurize_one  # noqa: E402

MODELS = ROOT / "models_rf"

TARGET_CLASSIFIERS = {
    "BBB": "Blood-brain barrier penetration",
    "AChE": "Acetylcholinesterase (Alzheimer's)",
    "BChE": "Butyrylcholinesterase (Alzheimer's)",
    "BACE1": "Beta-secretase 1 / amyloid (Alzheimer's)",
    "GSK3B": "GSK-3 beta / tau (Alzheimer's)",
    "MAO_A": "Monoamine oxidase A (mood)",
    "MAO_B": "Monoamine oxidase B (Parkinson's)",
    "hERG": "hERG cardiac-safety liability",
}
RECEPTOR_REGRESSORS = {
    "D2": "Dopamine D2 receptor", "A2A": "Adenosine A2A receptor",
    "HT2A": "Serotonin 5-HT2A receptor", "SERT": "Serotonin transporter",
}
ADME = {
    "kpuu": ("Unbound brain/plasma (Kp,uu)", "ratio", "10^"),
    "logbb": ("Brain penetration (logBB, total)", "log ratio", None),
    "caco2_permeability": ("Passive permeability (Caco-2)", "log Papp cm/s", None),
    "pgp_substrate": ("P-gp substrate (efflux)", "probability", "proba"),
    "pgp_inhibition": ("P-gp inhibition", "probability", "proba"),
    "solubility": ("Aqueous solubility (logS)", "log mol/L", None),
    "lipophilicity": ("Lipophilicity (logD)", "logD", None),
    "plasma_protein_binding": ("Plasma-protein binding", "% bound", None),
    "clearance_hepatocyte": ("Hepatocyte clearance", "uL/min/1e6", None),
}


@st.cache_resource
def load_models():
    m = {}
    for ep in TARGET_CLASSIFIERS:
        cal = MODELS / f"{ep}_calibrated.joblib"
        m[ep] = joblib.load(cal if cal.exists() else MODELS / f"{ep}.joblib")
    for ep in RECEPTOR_REGRESSORS:
        m[ep] = joblib.load(MODELS / f"{ep}.joblib")
    m["antioxidant_DPPH"] = joblib.load(MODELS / "antioxidant_DPPH.joblib")
    for ep in ADME:
        m[f"adme_{ep}"] = joblib.load(MODELS / "adme" / f"{ep}.joblib")
    return m


def predict_all(smiles, m):
    x = featurize_one(smiles)
    if x is None:
        return None
    x = x.reshape(1, -1)
    out = {"targets": {}, "receptors": {}, "adme": {}}
    for ep in TARGET_CLASSIFIERS:
        out["targets"][ep] = float(m[ep].predict_proba(x)[0, 1])
    for ep in RECEPTOR_REGRESSORS:
        out["receptors"][ep] = float(m[ep].predict(x)[0])
    out["antioxidant"] = float(m["antioxidant_DPPH"].predict(x)[0])
    for ep, (_, _, tf) in ADME.items():
        v = float(m[f"adme_{ep}"].predict(x)[0]) if tf != "proba" else float(m[f"adme_{ep}"].predict_proba(x)[0, 1])
        out["adme"][ep] = 10 ** v if tf == "10^" else v
    return out


def verdict(kpuu, bbb):
    if kpuu >= 0.3:
        return "Favourable", "#009E73", "Predicted to reach meaningful free concentration in the brain."
    if kpuu >= 0.1:
        return "Borderline", "#E69F00", "Some free brain exposure predicted; interpret with caution."
    return "Limited", "#D55E00", "Low predicted free brain exposure (poor penetration or active efflux)."


def main():
    st.set_page_config(page_title="BrainSafe AI", page_icon="🧠", layout="wide")
    st.title("🧠 BrainSafe AI")
    st.caption("Structure-based prediction of brain-relevant properties from measured public data "
               "(ChEMBL, BindingDB, B3DB). Random-forest models, scaffold-validated and calibrated. "
               "**Research triage tool, not a clinical or diagnostic device.**")

    with st.sidebar:
        st.header("Input")
        examples = {"Donepezil": "COc1cc2c(cc1OC)C(=O)C(CC1CCN(Cc3ccccc3)CC1)C2",
                    "Diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21",
                    "Atenolol (peripheral)": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
                    "Loperamide (effluxed)": "CN(C)C(=O)C(CCN1CCC(O)(c2ccc(Cl)cc2)CC1)(c1ccccc1)c1ccccc1"}
        pick = st.selectbox("Example compound", ["(type your own)"] + list(examples))
        default = examples.get(pick, "")
        smiles = st.text_area("SMILES", value=default, height=80)
        go = st.button("Predict", type="primary", use_container_width=True)

    if not (go and smiles.strip()):
        st.info("Enter a SMILES string (or pick an example) and press **Predict**.")
        return

    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        st.error("Could not parse that SMILES.")
        return
    m = load_models()
    r = predict_all(smiles.strip(), m)
    if r is None:
        st.error("Could not featurize that structure.")
        return

    c1, c2 = st.columns([1, 2])
    with c1:
        st.image(Draw.MolToImage(mol, size=(320, 240)))
    with c2:
        kpuu = r["adme"]["kpuu"]; bbb = r["targets"]["BBB"]
        label, colour, note = verdict(kpuu, bbb)
        st.markdown(f"### Free brain exposure: <span style='color:{colour}'>**{label}**</span>",
                    unsafe_allow_html=True)
        cc = st.columns(3)
        cc[0].metric("Predicted Kp,uu", f"{kpuu:.2f}")
        cc[1].metric("BBB penetration", f"{bbb:.0%}")
        cc[2].metric("hERG safety flag", f"{r['targets']['hERG']:.0%}")
        st.caption(note)

    st.divider()
    left, right = st.columns(2)
    with left:
        st.subheader("Target engagement (calibrated probability)")
        df = pd.DataFrame([{"Endpoint": ep, "Context": TARGET_CLASSIFIERS[ep],
                            "P(active)": f"{r['targets'][ep]:.0%}",
                            "Call": "active" if r["targets"][ep] >= 0.5 else "inactive"}
                           for ep in TARGET_CLASSIFIERS])
        st.dataframe(df, hide_index=True, use_container_width=True)

        st.subheader("Receptor potency (predicted pKi/pChEMBL)")
        dr = pd.DataFrame([{"Receptor": RECEPTOR_REGRESSORS[ep], "Predicted": f"{r['receptors'][ep]:.2f}"}
                          for ep in RECEPTOR_REGRESSORS] +
                          [{"Receptor": "Antioxidant (DPPH pIC50)", "Predicted": f"{r['antioxidant']:.2f}"}])
        st.dataframe(dr, hide_index=True, use_container_width=True)
    with right:
        st.subheader("ADME / exposure")
        da = pd.DataFrame([{"Property": ADME[ep][0], "Value": (f"{r['adme'][ep]:.2f}"),
                            "Units": ADME[ep][1]} for ep in ADME])
        st.dataframe(da, hide_index=True, use_container_width=True)

    st.divider()
    st.caption("Predictions are calibrated random-forest outputs trained on measured public bioactivity "
               "and ADME data; scaffold-split AUROC ~0.92 (target panel). Values are least reliable for "
               "compounds far from the training chemistry (applicability domain) and for the weakest "
               "endpoints (clearance, plasma-protein binding). Not for medical, diagnostic or treatment "
               "use. See docs/METHODS.md and docs/ADME_RESULTS.md.")


if __name__ == "__main__":
    main()
