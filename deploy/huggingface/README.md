---
title: BrainSafe AI
emoji: 🧠
colorFrom: yellow
colorTo: indigo
sdk: streamlit
sdk_version: 1.50.0
app_file: app.py
pinned: false
license: mit
short_description: Calibrated, exposure-gated prediction of small-molecule action in the human brain
---

# BrainSafe AI

Predicts, from a SMILES string or a compound name, whether a small molecule reaches the human brain,
what it engages there, what conditions that mechanism touches, and what would stop it being a drug.

Every endpoint is trained on measured experimental values, never on qualitative annotation. Every
prediction carries a calibrated probability, a conformal interval and the distance to the nearest
measured analogue, so a user can tell interpolation from extrapolation. Target scores are admitted
only in proportion to predicted blood-brain-barrier penetration, so potency at a target the compound
cannot reach contributes nothing.

The server reports where it does not work: five endpoints were trained, tested and withdrawn, one
adversarial check fails and is shown failing, and compounds outside the applicability domain are
flagged rather than guessed at.

- **Source, validation artefacts and reproduction:** https://github.com/krishna-g-999/brainsafe-ai
- **Licence:** MIT for the code. Underlying data retain their sources' licences (ChEMBL CC BY-SA 3.0,
  BindingDB, B3DB, Therapeutics Data Commons, MoleculeNet; pathway annotations from KEGG, Reactome
  and IUPHAR/BPS, whose terms restrict some redistribution).

Research preview, pending peer review. Not for medical, diagnostic or treatment use.
