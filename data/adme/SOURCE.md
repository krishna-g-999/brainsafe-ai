# Source: measured ADME / exposure datasets (ladder B)

Downloaded 2026-07-21 from TDC / Harvard Dataverse and MoleculeNet; standardised by
`src/brainsafe/adme/fetch_adme.py`. All values are measured.

| Endpoint | Task | n (standardised) | Description |
|---|---|---|---|
| solubility | regression | 9573 | aqueous solubility logS (mol/L), AqSolDB |
| caco2_permeability | regression | 897 | log Papp (cm/s), Wang et al. |
| pgp_inhibition | classification | 1212 | P-glycoprotein inhibitor 0/1, Broccatelli |
| plasma_protein_binding | regression | 1797 | percent bound, AstraZeneca |
| clearance_hepatocyte | regression | 1020 | uL/min/1e6 cells, AstraZeneca |
| lipophilicity | regression | 4200 | logD7.4, AstraZeneca (MoleculeNet) |