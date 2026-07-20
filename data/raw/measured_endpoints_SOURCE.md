# Source: measured bioactivity endpoints (labels)

These are the measured-data sources that provide the supervised training labels. The per-target
files currently live in `data/endpoints/` and `data/endpoints_reg/` (the current pipeline location)
and are consolidated into `data/processed/compound_library.csv` by
`src/brainsafe/data/build_compound_library.py`.

## Target activity — ChEMBL
- Database / release: ChEMBL version 37 (release 2026-05-01)
- URL / API: https://www.ebi.ac.uk/chembl/ (REST API)
- Retrieved by: `BS_fetch_endpoints.py`
- Records kept: activities with a defined pChEMBL value; standard types IC50, Ki, Kd, EC50, Potency
- Targets (ChEMBL target IDs): AChE CHEMBL220, BChE CHEMBL1914, BACE1 CHEMBL4822, GSK-3beta CHEMBL262,
  MAO-A CHEMBL1951, MAO-B CHEMBL2039, hERG CHEMBL240, D2 CHEMBL217, A2A CHEMBL251, 5-HT2A CHEMBL224,
  SERT CHEMBL228
- Labelling: per-compound median pChEMBL; active >= 6 (<= 1 uM), inactive < 5 (> 10 uM); 5-6 grey zone dropped
- Licence: ChEMBL data are released under CC BY-SA 3.0

## Blood-brain barrier — B3DB
- Dataset: B3DB (curated blood-brain-barrier permeability database)
- Reference: Meng F, Xi Y, Huang J, Ayers PW. Sci Data 2021;8:289
- Records kept: measured BBB+ / BBB- classification labels
- Retrieved by: `BS_fetch_endpoints.py`

## Antioxidant — ChEMBL DPPH assays
- Source: ChEMBL DPPH radical-scavenging assays; IC50/EC50 converted to pIC50 (-log10 M)
- Retrieved by: `BS_fetch_antioxidant.py`

## Clinical / translational reference (no training label)
- Source: ChEMBL ATC level-1 "N" (nervous system) molecules with a recorded clinical phase
- File: `data/clinical_cns_reference.csv` (504 compounds: name, SMILES, max phase, ATC, disease)
- Role: external reference only; not used as a supervised training label

## Note on standardisation
All SMILES are standardised on load (largest organic fragment, sanitisation, canonical SMILES,
InChIKey) before deduplication; the transformation is applied by the library-builder script and is
fully reproducible.
