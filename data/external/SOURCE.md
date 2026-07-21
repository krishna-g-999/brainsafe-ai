# Source: external compound libraries (structures, not training labels)

These libraries were transferred from the group HPC cluster
(`/home/krishnasalini-rs/...`, host 10.110.90.30) on 2026-07-20 and standardised by
`src/brainsafe/data/integrate_external.py`. They provide chemical structures only; none carries a
measured value for the twelve modelled endpoints, so they are **never used as training labels**.
Their roles are `evaluation` (held-out external testing) and `coverage` (applicability-domain and
the flavonoid / natural-product panel). The one measured field present anywhere in this set is the
blood-brain-barrier status in the FDA-curated BBB file, which is used as an external BBB test.

## Files and provenance

| File (raw) | Rows in | Unique after standardisation | Role | Notes |
|---|---|---|---|---|
| `drugbank_smiles.tsv` | 12,313 | 11,723 | evaluation | DrugBank approved + experimental drugs; 346 biologics / > 900 Da removed. Columns src,id,smiles. |
| `BBB_Final_GNN_Dataset_2026_CLEAN_v5.csv` | 1,690 | 1,683 | evaluation | FDA-curated BBB set (source tag FDA_Curated_KS). Measured `bbb_status` (1 permeable / 0 non-permeable); 306 not present in the B3DB training set. |
| `coconut_smiles.tsv` | 695,119 | 37,647 flavonoid hits | coverage | COCONUT natural-product database; flavonoid-core compounds extracted by substructure. |
| `natural_all_CNS.tsv` | 217,480 | 214,740 | coverage | COCONUT CNS-like natural products (carries MW/cLogP/TPSA in source). |
| `pubchem_natural_smiles_clean.tsv` | 2,450 | 21 flavonoid hits | coverage | PubChem natural products; folded into the flavonoid panel. |
| `chembl34_smiles_clean.tsv` | 1,673,554 | 1,673,554 (counted) | coverage_pool | ChEMBL 34 structures; recorded as a broad chemical-space pool, standardised on demand. |

Not used (inspected and rejected): `approved_drugs_targets.tsv` (5 biologics, SMILES = "Biologic"),
`COMPREHENSIVE_DRUG_DATABASE.tsv` (62 rheumatology biologics, no SMILES), `drug_smiles.tsv`
(single compound).

## Standardised outputs (`data/external/processed/`)

- `external_drugs.csv` - 11,723 DrugBank small molecules (canonical SMILES, InChIKey, descriptors).
- `external_bbb_test.csv` - 1,683 FDA-curated BBB compounds with measured label and an
  `in_b3db_training` flag so the 306 novel compounds can be used as a clean external test.
- `flavonoid_panel.csv` - 37,647 flavonoid-core compounds (COCONUT + PubChem natural).
- `natural_products_coverage.csv` - 214,740 CNS-like natural products.
- `external_summary.csv` - the counts above.

## Standardisation

Identical to the measured pipeline: largest organic fragment kept, sanitised, canonical SMILES and
InChIKey computed, deduplicated by InChIKey (`src/brainsafe/data/build_compound_library.py`).

## Licences

DrugBank: Creative Common­s Attribution-NonCommercial 4.0 (academic use). COCONUT: CC BY 4.0.
ChEMBL 34: CC BY-SA 3.0. PubChem: public domain. The FDA-curated BBB file is an internal curation
of publicly reported blood-brain-barrier outcomes for approved drugs.
