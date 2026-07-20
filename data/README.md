# Data

Layout and provenance rules for all data in this project.

## Folders

- **raw/** — data exactly as downloaded from a source, never edited by hand. Each source lives in its
  own subfolder and carries a `SOURCE.md` (see template below). Raw files are the single point of
  truth for provenance.
- **interim/** — standardised intermediates produced by scripts: canonical SMILES, InChIKey
  deduplication, unit harmonisation. Regenerable from `raw/` by the curation scripts.
- **processed/** — the final modelling tables, one CSV per endpoint, containing the exact rows and
  columns used to train and test a model (SMILES, encoded features or a pointer to them, label,
  provenance flag, cross-validation group). These are what the models consume.
- **external/** — reference and evaluation sets held apart from training (clinical-precedent
  compounds, flavonoid panels, approved-drug external test sets).

## Original compounds versus training variants

Two records are kept for every compound so the transformation is fully traceable:
1. the **original** entry as it arrived from the source (in `raw/`), and
2. the **standardised training variant** actually used (in `interim/` and `processed/`), i.e. the
   canonical SMILES after salt stripping and deduplication, with its numeric label.
A per-compound provenance table links the two and records the source and role of each compound.

## Role of each compound

Every compound is tagged with a role so annotation-only sources are never used as training labels:
- `label` — contributes a measured endpoint label (ChEMBL pChEMBL, B3DB permeability, DPPH pIC50);
- `evaluation` — used only for external testing (for example DrugBank/FDA approved drugs);
- `coverage` — used only to widen chemical diversity / the applicability domain.

## SOURCE.md template (copy into each raw/<source>/ folder)

```
# Source: <database or dataset name>
- Version / release: <e.g. ChEMBL 37, release 2026-05-01>
- URL: <download URL or API endpoint>
- Access date: <YYYY-MM-DD>
- Retrieved by: <script name, e.g. src/brainsafe/data/fetch_chembl.py>
- Records: <count>
- Fields kept: <list>
- Licence / terms: <e.g. CC BY-SA 3.0 for ChEMBL>
- Notes: <query filters, standard types retained, any exclusions>
```
