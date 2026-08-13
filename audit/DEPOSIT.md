# Depositing the models and the source caches

Two archives accompany this repository, because neither belongs in git: the trained estimators are
0.78 GB, and the raw API responses are the only record of what ChEMBL, BindingDB and PubChem returned
on the retrieval dates.

**I cannot do the upload.** It needs your Zenodo credentials and it publishes under your name, so the
last step is yours. Everything up to it is automated, and both manifests are written so that a
download can be verified rather than trusted.

---

## Why the existing DOI must not be reused

`10.5281/zenodo.21858576` holds the **pre-audit** models. Those are not the models in this repository
any more: the panel was retrained with duplicates collapsed, actives withheld by scaffold, thresholds
set from withheld inactives, and the measured negative class recovered. Regeneration changed the
majority of the files the old archive contains.

`models_manifest.json` therefore records `version: 1.1`, an empty `doi`, an empty `urls`, and the old
identifier under `supersedes_doi`. That is deliberate. A manifest pointing at a record which serves
different bytes than it describes is worse than one pointing nowhere: the first fails silently or
downloads the wrong models over the right ones, the second fails honestly.

**If the old record is reused, a fresh deployment will download the pre-audit models and every number
in the manuscript will describe something the server is not running.**

---

## Step 1 — rebuild both archives from what is now on disk

Run these *after* the regeneration is complete, so the checksums describe the final models:

```bash
python src/brainsafe/models/package_models.py
python src/brainsafe/data/package_caches.py
```

They write:

| Archive | Contents | Manifest (committed) |
|---|---|---|
| `dist/brainsafe_models_v1.1.tar.gz` | every deployed estimator and supporting artefact | `models_manifest.json` |
| `dist/brainsafe_source_caches_v1.0.tar.gz` | the raw ChEMBL, BindingDB and PubChem responses | `caches_manifest.json` |

The archives stay out of git (`dist/` is ignored). Only the manifests are committed, which is the
point: they are small, they are versioned with the code, and they carry the SHA-256 of the archive
and of every file inside it.

---

## Step 2 — verify locally before uploading anything

```bash
python src/brainsafe/data/package_caches.py --verify
```

This compares the working tree against the published checksums and exits non-zero on any drift. It
detects a single flipped bit. Do this first: uploading an archive you have not verified means the DOI
certifies whatever happened to be on disk.

---

## Step 3 — deposit

Create a **new version** of the Zenodo record rather than a new record, so the DOI lineage is
preserved and the old version stays citable for anyone reproducing the pre-audit numbers.

1. Open the existing record, choose "New version".
2. Upload `dist/brainsafe_models_v1.1.tar.gz` and `dist/brainsafe_source_caches_v1.0.tar.gz`.
3. In the description, state plainly what changed, so nobody reproduces the old numbers by accident:

   > Version 1.1. Models retrained after an internal audit: duplicate compounds collapsed before
   > splitting, actives withheld by scaffold, decision thresholds set from withheld inactives, and the
   > measured negative class recovered from ChEMBL. Reported performance differs from version 1.0.
   > Version 1.0 remains available for reproducing the earlier figures.

4. Publish, then copy the new DOI and the two file URLs.

---

## Step 4 — record the identifiers, then re-verify

Fill in both manifests:

```jsonc
// models_manifest.json
"doi":  "10.5281/zenodo.XXXXXXXX",
"urls": ["https://zenodo.org/records/XXXXXXXX/files/brainsafe_models_v1.1.tar.gz?download=1",
         "https://zenodo.org/api/records/XXXXXXXX/files/brainsafe_models_v1.1.tar.gz/content"],
```

```jsonc
// caches_manifest.json
"doi":  "10.5281/zenodo.XXXXXXXX",
"urls": ["https://zenodo.org/records/XXXXXXXX/files/brainsafe_source_caches_v1.0.tar.gz?download=1"],
```

Then prove the round trip actually works, from a clean directory:

```bash
mv models_rf models_rf.bak
python model_fetch.py          # downloads, verifies every file, extracts
python src/brainsafe/evaluation/app_health.py
```

`model_fetch.py` checks the archive SHA-256 and then all extracted files, and is fatal on mismatch.
If it passes and the health gate exits 0, the deposit is correct. Restore `models_rf.bak` if not.

---

## Step 5 — cite it where a reader will look

The DOI currently appears in `REPOSITORY_MAP.md`, `AUDIT_PACKAGE/00_READ_ME_FIRST.md` and three
reviewer-package files, and is **absent from `README.md` and from the manuscript**. The manuscript's
data-availability statement points only at GitHub, which does not contain the models. Add the new DOI
to both, or the statement promises something the repository does not hold.

---

## What is still outstanding after this

The 49 expansion endpoints have no `_inactive.json` cache, so the negative-class recovery has been
applied to the eleven core targets only. Until the same query is run for the rest, 37 of 60 endpoints
remain above 90 per cent active and their class balance is an artefact of the query rather than of
the chemistry:

```bash
BRAINSAFE_ALLOW_NONSTRICT_TLS=1 python src/brainsafe/data/fetch_endpoints.py   # core eleven
# the expansion fetchers need the same standard_relation=">" query added
```
