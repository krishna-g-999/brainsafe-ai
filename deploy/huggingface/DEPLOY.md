# Deploying BrainSafe AI to a free public URL

The constraint is the panel. It is roughly 0.85 GB of fitted estimators on disk after compression,
and a warmed process holds **2.57 GB resident**: 2.31 GB once every model is loaded, rising to 2.57 GB
after the applicability-domain reference and the read-across index are also in memory. That measured
figure, not the size on disk, decides the options, and it is the number to check against any host's
limit. Measure it again after any retrain with `psutil` rather than trusting this line.

| Host | Memory | Disk | Verdict |
|---|---|---|---|
| **Hugging Face Spaces** | 16 GB | 50 GB | **works**, with 6x headroom over the measured 2.57 GB. Free, permanent URL, models via git-LFS at no cost |
| Google Cloud Run | configurable | image | works, scales to zero, but needs a billing account attached even to stay inside the free tier |
| Oracle Cloud Always Free | 24 GB | 200 GB | works, genuinely free with no time limit, but you administer the machine |
| Streamlit Community Cloud | 1 GB | small | **does not fit.** The panel exceeds the memory limit once loaded |
| Render / Fly.io free tiers | 256 MB - 1 GB | small | do not fit |

Hugging Face Spaces is the recommendation: it is the only one that is free without qualification,
needs no card, and is designed for exactly this shape of application.

## Before deploying: shrink the panel

Do this first, or the upload is three times larger than it needs to be.

```bash
python tools/compress_models.py            # report what would be saved, change nothing
python tools/compress_models.py --apply    # recompress in place, verifying each file
python src/brainsafe/models/package_models.py 1.1   # regenerate the manifest
```

joblib writes an uncompressed pickle unless asked otherwise. The binder models were saved with
`compress=3` and average 3 MB; the calibrated classifiers were saved without and average 117 MB.
Recompressing takes the panel from 2.05 GB to 0.85 GB, and the tool verifies every file it rewrites:
each is reloaded and scored against the original, and replaced only if predictions agree to 1e-12.
Anything that cannot be verified is left exactly as it was.

`models_rf/holdout/` is another 155 MB and is not needed to serve predictions. It holds the
scaffold-split twins used for validation. Exclude it.

## Creating the Space under your institute organisation

You have already joined the organisation, which is the part people usually get stuck on. What
follows assumes the Space belongs to the organisation rather than to you personally, so that it
survives you changing accounts and carries the institute's name in the URL.

**1. Check you can write to the organisation.** Open
`https://huggingface.co/organizations/<org>/settings/members` and confirm your role is `write` or
`admin`. With `read` you can see the organisation but the Owner dropdown in step 2 will not offer
it, and nothing later will work. If you only have `read`, an admin has to raise it.

**2. Create the Space.** Go to https://huggingface.co/new-space and set:

| Field | Value |
|---|---|
| Owner | **the organisation**, not your username. This is the dropdown people miss |
| Space name | `brainsafe-ai` |
| Licence | MIT |
| SDK | **Streamlit** |
| Hardware | CPU basic, 2 vCPU, 16 GB. Free |
| Visibility | Public |

The URL is then `https://huggingface.co/spaces/<org>/brainsafe-ai`. That is the address for the
manuscript.

**3. Get a token that can write to the organisation.** A personal token is not automatically an
organisational one. At https://huggingface.co/settings/tokens create a token with **Write** access,
and if you choose a fine-grained token, tick the organisation and give it write permission on
repositories. Copy it; it is shown once.

**4. Clone the Space and fill it.**

```bash
git clone https://huggingface.co/spaces/<org>/brainsafe-ai
python deploy/huggingface/prepare_space.py --out brainsafe-ai
```

`prepare_space.py` copies only what answers a query: the app, `src/`, `assets/`, `results/`,
`docs/`, the four `data/` subdirectories the server reads, and `models_rf/` without `holdout/`. It
also writes the Space card and the git-LFS rules. That is 0.81 GB. Copying the repository instead
would be 1.36 GB, most of it raw pulls and API caches that are never opened to serve a prediction.

It ships `deploy/huggingface/requirements.txt`, not the repository's. The root file is the
environment that trains and validates, and carries matplotlib, python-docx, `pypandoc_binary` and
xgboost, none of which is imported to answer a query; `pypandoc_binary` alone pulls a pandoc binary
of over a hundred megabytes. The runtime set is nine packages, pinned to the versions the estimators
were fitted under. streamlit is deliberately not among them, because the Space card's `sdk_version`
installs it and naming it twice invites a conflict. Keep `sdk_version` in the card equal to the
streamlit version this project validates against; it is 1.58.0 today.

**Verify the assembled directory before pushing.**

```bash
python deploy/huggingface/verify_space.py brainsafe-ai
```

This runs the app with the Space as the working directory and the repository's source roots removed
from the path, loads every model, and predicts for three compounds chosen to exercise different
paths. It exists because a file that failed to travel is invisible in a directory listing and
obvious to the first visitor. It also checks that `results/tables/external_novelty_strata.csv`
arrived, since the interface quotes an expected recall from it and would drop that row without
comment if it were missing.

**5. Track LFS before adding the models.** This ordering matters more than anything else here. Git
refuses single files over 10 MB on the Hub, and once a large file is in a commit, adding LFS
afterwards does not fix that commit.

```bash
cd brainsafe-ai
git lfs install
git add .gitattributes && git commit -m "Track model files with LFS"
git add -A && git commit -m "BrainSafe AI"
git push
```

When prompted, the username is your Hugging Face username and the **password is the token** from
step 3, not your account password. Alternatively authenticate once and let git use the stored
credential, which avoids pasting the token into a terminal prompt where it may be logged:

```bash
brainsafe_env/Scripts/hf.exe auth login
```

The token is entered into the Hugging Face client, which stores it under `~/.cache/huggingface`.

The push moves 0.8 GB and takes a while. The first build then takes several minutes because the
scientific stack is large.

**6. Silence the model-fetch warning.** `model_fetch.py` compares what is on disk against
`models_manifest.json`, which lists the hold-out files that were deliberately not shipped, so the
log opens with "102 of 252 model files missing". The app runs correctly regardless. To keep the log
clean, add a Space variable under Settings, Variables and secrets:

    BRAINSAFE_SKIP_MODEL_FETCH = 1

## Two things to check after it is live

**The models load.** The first prediction is slow because the panel is read into memory once; every
prediction after it is fast. If the Space restarts on the first query, the memory limit was hit and
the panel needs trimming further, though at 2.57 GB measured against 16 GB available that should not
happen.

**The expected-recall row appears.** Submit a compound far from the training chemistry, a steroidal
natural product will do, and confirm the result carries both an applicability-domain distance and
the measured recall at that distance. If the recall row is missing, the novelty-strata table did not
travel and `verify_space.py` was skipped.

**The URL goes in the manuscript.** NAR requires the server address in the abstract. It is currently
`[SERVER URL TO BE SUPPLIED]` in `manuscript/NAR_WebServer_BrainSafe_draft.md`, and the manuscript
must be rebuilt after it is filled in.

## The alternative, if the Space memory is ever a problem

`model_fetch.py` already downloads the panel from a URL recorded in `models_manifest.json` and
verifies it against the committed SHA-256 of the archive and of every file inside it. Publishing the
archive as a Hugging Face model repository and recording that URL keeps the Space small and the
models versioned separately, which is also what the manuscript's data-availability statement
describes. That is the tidier arrangement long term; putting the models in the Space is the faster
one to get a working URL today.
