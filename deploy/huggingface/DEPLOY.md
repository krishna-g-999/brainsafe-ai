# Deploying BrainSafe AI to a free public URL

The constraint is the panel. It is roughly 0.85 GB of fitted estimators after compression, and about
0.70 GB of that is needed at runtime. Any host with less memory than that cannot serve the tool, and
that single number decides the options.

| Host | Memory | Disk | Verdict |
|---|---|---|---|
| **Hugging Face Spaces** | 16 GB | 50 GB | **works.** Free, permanent URL, models via git-LFS at no cost |
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

## Creating the Space

1. Create a Space at https://huggingface.co/new-space, SDK **Streamlit**, hardware **CPU basic**
   (free). This gives the URL `https://huggingface.co/spaces/<user>/brainsafe-ai`.

2. Clone it and copy this directory's `README.md` to the Space root. The YAML block at the top is
   what tells Hugging Face which SDK to run and which file is the entry point; without it the Space
   will not start.

3. Copy in what the server needs, and nothing else:

   ```
   app.py  api.py  serve.py  model_fetch.py  models_manifest.json  requirements.txt
   src/  assets/  results/  docs/  data/
   models_rf/            (excluding holdout/)
   ```

4. Track the models with git-LFS before adding them, or the push is rejected for file size:

   ```bash
   git lfs install
   git lfs track "*.joblib" "*.pkl"
   git add .gitattributes
   ```

5. Push. The first build takes several minutes because the scientific stack is large.

## Two things to check after it is live

**The models load.** The first prediction is slow because the panel is read into memory once; every
prediction after it is fast. If the Space restarts on the first query, the memory limit was hit and
the panel needs trimming further.

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
