# Deploying BrainSafe AI

The server is a single container running two things that share one copy of the loaded models: the
Streamlit interface on `$PORT` (default 8501) and the REST API on `$API_PORT` (default 8000).

Everything below is prepared and tested locally. The one step that cannot be done from here is
creating an account with a host and pressing deploy, because that needs your credentials. Each option
is written so that it is a single command or a single form.

## What must be true before you deploy

```bash
brainsafe_env/Scripts/python.exe src/brainsafe/evaluation/app_health.py
```

Twenty checks, exit code 0 required. It verifies that every declared dependency resolves at its
pinned version, that all model artefacts load, that the knowledge graph is internally consistent,
that chemically unrelated compounds give distinct and directionally correct profiles, that every
export format is well formed and self-contained, that the coverage panel matches the deployed panel,
and that the API answers. Do not deploy a build that fails it.

## Option 1: Hugging Face Spaces (recommended, free, no card)

Best fit for a journal submission: free, no payment method, gives a permanent public URL, and the
model files are handled by git-lfs rather than needing separate hosting.

1. Create a Space at https://huggingface.co/new-space, SDK **Docker**, visibility **Public**.
2. Add the Space as a remote and push:

```bash
git remote add space https://huggingface.co/spaces/<your-user>/brainsafe-ai
git push space main
```

The `Dockerfile` is used as-is. Spaces sets `$PORT`, which `serve.py` reads. The URL becomes
`https://<your-user>-brainsafe-ai.hf.space` and goes in the manuscript.

### The model files, which are the real obstacle

`.gitignore` excludes `*.joblib`, so the repository tracks 64 metadata files and **no model
binaries**. A clone cannot run the server. This is a larger problem than not having a URL and must be
solved before any deployment succeeds.

The binaries are now 0.78 GB, down from 1.8 GB: `src/brainsafe/models/compress_models.py`
recompressed every model and verified per model that predictions are unchanged, checking the round
trip against each estimator's own run-to-run variation rather than demanding bit equality, because
these forests are not bit-deterministic between calls.

Two ways to get them to the host, in order of preference:

1. **Track them with git-lfs on the Space.** 0.78 GB fits comfortably within a Hugging Face Space.
   Do this in the Space clone, not in the GitHub repository, whose free LFS quota is 1 GB total:

   ```bash
   git lfs install
   git lfs track "models_rf/**/*.joblib" "models_rf/*.pkl"
   git add .gitattributes models_rf -f      # -f overrides .gitignore
   git commit -m "Add model binaries via LFS for deployment"
   git push space main
   ```

2. **Publish a release archive and fetch at container start.** Preferable if you want the models
   citable, and it is what `README.md` already anticipates by naming Zenodo. Upload
   `models_rf/` as one archive, then add a fetch step to the Dockerfile before `COPY models_rf/`.
   A Zenodo DOI for the model set is worth having in the manuscript regardless of which route you
   take.

## Option 2: Streamlit Community Cloud (free, simplest, interface only)

Fastest route to a URL, but it runs `app.py` directly and does **not** run the container, so the REST
API is not served. Acceptable if you only need the interface reviewed.

1. https://share.streamlit.io, connect the GitHub repository.
2. Main file `app.py`, Python 3.13, `requirements.txt` is picked up automatically.

## Option 3: any container host (Render, Railway, Fly.io, Cloud Run)

```bash
docker build -t brainsafe-ai .
docker run -p 8501:8501 -p 8000:8000 brainsafe-ai
```

Then push the image to the host's registry. All four read `$PORT`. Give the instance at least 2 GB of
memory: the models are about 1 GB resident and the applicability-domain index adds several hundred
megabytes.

## After deploying

Confirm both halves are live:

```bash
curl -fsS https://<your-url>/_stcore/health
curl -fsS https://<your-url>:8000/health
```

Then put the URL in three places, all of which currently read `[URL]`:

- `manuscript/NAR_WebServer_BrainSafe_draft.md`, the Abstract and the Data availability section
- `README.md`
- the `source` field returned by `GET /` on the API

## The REST API

```
GET  /                 documentation, as JSON
GET  /health           liveness
GET  /version          model counts and the knowledge-graph fingerprint
GET  /targets          every deployed endpoint with threshold and measured sensitivity
GET  /predict?q=...    full profile for one compound, by name or SMILES
POST /batch            {"compounds": [...]} up to 300
```

```bash
curl "https://<your-url>:8000/predict?q=rimegepant"
curl -X POST https://<your-url>:8000/batch \
     -H "Content-Type: application/json" \
     -d '{"compounds":["donepezil","fluoxetine","CC(=O)Nc1ccc(O)cc1"]}'
```

The API is built on the Python standard library rather than a web framework. That was forced by this
environment, where the network appliance blocks `pydantic_core` as an unscannable binary, but it is
also the better choice for a container whose value rests on a pinned, reproducible scientific stack:
a read-only JSON service over already-loaded models gains nothing from a framework and would gain a
dependency to pin.

## What is deliberately not here

No authentication, no rate limiting, no request logging of submitted structures. The service is
read-only and stateless, and a submitted structure is not written to disk. If you deploy behind an
institutional gateway, add rate limiting there rather than in the application.
