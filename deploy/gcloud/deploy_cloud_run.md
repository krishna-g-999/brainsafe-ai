# Deploying BrainSafe AI to Google Cloud Run

Cloud Run is the right home for this server. Streamlit Community Cloud refused it on memory, and
Hugging Face now requires a paid plan for Docker Spaces on new accounts. Cloud Run lets you set the
memory the panel actually needs, and it scales to zero so an idle server costs nothing.

You do not need Docker installed. Cloud Build builds the image in Google's infrastructure from your
source directory.

## What this will cost

Far less than your trial credit. Cloud Run bills only while a request is being served, and scales to
zero in between. A server this size, used by reviewers rather than by traffic, typically costs a few
rupees a day, and the first 2 million requests a month sit in the always-free tier. Your credit is
₹28,694 over 90 days; this will not come close.

Set a budget alert anyway, at step 2. It costs nothing and removes the worry.

## 1. Install the Google Cloud CLI

Download and run the installer: https://cloud.google.com/sdk/docs/install-sdk#windows

Then, in a **new** PowerShell window:

```bash
gcloud init
```

Sign in with the same account, and choose the project `project-47d4c9c7-5bfb-4527-b67`.

## 2. Set a budget alert, so nothing can surprise you

https://console.cloud.google.com/billing → Budgets & alerts → Create budget. Set it to ₹1,000 a
month with alerts at 50 and 100 per cent. This does not cap spending, it emails you.

## 3. Enable the two services this needs

```bash
gcloud services enable run.googleapis.com cloudbuild.googleapis.com artifactregistry.googleapis.com
```

## 4. Deploy

From the project directory:

```bash
cd D:\BRAINSAFE_AI
```

```bash
gcloud run deploy brainsafe-ai --source . --region asia-south1 --allow-unauthenticated --memory 4Gi --cpu 2 --timeout 300 --port 8501 --max-instances 3
```

Answer `y` if it offers to create an Artifact Registry repository.

What each flag is for, since these are the ones that matter:

| Flag | Why |
|---|---|
| `--source .` | Cloud Build builds the image remotely; no local Docker needed |
| `--region asia-south1` | Mumbai, closest to you. Any region works |
| `--allow-unauthenticated` | the server must be public for a journal reviewer to use it |
| `--memory 4Gi` | the panel is roughly 1.8 GB resident; 4 GB leaves headroom. **This is the setting Streamlit could not give you** |
| `--cpu 2` | model loading is the slow part and benefits from two cores |
| `--timeout 300` | allows five minutes for the first request while models load |
| `--port 8501` | the port `serve.py` listens on when `PORT` is unset |
| `--max-instances 3` | a hard ceiling on cost. Three instances is ample for review traffic |

The first build takes 10 to 20 minutes: it installs RDKit and scikit-learn, then downloads and
verifies the 0.79 GB model archive from Zenodo. Later deploys reuse cached layers and are much faster.

## 5. Check it

The deploy prints a URL like `https://brainsafe-ai-<hash>-el.a.run.app`.

```bash
curl -fsS https://<your-url>/_stcore/health
```

Then open the URL and try `donepezil`.

## If the build fails

Read the build log link the command prints. The two likely causes:

- **Cloud Build timeout on a slow layer.** Raise it: `gcloud config set builds/timeout 3600`
- **Out of memory during the build.** Use a larger builder:
  `gcloud run deploy ... --memory 4Gi` unchanged, but add
  `gcloud config set builds/machine_type e2-highcpu-8`

## Serving the REST API as well

Cloud Run exposes one port per service, like the other hosts. Deploy a second service for the API:

```bash
gcloud run deploy brainsafe-api --source . --region asia-south1 --allow-unauthenticated --memory 4Gi --cpu 2 --timeout 300 --port 8501 --max-instances 2 --set-env-vars BRAINSAFE_API_ONLY=1
```

`BRAINSAFE_API_ONLY=1` makes `serve.py` hand the published port to the API instead of the interface.

## Why the models are fetched during the build

`.gitignore` excludes the binaries, so no build context contains them, and `.gcloudignore` keeps them
out of the upload deliberately: sending 0.79 GB from a home connection on every deploy would be slow
for no gain. The Dockerfile runs `model_fetch.py` during the build, which downloads the archive from
the Zenodo DOI, verifies its checksum and the checksum of all 195 extracted files, and fails the
build if either disagrees. The models are then baked into the image, so cold starts are fast and the
transfer happens once per build rather than once per instance.
