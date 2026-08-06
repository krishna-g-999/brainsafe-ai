# BrainSafe AI web server.
#
# Pinned to python:3.13-slim because the deployed estimators were fitted under scikit-learn 1.8.0 on
# that line; see requirements.txt for why the pin matters. The image installs dependencies in a
# separate layer from the application so that a code change does not rebuild the scientific stack,
# which is the slow part.
FROM python:3.13-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false

# libxrender and libxext are needed by RDKit's molecule drawing; without them the structure image
# fails at runtime rather than at build time, which is the worst place to discover it.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 libxext6 libsm6 curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY requirements.txt .
RUN pip install --upgrade pip && pip install -r requirements.txt

COPY app.py api.py serve.py ./
COPY src/ ./src/
COPY models_rf/ ./models_rf/
COPY assets/ ./assets/
COPY results/ ./results/
COPY docs/ ./docs/

# Run as an unprivileged user. The application only ever reads from disk.
RUN useradd --create-home --uid 10001 brainsafe && chown -R brainsafe:brainsafe /app
USER brainsafe

# 8501 serves the interface, 8000 the REST API. Both are started by serve.py in one container
# because they share a single copy of the loaded models; running them as separate containers would
# double the memory for no benefit, the models being read-only.
EXPOSE 8501 8000

# Streamlit's own health endpoint, so an orchestrator restarts a wedged container rather than
# leaving it serving errors. The API's health endpoint is checked too, since the interface can be
# healthy while the API thread has died.
HEALTHCHECK --interval=30s --timeout=8s --start-period=90s --retries=3 \
    CMD curl -fsS http://localhost:8501/_stcore/health && curl -fsS http://localhost:8000/health || exit 1

CMD ["python", "serve.py"]
