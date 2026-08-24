# Container for the Hugging Face Space.
#
# This exists because the Streamlit SDK was retired: the Hub now accepts only gradio, docker or
# static, so Streamlit runs inside a container we define rather than one the platform provides. The
# application is unchanged. Everything below reproduces what the SDK used to do implicitly, which is
# also why each line is worth stating.
#
# Python 3.13 matches the interpreter the deployed estimators were fitted under. scikit-learn records
# its version inside each pickle and warns on a mismatch, and while that warning was checked and
# found to cost one value in 630 a difference of 1e-12, the point of a container is that it does not
# have to be checked again.
FROM python:3.13-slim

# RDKit's wheels link against these at import time. Without them the import fails with a bare
# "libXrender.so.1: cannot open shared object file", which reads like a broken install rather than a
# missing system package.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libxrender1 libxext6 libsm6 libglib2.0-0 \
    && rm -rf /var/lib/apt/lists/*

# Spaces run as uid 1000. Building as root and copying without --chown leaves the application files
# owned by root and unreadable at runtime, which surfaces much later as a permission error inside
# Streamlit rather than at build time.
RUN useradd -m -u 1000 user
USER user
ENV PATH="/home/user/.local/bin:$PATH" \
    PYTHONUNBUFFERED=1 \
    STREAMLIT_SERVER_HEADLESS=true \
    STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
WORKDIR /home/user/app

# Requirements first, so a change to the application does not invalidate the layer that installs a
# large scientific stack.
COPY --chown=user:user requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

COPY --chown=user:user . .

# 7860 is the port a Space exposes, and it must be bound on 0.0.0.0: bound to localhost the container
# starts cleanly and the Space stays blank, which is a slow thing to diagnose from outside.
EXPOSE 7860
CMD ["streamlit", "run", "app.py", \
     "--server.port=7860", \
     "--server.address=0.0.0.0", \
     "--server.enableXsrfProtection=false", \
     "--server.enableCORS=false"]
