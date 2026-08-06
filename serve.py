"""Run the interface and the REST API together in one container.

They are deliberately not separate containers. Both load the same read-only model set, roughly a
gigabyte resident, so splitting them would double the memory for no isolation benefit. The API runs
in a background thread of this process and the Streamlit server in a child process, because Streamlit
insists on owning its own process and does not embed.

Ports are read from the environment so a host that assigns them can do so: PORT for the interface,
which is the convention most platforms use, and API_PORT for the API.

If either dies the other is stopped, so an orchestrator sees a failed container and restarts it,
rather than a half-working one that answers health checks on the surviving half.
"""
from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
UI_PORT = os.environ.get("PORT", "8501")
API_PORT = os.environ.get("API_PORT", "8000")


def run_api():
    sys.argv = ["api.py", API_PORT]
    sys.path.insert(0, str(ROOT))
    import api
    api.main()


def main():
    print(f"starting API on {API_PORT} and interface on {UI_PORT}", flush=True)
    t = threading.Thread(target=run_api, daemon=True, name="brainsafe-api")
    t.start()

    ui = subprocess.Popen([sys.executable, "-m", "streamlit", "run", str(ROOT / "app.py"),
                           f"--server.port={UI_PORT}", "--server.address=0.0.0.0",
                           "--server.headless=true", "--browser.gatherUsageStats=false"])

    def stop(*_):
        ui.terminate()
        sys.exit(0)

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)

    while True:
        code = ui.poll()
        if code is not None:
            print(f"interface exited with {code}; stopping container", flush=True)
            sys.exit(code or 1)
        if not t.is_alive():
            print("API thread died; stopping container", flush=True)
            ui.terminate()
            sys.exit(1)
        time.sleep(5)


if __name__ == "__main__":
    main()
