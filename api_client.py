import os
import time
import requests
import streamlit as st


def _get_backend_url():
    """Reads the docking API base URL from Streamlit secrets or an environment variable."""
    url = None
    try:
        url = st.secrets.get("DOCKING_API_URL")
    except Exception:
        url = None
    return (url or os.environ.get("DOCKING_API_URL", "")).rstrip("/")


def submit_cloud_docking(smiles, target_variant, max_wait_seconds=600, poll_every=3):
    """
    Submits a docking job to the remote backend and waits for the result.

    Expected backend contract:
      POST {base}/dock            body: {"smiles": "...", "target": "..."}  ->  {"job_id": "..."}
      GET  {base}/jobs/{job_id}   -> {"status": "queued|running|done|failed",
                                      "results": [ {"Ligand": ..., "Target": ...,
                                                    "Affinity (kcal/mol)": float, ...}, ... ],
                                      "error": "..."}

    Returns a list of result rows (dicts), or None on any failure (after showing an error).
    """
    base = _get_backend_url()
    if not base:
        st.error("Docking backend is not configured. Set DOCKING_API_URL in Streamlit secrets "
                 "(Manage app > Settings > Secrets) or as an environment variable.")
        return None

    try:
        resp = requests.post(
            f"{base}/dock",
            json={"smiles": smiles, "target": target_variant},
            timeout=30,
        )
        resp.raise_for_status()
        job_id = resp.json().get("job_id")
        if not job_id:
            st.error("Backend did not return a job_id.")
            return None

        waited = 0
        while waited <= max_wait_seconds:
            status_resp = requests.get(f"{base}/jobs/{job_id}", timeout=30)
            status_resp.raise_for_status()
            data = status_resp.json()
            status = data.get("status")

            if status == "done":
                return data.get("results") or None
            if status == "failed":
                st.error(f"Docking job failed: {data.get('error', 'unknown error')}")
                return None

            time.sleep(poll_every)
            waited += poll_every

        st.error("Docking job timed out. Try again later.")
        return None

    except requests.exceptions.RequestException as e:
        st.error(f"Could not reach the docking backend: {e}")
        return None
    except ValueError:
        st.error("The docking backend returned an invalid response.")
        return None  
