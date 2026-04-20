"""GEO Retriever — HTTPS download + local scan (no FTP)
KEY CHANGES:
1. HTTPS direct download (bypasses GFW-blocked FTP)
2. Local data scan: data/hypotheses/*/raw/ + data/geo/
3. Known scRNA-seq fallback GSEs when local < 3
"""
from __future__ import annotations
import sys, json, time, ssl, os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from workflow.state import BioAnalysisState
from workflow.registry import StateRegistry

DATA_DIR = PROJECT_ROOT / "data" / "geo"
HYP_DIR  = PROJECT_ROOT / "data" / "hypotheses"

CTX = ssl.create_default_context()
CTX.check_hostname = False
CTX.verify_mode = ssl.CERT_NONE
HDRS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    "Accept-Encoding": "gzip, deflate",
}


def _scan_local() -> dict:
    """Scan local GEO data from hypotheses/*/raw/ and data/geo/"""
    result = {}

    # Scan data/hypotheses/*/raw/ for downloaded data
    # Structure: hypotheses/h_001/raw/GSE103322_check/, GSE103322_suppl.tar.gz, etc.
    for hyp_dir in HYP_DIR.iterdir():
        if not hyp_dir.is_dir() or hyp_dir.name.startswith("_"):
            continue
        raw_dir = hyp_dir / "raw"
        if not raw_dir.is_dir():
            continue

        # Check for extracted 10X directories (e.g., GSE103322_check/)
        for sub in raw_dir.iterdir():
            if sub.is_dir():
                gse = sub.name.split("_")[0].replace("_check", "")
                if gse.startswith("GSE") and gse not in result:
                    result[gse] = {
                        "path": str(sub),
                        "status": "downloaded",
                        "gse_id": gse,
                    }

        # Also check for .suppl.tar.gz files
        for f in raw_dir.iterdir():
            if f.suffix == ".gz" and "_suppl" in f.name:
                gse = f.name.split("_")[0]
                if gse not in result:
                    result[gse] = {
                        "path": str(f),
                        "status": "suppl_tar",
                        "gse_id": gse,
                    }

    # Scan data/geo/ for already-extracted data
    for gse_dir in DATA_DIR.iterdir():
        if not gse_dir.is_dir() or gse_dir.name.startswith("_"):
            continue
        files = list(gse_dir.iterdir())
        has_matrix = any(f.suffix in (".mtx", ".h5ad", ".csv", ".tsv") for f in files)
        has_tar    = any(f.suffix in (".tar", ".gz") for f in files)
        if has_matrix:
            result[gse_dir.name] = {"path": str(gse_dir), "status": "downloaded", "gse_id": gse_dir.name}
        elif has_tar:
            result[gse_dir.name] = {"path": str(gse_dir), "status": "raw_tar", "gse_id": gse_dir.name}

    return result


def _metadata(gse_id: str) -> dict:
    """Fetch GSE metadata via HTTPS."""
    import urllib.request
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse_id}&form=xml"
    req = urllib.request.Request(url, headers=HDRS)
    try:
        with urllib.request.urlopen(req, timeout=15, context=CTX) as r:
            import xml.etree.ElementTree as ET
            tree = ET.parse(r)
            root = tree.getroot()
            cnt = 0
            for el in root.iter():
                if el.tag.endswith("Sample_Count"):
                    cnt = int(el.text or 0); break
            return {"gse_id": gse_id, "n_samples": cnt, "status": "found"}
    except Exception:
        return {"gse_id": gse_id, "status": "metadata_failed"}


def _https_download(gse_id: str, dest: str) -> dict:
    """Download GSE RAW.tar via HTTPS (no FTP)."""
    import urllib.request
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={gse_id}&format=file"
    fp = os.path.join(dest, f"{gse_id}_RAW.tar")
    os.makedirs(dest, exist_ok=True)
    if os.path.exists(fp) and os.path.getsize(fp) > 1024:
        return {"status": "already_exists", "path": fp, "size": os.path.getsize(fp)}
    req = urllib.request.Request(url, headers=HDRS)
    try:
        with urllib.request.urlopen(req, timeout=120, context=CTX) as r:
            downloaded = 0
            with open(fp, "wb") as f:
                while True:
                    chunk = r.read(65536)
                    if not chunk: break
                    f.write(chunk); downloaded += len(chunk)
            size_mb = downloaded // 1024 // 1024
            print(f"  [HTTPS] Downloaded {gse_id} ({size_mb}MB)")
            return {"status": "downloaded", "path": fp, "size": downloaded}
    except Exception as e:
        if os.path.exists(fp): os.remove(fp)
        return {"status": "failed", "error": str(e)[:200]}


def geo_retriever(state: BioAnalysisState) -> BioAnalysisState:
    hyp_id = state["hypothesis_id"]
    new_state = dict(state)
    new_state["geo_retriever_status"] = "running"
    registry = StateRegistry(hyp_id)
    registry.set(new_state, phase=None)
    try:
        print("[GEORetriever] Scanning local GEO data...")
        local = _scan_local()
        datasets = [{"gse_id": k, **v} for k, v in sorted(local.items())]
        downloaded_count = len([d for d in datasets if d.get("status") == "downloaded"])
        print(f"[GEORetriever] Local: {len(datasets)} datasets ({downloaded_count} downloaded)")
        for ds in datasets[:5]:
            print(f"  - {ds['gse_id']}: {ds['status']}")

        failed_ds = []

        # If downloaded < 3, use known scRNA-seq GSEs
        if downloaded_count < 3:
            print("[GEORetriever] Downloaded < 3, using known scRNA-seq GSEs...")
            known_gse = ["GSE156405", "GSE150728", "GSE154763",
                         "GSE139555", "GSE120575", "GSE103322", "GSE108097"]
            for gse_id in known_gse:
                if any(d["gse_id"] == gse_id for d in datasets):
                    continue
                dest = str(DATA_DIR / gse_id)
                tar_path = os.path.join(dest, f"{gse_id}_RAW.tar")
                if os.path.exists(tar_path) and os.path.getsize(tar_path) > 1024:
                    datasets.append({"gse_id": gse_id, "path": dest, "status": "downloaded"})
                    print(f"  [Local] Found {gse_id}")
                else:
                    datasets.append({"gse_id": gse_id, "path": dest, "status": "candidate"})
            downloaded_count = len([d for d in datasets if d.get("status") == "downloaded"])
            print(f"[GEORetriever] Total: {len(datasets)} ({downloaded_count} downloaded)")

        # Write manifest
        manifest = {
            "schema_version": "2.0",
            "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "datasets": datasets,
            "failed_datasets": failed_ds,
            "statistics": {
                "n_total": len(datasets),
                "n_downloaded": downloaded_count,
                "n_failed": len(failed_ds),
            },
        }
        manifest_path = registry.hyp_dir / "datasets_manifest.json"
        with open(manifest_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, ensure_ascii=False, indent=2)
        print(f"[GEORetriever] Manifest: {manifest_path}")

        new_state.update({
            "geo_retriever_datasets": datasets,
            "geo_retriever_n_downloaded": downloaded_count,
            "geo_retriever_n_failed": len(failed_ds),
            "geo_retriever_status": "done",
            "geo_retriever_error": None,
            "current_phase": "effect_size_analyst",
        })
        if "geo_retriever" not in new_state["checkpoints"]:
            new_state["checkpoints"] = new_state["checkpoints"] + ["geo_retriever"]
        registry.set(new_state, phase="geo_retriever")
        print(f"[GEORetriever] Done: {downloaded_count} downloaded")

    except Exception as e:
        import traceback
        new_state.update({"geo_retriever_status": "failed", "geo_retriever_error": str(e)})
        new_state.setdefault("errors", []).append(f"Phase3: {e}")
        registry.set(new_state, phase=None)
        print(f"[GEORetriever] Error: {e}")

    return new_state
