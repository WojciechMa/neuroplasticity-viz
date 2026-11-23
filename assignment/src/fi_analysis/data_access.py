"""Allen Patch-seq metadata and the NWB recordings on DANDI."""
import json
import re
import urllib.request
from pathlib import Path

import pandas as pd

from . import config


def download(url: str, target: Path) -> Path:
    """Download url to target unless it is already there."""
    if not url.startswith("https://"):
        raise ValueError(f"refusing non-https URL: {url}")
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        partial = target.with_name(target.name + ".part")
        urllib.request.urlretrieve(url, partial)
        partial.replace(target)
    return target


def fetch_metadata(dest: Path = config.METADATA_DIR) -> None:
    for name, url in config.METADATA_URLS.items():
        download(url, dest / name)


def load_cells(metadata_dir: Path = config.METADATA_DIR) -> pd.DataFrame:
    """One row per Patch-seq cell with a t-type: specimen_id, ephys_session_id, t_type, source, nwb_url.

    nwb_url is only known for the excitatory study; fill the rest with dandi_nwb_urls().
    """
    exc_types = pd.read_csv(metadata_dir / "exc_inferred_met_types.csv", index_col=0)
    manifest = pd.read_csv(metadata_dir / "exc_patchseq_file_manifest.csv")
    ephys = manifest[manifest["technique"] == "electrophysiology"]
    exc = (
        exc_types[["t_type"]]
        .rename_axis("specimen_id")
        .reset_index()
        .merge(ephys[["specimen_id", "file_name", "archive_uri"]], on="specimen_id")
    )
    exc = pd.DataFrame({
        "specimen_id": exc["specimen_id"],
        "ephys_session_id": exc["file_name"].str.extract(r"ses-(\d+)")[0].astype(int),
        "t_type": exc["t_type"],
        "source": "sorensen",
        "nwb_url": exc["archive_uri"],
    })

    inh = pd.read_csv(metadata_dir / "patchseq_metadata_mouse.csv")
    inh = pd.DataFrame({
        "specimen_id": inh["cell_specimen_id"],
        "ephys_session_id": inh["ephys_session_id"],
        "t_type": inh["T-type Label"],
        "source": "gouwens2020",
        "nwb_url": None,
    })

    # 75 L2/3 IT cells appear in both tables; only the excitatory study gives them a t-type.
    cells = pd.concat([exc, inh.dropna(subset=["t_type"])], ignore_index=True)
    return cells.drop_duplicates("specimen_id")


def dandi_nwb_urls(dandiset_id: str = config.INHIBITORY_DANDISET) -> dict[int, str]:
    """Map ephys session ID -> download URL for every NWB asset of a dandiset."""
    urls = {}
    page = f"{config.DANDI_API}/dandisets/{dandiset_id}/versions/draft/assets/?page_size=1000"
    while page:
        with urllib.request.urlopen(page) as response:
            data = json.load(response)
        for asset in data["results"]:
            match = re.search(r"ses-(\d+)", asset["path"])
            if match and asset["path"].endswith(".nwb"):
                urls[int(match.group(1))] = f"{config.DANDI_API}/assets/{asset['asset_id']}/download/"
        page = data["next"]
    return urls
