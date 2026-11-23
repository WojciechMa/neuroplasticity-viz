"""Step 1: pick cells with both ephys and a t-type label, and download their NWB files from DANDI."""
from concurrent.futures import ThreadPoolExecutor

from fi_analysis import config, data_access
from fi_analysis.cell_class import classify_ttype

CELLS_PER_CLASS = 50
DOWNLOAD_WORKERS = 4


def main() -> None:
    data_access.fetch_metadata()
    cells = data_access.load_cells()
    cells["cell_class"] = cells["t_type"].map(classify_ttype)
    cells = cells.dropna(subset=["cell_class"])
    sample = cells.sample(frac=1, random_state=0).groupby("cell_class").head(CELLS_PER_CLASS).copy()

    missing = sample["nwb_url"].isna()
    sample.loc[missing, "nwb_url"] = sample.loc[missing, "ephys_session_id"].map(data_access.dandi_nwb_urls())
    sample = sample.dropna(subset=["nwb_url"])

    sample["nwb_file"] = sample["ephys_session_id"].astype(str) + ".nwb"
    with ThreadPoolExecutor(DOWNLOAD_WORKERS) as pool:
        ok = list(pool.map(_fetch, sample["nwb_url"], sample["nwb_file"]))
    sample = sample[ok]
    print(sample["cell_class"].value_counts().to_string())

    config.PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    sample.to_csv(config.PROCESSED_DIR / "selected_cells.csv", index=False)


def _fetch(url: str, name: str) -> bool:
    try:
        data_access.download(url, config.RAW_DIR / name)
    except OSError as err:
        print(f"skipped {name}: {err}", flush=True)
        return False
    print(f"have {name}", flush=True)
    return True


if __name__ == "__main__":
    main()
