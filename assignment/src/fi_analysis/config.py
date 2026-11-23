from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
METADATA_DIR = DATA_DIR / "metadata"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"

# Mouse visual cortex Patch-seq. Inhibitory cells: Gouwens et al. 2020; excitatory cells: Sorensen et al.
_ALLEN_CDN = "https://cdn.prod.website-files.com/689cfbd308fa7373b604d290"
_EXC_REPO = "https://raw.githubusercontent.com/AllenInstitute/exc_vis_manuscript/main"
METADATA_URLS = {
    "patchseq_metadata_mouse.csv": f"{_ALLEN_CDN}/68d61b068b5cfcfd77c92050_20200711_patchseq_metadata_mouse.csv",
    "exc_inferred_met_types.csv": f"{_EXC_REPO}/derived_data/inferred_met_types.csv",
    "exc_patchseq_file_manifest.csv": f"{_EXC_REPO}/exc_patchseq_file_manifest.csv",
}
DANDI_API = "https://api.dandiarchive.org/api"
INHIBITORY_DANDISET = "000020"

# 1-s current steps. The 0.5-s search steps (X2LP, Y2LP, Rheobase_finder) are left out.
OLDER_PROTOCOL_STIMULI = ("C1LSCOARSE", "C1LSFINEST", "subthreshold", "Rheo_Supra")
NEWER_PROTOCOL_STIMULI = ("X1PS_SubThresh", "X3LP_Rheo", "X4PS_SupraThresh", "Y3LP_Rheo", "Y4PS_SupraThresh")
LONG_SQUARE_STIMULI = OLDER_PROTOCOL_STIMULI + NEWER_PROTOCOL_STIMULI
LONG_SQUARE_DURATION_S = 1.0

# Every cell has 1-s steps up to at least this far above its rheobase; E vs I statistics use only this range.
SHARED_RANGE_ABOVE_RHEOBASE_PA = 60.0

SPIKE_DVDT_THRESHOLD_MV_PER_MS = 20.0
SPIKE_MIN_PEAK_MV = -20.0
SPIKE_REFRACTORY_S = 0.002
BASELINE_WINDOW_S = 0.1
STEADY_STATE_WINDOW_S = 0.1

# Sweep QC on the pre-step baseline
BASELINE_MAX_DEVIATION_MV = 5.0
BASELINE_MAX_NOISE_MV = 1.0
