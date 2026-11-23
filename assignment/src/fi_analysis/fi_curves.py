"""Per-cell f-I curves and features, group f-I curves, and choice of representative cells."""
import numpy as np
import pandas as pd

from . import config, ephys
from .ephys import Sweep


def cell_fi_table(sweeps: list[Sweep]) -> pd.DataFrame:
    """Firing rate per non-negative step amplitude, repeated amplitudes averaged."""
    rows = [(s.amplitude_pa, ephys.firing_rate_hz(s)) for s in sweeps if s.amplitude_pa >= 0]
    fi = pd.DataFrame(rows, columns=["current_pa", "rate_hz"])
    return fi.groupby("current_pa", as_index=False)["rate_hz"].mean().sort_values("current_pa")


def rheobase_pa(fi: pd.DataFrame) -> float:
    return float(fi.loc[fi["rate_hz"] > 0, "current_pa"].min())


def fi_slope_hz_per_pa(fi: pd.DataFrame, max_above_rheobase_pa: float = np.inf) -> float:
    rheo = rheobase_pa(fi)
    supra = fi[(fi["current_pa"] >= rheo) & (fi["current_pa"] <= rheo + max_above_rheobase_pa)]
    if len(supra) < 2:
        return np.nan
    return float(np.polyfit(supra["current_pa"], supra["rate_hz"], 1)[0])


def protocol(sweeps: list[Sweep]) -> str:
    older = any(s.stimulus.startswith(config.OLDER_PROTOCOL_STIMULI) for s in sweeps)
    newer = any(s.stimulus.startswith(config.NEWER_PROTOCOL_STIMULI) for s in sweeps)
    return "both" if older and newer else "older" if older else "newer"


def summarize_cell(sweeps: list[Sweep]) -> dict:
    fi = cell_fi_table(sweeps)
    rheo = rheobase_pa(fi)
    shared = config.SHARED_RANGE_ABOVE_RHEOBASE_PA
    in_shared = fi[(fi["current_pa"] >= rheo) & (fi["current_pa"] <= rheo + shared)]
    hyper = [s for s in sweeps if s.amplitude_pa < 0]
    rheo_sweep = next((s for s in sweeps if s.amplitude_pa == rheo and ephys.firing_rate_hz(s) > 0), None)
    v_thr, v_trough = ephys.spike_threshold_and_trough_mv(rheo_sweep) if rheo_sweep else (np.nan, np.nan)
    return {
        "protocol": protocol(sweeps),
        "n_sweeps": len(sweeps),
        "n_hyperpolarising": len(hyper),
        "n_points_shared_range": len(in_shared),
        "rheobase_pa": rheo,
        "fi_slope_hz_per_pa": fi_slope_hz_per_pa(fi),
        "fi_slope_shared_hz_per_pa": fi_slope_hz_per_pa(fi, shared),
        "max_rate_shared_hz": float(in_shared["rate_hz"].max()) if len(in_shared) else np.nan,
        "max_rate_hz": float(fi["rate_hz"].max()) if len(fi) else np.nan,
        "max_current_above_rheobase_pa": float(fi["current_pa"].max() - rheo) if len(fi) else np.nan,
        "input_resistance_mohm": ephys.input_resistance_mohm(sweeps),
        "tau_ms": float(np.nanmedian([ephys.membrane_tau_ms(s) for s in hyper])) if hyper else np.nan,
        "v_rest_mv": float(np.median([ephys.baseline_mv(s) for s in sweeps])),
        "v_threshold_mv": v_thr,
        "v_reset_mv": v_trough,
    }


def cell_qc_failure(summary: dict) -> str:
    """Reason the cell cannot be used, or an empty string."""
    if summary["n_hyperpolarising"] == 0:
        return "no hyperpolarising steps"
    if np.isnan(summary["rheobase_pa"]):
        return "no spikes"
    if summary["n_points_shared_range"] < 2:
        return "fewer than 2 points in shared range"
    if not 0 < summary["input_resistance_mohm"] < 5000 or np.isnan(summary["tau_ms"]):
        return "passive properties not measurable"
    return ""


def group_fi_curves(points: pd.DataFrame, x: str = "current_pa", bin_width_pa: float = 20.0) -> pd.DataFrame:
    """Mean, SEM and n of firing rate per cell class and current bin; cells are averaged first."""
    binned = points.assign(bin=np.floor(points[x] / bin_width_pa + 0.5) * bin_width_pa)
    per_cell = binned.groupby(["cell_class", "cell_id", "bin"], as_index=False)["rate_hz"].mean()
    return (
        per_cell.groupby(["cell_class", "bin"])["rate_hz"]
        .agg(mean="mean", sem="sem", n="count")
        .reset_index()
    )


def pick_representative(cells: pd.DataFrame, features: tuple[str, ...] = ("rheobase_pa", "fi_slope_shared_hz_per_pa")):
    """Index of the cell closest to the group median in z-scored feature space."""
    x = cells[list(features)].dropna()
    z = (x - x.median()) / x.std(ddof=0)
    return z.pow(2).sum(axis=1).idxmin()


def representative_cells(summaries: pd.DataFrame) -> pd.DataFrame:
    """One row per cell class: the QC-passed cell closest to its class median."""
    good = summaries[summaries["qc_failure"].fillna("") == ""]
    return pd.DataFrame([good.loc[pick_representative(group)] for _, group in good.groupby("cell_class")])


def curve_metrics(current_pa, rate_hz) -> dict[str, float]:
    """Rheobase and shared-range f-I slope of a (model) f-I curve."""
    fi = pd.DataFrame({"current_pa": current_pa, "rate_hz": rate_hz})
    return {
        "rheobase_pa": rheobase_pa(fi),
        "fi_slope_shared_hz_per_pa": fi_slope_hz_per_pa(fi, config.SHARED_RANGE_ABOVE_RHEOBASE_PA),
    }
