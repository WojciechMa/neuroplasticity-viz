"""Step 2: per-cell f-I curves and features, group f-I curves and E vs I statistics."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from fi_analysis import config, ephys, fi_curves

MIN_CELLS_PER_BIN = 5
CLASS_COLORS = {"E": "tab:red", "I": "tab:blue"}
FEATURES = {
    "rheobase_pa": "rheobase (pA)",
    "fi_slope_shared_hz_per_pa": "f-I slope, shared range (Hz/pA)",
    "max_rate_shared_hz": "rate at rheobase + 60 pA (Hz)",
    "input_resistance_mohm": "input resistance (MOhm)",
    "tau_ms": "membrane time constant (ms)",
    "v_rest_mv": "resting potential (mV)",
    "v_threshold_mv": "spike threshold (mV)",
    "v_reset_mv": "post-spike trough (mV)",
}


def main() -> None:
    cells = pd.read_csv(config.PROCESSED_DIR / "selected_cells.csv")
    points, summaries = [], []
    for cell in cells.itertuples():
        loaded = ephys.load_long_square_sweeps(config.RAW_DIR / cell.nwb_file)
        sweeps = ephys.qc_sweeps(loaded)
        ids = {"cell_id": cell.ephys_session_id, "cell_class": cell.cell_class, "t_type": cell.t_type}
        summary = fi_curves.summarize_cell(sweeps)
        summary["n_sweeps_failed_qc"] = len(loaded) - len(sweeps)
        summary["qc_failure"] = fi_curves.cell_qc_failure(summary)
        summaries.append({**ids, **summary})
        points.append(fi_curves.cell_fi_table(sweeps).assign(cell_id=cell.ephys_session_id, cell_class=cell.cell_class))

    summaries = pd.DataFrame(summaries)
    summaries.to_csv(config.PROCESSED_DIR / "cell_summaries.csv", index=False)
    failed = summaries[summaries["qc_failure"] != ""]
    print(f"cells failing QC: {len(failed)}")
    print(failed.groupby(["cell_class", "qc_failure"]).size().to_string())

    good = summaries[summaries["qc_failure"] == ""]
    print("\ncells used:", good.groupby(["cell_class", "protocol"]).size().to_dict())
    points = pd.concat(points, ignore_index=True)
    points = points.merge(good[["cell_id", "rheobase_pa"]], on="cell_id")
    points["current_above_rheobase_pa"] = points["current_pa"] - points["rheobase_pa"]
    points.to_csv(config.PROCESSED_DIR / "fi_points.csv", index=False)

    config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    plot_group_curves(points)
    plot_features(good)

    stats = pd.concat([compare(good, "all"), compare(good[good["protocol"] == "newer"], "newer protocol only")])
    stats.to_csv(config.RESULTS_DIR / "stats_e_vs_i.csv", index=False)
    print()
    print(stats.to_string(index=False, float_format=lambda x: f"{x:.3g}"))


def compare(cells: pd.DataFrame, subset: str) -> pd.DataFrame:
    """Two-sided Mann-Whitney U test per feature, E vs I."""
    rows = []
    for feature in FEATURES:
        e = cells.loc[cells["cell_class"] == "E", feature].dropna()
        i = cells.loc[cells["cell_class"] == "I", feature].dropna()
        u, p = mannwhitneyu(e, i, alternative="two-sided")
        rows.append({
            "subset": subset, "feature": feature, "n_E": len(e), "n_I": len(i),
            "median_E": e.median(), "median_I": i.median(), "U": u, "p": p,
        })
    return pd.DataFrame(rows)


def plot_group_curves(points: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, x in zip(axes, ("current_pa", "current_above_rheobase_pa")):
        for (cell_class, _), cell in points.groupby(["cell_class", "cell_id"]):
            ax.plot(cell[x], cell["rate_hz"], color=CLASS_COLORS[cell_class], alpha=0.12, lw=0.8)
        curves = fi_curves.group_fi_curves(points, x=x)
        curves = curves[curves["n"] >= MIN_CELLS_PER_BIN]
        for cell_class, curve in curves.groupby("cell_class"):
            n = points.loc[points["cell_class"] == cell_class, "cell_id"].nunique()
            ax.errorbar(
                curve["bin"], curve["mean"], yerr=curve["sem"], color=CLASS_COLORS[cell_class],
                lw=2, capsize=2, label=f"{cell_class} (n = {n})",
            )
    axes[1].axvspan(0, config.SHARED_RANGE_ABOVE_RHEOBASE_PA, color="0.9", zorder=0, label="shared range")
    axes[0].set_xlabel("injected current (pA)")
    axes[1].set_xlabel("current above rheobase (pA)")
    axes[0].set_ylabel("mean firing rate over 1-s step (Hz)")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.set_xlim(left=-10 if ax is axes[1] else 0)
    fig.suptitle(f"Empirical f-I curves (mean ± SEM; bins with fewer than {MIN_CELLS_PER_BIN} cells hidden)")
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "empirical_fi.png", dpi=150)
    plt.close(fig)


def plot_features(cells: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(13, 6))
    rng = np.random.default_rng(0)
    for ax, (feature, label) in zip(axes.flat, FEATURES.items()):
        groups = [cells.loc[cells["cell_class"] == c, feature].dropna() for c in CLASS_COLORS]
        ax.boxplot(groups, tick_labels=list(CLASS_COLORS), widths=0.5, showfliers=False)
        for k, (values, color) in enumerate(zip(groups, CLASS_COLORS.values()), start=1):
            ax.scatter(k + rng.uniform(-0.15, 0.15, len(values)), values, s=8, color=color, alpha=0.6)
        ax.set_title(label, fontsize=9)
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "features_e_vs_i.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
