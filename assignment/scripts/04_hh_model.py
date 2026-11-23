"""Step 4: the same perturbation in a Hodgkin-Huxley-type model (Na, K and M channels), plus g_Na and g_K."""
from dataclasses import asdict, replace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fi_analysis import config, fi_curves
from fi_analysis.model import HHParams, fit_hh, hh_rate_hz, passive_parameters

SCALES = (0.5, 0.75, 1.0, 1.5, 2.0)
CONDUCTANCES = {"g_L": "leak", "g_Na": "Na", "g_K": "delayed-rectifier K"}
TEMPLATE = {"E": "RS", "I": "FS"}
V_T_GRID = np.arange(-75.0, -34.5, 1.0)
G_M_DENSITY_GRID = np.array([0.0, 2.5e-5, 5e-5, 1e-4, 2e-4, 4e-4])  # S/cm2; the RS template has 7e-5
CURRENTS_PA = np.arange(0.0, 400.5, 5.0)
CLASS_COLORS = {"E": "tab:red", "I": "tab:blue"}


def main() -> None:
    summaries = pd.read_csv(config.PROCESSED_DIR / "cell_summaries.csv")
    points = pd.read_csv(config.PROCESSED_DIR / "fi_points.csv")
    reps = fi_curves.representative_cells(summaries)

    fits, scan = [], []
    fig, axes = plt.subplots(2, len(CONDUCTANCES), figsize=(14, 8), sharey="row")
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(SCALES)))
    for row, rep in zip(axes, reps.itertuples()):
        recorded = points[points["cell_id"] == rep.cell_id]
        g_L, C = passive_parameters(rep.input_resistance_mohm, rep.tau_ms)
        base = HHParams.from_template(TEMPLATE[rep.cell_class], C, g_L, rep.v_rest_mv)
        fitted = fit_hh(recorded["current_pa"], recorded["rate_hz"], base, V_T_GRID, G_M_DENSITY_GRID * C * 1e3)
        rmse = np.sqrt(np.mean((hh_rate_hz(recorded["current_pa"].to_numpy(), fitted) - recorded["rate_hz"]) ** 2))
        fits.append({
            "cell_class": rep.cell_class, "cell_id": rep.cell_id, "t_type": rep.t_type,
            "template": TEMPLATE[rep.cell_class], **asdict(fitted), "rmse_hz": rmse,
        })
        print(f"{rep.cell_class}: V_T = {fitted.V_T:g} mV, g_M = {fitted.g_M:.1f} nS, RMSE = {rmse:.1f} Hz", flush=True)

        for ax, (conductance, name) in zip(row, CONDUCTANCES.items()):
            scales = np.asarray(SCALES)[:, None]
            rates = hh_rate_hz(CURRENTS_PA[None, :], replace(fitted, **{conductance: getattr(fitted, conductance) * scales}))
            for scale, r, color in zip(SCALES, rates, colors):
                scan.append({
                    "cell_class": rep.cell_class, "conductance": conductance, "scale": scale,
                    **fi_curves.curve_metrics(CURRENTS_PA, r),
                })
                ax.plot(CURRENTS_PA, r, color=color, lw=2.5 if scale == 1 else 1.2, label=f"{scale:g} × {conductance}")
            ax.plot(recorded["current_pa"], recorded["rate_hz"], "ko", ms=4, label="recorded")
            ax.set_title(f"{rep.cell_class} ({rep.t_type}): {name} scaled", fontsize=9)
            ax.set_xlabel("injected current (pA)")
            ax.legend(frameon=False, fontsize=7)
        row[0].set_ylabel("firing rate over 1-s step (Hz)")
    fig.suptitle("Hodgkin-Huxley-type model (Pospischil et al. 2008 templates, passive properties from the cell)")
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "hh_perturbation.png", dpi=150)
    plt.close(fig)

    fits, scan = pd.DataFrame(fits), pd.DataFrame(scan)
    fits.to_csv(config.RESULTS_DIR / "hh_fits.csv", index=False)
    scan.to_csv(config.RESULTS_DIR / "hh_scan.csv", index=False)
    print(scan.round(3).to_string(index=False))

    lif = pd.read_csv(config.RESULTS_DIR / "lif_leak_scan.csv")
    plot_leak_comparison({"LIF": lif, "HH": scan[scan["conductance"] == "g_L"]})


def plot_leak_comparison(scans: dict[str, pd.DataFrame]) -> None:
    """Rheobase and shared-range slope relative to the unperturbed model, against the g_L scale."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    styles = {"LIF": "--s", "HH": "-o"}
    titles = {
        "rheobase_pa": "rheobase",
        "fi_slope_shared_hz_per_pa": f"f-I slope, rheobase to +{config.SHARED_RANGE_ABOVE_RHEOBASE_PA:g} pA",
    }
    for ax, (metric, title) in zip(axes, titles.items()):
        for model, scan in scans.items():
            for cell_class, group in scan.groupby("cell_class"):
                ref = group.loc[group["scale"] == 1.0, metric].iloc[0]
                ax.plot(group["scale"], group[metric] / ref, styles[model], color=CLASS_COLORS[cell_class],
                        label=f"{model}, {cell_class}")
        ax.axhline(1, color="0.7", lw=0.8)
        ax.set_xscale("log", base=2)
        ax.set_xticks(SCALES, [f"{s:g}" for s in SCALES])
        ax.set_xlabel("g_L scale")
        ax.set_ylabel("relative to 1 × g_L")
        ax.set_title(title, fontsize=10)
        ax.legend(frameon=False, fontsize=8)
    fig.suptitle("Effect of the leak conductance in both models")
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "leak_effect_lif_vs_hh.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
