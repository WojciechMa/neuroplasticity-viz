"""Step 3: vary g_L in a LIF model fitted to a representative E and I cell (a model perturbation only)."""
from dataclasses import asdict, replace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from fi_analysis import config, fi_curves
from fi_analysis.model import fit_lif, lif_rate_hz, passive_parameters

LEAK_SCALES = (0.5, 0.75, 1.0, 1.5, 2.0)
CURRENTS_PA = np.arange(0.0, 400.5, 1.0)


def main() -> None:
    summaries = pd.read_csv(config.PROCESSED_DIR / "cell_summaries.csv")
    points = pd.read_csv(config.PROCESSED_DIR / "fi_points.csv")
    reps = fi_curves.representative_cells(summaries)

    fits, scan = [], []
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(LEAK_SCALES)))
    for ax, rep in zip(axes, reps.itertuples()):
        recorded = points[points["cell_id"] == rep.cell_id]
        g_L, C = passive_parameters(rep.input_resistance_mohm, rep.tau_ms)
        base = fit_lif(recorded["current_pa"], recorded["rate_hz"], C=C, g_L=g_L, E_L=rep.v_rest_mv)
        rmse = np.sqrt(np.mean((lif_rate_hz(recorded["current_pa"], base) - recorded["rate_hz"]) ** 2))
        fits.append({"cell_class": rep.cell_class, "cell_id": rep.cell_id, "t_type": rep.t_type, **asdict(base), "rmse_hz": rmse})

        for scale, color in zip(LEAK_SCALES, colors):
            rates = lif_rate_hz(CURRENTS_PA, replace(base, g_L=base.g_L * scale))
            scan.append({"cell_class": rep.cell_class, "scale": scale, **fi_curves.curve_metrics(CURRENTS_PA, rates)})
            ax.plot(CURRENTS_PA, rates, color=color, lw=2.5 if scale == 1 else 1.2, label=f"{scale:g} × g_L")
        ax.plot(recorded["current_pa"], recorded["rate_hz"], "ko", ms=5, label="recorded")
        ax.set_title(f"{rep.cell_class}: {rep.t_type} (cell {rep.cell_id})", fontsize=9)
        ax.set_xlabel("injected current (pA)")
        ax.set_ylabel("firing rate (Hz)")
        ax.set_ylim(0, max(40, recorded["rate_hz"].max() * 1.5))
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("LIF model: g_L scaled; C, E_L, V_th, V_reset and t_ref fixed")
    fig.tight_layout()
    config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(config.FIGURES_DIR / "lif_leak_perturbation.png", dpi=150)

    fits, scan = pd.DataFrame(fits), pd.DataFrame(scan)
    fits.to_csv(config.RESULTS_DIR / "lif_fits.csv", index=False)
    scan.to_csv(config.RESULTS_DIR / "lif_leak_scan.csv", index=False)
    print(fits.round(2).to_string(index=False))
    print(scan.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
