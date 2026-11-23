import numpy as np
import pandas as pd

from fi_analysis.ephys import Sweep, _longest_step, firing_rate_hz, qc_sweeps
from fi_analysis.fi_curves import curve_metrics, fi_slope_hz_per_pa, rheobase_pa

RATE = 50_000.0


def make_sweep(baseline_mv=-70.0, spike_times_s=(), noise_mv=0.0, seed=0) -> Sweep:
    """1.5-s sweep with a step from 0.25 to 1.25 s and 0.5-ms triangular spikes to +30 mV."""
    t = np.arange(int(1.5 * RATE)) / RATE
    v = np.full_like(t, baseline_mv) + np.random.default_rng(seed).normal(0, noise_mv, t.size)
    for ts in spike_times_s:
        rise = (t >= ts) & (t < ts + 5e-4)
        v[rise] = baseline_mv + (t[rise] - ts) / 5e-4 * (30.0 - baseline_mv)
    return Sweep(t=t, v=v, amplitude_pa=100.0, stim_start=0.25, stim_end=1.25)


def test_longest_step_skips_test_pulse():
    i = np.zeros(1000)
    i[10:20] = -10.0  # test pulse
    i[300:800] = 50.0  # step
    assert _longest_step(i) == slice(300, 800)
    assert _longest_step(np.zeros(10)) is None


def test_firing_rate_counts_only_spikes_in_step():
    sweep = make_sweep(spike_times_s=(0.1, 0.3, 0.5, 0.7, 1.4))
    assert firing_rate_hz(sweep) == 3.0


def test_qc_drops_shifted_and_noisy_sweeps():
    good = [make_sweep(seed=k) for k in range(3)]
    shifted = make_sweep(baseline_mv=-55.0)
    noisy = make_sweep(noise_mv=3.0)
    assert qc_sweeps([*good, shifted, noisy]) == good


def test_shared_range_slope_ignores_currents_beyond_limit():
    fi = pd.DataFrame({"current_pa": [50, 100, 120, 140, 160, 300], "rate_hz": [0, 2, 4, 6, 8, 100]})
    assert rheobase_pa(fi) == 100
    assert np.isclose(fi_slope_hz_per_pa(fi, max_above_rheobase_pa=60), 0.1)
    assert fi_slope_hz_per_pa(fi) > 0.1


def test_curve_metrics_on_threshold_linear_curve():
    current = np.arange(0.0, 300.0, 1.0)
    rate = np.clip(0.2 * (current - 100.0), 0, None)
    metrics = curve_metrics(current, rate)
    assert metrics["rheobase_pa"] == 101.0
    assert np.isclose(metrics["fi_slope_shared_hz_per_pa"], 0.2)
