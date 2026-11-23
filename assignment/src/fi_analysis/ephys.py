"""Single-cell current-clamp analysis. Units: time s, voltage mV, current pA."""
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import curve_fit

from . import config


@dataclass
class Sweep:
    t: np.ndarray
    v: np.ndarray
    amplitude_pa: float
    stim_start: float
    stim_end: float
    stimulus: str = ""

    @property
    def duration(self) -> float:
        return self.stim_end - self.stim_start


def load_long_square_sweeps(nwb_path: Path) -> list[Sweep]:
    """All 1-s current-step sweeps of an Allen Patch-seq NWB file."""
    sweeps = []
    with h5py.File(nwb_path, "r") as f:
        stimuli = {int(s.attrs["sweep_number"]): s for s in f["stimulus/presentation"].values()}
        for series in f["acquisition"].values():
            if series.attrs.get("neurodata_type") != "CurrentClampSeries":
                continue
            description = str(series.attrs.get("stimulus_description", ""))
            if not description.startswith(config.LONG_SQUARE_STIMULI):
                continue
            stim = stimuli[int(series.attrs["sweep_number"])]
            rate = float(series["starting_time"].attrs["rate"])
            v = series["data"][:] * series["data"].attrs["conversion"] * 1e3
            i = stim["data"][:] * stim["data"].attrs["conversion"] * 1e12
            step = _longest_step(i)
            # some sweeps hold only the test pulse; aborted sweeps are padded with NaN or exact zeros
            if step is None or abs((step.stop - step.start) / rate - config.LONG_SQUARE_DURATION_S) > 0.05:
                continue
            if np.isnan(v[: step.stop]).any() or np.count_nonzero(v[: step.stop] == 0) > 100:
                continue
            sweeps.append(
                Sweep(
                    t=np.arange(len(v)) / rate,
                    v=v,
                    amplitude_pa=round(float(np.median(i[step])), 1),
                    stim_start=step.start / rate,
                    stim_end=step.stop / rate,
                    stimulus=description,
                )
            )
    return sweeps


def _longest_step(i: np.ndarray) -> slice | None:
    """Longest run of non-zero command current, which skips the short test pulse at sweep start."""
    on = np.flatnonzero(np.abs(i) > 1e-3)
    if len(on) == 0:
        return None
    runs = np.split(on, np.flatnonzero(np.diff(on) > 1) + 1)
    run = max(runs, key=len)
    return slice(int(run[0]), int(run[-1]) + 1)


def detect_spikes(
    t: np.ndarray,
    v: np.ndarray,
    dvdt_threshold: float = config.SPIKE_DVDT_THRESHOLD_MV_PER_MS,
    min_peak_mv: float = config.SPIKE_MIN_PEAK_MV,
    refractory_s: float = config.SPIKE_REFRACTORY_S,
) -> np.ndarray:
    """Return sample indices of spike onsets (upward crossings of the dV/dt threshold)."""
    dvdt = np.diff(v) / (np.diff(t) * 1e3)
    above = dvdt >= dvdt_threshold
    crossings = np.flatnonzero(~above[:-1] & above[1:]) + 1
    peak_window = max(1, int(round(refractory_s / (t[1] - t[0]))))
    onsets: list[int] = []
    for i in crossings:
        if onsets and t[i] - t[onsets[-1]] < refractory_s:
            continue
        if v[i : i + peak_window].max() >= min_peak_mv:
            onsets.append(i)
    return np.asarray(onsets, dtype=int)


def step_spikes(sweep: Sweep) -> np.ndarray:
    idx = detect_spikes(sweep.t, sweep.v)
    in_step = (sweep.t[idx] >= sweep.stim_start) & (sweep.t[idx] < sweep.stim_end)
    return idx[in_step]


def firing_rate_hz(sweep: Sweep) -> float:
    return len(step_spikes(sweep)) / sweep.duration


def baseline_mv(sweep: Sweep, window_s: float = config.BASELINE_WINDOW_S) -> float:
    mask = (sweep.t >= sweep.stim_start - window_s) & (sweep.t < sweep.stim_start)
    return float(sweep.v[mask].mean())


def baseline_noise_mv(sweep: Sweep, window_s: float = config.BASELINE_WINDOW_S) -> float:
    mask = (sweep.t >= sweep.stim_start - window_s) & (sweep.t < sweep.stim_start)
    return float(sweep.v[mask].std())


def qc_sweeps(
    sweeps: list[Sweep],
    max_deviation_mv: float = config.BASELINE_MAX_DEVIATION_MV,
    max_noise_mv: float = config.BASELINE_MAX_NOISE_MV,
) -> list[Sweep]:
    """Drop sweeps whose pre-step baseline is noisy or far from the cell's median baseline."""
    if not sweeps:
        return []
    base = np.array([baseline_mv(s) for s in sweeps])
    median = np.median(base)
    return [
        s for s, b in zip(sweeps, base)
        if abs(b - median) <= max_deviation_mv and baseline_noise_mv(s) <= max_noise_mv
    ]


def input_resistance_mohm(sweeps: list[Sweep], window_s: float = config.STEADY_STATE_WINDOW_S) -> float:
    """Slope through the origin of steady-state deflection vs current, over hyperpolarising steps."""
    currents, deflections = [], []
    for s in sweeps:
        if s.amplitude_pa >= 0:
            continue
        mask = (s.t >= s.stim_end - window_s) & (s.t < s.stim_end)
        currents.append(s.amplitude_pa)
        deflections.append(s.v[mask].mean() - baseline_mv(s))
    if not currents:
        return np.nan
    i, dv = np.asarray(currents), np.asarray(deflections)
    return float(dv @ i / (i @ i) * 1e3)  # mV/pA = GOhm


def membrane_tau_ms(sweep: Sweep) -> float:
    """Exponential fit from step onset to 90 % of the peak hyperpolarisation, as in the Allen pipeline."""
    v0 = baseline_mv(sweep)
    mask = (sweep.t >= sweep.stim_start) & (sweep.t < sweep.stim_end)
    t = sweep.t[mask] - sweep.stim_start
    v = sweep.v[mask]
    peak = int(np.argmin(v))
    end = int(np.argmax(v[: peak + 1] <= v0 + 0.9 * (v[peak] - v0)))
    t, v = t[: end + 1], v[: end + 1]
    if len(t) < 10:
        return np.nan

    def exp_decay(t, v_inf, a, tau):
        return v_inf + a * np.exp(-t / tau)

    try:
        (_, _, tau), _ = curve_fit(exp_decay, t, v, p0=(v[-1], v[0] - v[-1], 0.01), maxfev=10_000)
    except RuntimeError:
        return np.nan
    return float(tau * 1e3)


def spike_threshold_and_trough_mv(sweep: Sweep, trough_window_s: float = 0.02) -> tuple[float, float]:
    """Mean voltage at spike onset and mean minimum after each spike, within the step."""
    idx = step_spikes(sweep)
    if len(idx) == 0:
        return np.nan, np.nan
    window = int(round(trough_window_s / (sweep.t[1] - sweep.t[0])))
    following = np.append(idx[1:], len(sweep.v))
    troughs = [sweep.v[i : min(j, i + window)].min() for i, j in zip(idx, following)]
    return float(sweep.v[idx].mean()), float(np.mean(troughs))
