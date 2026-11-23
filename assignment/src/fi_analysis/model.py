"""Neuron models: leaky integrate-and-fire and a single-compartment Hodgkin-Huxley type.

Units: pA, nS, pF, mV, ms.
"""
from dataclasses import dataclass, replace

import numpy as np
from scipy.optimize import least_squares


def passive_parameters(input_resistance_mohm: float, tau_ms: float) -> tuple[float, float]:
    """Leak conductance (nS) and capacitance (pF) from input resistance and membrane time constant."""
    g_L = 1e3 / input_resistance_mohm
    return g_L, tau_ms * g_L


@dataclass(frozen=True)
class LIFParams:
    C: float
    g_L: float
    E_L: float
    V_th: float
    V_reset: float
    t_ref: float

    @property
    def tau(self) -> float:
        return self.C / self.g_L

    @property
    def rheobase(self) -> float:
        return self.g_L * (self.V_th - self.E_L)


def lif_rate_hz(current_pa, p: LIFParams) -> np.ndarray:
    """Analytic steady firing rate under constant current."""
    i = np.atleast_1d(np.asarray(current_pa, dtype=float))
    rate = np.zeros_like(i)
    supra = i > p.rheobase
    isi_ms = p.t_ref + p.tau * np.log((i[supra] - p.g_L * (p.V_reset - p.E_L)) / (i[supra] - p.rheobase))
    rate[supra] = 1e3 / isi_ms
    return rate


def fit_lif(current_pa, rate_hz, C: float, g_L: float, E_L: float) -> LIFParams:
    """LIF matched to recorded f-I points; C, g_L and E_L stay as measured.

    V_th puts the model rheobase midway between the last silent and the first spiking step;
    V_reset (at most 20 mV below rest) and t_ref are then fitted by least squares. Without
    adaptation, the fitted t_ref acts as an effective refractory period rather than a measured one.
    """
    current_pa = np.asarray(current_pa, dtype=float)
    rate_hz = np.asarray(rate_hz, dtype=float)
    rheobase = current_pa[rate_hz > 0].min()
    silent = current_pa[(rate_hz == 0) & (current_pa < rheobase)]
    V_th = E_L + (0.5 * (silent.max() + rheobase) if len(silent) else rheobase) / g_L

    def params(x) -> LIFParams:
        reset_below_th, t_ref = x
        return LIFParams(C=C, g_L=g_L, E_L=E_L, V_th=V_th, V_reset=V_th - reset_below_th, t_ref=t_ref)

    fit = least_squares(
        lambda x: lif_rate_hz(current_pa, params(x)) - rate_hz,
        x0=(V_th - E_L, 2.0),
        bounds=([0.5, 0.0], [V_th - E_L + 20.0, 200.0]),
    )
    return params(fit.x)


def simulate_lif(current_pa: float, p: LIFParams, duration_ms: float, dt_ms: float = 0.01) -> np.ndarray:
    """Forward-Euler simulation under constant current; returns spike times in ms."""
    ref_steps = int(round(p.t_ref / dt_ms))
    v = p.E_L
    hold = 0
    spikes = []
    for k in range(int(duration_ms / dt_ms)):
        if hold > 0:
            hold -= 1
            continue
        v += dt_ms * (-p.g_L * (v - p.E_L) + current_pa) / p.C
        if v >= p.V_th:
            spikes.append((k + 1) * dt_ms)
            v = p.V_reset
            hold = ref_steps
    return np.asarray(spikes)


# Channel densities (S/cm2) of the regular-spiking and fast-spiking cells of Pospischil et al. 2008,
# taken from their ModelDB code (accession 123623).
POSPISCHIL_DENSITIES = {
    "RS": {"g_Na": 0.05, "g_K": 0.005, "g_M": 7e-5},
    "FS": {"g_Na": 0.05, "g_K": 0.01, "g_M": 0.0},
}


@dataclass(frozen=True)
class HHParams:
    """Na, delayed-rectifier K and slow M-type K channels plus leak (Pospischil et al. 2008).

    Any field may be a numpy array; hh_rate_hz broadcasts them against the currents.
    """

    C: float
    g_L: float
    E_L: float
    g_Na: float
    g_K: float
    g_M: float
    V_T: float = -55.0
    tau_max: float = 1000.0
    E_Na: float = 50.0
    E_K: float = -100.0

    @classmethod
    def from_template(cls, template: str, C: float, g_L: float, E_L: float) -> "HHParams":
        """Template channel densities scaled to the cell's capacitance, assuming 1 uF/cm2."""
        nS_per_S_cm2 = C * 1e3
        densities = POSPISCHIL_DENSITIES[template]
        return cls(C=C, g_L=g_L, E_L=E_L, **{k: d * nS_per_S_cm2 for k, d in densities.items()})


def _vtrap(x, y):
    """x / (exp(x/y) - 1) without the 0/0 at x = 0."""
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = x / np.expm1(x / y)
    return np.where(np.abs(x / y) < 1e-6, y * (1 - x / y / 2), ratio)


def _gate_rates(v, V_T):
    """Opening and closing rates (1/ms) of m, h and n (Traub & Miles kinetics, as in Pospischil et al.)."""
    u = v - V_T
    return (
        (0.32 * _vtrap(13 - u, 4), 0.28 * _vtrap(u - 40, 5)),
        (0.128 * np.exp((17 - u) / 18), 4 / (1 + np.exp((40 - u) / 5))),
        (0.032 * _vtrap(15 - u, 5), 0.5 * np.exp((10 - u) / 40)),
    )


def _m_current_gate(v, tau_max):
    x = (v + 35) / 20
    return 1 / (1 + np.exp(-(v + 35) / 10)), tau_max / (3.3 * np.exp(x) + np.exp(-x))


def hh_rate_hz(
    current_pa,
    p: HHParams,
    duration_ms: float = 1000.0,
    settle_ms: float = 200.0,
    dt_ms: float = 0.05,  # spike counts match dt = 0.005 ms
    spike_mv: float = 0.0,
) -> np.ndarray:
    """Mean firing rate during a current step that follows a settling period at zero current."""
    current = np.asarray(current_pa, dtype=float)
    fields = {k: np.asarray(v, dtype=float) for k, v in vars(p).items()}
    shape = np.broadcast_shapes(current.shape, *(f.shape for f in fields.values()))
    f = {k: np.broadcast_to(v, shape) for k, v in fields.items()}

    v = f["E_L"].copy()
    (am, bm), (ah, bh), (an, bn) = _gate_rates(v, f["V_T"])
    m, h, n = am / (am + bm), ah / (ah + bh), an / (an + bn)
    w = _m_current_gate(v, f["tau_max"])[0]

    n_settle = int(round(settle_ms / dt_ms))
    count = np.zeros(shape)
    above = np.zeros(shape, dtype=bool)
    for k in range(n_settle + int(round(duration_ms / dt_ms))):
        i_ext = current if k >= n_settle else 0.0
        (am, bm), (ah, bh), (an, bn) = _gate_rates(v, f["V_T"])
        m = am / (am + bm) + (m - am / (am + bm)) * np.exp(-dt_ms * (am + bm))
        h = ah / (ah + bh) + (h - ah / (ah + bh)) * np.exp(-dt_ms * (ah + bh))
        n = an / (an + bn) + (n - an / (an + bn)) * np.exp(-dt_ms * (an + bn))
        w_inf, tau_w = _m_current_gate(v, f["tau_max"])
        w = w_inf + (w - w_inf) * np.exp(-dt_ms / tau_w)

        g_na = f["g_Na"] * m**3 * h
        g_k = f["g_K"] * n**4 + f["g_M"] * w
        g_total = f["g_L"] + g_na + g_k
        v_inf = (f["g_L"] * f["E_L"] + g_na * f["E_Na"] + g_k * f["E_K"] + i_ext) / g_total
        v = v_inf + (v - v_inf) * np.exp(-dt_ms * g_total / f["C"])

        now_above = v > spike_mv
        if k >= n_settle:
            count += now_above & ~above
        above = now_above
    return count / (duration_ms / 1e3)


def fit_hh(current_pa, rate_hz, base: HHParams, v_t_grid, g_m_grid) -> HHParams:
    """Grid search over V_T and g_M for the smallest squared error to the recorded f-I points."""
    v_t_grid, g_m_grid = np.asarray(v_t_grid, dtype=float), np.asarray(g_m_grid, dtype=float)
    grid = replace(base, V_T=v_t_grid[:, None, None], g_M=g_m_grid[None, :, None])
    rates = hh_rate_hz(np.asarray(current_pa, dtype=float)[None, None, :], grid)
    error = ((rates - np.asarray(rate_hz, dtype=float)) ** 2).sum(axis=-1)
    k, j = np.unravel_index(np.argmin(error), error.shape)
    return replace(base, V_T=float(v_t_grid[k]), g_M=float(g_m_grid[j]))
