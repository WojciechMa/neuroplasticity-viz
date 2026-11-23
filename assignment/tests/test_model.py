from dataclasses import replace

import numpy as np

from fi_analysis.model import HHParams, LIFParams, fit_lif, hh_rate_hz, lif_rate_hz, simulate_lif

P = LIFParams(C=200.0, g_L=10.0, E_L=-70.0, V_th=-50.0, V_reset=-65.0, t_ref=2.0)
HH = HHParams.from_template("RS", C=100.0, g_L=5.0, E_L=-70.0)


def test_silent_up_to_rheobase():
    assert P.rheobase == 200.0
    assert np.all(lif_rate_hz([0.0, 100.0, 200.0], P) == 0.0)


def test_analytic_rate_matches_simulation():
    for current in (250.0, 400.0, 800.0):
        spikes = simulate_lif(current, P, duration_ms=1000.0)
        simulated = 1e3 / np.mean(np.diff(spikes))
        assert abs(simulated - lif_rate_hz(current, P)[0]) / simulated < 0.01


def test_more_leak_raises_rheobase():
    assert replace(P, g_L=2 * P.g_L).rheobase == 2 * P.rheobase


def test_fit_lif_recovers_rheobase_and_rates():
    current = np.arange(150.0, 450.0, 20.0)
    fitted = fit_lif(current, lif_rate_hz(current, P), C=P.C, g_L=P.g_L, E_L=P.E_L)
    assert 190.0 <= fitted.rheobase <= 210.0
    assert np.allclose(lif_rate_hz(current, fitted), lif_rate_hz(current, P), atol=1.0)


def test_hh_silent_at_rest_and_fires_when_driven():
    rates = hh_rate_hz(np.array([0.0, 400.0]), HH, duration_ms=300.0)
    assert rates[0] == 0.0
    assert rates[1] > 10.0


def test_hh_more_leak_raises_rheobase():
    current = np.arange(0.0, 400.0, 10.0)
    rates = hh_rate_hz(current[None, :], replace(HH, g_L=HH.g_L * np.array([[1.0], [2.0]])), duration_ms=300.0)
    rheobase = [current[r > 0].min() for r in rates]
    assert rheobase[1] > rheobase[0]
