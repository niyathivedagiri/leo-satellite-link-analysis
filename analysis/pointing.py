"""Pointing-error scenarios and link-quality metrics (used by Phases 3 and 4).

A scenario turns one pass into a time series of pointing error epsilon(t) [deg].
The link model then subtracts the pointing loss from the perfectly pointed link (Phase 2).
"""
import numpy as np
import pandas as pd

import ntn_model as M

CN_THRESHOLD_DB = -2.0   # [A] link counted as usable at C/N >= -2 dB: the simulation's own "connect" threshold,
                         #     roughly where the most robust 5G modulation and coding still work


# ---------------------------------------------------------------------------
# Error scenarios: each returns epsilon(t) in degrees for the samples of one pass
# ---------------------------------------------------------------------------
def constant_offset(out, offset_deg):
    """The dish is always offset by the same angle (e.g. a misaligned mount or a bent feed)."""
    return np.full(len(out), float(offset_deg))


def stale_orbit_data(out, s, delay_s):
    """The terminal points where the satellite WAS delay_s seconds ago (out-of-date orbit data or a slow loop).
    Error = angle between the satellite's direction at t and at t - delay_s. It changes during the pass:
    largest near the top of the pass, where the satellite moves fastest across the sky."""
    t = out["time_s"].to_numpy()
    az1, el1, _ = M.look_angles(t, s)
    az0, el0, _ = M.look_angles(t - delay_s, s)
    return M.angle_between_deg(az1, el1, az0, el0)


def random_jitter(out, sigma_deg, seed=0):
    """Random wobble (servo noise, wind buffeting): independent Gaussian errors on both axes, std sigma each.
    The total error then follows a Rayleigh distribution with mean sigma*sqrt(pi/2)."""
    rng = np.random.default_rng(seed)
    x, y = rng.normal(0, sigma_deg, len(out)), rng.normal(0, sigma_deg, len(out))
    return np.hypot(x, y)


# ---------------------------------------------------------------------------
# Apply an error to the ideal link and measure the result
# ---------------------------------------------------------------------------
def apply_error(out, eps_deg, beamwidth_deg):
    """Add pointing loss to the Phase 2 (perfect pointing) link."""
    r = out[["time_s", "t_rel_s", "el_deg", "range_km", "doppler_khz", "prx_dbm", "cn_db", "noise_dbm"]].copy()
    r["pointing_error_deg"] = eps_deg
    r["pointing_loss_db"] = M.pointing_loss_db(eps_deg, beamwidth_deg)
    r["prx_dbm"] = out["prx_dbm"].to_numpy() - r["pointing_loss_db"].to_numpy()
    r["cn_db"] = r["prx_dbm"] - r["noise_dbm"]
    return r


def metrics(r, s, bandwidth_hz, cn_threshold_db=CN_THRESHOLD_DB):
    """Link-quality numbers over the part of the pass above the elevation mask."""
    t = r["time_s"].to_numpy()
    dt = np.gradient(t)
    vis = r["el_deg"].to_numpy() >= s["minElevation"]
    cn = r["cn_db"].to_numpy()
    up = vis & (cn >= cn_threshold_db)
    se = M.spectral_efficiency(cn)
    return {
        "mean_error_deg": r.loc[vis, "pointing_error_deg"].mean(),
        "max_error_deg": r.loc[vis, "pointing_error_deg"].max(),
        "mean_pointing_loss_db": r.loc[vis, "pointing_loss_db"].mean(),
        "max_pointing_loss_db": r.loc[vis, "pointing_loss_db"].max(),
        "min_cn_db": cn[vis].min(),
        "mean_cn_db": cn[vis].mean(),
        "visible_s": dt[vis].sum(),
        "usable_s": dt[up].sum(),
        "availability_pct": 100 * dt[up].sum() / dt[vis].sum(),
        "data_gb": (se[up] * bandwidth_hz * dt[up]).sum() / 8e9,
    }


def run_scenarios(out, s, p, scenarios):
    """scenarios: list of (label, eps array). Returns a table with data relative to perfect pointing."""
    rows = []
    for label, eps in scenarios:
        m = metrics(apply_error(out, eps, s["beamwidthDeg"]), s, p.bandwidth_mhz * 1e6)
        rows.append({"scenario": label, **m})
    df = pd.DataFrame(rows)
    ideal = df.loc[0, "data_gb"]
    df["data_vs_perfect_pct"] = 100 * df["data_gb"] / ideal
    return df


def allowed_error_deg(margin_db, beamwidth_deg):
    """Largest pointing error that costs no more than margin_db: invert L = 12.04 (e/theta3)^2."""
    return beamwidth_deg * np.sqrt(np.maximum(margin_db, 0) / M.GAUSS_COEFF_DB)
