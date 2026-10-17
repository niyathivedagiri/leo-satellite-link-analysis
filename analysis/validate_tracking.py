"""PHASE 1 - Validate satellite tracking.

Checks the CoppeliaSim run against independent calculations and prints PASS / FAIL:

  V1  orbit period and speed (Kepler's third law)                 - information
  V2  az / el / range in the log  =  an independent Python orbit   - the orbit + geometry code
  V3  logged range  =  closed-form slant range from elevation      - the geometry, a second way
  V4  range rate crosses zero exactly when elevation peaks         - the range-rate (Doppler) input
  V5  the 3D dish in the scene points where the maths says         - joint signs and axes
  V6  pointing error is ~0 whenever the motors are fast enough     - the tracking controller
  V7  the joints never move faster than the motor limit            - the actuator model
  V8  the run used zero bias and no step-track                     - the run is a clean baseline

Usage (from the project folder):
    python analysis/validate_tracking.py                       # reads ~/ntn_log.csv
    python analysis/validate_tracking.py --log path/to/ntn_log.csv --pass 2
Outputs go to results/phase1/.
"""
import argparse
import os

import numpy as np
import pandas as pd

import ntn_model as M
import plot_style as S
from runlog import ON_SKY, DEFAULT_LOG, load_log, pick_pass, rate_deg_s, time_since_rise

plt = S.plt


def run_checks(df, s):
    rows = []

    def check(cid, name, value, limit, ok, note=""):
        rows.append({"id": cid, "check": name, "value": value, "limit": limit,
                     "result": "PASS" if ok else "FAIL", "note": note})

    # V1 - orbit period and speed (no log needed; it is what the log must be consistent with)
    T, v = M.orbit_period_s(s), M.orbit_speed_kms(s)
    rows.append({"id": "V1", "check": "Orbit period / speed", "value": f"{T:.1f} s ({T / 60:.2f} min) / {v:.3f} km/s",
                 "limit": "-", "result": "INFO", "note": f"{s['altitudeKm']} km circular orbit"})

    # V2 - independent re-computation of the satellite geometry at every logged time
    az_py, el_py, r_py = M.look_angles(df["time_s"].to_numpy(), s)
    d_el = np.abs(el_py - df["sat_el_deg"]).max()
    vis = df["sat_el_deg"].to_numpy() > 1
    d_az = np.abs(M.wrap180(az_py[vis] - df["sat_az_deg"].to_numpy()[vis])).max() if vis.any() else 0
    d_r = np.abs(r_py - df["range_km"]).max()
    check("V2a", "Elevation: Lua vs Python (max diff)", f"{d_el:.2e} deg", "< 1e-3 deg", d_el < 1e-3)
    check("V2b", "Azimuth: Lua vs Python (max diff, el > 1 deg)", f"{d_az:.2e} deg", "< 1e-3 deg", d_az < 1e-3)
    check("V2c", "Range: Lua vs Python (max diff)", f"{d_r:.2e} km", "< 0.01 km", d_r < 0.01)

    # V3 - closed-form slant range from elevation (independent of the orbit code)
    up = df["sat_el_deg"] > 0
    d_cf = np.abs(M.slant_range_from_elevation_km(df.loc[up, "sat_el_deg"], s) - df.loc[up, "range_km"]).max()
    check("V3", "Range vs closed form R(el) (max diff)", f"{d_cf:.2e} km", "< 0.01 km", d_cf < 0.01)

    # V4 - range rate is zero at maximum elevation, in every complete pass
    worst = 0.0
    for _, g in df[df["sat_el_deg"] > 0].groupby("pass_id"):
        if len(g) < 10:
            continue
        t = g["time_s"].to_numpy()
        t_top = t[np.argmax(g["sat_el_deg"].to_numpy())]
        rr = g["range_rate_kms"].to_numpy()
        k = np.where(np.diff(np.sign(rr)) > 0)[0]          # - to + : closest approach
        if len(k):
            i = k[0]
            t_zero = t[i] - rr[i] * (t[i + 1] - t[i]) / (rr[i + 1] - rr[i])
            worst = max(worst, abs(t_zero - t_top))
    step = np.median(np.diff(df["time_s"]))
    check("V4", "Range rate = 0 at max elevation (time gap)", f"{worst:.3f} s", f"< {2 * step:.2f} s (2 steps)",
          worst < 2 * step)

    # V5 - the 3D scene: angle between the dish marker's axis and the computed dish direction
    sc = df["scene_check_deg"].max()
    check("V5", "3D dish direction vs maths (max angle)", f"{sc:.2e} deg", "< 0.01 deg", sc < 0.01,
          "fails if a joint turns the wrong way (panSign / tiltSign)")

    # V6 - pointing error while tracking above the mask, when the motors can keep up
    trk = df["state"].isin(ON_SKY) & (df["sat_el_deg"] >= s["minElevation"])
    limited = np.zeros(len(df), dtype=bool)
    for _, g in df.groupby("pass_id"):
        t = g["time_s"].to_numpy()
        need = np.maximum(np.abs(rate_deg_s(g["sat_az_deg"], t, wrap=True)), np.abs(rate_deg_s(g["sat_el_deg"], t)))
        lim = need > s["slewDegPerSec"]
        lim = np.convolve(lim, np.ones(41), mode="same") > 0     # +/- 20 samples to settle
        limited[df.index.get_indexer(g.index)] = lim
    ok_rows = trk.to_numpy() & ~limited
    err = df.loc[ok_rows, "pointing_error_deg"]
    frac = (err < 0.01).mean() if len(err) else np.nan
    check("V6", "Pointing error < 0.01 deg while tracking (motors not saturated)",
          f"{100 * frac:.2f} % of {len(err)} samples", ">= 99 %", frac >= 0.99,
          f"median {err.median():.1e} deg, max {err.max():.1e} deg" if len(err) else "no tracking samples")
    n_lim = int((trk.to_numpy() & limited).sum())
    if n_lim:
        e_lim = df.loc[trk.to_numpy() & limited, "pointing_error_deg"].max()
        rows.append({"id": "V6b", "check": "Samples where the satellite moves faster than the motors",
                     "value": f"{n_lim} samples, max error {e_lim:.2f} deg", "limit": "-", "result": "INFO",
                     "note": "expected near zenith (keyhole); real tracking error, not a bug"})

    # V7 - joint speed never above the motor limit
    worst_rate = 0.0
    for _, g in df.groupby("pass_id"):
        t = g["time_s"].to_numpy()
        dt = np.diff(t)
        good = dt > 0
        r_az = np.abs(M.wrap180(np.diff(g["dish_az_deg"].to_numpy())))[good] / dt[good]
        r_el = np.abs(np.diff(g["dish_el_deg"].to_numpy()))[good] / dt[good]
        worst_rate = max(worst_rate, r_az.max(initial=0), r_el.max(initial=0))
    check("V7", "Fastest joint motion in the log", f"{worst_rate:.3f} deg/s",
          f"<= {s['slewDegPerSec']} deg/s", worst_rate <= s["slewDegPerSec"] * 1.001)

    # V8 - clean baseline
    clean = s["biasAzDeg"] == 0 and s["biasElDeg"] == 0 and not s["linkAssist"]
    check("V8", "Run used zero bias, no step-track", f"bias {s['biasAzDeg']}/{s['biasElDeg']} deg, "
          f"assist {s['linkAssist']}", "0 / 0, false", clean, "wind button must not be pressed either")
    return pd.DataFrame(rows)


def plot_pass(g, s, out_dir):
    t = time_since_rise(g)
    pid = int(g["pass_id"].iloc[0])

    # 1. sky plot: the satellite's path and the dish's path, seen from the station
    fig = plt.figure(figsize=(5.2, 5.2))
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)                       # compass: clockwise from North
    up = g["sat_el_deg"] > 0
    ax.plot(np.radians(g.loc[up, "sat_az_deg"]), 90 - g.loc[up, "sat_el_deg"], color=S.SAT, label="Satellite")
    on = g["state"].isin(ON_SKY)
    ax.plot(np.radians(g.loc[on, "dish_az_deg"]), 90 - g.loc[on, "dish_el_deg"], color=S.DISH, lw=1.2,
            ls=(0, (4, 3)), label="Dish (while tracking)")
    ring = np.linspace(0, 2 * np.pi, 200)
    ax.plot(ring, np.full_like(ring, 90 - s["minElevation"]), color=S.INK2, lw=0.8, ls=":")
    ax.set_rlim(0, 90)
    ax.set_yticks([0, 30, 60, 80])
    ax.set_yticklabels(["90", "60", "30", "10"], color=S.INK2, fontsize=8)
    ax.set_title(f"Pass {pid}: sky track (max elevation {g['sat_el_deg'].max():.1f} deg)", pad=14)
    ax.legend(loc="lower left", bbox_to_anchor=(-0.08, -0.12), fontsize=8)
    fig.savefig(os.path.join(out_dir, f"phase1_sky_pass{pid}.png"))
    plt.close(fig)

    # 2. az, el and pointing error over time
    fig, axs = plt.subplots(3, 1, figsize=(8, 7.5), sharex=True)
    mask_lo = g["sat_el_deg"].to_numpy() < s["minElevation"]
    for a in axs:
        a.fill_between(t, 0, 1, where=mask_lo, transform=a.get_xaxis_transform(), color=S.BAND, lw=0)
    sat_az = np.degrees(np.unwrap(np.radians(g["sat_az_deg"])))           # continuous, no 360 -> 0 jump
    dish_az = sat_az + M.wrap180(g["dish_az_deg"].to_numpy() - g["sat_az_deg"].to_numpy())   # same branch
    axs[0].plot(t, sat_az, color=S.SAT, label="Satellite")
    axs[0].plot(t, dish_az, color=S.DISH, lw=1.2, ls=(0, (4, 3)), label="Dish")
    axs[0].set_ylabel("Azimuth (deg)")
    axs[0].legend(loc="best")
    axs[0].set_title(f"Pass {pid}: the dish follows the satellite (shaded = below {s['minElevation']:g} deg)")
    axs[1].plot(t, g["sat_el_deg"], color=S.SAT, label="Satellite")
    axs[1].plot(t, g["dish_el_deg"], color=S.DISH, lw=1.2, ls=(0, (4, 3)), label="Dish")
    axs[1].set_ylabel("Elevation (deg)")
    err = np.clip(g["pointing_error_deg"].to_numpy(), 1e-6, None)
    axs[2].semilogy(t, err, color=S.DISH, lw=1.2, label="Pointing error")
    axs[2].axhline(s["beamwidthDeg"] / 2, color=S.INK2, lw=0.8, ls=":")
    axs[2].text(t[-1], s["beamwidthDeg"] / 2, " half beamwidth (3 dB loss)", va="bottom", ha="right",
                fontsize=8, color=S.INK2)
    axs[2].set_ylabel("Pointing error (deg)\n(floor 1e-6)")
    axs[2].set_xlabel("Time since the satellite rose (s)")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"phase1_tracking_pass{pid}.png"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=DEFAULT_LOG)
    ap.add_argument("--settings", default=None, help="default: ntn_settings.json next to the log")
    ap.add_argument("--pass", dest="pass_id", type=int, default=None, help="pass to plot (default: highest)")
    ap.add_argument("--out", default="results/phase1")
    a = ap.parse_args()

    settings_path = a.settings or os.path.join(os.path.dirname(os.path.abspath(a.log)), "ntn_settings.json")
    if not os.path.exists(settings_path):
        print(f"(no {settings_path}; using the default settings)")
    s = M.load_settings(settings_path)
    df = load_log(a.log)
    os.makedirs(a.out, exist_ok=True)

    print(f"Log: {a.log}  ({len(df)} rows, {df['pass_id'].nunique()} pass(es))\n")
    res = run_checks(df, s)
    with pd.option_context("display.max_colwidth", 60, "display.width", 200):
        print(res[["id", "check", "value", "limit", "result"]].to_string(index=False))
    for _, r in res.iterrows():
        if r["note"]:
            print(f"  {r['id']}: {r['note']}")
    res.to_csv(os.path.join(a.out, "phase1_checks.csv"), index=False)

    g = pick_pass(df, a.pass_id)
    plot_pass(g, s, a.out)
    n_fail = (res["result"] == "FAIL").sum()
    print(f"\n{'ALL CHECKS PASSED' if n_fail == 0 else f'{n_fail} CHECK(S) FAILED'} - "
          f"plots and phase1_checks.csv saved in {a.out}/")


if __name__ == "__main__":
    main()
