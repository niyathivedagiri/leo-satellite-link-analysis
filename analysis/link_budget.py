"""PHASE 2 - Model the communication link (perfectly pointed dish).

For one satellite pass, computes in Python:
  slant range, one-way delay, range rate, Doppler shift, free-space path loss,
  atmospheric loss, received power (link budget), noise power and C/N.

It also cross-checks Python against the numbers CoppeliaSim logged (Doppler, received power).

Usage (from the project folder):
    python analysis/link_budget.py                      # the pass in ~/ntn_log.csv
    python analysis/link_budget.py --log path/to/ntn_log.csv --pass 2
    python analysis/link_budget.py --synthetic          # no CoppeliaSim needed: Python makes a pass
Outputs go to results/phase2/.
"""
import argparse
import os

import numpy as np
import pandas as pd

import ntn_model as M
import plot_style as S
from runlog import DEFAULT_LOG, load_log, pick_pass

plt = S.plt


def synthetic_pass(s, min_max_el=30.0, step_s=0.5):
    """Find the first pass with max elevation above min_max_el using the Python orbit alone."""
    t = np.arange(0, 2 * 86400, 10.0)
    _, el, _ = M.look_angles(t, s)
    up = el > 0
    k = 0
    while k < len(t):
        if up[k]:
            j = k
            while j < len(t) and up[j]:
                j += 1
            if el[k:j].max() >= min_max_el:
                tt = np.arange(t[max(k - 1, 0)], t[min(j, len(t) - 1)] + step_s, step_s)
                az, e, r = M.look_angles(tt, s)
                keep = e > 0
                return pd.DataFrame({"time_s": tt[keep], "sat_az_deg": az[keep], "sat_el_deg": e[keep],
                                     "range_km": r[keep]})
            k = j
        k += 1
    raise SystemExit("no pass found in 2 days")


def analyse(g, s, p):
    g = g[g["sat_el_deg"] > 0].copy()                       # the link only exists above the horizon
    t = g["time_s"].to_numpy()
    rng, el = g["range_km"].to_numpy(), g["sat_el_deg"].to_numpy()
    out = pd.DataFrame({"time_s": t, "t_rel_s": t - t[0], "el_deg": el, "range_km": rng})
    out["range_rate_kms"] = np.gradient(rng, t)             # Python's own central differences
    out["doppler_khz"] = M.doppler_hz(out["range_rate_kms"], p.f_ghz) / 1e3
    out["doppler_rate_hz_s"] = np.gradient(out["doppler_khz"].to_numpy() * 1e3, t)
    out["delay_ms"] = M.delay_ms(rng)
    lb = M.link_budget(rng, el, p, pointing_loss_db=0.0)    # Phase 2: perfect pointing
    for k, v in lb.items():
        out[k] = v
    # cross-check columns (only when the data came from CoppeliaSim)
    if "doppler_khz" in g:
        out["lua_doppler_khz"] = g["doppler_khz"].to_numpy()
    if "prx_dbm" in g:
        out["lua_prx_dbm"] = g["prx_dbm"].to_numpy()
        out["lua_pointing_error_deg"] = g["pointing_error_deg"].to_numpy()
    return out


def budget_table(out, s):
    """The link budget, line by line, at three moments of the pass."""
    top = int(np.argmax(out["el_deg"]))
    rising = out.iloc[:top + 1]
    picks = {}
    for target in (s["minElevation"], 30.0):
        k = np.where(rising["el_deg"] >= target)[0]
        if len(k):
            picks[f"Rising, el = {target:g} deg"] = k[0]
    picks[f"Highest point (el = {out['el_deg'].iloc[top]:.1f} deg)"] = top
    lines = [("Elevation", "deg", "el_deg", "{:.1f}"), ("Slant range", "km", "range_km", "{:.0f}"),
             ("One-way delay", "ms", "delay_ms", "{:.2f}"), ("Doppler shift", "kHz", "doppler_khz", "{:+.1f}"),
             ("Satellite EIRP", "dBm", "eirp_dbm", "{:.1f}"), ("+ Receive antenna gain", "dBi", "rx_gain_dbi", "{:.1f}"),
             ("- Free-space path loss", "dB", "fspl_db", "{:.2f}"), ("- Atmospheric loss", "dB", "atm_loss_db", "{:.2f}"),
             ("- Pointing loss", "dB", "pointing_loss_db", "{:.2f}"), ("- Other losses", "dB", "misc_loss_db", "{:.2f}"),
             ("= Received power P_rx", "dBm", "prx_dbm", "{:.2f}"), ("RSRP (per resource element)", "dBm", "rsrp_dbm", "{:.2f}"),
             ("Noise power N = kTB", "dBm", "noise_dbm", "{:.2f}"), ("C/N = P_rx - N", "dB", "cn_db", "{:.2f}")]
    rows = []
    for name, unit, col, fmt in lines:
        row = {"Quantity": name, "Unit": unit}
        for label, i in picks.items():
            row[label] = fmt.format(out[col].iloc[i])
        rows.append(row)
    return pd.DataFrame(rows)


def sanity_checks(out, s, p):
    rows = []

    def check(name, value, expect, ok):
        rows.append({"check": name, "value": value, "expected": expect, "result": "PASS" if ok else "FAIL"})

    check("FSPL at 600 km, 20 GHz", f"{M.fspl_db(600, 20):.2f} dB", "174.03 dB", abs(M.fspl_db(600, 20) - 174.03) < 0.01)
    approx = 92.45 + 20 * np.log10(1932) + 20 * np.log10(20)
    check("Exact FSPL vs 92.45 + 20log R + 20log f", f"{abs(M.fspl_db(1932, 20) - approx):.4f} dB", "< 0.01 dB",
          abs(M.fspl_db(1932, 20) - approx) < 0.01)
    check("Delay at 600 km", f"{M.delay_ms(600):.3f} ms", "2.001 ms", abs(M.delay_ms(600) - 2.001) < 0.001)
    n = M.noise_power_dbm(p.t_sys_k, p.bandwidth_mhz * 1e6)
    check("Noise for T_sys from G/T, 400 MHz", f"{n:.2f} dBm (T_sys = {p.t_sys_k:.0f} K)", "about -88.8 dBm",
          abs(n + 88.8) < 0.1)
    g_calc = M.parabolic_gain_dbi(0.6, p.f_ghz, 0.6)
    check("0.6 m dish gain (eta = 0.6) vs the 39.7 dBi used", f"{g_calc:.2f} dBi", "within 0.2 dB of 39.7",
          abs(g_calc - p.rx_gain_dbi) < 0.2)
    check("0.6 m dish beamwidth", f"{M.beamwidth_deg(0.6, p.f_ghz):.2f} deg", "about 1.75 deg",
          abs(M.beamwidth_deg(0.6, p.f_ghz) - 1.75) < 0.02)
    top = int(np.argmax(out["el_deg"]))
    sign_ok = (out["doppler_khz"].iloc[:max(top - 5, 1)] > 0).all() and (out["doppler_khz"].iloc[top + 5:] < 0).all()
    check("Doppler positive while rising, negative while setting", "yes" if sign_ok else "no", "yes", sign_ok)
    vmax = M.orbit_speed_kms(s)
    rr = out["range_rate_kms"].abs().max()
    check("|range rate| below orbital speed", f"{rr:.3f} km/s", f"< {vmax:.3f} km/s", rr < vmax)
    if "lua_doppler_khz" in out:
        d = (out["doppler_khz"] - out["lua_doppler_khz"]).abs().iloc[2:-2].max()
        check("Doppler: Python vs CoppeliaSim log (max diff)", f"{d:.3f} kHz", "< 0.5 kHz", d < 0.5)
    if "lua_prx_dbm" in out:
        ok_rows = out["lua_pointing_error_deg"] < 0.01             # Lua adds pointing loss; compare where it is ~0
        if ok_rows.any():
            d = (out.loc[ok_rows, "prx_dbm"] - out.loc[ok_rows, "lua_prx_dbm"]).abs().max()
            check("Received power: Python vs CoppeliaSim log (max diff)", f"{d:.4f} dB", "< 0.01 dB", d < 0.01)
        else:
            print("(run had a pointing error, so the perfect-pointing power check is skipped; see evaluate.py)")
    return pd.DataFrame(rows)


def plot_link(out, s, out_dir, tag):
    t = out["t_rel_s"].to_numpy()
    fig, axs = plt.subplots(5, 1, figsize=(8, 11), sharex=True)
    low = out["el_deg"].to_numpy() < s["minElevation"]
    for a in axs:
        a.fill_between(t, 0, 1, where=low, transform=a.get_xaxis_transform(), color=S.BAND, lw=0)
    panels = [("el_deg", "Elevation (deg)", "The satellite rises, peaks and sets"),
              ("range_km", "Slant range (km)", "Distance is smallest at the highest point"),
              ("fspl_db", "Free-space path loss (dB)", "Path loss follows distance: 20 log10(R)"),
              ("doppler_khz", "Doppler shift (kHz)", "Doppler: + approaching, 0 at closest point, - leaving"),
              ("prx_dbm", "Received power (dBm)", "Received power, perfectly pointed dish")]
    for a, (col, ylab, title) in zip(axs, panels):
        a.plot(t, out[col], color=S.SAT)
        a.set_ylabel(ylab)
        a.set_title(title, loc="left", fontsize=10)
    axs[3].axhline(0, color=S.INK2, lw=0.8)
    axs[-1].set_xlabel(f"Time since the satellite rose (s)   [shaded: elevation below {s['minElevation']:g} deg]")
    fig.suptitle(f"Phase 2: link over one pass ({tag})", fontweight="bold", y=0.995)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase2_link_over_pass.png"))
    plt.close(fig)

    # received power against elevation: the single most useful link-budget picture
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.plot(out["el_deg"], out["prx_dbm"], color=S.SAT)
    ax.axvline(s["minElevation"], color=S.INK2, lw=0.8, ls=":")
    ax.set_xlabel("Elevation (deg)")
    ax.set_ylabel("Received power (dBm)")
    ax.set_title("Higher satellite = shorter path = stronger signal")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase2_power_vs_elevation.png"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=DEFAULT_LOG)
    ap.add_argument("--settings", default=None)
    ap.add_argument("--pass", dest="pass_id", type=int, default=None)
    ap.add_argument("--synthetic", action="store_true", help="make a pass in Python instead of reading the log")
    ap.add_argument("--out", default="results/phase2")
    a = ap.parse_args()

    settings_path = a.settings or os.path.join(os.path.dirname(os.path.abspath(a.log)), "ntn_settings.json")
    s = M.load_settings(None if a.synthetic else settings_path)
    p = M.LinkParams(f_ghz=s["freqGHz"], bandwidth_mhz=s["bandwidthMHz"], eirp_dbw=s["eirpDbm"] - 30,
                     rx_gain_dbi=s["dishGainDbi"])
    if a.synthetic:
        g, tag = synthetic_pass(s), "Python orbit"
    else:
        g = pick_pass(load_log(a.log), a.pass_id)
        tag = f"CoppeliaSim log, pass {int(g['pass_id'].iloc[0])}"
    os.makedirs(a.out, exist_ok=True)

    print("Link parameters ([A] = assumed reference value, see ntn_model.LinkParams):")
    for k, v in M.params_table(p).items():
        print(f"  {k:22s} {v}")
    out = analyse(g, s, p)
    print(f"\nLink budget ({tag}):")
    tbl = budget_table(out, s)
    print(tbl.to_string(index=False))
    print("\nSanity checks:")
    chk = sanity_checks(out, s, p)
    print(chk.to_string(index=False))

    out.to_csv(os.path.join(a.out, "phase2_link_timeseries.csv"), index=False)
    tbl.to_csv(os.path.join(a.out, "phase2_link_budget.csv"), index=False)
    chk.to_csv(os.path.join(a.out, "phase2_checks.csv"), index=False)
    plot_link(out, s, a.out, tag)
    n_fail = (chk["result"] == "FAIL").sum()
    print(f"\n{'ALL CHECKS PASSED' if n_fail == 0 else f'{n_fail} CHECK(S) FAILED'} - outputs in {a.out}/")


if __name__ == "__main__":
    main()


def load_passes(log=DEFAULT_LOG, settings=None, synthetic=False, pass_id=None, all_complete=False):
    """Shared by Phases 3 and 4: returns settings, link parameters and [(label, Phase 2 table), ...]."""
    from runlog import complete_passes
    settings_path = settings or os.path.join(os.path.dirname(os.path.abspath(log)), "ntn_settings.json")
    s = M.load_settings(None if synthetic else settings_path)
    p = M.LinkParams(f_ghz=s["freqGHz"], bandwidth_mhz=s["bandwidthMHz"], eirp_dbw=s["eirpDbm"] - 30,
                     rx_gain_dbi=s["dishGainDbi"])
    if synthetic:
        out = analyse(synthetic_pass(s), s, p)
        return s, p, [(f"Python pass (max el {out['el_deg'].max():.1f} deg)", out)]
    df = load_log(log)
    if all_complete:
        ids = complete_passes(df)
    else:
        ids = [int(pick_pass(df, pass_id)["pass_id"].iloc[0])]
    res = []
    for i in ids:
        out = analyse(df[df["pass_id"] == i], s, p)
        if all_complete and out["el_deg"].max() < s["minElevation"] + 5:
            print(f"(skipping pass {i}: max elevation {out['el_deg'].max():.1f} deg is too close to the "
                  f"{s['minElevation']:g} deg mask to be useful)")
            continue
        res.append((f"Pass {i} (max el {out['el_deg'].max():.1f} deg)", out))
    return s, p, res
