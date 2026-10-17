"""PHASE 3 - How does pointing error change the link?

1. Antenna pattern: the Gaussian main-beam model (used everywhere) vs the textbook Airy pattern
   -> shows where the simple model is valid.
2. One pass, four cases (perfect, 0.5 deg offset, 1.0 deg offset, orbit data 1 s out of date):
   pointing error, received power, C/N and Doppler over time.
3. A scenario table for that pass, and the pointing-error budget the link can afford.

Usage (from the project folder):
    python analysis/pointing_analysis.py                  # pass with the highest elevation in ~/ntn_log.csv
    python analysis/pointing_analysis.py --synthetic      # no CoppeliaSim needed
Outputs go to results/phase3/.
"""
import argparse
import os

import numpy as np
import pandas as pd

import ntn_model as M
import plot_style as S
import pointing as P
from link_budget import load_passes
from runlog import DEFAULT_LOG

plt = S.plt
COLORS = [S.SAT, S.DISH, S.THIRD, "#eda100"]


def plot_pattern(theta3, out_dir):
    e = np.linspace(0, 1.3 * theta3, 400)
    g, a = M.pointing_loss_db(e, theta3), M.airy_loss_db(e, theta3)
    valid = e[np.abs(a - g) <= 1.0].max()            # where the two agree within 1 dB
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.plot(e, g, color=S.SAT, label="Gaussian model: 12.04 (e / theta3)^2  (used)")
    ax.plot(e, a, color=S.DISH, lw=1.4, ls=(0, (4, 3)), label="Airy pattern, uniform dish (check)")
    ax.axvspan(valid, e.max(), color=S.BAND, lw=0)
    ax.text(valid + 0.02, 27, "models differ\nby > 1 dB", fontsize=8, color=S.INK2, va="top")
    for x, y, txt in [(theta3 / 2, 3.01, " half beamwidth: 3 dB"), (theta3, 12.04, " full beamwidth: 12 dB")]:
        ax.plot([x], [y], "o", color=S.INK, ms=4)
        ax.annotate(txt, (x, y), fontsize=8, color=S.INK, xytext=(4, -10), textcoords="offset points")
    ax.set_ylim(0, 30)
    ax.set_xlabel("Pointing error (deg)")
    ax.set_ylabel("Loss in received power (dB)")
    ax.set_title(f"How much signal a {theta3:.2f} deg-wide beam loses when it misses")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase3_antenna_pattern.png"))
    plt.close(fig)
    return valid


def plot_pass(out, cases, s, label, out_dir):
    t = out["t_rel_s"].to_numpy()
    fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
    low = out["el_deg"].to_numpy() < s["minElevation"]
    for a in axs:
        a.fill_between(t, 0, 1, where=low, transform=a.get_xaxis_transform(), color=S.BAND, lw=0)
    for (name, eps), c in zip(cases, COLORS):
        r = P.apply_error(out, eps, s["beamwidthDeg"])
        axs[0].plot(t, eps, color=c, label=name, lw=1.6)
        axs[1].plot(t, r["prx_dbm"], color=c, lw=1.6)
        axs[2].plot(t, r["cn_db"], color=c, lw=1.6)
    axs[2].axhline(P.CN_THRESHOLD_DB, color=S.INK2, lw=0.8, ls=":")
    axs[2].text(t[-1], P.CN_THRESHOLD_DB, f"usable above {P.CN_THRESHOLD_DB:g} dB ", fontsize=8, color=S.INK2,
                ha="right", va="bottom")
    axs[3].plot(t, out["doppler_khz"], color=S.SAT, lw=1.6)
    axs[3].axhline(0, color=S.INK2, lw=0.8)
    axs[0].set_ylabel("Pointing error (deg)")
    axs[0].legend(loc="upper left", fontsize=8, ncol=2)
    axs[0].set_title(f"{label}: same pass, four pointing cases", loc="left")
    axs[1].set_ylabel("Received power (dBm)")
    axs[2].set_ylabel("C/N (dB)")
    axs[3].set_ylabel("Doppler (kHz)")
    axs[3].set_title("Doppler is the same in every case: it depends on motion, not on pointing", loc="left",
                     fontsize=9)
    axs[3].set_xlabel(f"Time since the satellite rose (s)   [shaded: below {s['minElevation']:g} deg]")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase3_pass_cases.png"))
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=DEFAULT_LOG)
    ap.add_argument("--settings", default=None)
    ap.add_argument("--pass", dest="pass_id", type=int, default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--out", default="results/phase3")
    a = ap.parse_args()
    s, p, passes = load_passes(a.log, a.settings, a.synthetic, a.pass_id)
    label, out = passes[0]
    os.makedirs(a.out, exist_ok=True)
    th = s["beamwidthDeg"]

    valid = plot_pattern(th, a.out)
    print(f"Beamwidth theta3 = {th} deg. Gaussian model within 1 dB of the Airy pattern up to {valid:.2f} deg "
          f"({valid / th:.2f} theta3).")

    cases = [("Perfect pointing", P.constant_offset(out, 0)),
             ("Constant 0.5 deg", P.constant_offset(out, 0.5)),
             ("Constant 1.0 deg", P.constant_offset(out, 1.0)),
             ("Orbit data 1 s old", P.stale_orbit_data(out, s, 1.0))]
    plot_pass(out, cases, s, label, a.out)

    scen = [("Perfect pointing", P.constant_offset(out, 0))]
    scen += [(f"Constant {e:g} deg", P.constant_offset(out, e)) for e in (0.25, 0.5, th / 2, 1.0, 1.5)]
    scen += [(f"Orbit data {d:g} s old", P.stale_orbit_data(out, s, d)) for d in (0.5, 1.0, 2.0)]
    scen += [(f"Random jitter sigma {j:g} deg", P.random_jitter(out, j, seed=1)) for j in (0.1, 0.2, 0.3)]
    tbl = P.run_scenarios(out, s, p, scen)
    show = tbl[["scenario", "mean_error_deg", "max_error_deg", "mean_pointing_loss_db", "min_cn_db",
                "availability_pct", "data_vs_perfect_pct"]]
    print(f"\n{label}: link quality above {s['minElevation']:g} deg elevation (usable = C/N >= {P.CN_THRESHOLD_DB:g} dB)")
    print(show.round(2).to_string(index=False))
    tbl.to_csv(os.path.join(a.out, "phase3_scenarios.csv"), index=False)

    # pointing-error budget: how much error can the link afford at its weakest point?
    vis = out[out["el_deg"] >= s["minElevation"]]
    margin = vis["cn_db"].min() - P.CN_THRESHOLD_DB
    print("\nPointing-error budget:")
    print(f"  Weakest point of the pass (el = {vis['el_deg'].iloc[vis['cn_db'].argmin()]:.1f} deg): "
          f"C/N = {vis['cn_db'].min():.2f} dB -> margin {margin:.2f} dB above {P.CN_THRESHOLD_DB:g} dB")
    for name, m in [("keep loss <= 1 dB", 1.0), ("keep loss <= 3 dB", 3.0),
                    ("never lose the link on this pass", margin)]:
        print(f"  {name:34s}: error <= {P.allowed_error_deg(m, th):.2f} deg")
    print(f"\nPlots and phase3_scenarios.csv saved in {a.out}/")


if __name__ == "__main__":
    main()
