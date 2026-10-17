"""PHASE 4 - Evaluate: many pointing errors, every pass, a results table.

1. Sweep a constant pointing error from 0 to 1.5 deg on every complete pass in the log:
   link availability and data per pass vs error.
2. Error TYPE matters: constant offset vs random jitter (Monte Carlo, 30 seeds) vs out-of-date orbit data,
   compared at the same average error.
3. A concise results table (CSV + Markdown) for the report.
4. If the CoppeliaSim run itself had a pointing error (biasElDeg / biasAzDeg set), checks that the received
   power CoppeliaSim logged equals the Python model -> the 3D demo and the analysis tell the same story.

Usage (from the project folder):
    python analysis/evaluate.py                 # all complete passes in ~/ntn_log.csv
    python analysis/evaluate.py --synthetic
Outputs go to results/phase4/.
"""
import argparse
import os

import numpy as np
import pandas as pd

import ntn_model as M
import plot_style as S
import pointing as P
from link_budget import load_passes
from runlog import DEFAULT_LOG, ON_SKY, load_log

plt = S.plt
PASS_COLORS = [S.SAT, S.DISH, S.THIRD, "#eda100"]
SEEDS = range(30)


def sweep(passes, s, p, errors):
    rows = []
    for label, out in passes:
        for e in errors:
            m = P.metrics(P.apply_error(out, P.constant_offset(out, e), s["beamwidthDeg"]), s, p.bandwidth_mhz * 1e6)
            rows.append({"pass": label, "error_deg": e, **m})
    df = pd.DataFrame(rows)
    df["data_vs_perfect_pct"] = df.groupby("pass")["data_gb"].transform(lambda x: 100 * x / x.iloc[0])
    return df


def plot_sweep(df, valid, out_dir):
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    for (label, g), c in zip(df.groupby("pass", sort=False), PASS_COLORS):
        axs[0].plot(g["error_deg"], g["availability_pct"], color=c, label=label)
        axs[1].plot(g["error_deg"], g["data_vs_perfect_pct"], color=c, label=label)
    for a, yl, title in [(axs[0], "Link availability (% of visible time)", "How long the link stays usable"),
                         (axs[1], "Data per pass (% of perfect pointing)", "How much data gets through")]:
        a.axvspan(valid, df["error_deg"].max(), color=S.BAND, lw=0)
        a.set_xlabel("Constant pointing error (deg)")
        a.set_ylabel(yl)
        a.set_title(title)
        a.set_ylim(0, 103)
    axs[0].legend(loc="lower left", fontsize=8)
    axs[1].text(valid + 0.02, 98, "beam model\noptimistic here", fontsize=8, color=S.INK2, va="top")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase4_sweep_all_passes.png"))
    plt.close(fig)


def error_types(out, s, p):
    """Data per pass vs MEAN pointing error, for three kinds of error."""
    bw, rows = p.bandwidth_mhz * 1e6, []
    ideal = P.metrics(P.apply_error(out, P.constant_offset(out, 0), s["beamwidthDeg"]), s, bw)["data_gb"]
    for e in np.arange(0, 1.01, 0.05):
        m = P.metrics(P.apply_error(out, P.constant_offset(out, e), s["beamwidthDeg"]), s, bw)
        rows.append({"type": "Constant offset", "param": e, "mean_err": m["mean_error_deg"],
                     "data_pct": 100 * m["data_gb"] / ideal, "data_sd": 0.0})
    for sig in np.arange(0.025, 0.501, 0.025):
        ms = [P.metrics(P.apply_error(out, P.random_jitter(out, sig, k), s["beamwidthDeg"]), s, bw) for k in SEEDS]
        d = np.array([100 * m["data_gb"] / ideal for m in ms])
        rows.append({"type": "Random jitter", "param": sig, "mean_err": np.mean([m["mean_error_deg"] for m in ms]),
                     "data_pct": d.mean(), "data_sd": d.std(ddof=1)})
    for dly in np.arange(0.25, 4.01, 0.25):
        m = P.metrics(P.apply_error(out, P.stale_orbit_data(out, s, dly), s["beamwidthDeg"]), s, bw)
        rows.append({"type": "Out-of-date orbit data", "param": dly, "mean_err": m["mean_error_deg"],
                     "data_pct": 100 * m["data_gb"] / ideal, "data_sd": 0.0})
    return pd.DataFrame(rows)


def plot_error_types(et, label, out_dir):
    fig, ax = plt.subplots(figsize=(7, 4.3))
    for (typ, g), c in zip(et.groupby("type", sort=False), [S.SAT, S.DISH, S.THIRD]):
        g = g[g["mean_err"] <= 1.0]
        ax.plot(g["mean_err"], g["data_pct"], color=c, label=typ)
        if g["data_sd"].max() > 0:
            ax.fill_between(g["mean_err"], g["data_pct"] - 2 * g["data_sd"], g["data_pct"] + 2 * g["data_sd"],
                            color=c, alpha=0.18, lw=0)
    ax.set_xlabel("Average pointing error over the pass (deg)")
    ax.set_ylabel("Data per pass (% of perfect pointing)")
    ax.set_title(f"Same average error, different damage ({label})", fontsize=10)
    ax.legend(loc="lower left", fontsize=8)
    ax.text(0.99, 0.02, "jitter band: mean +/- 2 SD over 30 seeds", transform=ax.transAxes, ha="right",
            fontsize=8, color=S.INK2)
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "phase4_error_types.png"))
    plt.close(fig)


def results_table(passes, s, p):
    bw, rows = p.bandwidth_mhz * 1e6, []
    for label, out in passes:
        ideal = None
        cases = [("Perfect pointing", [P.constant_offset(out, 0)])]
        cases += [(f"Constant {e:g} deg", [P.constant_offset(out, e)]) for e in (0.25, 0.5, s["beamwidthDeg"] / 2, 1.0)]
        cases += [(f"Orbit data {d:g} s old", [P.stale_orbit_data(out, s, d)]) for d in (1.0, 2.0)]
        cases += [(f"Jitter sigma {j:g} deg (30 seeds)", [P.random_jitter(out, j, k) for k in SEEDS]) for j in (0.2, 0.3)]
        for name, eps_list in cases:
            ms = pd.DataFrame([P.metrics(P.apply_error(out, e, s["beamwidthDeg"]), s, bw) for e in eps_list])
            if ideal is None:
                ideal = ms["data_gb"].mean()
            rows.append({"Pass": label, "Pointing case": name,
                         "Mean error (deg)": round(ms["mean_error_deg"].mean(), 2),
                         "Min C/N (dB)": round(ms["min_cn_db"].mean(), 2),
                         "Availability (%)": round(ms["availability_pct"].mean(), 1),
                         "Data (GB)": round(ms["data_gb"].mean(), 1),
                         "Data vs perfect (%)": round(100 * ms["data_gb"].mean() / ideal, 1)})
    return pd.DataFrame(rows)


def coppelia_crosscheck(log_path, s, p):
    """If the CoppeliaSim run had a pointing error, compare its logged received power with the Python model."""
    df = load_log(log_path)
    trk = df[df["state"].isin(ON_SKY) & (df["sat_el_deg"] > 0)]
    if trk.empty or trk["pointing_error_deg"].median() < 0.05:
        return None
    lb = M.link_budget(trk["range_km"], trk["sat_el_deg"], p,
                       pointing_loss_db=M.pointing_loss_db(trk["pointing_error_deg"], s["beamwidthDeg"]))
    diff = np.abs(lb["prx_dbm"] - trk["prx_dbm"].to_numpy())
    return {"samples": len(trk), "median_error_deg": trk["pointing_error_deg"].median(),
            "max_prx_diff_db": float(np.nanmax(diff))}


def to_markdown(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--log", default=DEFAULT_LOG)
    ap.add_argument("--settings", default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--out", default="results/phase4")
    a = ap.parse_args()
    s, p, passes = load_passes(a.log, a.settings, a.synthetic, all_complete=True)
    os.makedirs(a.out, exist_ok=True)
    th = s["beamwidthDeg"]
    e_grid = np.linspace(0, 1.5, 61)
    valid = e_grid[np.abs(M.airy_loss_db(e_grid, th) - M.pointing_loss_db(e_grid, th)) <= 1.0].max()

    print(f"Passes analysed: {', '.join(l for l, _ in passes)}")
    sw = sweep(passes, s, p, e_grid)
    sw.to_csv(os.path.join(a.out, "phase4_sweep.csv"), index=False)
    plot_sweep(sw, valid, a.out)

    print("\nLargest constant error that keeps the link usable for the WHOLE visible pass (availability 100 %):")
    for label, g in sw.groupby("pass", sort=False):
        ok = g[g["availability_pct"] >= 99.999]["error_deg"].max()
        half = g[g["data_vs_perfect_pct"] >= 50]["error_deg"].max()
        print(f"  {label}: {ok:.3f} deg   (data falls to half at about {half:.2f} deg)")

    best_label, best = max(passes, key=lambda x: x[1]["el_deg"].max())
    et = error_types(best, s, p)
    et.to_csv(os.path.join(a.out, "phase4_error_types.csv"), index=False)
    plot_error_types(et, best_label, a.out)
    at = lambda typ, e: np.interp(e, *et[et["type"] == typ][["mean_err", "data_pct"]].to_numpy().T)  # noqa: E731
    print(f"\nAt the same AVERAGE error of 0.4 deg ({best_label}), data per pass is:")
    for typ in ("Constant offset", "Random jitter", "Out-of-date orbit data"):
        print(f"  {typ:24s} {at(typ, 0.4):5.1f} % of perfect")

    tbl = results_table(passes, s, p)
    tbl.to_csv(os.path.join(a.out, "phase4_results_table.csv"), index=False)
    with open(os.path.join(a.out, "phase4_results_table.md"), "w") as f:
        f.write(to_markdown(tbl) + "\n")
    print("\nResults table (usable = C/N >= %g dB, above %g deg elevation):" % (P.CN_THRESHOLD_DB, s["minElevation"]))
    print(tbl.to_string(index=False))

    if not a.synthetic:
        cc = coppelia_crosscheck(a.log, s, p)
        if cc:
            print(f"\nCoppeliaSim run had pointing error (median {cc['median_error_deg']:.2f} deg): logged received "
                  f"power vs Python model, max difference {cc['max_prx_diff_db']:.3f} dB over {cc['samples']} samples "
                  f"-> {'PASS' if cc['max_prx_diff_db'] < 0.01 else 'CHECK'}")
    print(f"\nOutputs in {a.out}/")


if __name__ == "__main__":
    main()
