"""Load the CoppeliaSim run log (ntn_log.csv) and split it into passes."""
import os

import numpy as np
import pandas as pd

ON_SKY = {"ACQUIRE", "TRACK", "CONNECTED", "DEGRADED", "LOST"}
DEFAULT_LOG = os.path.expanduser("~/ntn_log.csv")


def load_log(path=DEFAULT_LOG):
    if not os.path.exists(path):
        raise SystemExit(f"Log not found: {path}\nRun the scene in CoppeliaSim first (it writes ~/ntn_log.csv), "
                         "or pass --log <path>.")
    df = pd.read_csv(path)
    needed = {"time_s", "state", "sat_az_deg", "sat_el_deg", "range_km", "dish_az_deg", "scene_check_deg"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"{path} is from the old script (missing {sorted(missing)}). "
                         "Paste the new ntn_ground_station.lua into CoppeliaSim and run again.")
    # a new pass starts every time the state machine enters PRE-POINT
    starts = (df["state"] == "PREPOINT") & (df["state"].shift() != "PREPOINT")
    df["pass_id"] = starts.cumsum()
    df = df[df["pass_id"] > 0].copy()
    return df


def complete_passes(df):
    """Pass ids that reached PARK (the satellite rose AND set during the run)."""
    return [p for p, g in df.groupby("pass_id") if (g["state"] == "PARK").any() and g["state"].isin(ON_SKY).any()]


def pick_pass(df, pass_id=None):
    """The requested pass, or the complete pass with the highest maximum elevation."""
    if pass_id is not None:
        return df[df["pass_id"] == pass_id].copy()
    done = complete_passes(df) or sorted(df["pass_id"].unique())
    best = max(done, key=lambda p: df.loc[df["pass_id"] == p, "sat_el_deg"].max())
    return df[df["pass_id"] == best].copy()


def time_since_rise(g):
    """Seconds since the satellite crossed 0 deg elevation in this pass."""
    above = g.index[g["sat_el_deg"] > 0]
    t0 = g.loc[above[0], "time_s"] if len(above) else g["time_s"].iloc[0]
    return g["time_s"].to_numpy() - t0


def rate_deg_s(angle_deg, t_s, wrap=False):
    """Rate of change of an angle along the log (central differences)."""
    a = np.unwrap(np.radians(angle_deg)) if wrap else np.radians(angle_deg)
    return np.degrees(np.gradient(a, t_s))
