# 5G NTN LEO Satellite Tracking Terminal (CoppeliaSim)

A CoppeliaSim simulation of a ground-station dish that tracks a low-Earth-orbit (LEO) satellite,
keeps a 5G non-terrestrial network (NTN) radio link to it, and uses the link quality to correct
its own pointing.

## What it does

- **Orbit:** a satellite at 600 km altitude (Kepler's 3rd law), with the Earth turning underneath,
  seen from a station in Bengaluru (12.97° N, 77.59° E).
- **Pass prediction:** finds when the satellite will rise, its path across the sky and its highest point.
- **Antenna control:** a pan/tilt dish pre-points to the rise position, then follows the satellite
  with a motor speed limit.
- **Radio link:** Ka-band (20 GHz, 400 MHz) link budget: free-space loss, atmospheric loss,
  pointing loss → RSRP, SINR, throughput, delay and Doppler.
- **Step-track correction:** the dish is installed slightly wrong and can be pushed by wind; the
  terminal probes small offsets and moves toward the strongest signal to remove the error.
- **State machine:** WAIT → PRE-POINT → ACQUIRE → TRACK → CONNECTED ⇄ DEGRADED / LOST → PARK.
- **Visuals:** sky dome, predicted path, trail, beam cone, data packets, status light,
  telemetry panel, two camera views.
- **Log:** every step is written to `ntn_log.csv` in your home folder.

## How to run

1. Open CoppeliaSim (4.6 or newer) with an empty scene.
2. *Add → Script → Simulation script → Non-threaded → Lua*.
3. Paste the contents of `ntn_ground_station.lua` and press **Play**.

The script builds the whole scene itself. To use a scene built by hand instead, name the objects
`PanJoint`, `TiltJoint` (child of PanJoint), `Satellite`, and optionally `StatusLED` and `Beacon`.

## Panel buttons

| Button | Effect |
|---|---|
| Link assist | Turns step-track correction on/off (off → pointing error grows, link degrades) |
| Wind gust | Pushes the dish 1.5° off target |
| Speed | Time speed during a pass: x2, x5, x10, x25 |

## Settings

All tuning values are in the `SET` table at the top of the script (orbit, station location,
radio parameters, thresholds, motor speed).

## Research project: pointing accuracy vs link quality

The project is being extended to answer one question: **how does antenna pointing accuracy affect
the quality of the satellite-to-ground link?**

| Phase | What | Run |
|---|---|---|
| 1 | Validate the tracking (orbit, geometry, joints, controller) against independent Python checks | `python analysis/validate_tracking.py` |
| 2 | Link model for a perfectly pointed dish: range, FSPL, Doppler, delay, received power, C/N | `python analysis/link_budget.py` |
| 3 | Antenna beam model; one pass with different pointing errors; pointing-error budget | `python analysis/pointing_analysis.py` |
| 4 | Error sweep on every pass, error types compared, results table, CoppeliaSim cross-check | `python analysis/evaluate.py` |

Both read `~/ntn_log.csv` and `~/ntn_settings.json`, which the script writes when you press Play
(press Stop before running Python). Results land in `results/phase1` ... `results/phase4`.
Theory and equations for Phases 1-2: [docs/PHASE_1_2.md](docs/PHASE_1_2.md).

**3D demo of a pointing error:** set `biasElDeg = 1.0` in `SET`, press Play: the beam visibly misses the
satellite and the link drops to LINK DEGRADED. `python analysis/evaluate.py` then confirms the power
CoppeliaSim logged equals the Python model.

```
pip install -r analysis/requirements.txt
```

## Offline test

`tests/` runs the script outside CoppeliaSim against a mock API (needs Python + `lupa`):

```
pip install lupa
python tests/run_mock.py
```
