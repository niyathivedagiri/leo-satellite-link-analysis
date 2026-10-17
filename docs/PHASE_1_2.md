# Phases 1 and 2: Validate the tracking, then model the link

**Research question:** how does antenna pointing accuracy affect the quality of a LEO satellite-to-ground link?

These two phases build the foundation for answering it:

- **Phase 1:** prove the tracking simulation is correct.
- **Phase 2:** compute what the radio link looks like when the dish points perfectly.

Phase 3 will then add pointing error and measure the damage against this perfect baseline.

```
 CoppeliaSim (Lua)                          Python (analysis/)
 ─────────────────                          ──────────────────
 orbit → az/el/range ─┐
 dish joints ─────────┼─► ~/ntn_log.csv ──► validate_tracking.py  (Phase 1: is it right?)
 3D dish direction ───┘   ~/ntn_settings.json
                                        └─► link_budget.py        (Phase 2: what is the link?)
                                              uses ntn_model.py  (all the equations, once)
```

---

## What I found in your existing script

| Part | Status | Action |
|---|---|---|
| Orbit (circular, Kepler, Earth rotation), az/el/range | Correct | Kept. Python re-checks it independently (V2, V3) |
| Pass prediction, pre-point, park | Works | Kept |
| Rate-limited pan/tilt control | Works | Kept. Python checks the limit (V7) |
| Pointing error = angle between dish and satellite directions | Correct | Kept. Python checks it is ~0 while tracking (V6) |
| **Hidden mount bias (0.6°, −0.4°) and step-track ON by default** | Hides the baseline | **Set to 0 / OFF** for Phase 1. Phase 3 will add error on purpose, in a controlled way |
| **Doppler used a forward difference over 1 s** | Off by half a second | **Changed to a central difference** (t − 0.5 s, t + 0.5 s) |
| **No check that the 3D dish really points where the maths says** | Missing | **Added a marker on the dish** whose direction is compared every step (V5) |
| **Log too coarse, and missing dish angles and range rate** | Can't validate | **New columns + more decimals + `ntn_settings.json`** |
| `?` before the state name on the panel | Encoding problem | Replaced with an HTML entity (`&#9679;` = ●) |

Nothing else changed. The scene, visuals, panel, buttons and state machine are the same.

---

## Phase 1: Validate satellite tracking

### Theory you need (minimum)

**1. Three coordinate frames.**

- **ECI (inertial):** fixed to the stars. The orbit is simplest here, as a circle tilted by inclination *i* and turned by RAAN *Ω*.
- **ECEF (Earth-fixed):** turns with the Earth, once every 86,164 s. The station is fixed in this frame.
- **ENU (East-North-Up):** a local frame centred on the station. Azimuth and elevation are read in this frame.

The chain is: satellite in ECI → rotate by the Earth's angle → ECEF → subtract the station → rotate into ENU.

**2. Orbit speed (Kepler's third law), for a circular orbit of radius a = R_e + h.**

- n = √(μ / a³) is the angular speed (rad/s), with μ = 398,600 km³/s².
- T = 2π / n is the period. At 600 km, T = 5,792 s = **96.5 min**.
- v = √(μ / a) = **7.56 km/s**.

**3. Look angles.** From the ENU vector (E, N, U) of the station-to-satellite line:

- Azimuth = atan2(E, N): the compass direction, clockwise from North.
- Elevation = asin(U / R): the angle above the horizon.
- Range R = √(E² + N² + U²): the straight-line (slant) distance.

**4. Slant range from elevation alone (closed form).** On a spherical Earth, the range at a given elevation does not depend on the orbit details:

> R(el) = √(a² − R_e² cos² el) − R_e sin el

This gives 600 km at 90°, 1,075 km at 30° and **1,932 km at 10°**. It's a second, independent way to check the geometry.

**5. Pointing error.** It is the angle between two unit vectors, where the dish points (**d**) and where the satellite is (**s**):

> ε = arccos(**d** · **s**)

Note that an azimuth error matters less near zenith (it is multiplied by cos el). That's why we use the vector angle, not "az error + el error".

**6. Why the motors can lose the satellite near zenith (the keyhole).** When the satellite passes almost overhead, its azimuth swings through 180° in a few seconds. The required azimuth rate is roughly (angular speed seen from the ground) ÷ (angular distance from zenith). At 600 km the satellite crosses the sky at about 0.72°/s near the top. With a 6°/s motor limit, any pass higher than about **83°** will be rate-limited. That is real tracking error, not a bug. V6b reports it when it happens.

### What the checks prove

| Check | Proves | How |
|---|---|---|
| V1 | The orbit speed is right | Prints T and v (expect 96.5 min, 7.56 km/s) |
| V2 | The Lua orbit and geometry code is right | Python re-computes az/el/range at every logged time with **different code** (rotation matrices). Max difference must be < 0.001° and < 0.01 km |
| V3 | The geometry is right a second way | Logged range vs the closed form R(el) |
| V4 | Range rate (the Doppler input) is right | Range rate must be exactly 0 at maximum elevation (closest approach) |
| V5 | **The 3D model points where the maths says** | A marker on the dish; its axis vs the computed direction. Fails if a joint turns the wrong way |
| V6 | The tracking controller is right | Pointing error < 0.01° in ≥ 99% of samples above 10° (except when the motors are saturated) |
| V7 | The motor limit is respected | Fastest joint motion ≤ 6°/s |
| V8 | The run is a clean baseline | Bias 0, step-track off |

**Why V2 is not circular:** it uses the same *equations*, written as different *code*. A typo in one (a sign, a sin/cos swap, degrees vs radians) would not be in the other. V3 and V4 don't use the orbit code at all.

**An honest point for the interview:** V6 will show a pointing error of essentially **zero**. That's because the simulated controller knows exactly where the satellite is: no ephemeris error, no bias, no wind. This is the *ideal* reference case. Real dishes always have some pointing error, and Phase 3 adds it deliberately so we can measure its effect.

---

## Phase 2: Model the communication link

### Theory you need (minimum)

**1. Decibels.** Losses and gains multiply, so in dB they add.

- dB = 10 log₁₀(power ratio).
- dBm = power relative to 1 mW, so 0 dBm = 1 mW and −80 dBm = 10 picowatts.
- dBW = power relative to 1 W, so dBm = dBW + 30.

**2. Free-space path loss (FSPL).** The signal spreads over a sphere, so its power density falls as 1/R². Then, how much of it a receiving antenna captures depends on wavelength.

> FSPL = 20 log₁₀(4π R f / c)  =  92.45 + 20 log₁₀(R in km) + 20 log₁₀(f in GHz)   [dB]

- Doubling the distance adds **6 dB** of loss.
- At 20 GHz: 174.0 dB at 600 km and 184.2 dB at 1,932 km. A low pass loses about **10 dB more** than an overhead one.

**3. Link budget.** This is the bookkeeping from transmitter to receiver:

> P_rx = EIRP + G_rx − FSPL − L_atm − L_point − L_other   [dBm]

| Term | Meaning | Value used | Source |
|---|---|---|---|
| EIRP | Power the satellite radiates towards us, including its antenna gain | 30 dBW = **60 dBm** (4 dBW/MHz × 400 MHz) | Assumed: 3GPP TR 38.821 Ka-band LEO-600 set |
| G_rx | Gain of our 0.6 m dish | **39.7 dBi** | Assumed, same set. Checked: η(πD/λ)² with η = 0.6 gives 39.8 dBi |
| FSPL | Spreading loss | 174–188 dB | Calculated from range |
| L_atm | Oxygen and water vapour absorption: 0.5 dB straight up, more when low (0.5 / sin el) | 0.5–2.9 dB | Simple model; clamped below 5° |
| L_point | Loss from not pointing exactly at the satellite | **0 in Phase 2** | Phase 3 models this |
| L_other | Polarisation, scintillation, implementation | 0.5 dB | Assumed |

**4. Noise and C/N.** The receiver adds thermal noise:

> N = k T B,  where k = 1.38×10⁻²³ J/K, T = system noise temperature, B = bandwidth

- **T** comes from the terminal's figure of merit: G/T = 15.9 dB/K, so T = 10^((39.7 − 15.9)/10) ≈ **240 K**.
- Over B = 400 MHz, N = **−88.8 dBm**.
- **C/N = P_rx − N** (in dB) is the link-quality number. It goes from about **+1 dB at 10°** to **+12 dB near the top** of a 63° pass.

**5. Doppler shift.** The satellite moves at 7.56 km/s, so the received frequency shifts:

> f_d = −f_c · Ṙ / c

- Ṙ is the range rate: negative while approaching, so f_d is **positive** while the satellite rises.
- At 20 GHz the shift reaches about **±420–460 kHz**. That's more than three times a 5G subcarrier (120 kHz), so a real terminal **must** correct it. It crosses 0 exactly at closest approach.

**6. Delay.** τ = R / c, which is 2.0 ms overhead and 6.4 ms at 10°.

**7. RSRP.** This is the 5G measure of power per resource element: RSRP = P_rx − 10 log₁₀(12 × 264) = P_rx − 35.0 dB.

### What Phase 2 produces

- **A link-budget table** at 10°, 30° and the highest point, line by line, like a textbook budget.
- **Five plots on one time axis:** elevation, range, FSPL, Doppler and received power. They show that everything is driven by geometry.
- **Received power against elevation.** This is the single most useful link-budget picture.
- **Cross-checks against CoppeliaSim:** Python's Doppler and received power must match the Lua log (< 0.5 kHz and < 0.01 dB).

**One thing to explain if asked:** the received-power curve has a small kink near 5° elevation. That's where the atmospheric model is clamped (the cosecant law blows up at the horizon). It's a model limitation, and only matters below the 10° mask, where we don't use the link.

---

## How to run it on your Mac

**1. Install the Python tools (one time):**

```bash
cd ~/path/to/5g-ntn-satellite-tracker
python3 -m venv .venv && source .venv/bin/activate
pip install -r analysis/requirements.txt
```

**2. Run the simulation:**

1. Open your scene (or a new one), open the simulation script, delete the old code, and paste the new `ntn_ground_station.lua`.
2. Press **Play**. Don't press any buttons; the wind and step-track would spoil the baseline.
3. Wait until the state goes PARK → WAIT at least once. That is one full pass, about 3–4 minutes at the default speed.
4. Press **Stop**. The log is only fully written when the simulation stops.

The console prints `Logging to /Users/<you>/ntn_log.csv`.

**3. Run Phase 1:**

```bash
python analysis/validate_tracking.py
```

**4. Run Phase 2:**

```bash
python analysis/link_budget.py
```

**No CoppeliaSim handy?** `python analysis/link_budget.py --synthetic` makes a pass in Python only.

Results land in `results/phase1/` and `results/phase2/`.

### What you should see

The model is deterministic, so your numbers should be very close to these, from the same script run outside CoppeliaSim:

- **Phase 1:** `ALL CHECKS PASSED`.
    - V1 is 5792.3 s (96.54 min) / 7.562 km/s.
    - V2 differences are below 1e-4°.
    - V4 gap is under 0.5 s.
    - V5 is about 1e-6°.
    - V6 is 100%.
    - V7 is 6.000°/s.
- **Phase 2:** first pass, max elevation **62.6°**:

| | el = 10° | el = 30° | Highest (62.6°) |
|---|---|---|---|
| Slant range | 1,932 km | 1,074 km | 668 km |
| Delay | 6.44 ms | 3.58 ms | 2.23 ms |
| Doppler | +420 kHz | +354 kHz | ~0 kHz |
| FSPL | 184.19 dB | 179.09 dB | 174.97 dB |
| Received power | −87.9 dBm | −80.9 dBm | −76.3 dBm |
| C/N | 0.9 dB | 7.9 dB | 12.5 dB |

### If a check fails

| Failure | Most likely cause | Fix |
|---|---|---|
| "Log not found" | Simulation not stopped, or the log is somewhere else | Press Stop; use `--log /path/ntn_log.csv` |
| "from the old script" | Old Lua still in the scene | Paste the new script again |
| V5 fails (large angle) | A hand-built joint turns the wrong way | Set `panSign` or `tiltSign = -1` in `SET` |
| V6 / V8 fail | A button was pressed, or bias/assist was changed | Re-run without touching the buttons |
| V2 fails | A change to the orbit code in Lua | Compare with `ntn_model.look_angles` |

---

## Interview questions you should be able to answer

1. Why do you need three coordinate frames, and what does each one hold fixed?
2. How did you prove the 3D antenna really points where your equations say? (V5: the marker on the dish.)
3. Why is the pointing error zero in Phase 1, and why is that not realistic? (Perfect knowledge, so it is the ideal reference.)
4. Why does a 10° pass have about 10 dB less signal than an overhead pass? (Range 1,932 vs 600 km gives 20 log₁₀(3.2) ≈ 10 dB.)
5. Why is Doppler largest near the horizon and zero at the top? (Range rate is greatest when moving along the line of sight; it is zero at closest approach.)
6. Which link-budget values are physics, and which are assumptions? (FSPL, delay and Doppler are physics; EIRP, G/T and the losses are assumed reference values.)
7. What is C/N, and why does 400 MHz of bandwidth give −88.8 dBm of noise?

---

## Next: Phase 3 (after both phases pass on your Mac)

We will:

- Model the dish's beam pattern, deriving the loss formula 12(ε/θ₃dB)² from the beam shape.
- Inject controlled pointing errors (0.1°–1.5°).
- Compare received power and C/N against the perfect baseline from Phase 2.
