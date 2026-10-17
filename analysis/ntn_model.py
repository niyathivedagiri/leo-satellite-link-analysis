"""ntn_model.py - the physics, written once and reused by every analysis script.

Phase 1 uses: orbit + geometry (to re-check the CoppeliaSim run independently).
Phase 2 uses: free-space path loss, Doppler, delay, noise and the link budget.

Units are in every name: _km, _deg, _db, _dbm, _hz, _ghz, _s ...
Nothing here imports CoppeliaSim; everything is plain NumPy.
"""
import json
import os
from dataclasses import dataclass, asdict

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
C_KMS = 299_792.458          # speed of light, km/s (exact, by definition)
K_BOLTZMANN_DBW = -228.6     # 10*log10(1.380649e-23 J/K), in dBW/K/Hz
MU_STANDARD = 398_600.4418   # Earth's gravitational parameter, km^3/s^2 (WGS-84)

# The values the Lua script uses if no ntn_settings.json is found (same as SET in the script)
DEFAULT_SETTINGS = {
    "altitudeKm": 600, "inclinationDeg": 25, "raanDeg": 300, "startDeg": 60,
    "stationLat": 12.97, "stationLon": 77.59, "freqGHz": 20, "bandwidthMHz": 400,
    "eirpDbm": 60, "dishGainDbi": 39.7, "beamwidthDeg": 1.75, "slewDegPerSec": 6,
    "biasAzDeg": 0, "biasElDeg": 0, "linkAssist": False, "minElevation": 10,
    "earthRadiusKm": 6371, "mu": 398600.4, "earthSpin": 7.2921e-5,
}


def load_settings(path=None):
    """Settings written by the Lua script at Play (ntn_settings.json), or the defaults."""
    s = dict(DEFAULT_SETTINGS)
    if path and os.path.exists(path):
        with open(path) as f:
            s.update(json.load(f))
    if s["mu"] == 398600:
        print("WARNING: ntn_settings.json has mu = 398600 (rounded by an older copy of the Lua script). "
              "The script itself uses 398600.4 - edit the file or re-paste the latest script.")
    return s


# ---------------------------------------------------------------------------
# PHASE 1 - orbit and geometry
# Same equations as the Lua script, written a different way (rotation matrices
# instead of expanded formulas), so a bug in one is unlikely to be in the other.
# ---------------------------------------------------------------------------
def _rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def orbit_radius_km(s):
    return s["earthRadiusKm"] + s["altitudeKm"]


def mean_motion_rad_s(s, mu=None):
    """n = sqrt(mu / a^3): how fast the satellite goes round its orbit (rad/s). Kepler's third law."""
    return np.sqrt((mu or s["mu"]) / orbit_radius_km(s) ** 3)


def orbit_period_s(s, mu=None):
    return 2 * np.pi / mean_motion_rad_s(s, mu)


def orbit_speed_kms(s, mu=None):
    """Circular orbit speed v = sqrt(mu / a)."""
    return np.sqrt((mu or s["mu"]) / orbit_radius_km(s))


def satellite_ecef_km(t_s, s, mu=None):
    """Satellite position in the Earth-fixed frame at orbit time t (array of times OK).

    1. position in the orbit plane:   [a cos u, a sin u, 0],  u = u0 + n t
    2. tilt by inclination i (about X), turn by RAAN (about Z)  -> inertial frame (ECI)
    3. the Earth has turned by theta = omega_E * t               -> Earth-fixed frame (ECEF)
    """
    t = np.atleast_1d(np.asarray(t_s, dtype=float))
    a = orbit_radius_km(s)
    u = np.radians(s["startDeg"]) + mean_motion_rad_s(s, mu) * t
    in_plane = np.stack([a * np.cos(u), a * np.sin(u), np.zeros_like(u)])          # 3 x N
    eci = _rot_z(np.radians(s["raanDeg"])) @ _rot_x(np.radians(s["inclinationDeg"])) @ in_plane
    theta = s["earthSpin"] * t
    ecef = np.empty_like(eci)
    for k in range(t.size):                     # R_z(-theta) for each time
        ecef[:, k] = _rot_z(-theta[k]) @ eci[:, k]
    return ecef


def station_ecef_km(s):
    """Ground station on a spherical Earth: R_e [cos(lat)cos(lon), cos(lat)sin(lon), sin(lat)]."""
    la, lo = np.radians(s["stationLat"]), np.radians(s["stationLon"])
    return s["earthRadiusKm"] * np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])


def enu_matrix(s):
    """Rows are the station's East, North and Up unit vectors in ECEF."""
    la, lo = np.radians(s["stationLat"]), np.radians(s["stationLon"])
    east = [-np.sin(lo), np.cos(lo), 0]
    north = [-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)]
    up = [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
    return np.array([east, north, up])


def look_angles(t_s, s, mu=None):
    """Azimuth (deg, from North towards East), elevation (deg) and slant range (km)."""
    d = satellite_ecef_km(t_s, s, mu) - station_ecef_km(s)[:, None]      # station -> satellite
    e, n, u = enu_matrix(s) @ d
    rng = np.sqrt(e ** 2 + n ** 2 + u ** 2)
    az = np.degrees(np.arctan2(e, n)) % 360
    el = np.degrees(np.arcsin(u / rng))
    return az, el, rng


def slant_range_from_elevation_km(el_deg, s):
    """Closed form on a spherical Earth: R = sqrt(a^2 - Re^2 cos^2(el)) - Re sin(el).
    Independent of the orbit code, so it is a good check of the geometry."""
    re, a = s["earthRadiusKm"], orbit_radius_km(s)
    el = np.radians(el_deg)
    return np.sqrt(a ** 2 - (re * np.cos(el)) ** 2) - re * np.sin(el)


def unit_vector(az_deg, el_deg):
    """Direction for an az/el in the local frame (x = East, y = North, z = Up)."""
    az, el = np.radians(az_deg), np.radians(el_deg)
    return np.stack([np.sin(az) * np.cos(el), np.cos(az) * np.cos(el), np.sin(el)])


def angle_between_deg(az1, el1, az2, el2):
    """Pointing error = angle between two directions (deg)."""
    d = np.sum(unit_vector(az1, el1) * unit_vector(az2, el2), axis=0)
    return np.degrees(np.arccos(np.clip(d, -1, 1)))


def wrap180(a_deg):
    return (np.asarray(a_deg) + 180) % 360 - 180


# ---------------------------------------------------------------------------
# PHASE 2 - communication link
# ---------------------------------------------------------------------------
def fspl_db(range_km, f_ghz):
    """Free-space path loss, exact form: 20 log10(4 pi R f / c). R in km, f in GHz.
    Equivalent to 92.45 + 20 log10(R_km) + 20 log10(f_GHz)."""
    r_m = np.asarray(range_km) * 1e3
    f_hz = f_ghz * 1e9
    return 20 * np.log10(4 * np.pi * r_m * f_hz / (C_KMS * 1e3))


def doppler_hz(range_rate_kms, f_ghz):
    """f_d = -f_c * Rdot / c. Positive while the satellite approaches (range shrinking)."""
    return -np.asarray(range_rate_kms) / C_KMS * f_ghz * 1e9


def delay_ms(range_km):
    """One-way propagation delay tau = R / c."""
    return np.asarray(range_km) / C_KMS * 1e3


def atmospheric_loss_db(el_deg, zenith_db=0.5, min_el_deg=5.0):
    """Clear-sky gas absorption grows with the path length through the air: A_zenith / sin(el).
    Simple 'cosecant' model; not valid near the horizon, so el is clamped at 5 deg."""
    return zenith_db / np.sin(np.radians(np.maximum(el_deg, min_el_deg)))


def noise_power_dbm(t_sys_k, bandwidth_hz):
    """Thermal noise N = k T B, in dBm."""
    return K_BOLTZMANN_DBW + 10 * np.log10(t_sys_k) + 10 * np.log10(bandwidth_hz) + 30


def parabolic_gain_dbi(diameter_m, f_ghz, efficiency):
    """Ideal dish gain G = eta (pi D / lambda)^2 (sanity check for the 39.7 dBi value)."""
    lam = C_KMS * 1e3 / (f_ghz * 1e9)
    return 10 * np.log10(efficiency * (np.pi * diameter_m / lam) ** 2)


def beamwidth_deg(diameter_m, f_ghz):
    """Half-power beamwidth of a dish, rule of thumb theta_3dB ~ 70 lambda / D (degrees)."""
    lam = C_KMS * 1e3 / (f_ghz * 1e9)
    return 70 * lam / diameter_m


@dataclass
class LinkParams:
    """Every number in the link budget, with where it comes from.

    [C] calculated from physics   [A] assumed / taken from a reference parameter set
    Reference set: 3GPP TR 38.821, Ka-band LEO-600 downlink to a VSAT-class terminal.
    """
    f_ghz: float = 20.0               # [A] Ka-band downlink carrier
    bandwidth_mhz: float = 400.0      # [A] 5G NR channel bandwidth at FR2
    eirp_dbw: float = 30.0            # [A] satellite EIRP: 4 dBW/MHz x 400 MHz = 4 + 26.0 dB
    rx_gain_dbi: float = 39.7         # [A] terminal antenna gain (0.6 m dish); check: eta=0.6 -> 39.8 dBi [C]
    g_over_t_dbk: float = 15.9        # [A] terminal figure of merit G/T
    zenith_atm_db: float = 0.5        # [A] clear-sky gas absorption straight up at 20 GHz
    misc_loss_db: float = 0.5         # [A] polarisation, scintillation, implementation
    n_rb: int = 264                   # [A] resource blocks in 400 MHz at 120 kHz subcarrier spacing

    @property
    def t_sys_k(self):
        """[C] system noise temperature from G/T: T = 10^((G - G/T)/10)."""
        return 10 ** ((self.rx_gain_dbi - self.g_over_t_dbk) / 10)


def link_budget(range_km, el_deg, p=LinkParams(), pointing_loss_db=0.0):
    """Received power and carrier-to-noise ratio, term by term (all arrays OK).

    P_rx = EIRP + G_rx - FSPL - L_atm - L_point - L_misc      [dBW, then +30 -> dBm]
    C/N  = P_rx - N,   N = k T_sys B
    Below the horizon the Earth blocks the path, so P_rx is set to NaN there.
    """
    el = np.asarray(el_deg, dtype=float)
    out = {
        "eirp_dbm": np.full_like(el, p.eirp_dbw + 30),
        "rx_gain_dbi": np.full_like(el, p.rx_gain_dbi),
        "fspl_db": fspl_db(range_km, p.f_ghz),
        "atm_loss_db": atmospheric_loss_db(el, p.zenith_atm_db),
        "pointing_loss_db": np.broadcast_to(np.asarray(pointing_loss_db, dtype=float), el.shape).copy(),
        "misc_loss_db": np.full_like(el, p.misc_loss_db),
    }
    prx = (out["eirp_dbm"] + out["rx_gain_dbi"] - out["fspl_db"] - out["atm_loss_db"]
           - out["pointing_loss_db"] - out["misc_loss_db"])
    prx = np.where(el > 0, prx, np.nan)
    out["prx_dbm"] = prx
    out["rsrp_dbm"] = prx - 10 * np.log10(12 * p.n_rb)            # power per resource element
    out["noise_dbm"] = np.full_like(el, noise_power_dbm(p.t_sys_k, p.bandwidth_mhz * 1e6))
    out["cn_db"] = prx - out["noise_dbm"]
    return out


def params_table(p=LinkParams()):
    """The parameter list as rows, for the report."""
    d = asdict(p)
    d["t_sys_k (derived)"] = round(p.t_sys_k, 1)
    return d


# ---------------------------------------------------------------------------
# PHASE 3 - antenna pattern and pointing loss
# ---------------------------------------------------------------------------
GAUSS_COEFF_DB = 40 * np.log10(np.e) * np.log(2)     # = 12.04 dB


def pointing_loss_db(error_deg, beamwidth_deg):
    """Gaussian main-beam model: G(e) = G0 * exp(-4 ln2 (e/theta3)^2)
    In dB this is a parabola: L = 12.04 * (e / theta3)^2.
    By construction L = 3.01 dB at e = theta3/2 (the half-power point).
    Valid inside the main lobe only (about e < 1.2 theta3); beyond that the real pattern has nulls/sidelobes."""
    return GAUSS_COEFF_DB * (np.asarray(error_deg) / beamwidth_deg) ** 2


def bessel_j1(x):
    """J1(x) from its integral form (1/pi) * int_0^pi cos(t - x sin t) dt - avoids needing SciPy."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    tau = np.linspace(0, np.pi, 801)
    trap = getattr(np, "trapezoid", None) or np.trapz          # NumPy 2 renamed trapz
    return trap(np.cos(tau[None, :] - x[:, None] * np.sin(tau[None, :])), tau, axis=1) / np.pi


AIRY_U_HALF = 1.6163                                   # [2 J1(u)/u]^2 = 0.5 here


def airy_loss_db(error_deg, beamwidth_deg):
    """Uniformly illuminated circular aperture, scaled to the same half-power beamwidth:
    G(e)/G0 = [2 J1(u)/u]^2 with u = 1.6163 * e / (theta3/2).
    This is the textbook diffraction pattern: it has a first null at e = 1.19 theta3 and sidelobes,
    which the Gaussian model ignores. Used only to check where the Gaussian model is valid."""
    u = 2 * AIRY_U_HALF * np.atleast_1d(np.asarray(error_deg, dtype=float)) / beamwidth_deg
    u = np.where(np.abs(u) < 1e-9, 1e-9, u)
    ratio = (2 * bessel_j1(u) / u) ** 2
    return -10 * np.log10(np.maximum(ratio, 1e-6))


def spectral_efficiency(snr_db, alpha=0.6, se_max=4.4, snr_min_db=-10.0):
    """Attenuated Shannon bound (bit/s/Hz): alpha*log2(1+SNR), capped; 0 below SNR_min.
    An estimate of what a practical modem achieves, not a measurement."""
    snr_db = np.asarray(snr_db, dtype=float)
    se = alpha * np.log2(1 + 10 ** (snr_db / 10))
    return np.where(np.isfinite(snr_db) & (snr_db >= snr_min_db), np.minimum(se, se_max), 0.0)
