"""Shared plot look: quiet axes, two fixed series colours, readable text."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SAT = "#2a78d6"      # series 1 (blue): the satellite / the ideal case
DISH = "#eb6834"     # series 2 (orange): the dish / the case with error
THIRD = "#1baf7a"    # series 3 (aqua)
INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
BAND = "#f1f0ec"     # shading for "below the 10 deg elevation mask"

plt.rcParams.update({
    "figure.dpi": 130, "savefig.bbox": "tight", "font.size": 10,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": INK2, "ytick.color": INK2, "lines.linewidth": 2, "legend.frameon": False,
})
