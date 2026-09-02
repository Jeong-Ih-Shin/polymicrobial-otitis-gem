"""Figure 2 - environment-dependent interaction switching.

Boost factors for each organism in each environment and organism pair.
The classification is derived in code from the same thresholds used
throughout the analysis, so changing LOW/HIGH updates the column labels,
the colour-bar annotations, the colour scale and the classification
together.

The pale midpoint of the colour scale is anchored to the neutral band
rather than to the middle of the plotted range, so that the colour itself
marks where an organism is unaffected. The colour bar describes what
happened to a single organism, because one cell holds one organism's
value; the ecological class above each column follows from the
combination of both organisms and cannot be read off the colour scale.

Output: figures/Figure2.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.bf'] = 'Arial:bold'

FIGS = Path(__file__).resolve().parent.parent / "figures"
FIGS.mkdir(exist_ok=True)

DATA = {("CSOM", "CP+PA"): (1.33, 1.59), ("CSOM", "CP+SA"): (1.76, 1.37),
        ("OE", "CP+PA"): (0.78, 2.54),   ("OE", "CP+SA"): (1.57, 2.47)}
ORDER = [("CSOM", "CP+PA"), ("CSOM", "CP+SA"),
         ("OE", "CP+PA"), ("OE", "CP+SA")]
LOW, HIGH = 0.90, 1.10
VMIN, VMAX = 0.5, 2.6

RED, CREAM, LBLUE, NAVY = "#E64B35", "#FDF3D8", "#7CA6C4", "#3C5488"
_f = lambda v: (v - VMIN) / (VMAX - VMIN)
CMAP = LinearSegmentedColormap.from_list(
    "npg_anchor",
    [(0.0, RED), (_f(LOW), CREAM), (_f(HIGH), CREAM),
     (0.62, LBLUE), (1.0, NAVY)])


def cls(f, b):
    s = lambda v: "+" if v > HIGH else ("-" if v < LOW else "o")
    t = {("+", "+"): "Mutualism", ("-", "-"): "Competition",
         ("o", "o"): "Neutralism", ("+", "-"): "Parasitism",
         ("-", "+"): "Parasitism", ("+", "o"): "Commensalism",
         ("o", "+"): "Commensalism", ("-", "o"): "Amensalism",
         ("o", "-"): "Amensalism"}
    return t[(s(f), s(b))]


fig, ax = plt.subplots(figsize=(5.6, 2.8))
M = np.array([[DATA[k][0] for k in ORDER], [DATA[k][1] for k in ORDER]])
im = ax.imshow(M, cmap=CMAP, vmin=VMIN, vmax=VMAX, aspect="auto")

rgba = im.cmap(im.norm(M))
for i in range(2):
    for j in range(4):
        r, g, b, _ = rgba[i, j]
        lum = 0.299 * r + 0.587 * g + 0.114 * b
        ax.text(j, i, f"{M[i, j]:.2f}\u00d7", ha="center", va="center",
                fontsize=8.5, fontweight="bold",
                color="white" if lum < 0.5 else "#1A1A1A")

ax.set_yticks([0, 1]); ax.set_yticklabels(["CP", "PA or SA"], fontsize=8.5)
ax.set_xticks(range(4))
ax.set_xticklabels([f"{e}\n{p.replace('CP+', 'CP + ')}" for e, p in ORDER],
                   fontsize=7)
for j, k in enumerate(ORDER):
    m = cls(*DATA[k])
    ax.text(j, -0.70, m, ha="center", fontsize=7.2, fontweight="bold",
            color="#1A1A1A")
ax.set_ylim(1.5, -1.05)
for s in ax.spines.values(): s.set_visible(False)
ax.tick_params(length=0)

cb = fig.colorbar(im, ax=ax, orientation="horizontal",
                  fraction=0.10, pad=0.32, aspect=34,
                  ticks=[VMIN, LOW, HIGH, 1.5, 2.0, VMAX])
cb.set_label("Boost factor (joint / solo growth)", fontsize=7, labelpad=4)
cb.ax.tick_params(labelsize=6)
for v in (LOW, HIGH):
    cb.ax.axvline(v, color="#1A1A1A", lw=1.1)
tr = cb.ax.get_xaxis_transform()

cb.ax.text(VMIN, 1.55, "growth suppressed (\u2212)", ha="left",
           va="bottom", fontsize=6.5, fontweight="bold", color="#8C2F1B",
           transform=tr, clip_on=False)
cb.ax.text(VMAX, 1.55, "growth increased (+)",
           ha="right", va="bottom", fontsize=6.5, fontweight="bold",
           color="#2F5488", transform=tr, clip_on=False)
cb.ax.annotate("no effect (o)", xy=(1.0, 1.0), xycoords=tr,
               xytext=(1.20, 1.55), textcoords=tr,
               ha="center", va="bottom", fontsize=6.5, fontweight="bold",
               color="#5A5A5A", clip_on=False,
               arrowprops=dict(arrowstyle="-", lw=0.6, color="#5A5A5A",
                               shrinkA=2, shrinkB=1))

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(FIGS / f"Figure2.{ext}", bbox_inches="tight", dpi=300)
print("Saved:", FIGS / "Figure2.pdf")
for k in ORDER:
    f, b = DATA[k]
    print(f"  {k[0]:5s} {k[1]:6s}  CP={f:.2f} bac={b:.2f}  -> {cls(f, b)}")
