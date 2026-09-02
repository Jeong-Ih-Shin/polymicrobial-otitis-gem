"""Figure 1 - shared-lumen genome-scale metabolic modelling framework.

Panel a  model architecture: two reconstructions exchanging metabolites
         through a shared lumen, with system-boundary exchanges set by the
         medium definition
Panel b  analysis workflow

Coauthor revisions applied: CP/PA/SA abbreviations inside the figure,
media boxes reduced to the metabolite count, explanatory italic lines
removed (they duplicate the main text), panel b titles only, box widths
aligned across columns.

Output: figures/Figure1.{pdf,png}
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.bf'] = 'Arial:bold'

FUNGUS, BACT = "#4C9F70", "#3B7EA1"
LUMEN, MEDIUM, GREY = "#E8E2D0", "#B0763F", "#5A5A5A"

FIGS = Path(__file__).resolve().parent.parent / "figures"
FIGS.mkdir(exist_ok=True)

fig, ax = plt.subplots(figsize=(7.09, 4.4))
ax.set_xlim(-3, 103); ax.set_ylim(0, 64); ax.axis("off")

# shared geometry so the columns line up
COL_W = 26
LEFT_X, MID_X, RIGHT_X = 6, 37, 68


def box(x, y, w, h, fc, ec, label, sub=None, fs=9, subfs=7.4, lw=1.1):
    ax.add_patch(FancyBboxPatch((x, y), w, h,
                 boxstyle="round,pad=0.6,rounding_size=1.6",
                 facecolor=fc, edgecolor=ec, linewidth=lw))
    if sub:
        ax.text(x + w / 2, y + h * 0.62, label, ha="center", va="center",
                fontsize=fs, fontweight="bold")
        ax.text(x + w / 2, y + h * 0.26, sub, ha="center", va="center",
                fontsize=subfs, linespacing=1.35)
    else:
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center",
                fontsize=fs, fontweight="bold")


def arrow(x1, y1, x2, y2, color=GREY, lw=1.2, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style,
                 mutation_scale=11, color=color, linewidth=lw,
                 shrinkA=2, shrinkB=2))


# ============================== panel a ====================================
ax.text(-2, 62, "a", fontsize=12, fontweight="bold", va="top")

box(LEFT_X, 44, COL_W, 11, "#EAF3EE", FUNGUS, "CP", "iDC1003")
box(RIGHT_X, 44, COL_W, 11, "#E7F0F5", BACT, "Bacterial partner\n(PA, SA)",
    "AGORA 1.03", fs=8.6, subfs=7.4)

ax.add_patch(FancyBboxPatch((MID_X, 44), COL_W, 11,
             boxstyle="round,pad=0.6,rounding_size=1.6",
             facecolor=LUMEN, edgecolor="#8A7F5F", linewidth=1.3))
ax.text(MID_X + COL_W / 2, 51.4, "Shared lumen", ha="center", va="center",
        fontsize=9, fontweight="bold")
ax.text(MID_X + COL_W / 2, 47.8, "KEGG-to-AGORA\nnamespace mapping",
        ha="center", va="center", fontsize=7.2, linespacing=1.35)

arrow(LEFT_X + COL_W + 0.6, 51.6, MID_X - 0.6, 51.6, FUNGUS, 1.4)
arrow(MID_X - 0.6, 47.6, LEFT_X + COL_W + 0.6, 47.6, FUNGUS, 1.4)
arrow(RIGHT_X - 0.6, 51.6, MID_X + COL_W + 0.6, 51.6, BACT, 1.4)
arrow(MID_X + COL_W + 0.6, 47.6, RIGHT_X - 0.6, 47.6, BACT, 1.4)

# boundary bar aligned to the lumen column
BAR_X, BAR_W = MID_X - 5, COL_W + 10
ax.add_patch(Rectangle((BAR_X, 33.0), BAR_W, 5.2, facecolor="#F3E7D8",
                       edgecolor=MEDIUM, linewidth=1.2))
ax.text(MID_X + COL_W / 2, 35.6, "System-boundary exchanges",
        ha="center", va="center", fontsize=7.4, color="#6B4423")
arrow(MID_X + COL_W / 2, 43.4, MID_X + COL_W / 2, 38.4, MEDIUM, 1.3)

MW = 21
box(MID_X - 12, 22, MW, 8, "#FFFFFF", MEDIUM, "OE medium",
    "155 metabolites", fs=8.6, subfs=7.4)
box(MID_X + 17, 22, MW, 8, "#FFFFFF", MEDIUM, "CSOM medium",
    "157 metabolites", fs=8.6, subfs=7.4)
arrow(MID_X - 12 + MW / 2, 30.4, MID_X + 3, 32.6, MEDIUM, 1.1)
arrow(MID_X + 17 + MW / 2, 30.4, MID_X + COL_W - 3, 32.6, MEDIUM, 1.1)

# ============================== panel b ====================================
ax.text(-2, 18, "b", fontsize=12, fontweight="bold", va="top")

# panel b reads as a sequence of questions, with the method named
# underneath, so the logic is followable without knowing the methods
STEPS = [("Growth\nalone",       "solo FBA"),
         ("Growth\ntogether",    "joint FBA"),
         ("Gain\nor loss",       "boost factor"),
         ("Interaction\ntype",   "classification"),
         ("Cross-\nfeeding",     "pFBA")]

w, gap = 16.4, 4.4
x0 = (100 - (len(STEPS) * w + (len(STEPS) - 1) * gap)) / 2
for i, (what, how) in enumerate(STEPS):
    x = x0 + i * (w + gap)
    ax.add_patch(FancyBboxPatch((x, 1.5), w, 12,
                 boxstyle="round,pad=0.6,rounding_size=1.6",
                 facecolor="#FFFFFF", edgecolor=GREY, linewidth=0.9))
    ax.text(x + w / 2, 9.6, what, ha="center", va="center",
            fontsize=7.6, fontweight="bold", linespacing=1.2)
    ax.text(x + w / 2, 4.6, how, ha="center", va="center",
            fontsize=5.9, color=GREY, style="italic", linespacing=1.2)
    if i < len(STEPS) - 1:
        arrow(x + w + 0.6, 7.5, x + w + gap - 0.6, 7.5, GREY, 1.0)

for ext in ("pdf", "png"):
    fig.savefig(FIGS / f"Figure1.{ext}", bbox_inches="tight", dpi=300)
print("Saved:", FIGS / "Figure1.pdf")
