"""Figure 4 - oxygen dependence of interaction mode.

Panel a-b : boost factor trajectories along the oxygen gradient (1D scan)
Panel c-d : interaction-class phase diagram over oxygen x mucin (2D scan)

Inputs : results/oxygen_scan.csv, results/phase_diagram.csv,
         results/phase_reference_points.csv
Output : figures/Figure4_oxygen_dependence.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 8, 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'axes.linewidth': 0.8,
})

ROOT    = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
FIGS    = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

CP_C, BAC_C = '#4C9F70', '#3B7EA1'
MODE_COLORS = {
    'Mutualism'   : '#3E7B54',
    'Commensalism': '#8FBF9F',
    'Neutralism'  : '#DEDEDE',
    'Amensalism'  : '#E2B07A',
    'Parasitism'  : '#B3472F',
    'Competition' : '#6B3A57',
    'infeasible'  : '#FFFFFF',
}
BASE_ENV = 'CSOM'          # which 1D scan base medium to plot
LOW, HIGH = 0.90, 1.10


def shade_modes(ax, x, modes):
    """Background bands showing the interaction class along the axis."""
    start = 0
    for i in range(1, len(modes) + 1):
        if i == len(modes) or modes[i] != modes[start]:
            ax.axvspan(x[start], x[min(i, len(x) - 1)],
                       color=MODE_COLORS.get(modes[start], '#EEEEEE'),
                       alpha=0.16, lw=0)
            start = i


def main():
    scan = pd.read_csv(RESULTS / 'oxygen_scan.csv')
    phase = pd.read_csv(RESULTS / 'phase_diagram.csv')
    try:
        ref = pd.read_csv(RESULTS / 'phase_reference_points.csv')
    except FileNotFoundError:
        ref = None

    pairs = ['CP+PA', 'CP+SA']
    titles = {'CP+PA': '$\\it{C.\\ parapsilosis}$ + $\\it{P.\\ aeruginosa}$',
              'CP+SA': '$\\it{C.\\ parapsilosis}$ + $\\it{S.\\ aureus}$'}

    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.9))

    # ------------------------------------------------ panels a, b : 1D scan
    for c, pair in enumerate(pairs):
        ax = axes[0, c]
        g = scan[(scan.base_env == BASE_ENV) & (scan.pair == pair)]
        g = g.sort_values('o2_lumen')
        x = g.o2_lumen.values

        shade_modes(ax, x, g['mode'].tolist())
        ax.axhline(1.0, color='#666666', lw=0.7, ls='-')
        ax.axhline(HIGH, color='#999999', lw=0.6, ls=':')
        ax.axhline(LOW, color='#999999', lw=0.6, ls=':')

        ax.plot(x, g.cp_boost, color=CP_C, lw=1.8)
        ax.plot(x, g.bac_boost, color=BAC_C, lw=1.8)

        # mark the transition
        m = g['mode'].tolist()
        for i in range(1, len(m)):
            if m[i] != m[i - 1]:
                ax.axvline(x[i], color='#B3472F', lw=0.9, ls='--')
                ax.annotate(f'{m[i-1]}\n\u2193\n{m[i]}',
                            xy=(x[i], ax.get_ylim()[1]),
                            xytext=(3, -3), textcoords='offset points',
                            fontsize=5.8, va='top', color='#B3472F')
                break

        ax.set_title(titles[pair], fontsize=8.5, fontweight='bold', pad=6)
        ax.set_xlabel('Oxygen uptake bound (mmol gDW$^{-1}$ h$^{-1}$)',
                      fontsize=7)
        if c == 0:
            ax.set_ylabel('Boost factor\n(joint / solo growth)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0.02, 1.06, 'ab'[c], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

    axes[0, 1].legend(handles=[
        Line2D([], [], color=CP_C, lw=1.8, label='$\\it{C.\\ parapsilosis}$'),
        Line2D([], [], color=BAC_C, lw=1.8, label='bacterial partner')],
        fontsize=6.2, frameon=False, loc='upper left')

    # ---------------------------------------------- panels c, d : 2D phase
    present = [m for m in MODE_COLORS if m in set(phase['mode'])]
    cmap = ListedColormap([MODE_COLORS[m] for m in present])
    norm = BoundaryNorm(range(len(present) + 1), cmap.N)
    idx = {m: i for i, m in enumerate(present)}

    for c, pair in enumerate(pairs):
        ax = axes[1, c]
        g = phase[phase.pair == pair]
        piv = g.pivot_table(index='mucin_frac', columns='o2_lumen',
                            values='mode', aggfunc='first')
        Z = np.vectorize(lambda m: idx.get(m, 0))(piv.values)

        ax.pcolormesh(piv.columns.values, piv.index.values, Z,
                      cmap=cmap, norm=norm, shading='nearest')

        if ref is not None:
            for _, r in ref.iterrows():
                ax.plot(r.o2_lumen, r.mucin_frac, marker='*',
                        ms=11, color='white', mec='black', mew=0.7, zorder=5)
                ax.annotate(r.env, (r.o2_lumen, r.mucin_frac),
                            xytext=(5, 5), textcoords='offset points',
                            fontsize=6.5, fontweight='bold')

        ax.set_xlabel('Oxygen uptake bound (mmol gDW$^{-1}$ h$^{-1}$)',
                      fontsize=7)
        if c == 0:
            ax.set_ylabel('Mucin substrate availability\n'
                          '(fraction of CSOM level)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.text(0.02, 1.06, 'cd'[c], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

    axes[1, 1].legend(
        handles=[Patch(facecolor=MODE_COLORS[m], label=m) for m in present
                 if m != 'infeasible'],
        fontsize=6.0, frameon=False, loc='center left',
        bbox_to_anchor=(1.01, 0.5))

    fig.tight_layout(w_pad=1.8, h_pad=2.2)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'Figure4_oxygen_dependence.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'Figure4_oxygen_dependence.pdf')


if __name__ == '__main__':
    main()
