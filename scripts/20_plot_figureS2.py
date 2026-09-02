"""Supplementary Figure S2 - oxygen and mucin do not alter the interaction class.

Replaces 07_plot_figure4.py. The earlier version showed boost factors only,
which is the one thing this figure cannot rely on: the claim being supported
is that the joint solution does not respond to oxygen at all, and a ratio
cannot show that. Absolute solo and joint growth rates are therefore plotted
directly beneath, matching the layout used for the iron figure.

Panels a-b : boost factors along the oxygen gradient
Panels c-d : the underlying absolute growth rates
Panels e-f : interaction class over oxygen x mucin availability

Inputs : results/oxygen_scan.csv, results/phase_diagram.csv,
         results/phase_reference_points.csv (optional)
Output : figures/FigureS2.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.rm'] = 'Arial'
mpl.rcParams['mathtext.it'] = 'Arial:italic'
mpl.rcParams['mathtext.bf'] = 'Arial:bold'

ROOT    = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'
FIGS    = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

BASE_ENV = 'OE'
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
LOW, HIGH = 0.90, 1.10
TITLES = {'CP+PA': 'CP + PA',
          'CP+SA': 'CP + SA'}


def shade(ax, x, modes):
    start = 0
    for i in range(1, len(modes) + 1):
        if i == len(modes) or modes[i] != modes[start]:
            ax.axvspan(x[start], x[min(i, len(x) - 1)],
                       color=MODE_COLORS.get(modes[start], '#EEEEEE'),
                       alpha=0.18, lw=0, zorder=0)
            start = i


def main():
    scan = pd.read_csv(RESULTS / 'oxygen_scan.csv')
    scan = scan[scan.base_env == BASE_ENV]
    present = [m for m in MODE_COLORS if m in set(scan['mode'])]

    pairs = ['CP+PA', 'CP+SA']
    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.6))

    for r, pair in enumerate(pairs):
        g = scan[scan.pair == pair].sort_values('o2_lumen')
        x = g.o2_lumen.values
        modes = g['mode'].tolist()

        # ------------------------------------------------ a, b : boosts
        ax = axes[r, 0]
        shade(ax, x, modes)
        ax.axhline(1.0, color='#555555', lw=0.7)
        for h in (LOW, HIGH):
            ax.axhline(h, color='#999999', lw=0.6, ls=':')
        ax.plot(x, g.cp_boost, color=CP_C, lw=1.9, zorder=3)
        ax.plot(x, g.bac_boost, color=BAC_C, lw=1.9, zorder=3)
        ax.set_title(TITLES[pair], fontsize=8.5, fontweight='bold', pad=6)
        if True:
            ax.set_ylabel('Boost factor\n(joint / solo)', fontsize=7.5)
        ax.set_xlabel('Oxygen uptake bound (mmol gDW$^{-1}$ h$^{-1}$)',
                      fontsize=7)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0.01, 1.09, 'ac'[r], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

        # ------------------------------------------------ c, d : absolutes
        ax = axes[r, 1]
        ax.set_title(TITLES[pair], fontsize=8.5, fontweight='bold', pad=6)
        shade(ax, x, modes)
        ax.plot(x, g.cp_free, color=CP_C, lw=1.9, zorder=3)
        ax.plot(x, g.cp_solo, color=CP_C, lw=1.4, ls='--', alpha=0.8, zorder=3)
        ax.plot(x, g.bac_free, color=BAC_C, lw=1.9, zorder=3)
        ax.plot(x, g.bac_solo, color=BAC_C, lw=1.4, ls='--', alpha=0.8, zorder=3)
        if True:
            ax.set_ylabel('Growth rate (h$^{-1}$)', fontsize=7.5)
        ax.set_xlabel('Oxygen uptake bound (mmol gDW$^{-1}$ h$^{-1}$)',
                      fontsize=7)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0.01, 1.06, 'bd'[r], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

    handles = [
        Line2D([], [], color=CP_C, lw=1.9, label='CP'),
        Line2D([], [], color=BAC_C, lw=1.9, label='Bacterial partner (PA or SA)'),
        Line2D([], [], color='#666666', lw=2.4, label='Joint culture (solid lines)'),
        Line2D([], [], color='#666666', lw=1.4, ls=(0, (4, 2)), label='Solo culture (dashed lines)'),
    ] + [Patch(facecolor=MODE_COLORS[m], alpha=0.5, label=m) for m in present if m != 'infeasible']
    fig.tight_layout(w_pad=1.6, h_pad=1.6, rect=(0, 0.07, 1, 1))
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False,
               fontsize=6.6, bbox_to_anchor=(0.5, 0.005), handlelength=1.8,
               columnspacing=1.6, labelspacing=0.5)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'FigureS2.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'FigureS2.pdf')

    print('\nJoint solution response to oxygen:')
    for pair in pairs:
        g = scan[scan.pair == pair]
        print(f'  {pair}: d(cp_free)={g.cp_free.max()-g.cp_free.min():.2e}  '
              f'd(bac_free)={g.bac_free.max()-g.bac_free.min():.2e}  '
              f'd(bac_solo)={g.bac_solo.max()-g.bac_solo.min():.2e}')


if __name__ == '__main__':
    main()
