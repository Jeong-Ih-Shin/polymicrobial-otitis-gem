"""Figure 4 - iron availability determines the interaction mode.

Panel a-b : boost factors along the iron gradient, with the interaction
            class shown as background bands
Panel c-d : the underlying absolute growth rates (solo and joint)

Panels c-d are not decoration. In the oxygen scan the boost factor rose
steadily while both joint growth rates stayed exactly constant, so the
apparent coupling was a denominator effect. Showing the absolute rates is
what distinguishes a real response from that artefact, and a reader cannot
check the claim without them.

Input  : results/axis_scan_iron.csv
Output : figures/Figure4_iron_dependence.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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

BASE_ENV = 'OE'            # base medium whose iron content is varied
CP_C, BAC_C = '#4C9F70', '#3B7EA1'
MODE_COLORS = {
    'Mutualism'   : '#3E7B54',
    'Commensalism': '#8FBF9F',
    'Neutralism'  : '#DEDEDE',
    'Amensalism'  : '#E2B07A',
    'Parasitism'  : '#B3472F',
    'Competition' : '#6B3A57',
}
LOW, HIGH = 0.90, 1.10
TITLES = {'CP+PA': '$\\it{C.\\ parapsilosis}$ + $\\it{P.\\ aeruginosa}$',
          'CP+SA': '$\\it{C.\\ parapsilosis}$ + $\\it{S.\\ aureus}$'}


def shade(ax, x, modes):
    start = 0
    for i in range(1, len(modes) + 1):
        if i == len(modes) or modes[i] != modes[start]:
            ax.axvspan(x[start], x[min(i, len(x) - 1)],
                       color=MODE_COLORS.get(modes[start], '#EEEEEE'),
                       alpha=0.18, lw=0, zorder=0)
            start = i


def main():
    df = pd.read_csv(RESULTS / 'axis_scan_iron.csv')
    df = df[df.base_env == BASE_ENV]
    pairs = ['CP+PA', 'CP+SA']

    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.4), sharex='col')

    for c, pair in enumerate(pairs):
        g = df[df.pair == pair].sort_values('axis_value_lumen')
        # plot uptake capacity as a positive magnitude
        x = -g.axis_value_lumen.values
        modes = g['mode'].tolist()

        # ------------------------------------------------ a, b : boosts
        ax = axes[0, c]
        shade(ax, x, modes)
        ax.axhline(1.0, color='#555555', lw=0.7)
        for h in (LOW, HIGH):
            ax.axhline(h, color='#999999', lw=0.6, ls=':')
        ax.plot(x, g.cp_boost, color=CP_C, lw=1.9, zorder=3)
        ax.plot(x, g.bac_boost, color=BAC_C, lw=1.9, zorder=3)

        for i in range(1, len(modes)):
            if modes[i] != modes[i - 1]:
                ax.axvline(x[i], color='#444444', lw=0.8, ls='--', zorder=2)
        for val, name in [(0.005, 'CSOM'), (0.010, 'OE')]:
            ax.axvline(val, color='#1A1A1A', lw=1.1, zorder=4)
            ax.annotate(name, xy=(val, ax.get_ylim()[1]),
                        xytext=(2, -2), textcoords='offset points',
                        fontsize=6.5, fontweight='bold', va='top')

        ax.set_title(TITLES[pair], fontsize=8.5, fontweight='bold', pad=6)
        if c == 0:
            ax.set_ylabel('Boost factor\n(joint / solo)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0.01, 1.10, 'ab'[c], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

        # ------------------------------------------------ c, d : absolutes
        ax = axes[1, c]
        shade(ax, x, modes)
        ax.plot(x, g.cp_free, color=CP_C, lw=1.9, zorder=3)
        ax.plot(x, g.cp_solo, color=CP_C, lw=1.4, ls='--', alpha=0.8, zorder=3)
        ax.plot(x, g.bac_free, color=BAC_C, lw=1.9, zorder=3)
        ax.plot(x, g.bac_solo, color=BAC_C, lw=1.4, ls='--', alpha=0.8, zorder=3)

        for val in (0.005, 0.010):
            ax.axvline(val, color='#1A1A1A', lw=1.1, zorder=4)
        ax.set_xlabel('Inorganic iron uptake capacity\n'
                      '(mmol gDW$^{-1}$ h$^{-1}$)', fontsize=7)
        if c == 0:
            ax.set_ylabel('Growth rate (h$^{-1}$)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(0.01, 1.06, 'cd'[c], transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

    axes[0, 0].legend(handles=[
        Line2D([], [], color=CP_C, lw=1.9, label='$\\it{C.\\ parapsilosis}$'),
        Line2D([], [], color=BAC_C, lw=1.9, label='bacterial partner')],
        fontsize=6.0, frameon=False, loc='upper left')
    axes[1, 0].legend(handles=[
        Line2D([], [], color='#666666', lw=1.9, label='joint'),
        Line2D([], [], color='#666666', lw=1.4, ls='--', label='solo')],
        fontsize=6.0, frameon=False, loc='upper left')

    present = [m for m in MODE_COLORS if m in set(df['mode'])]
    axes[0, 1].legend(
        handles=[Patch(facecolor=MODE_COLORS[m], alpha=0.5, label=m)
                 for m in present],
        fontsize=6.0, frameon=False, loc='center left',
        bbox_to_anchor=(1.01, 0.5))

    fig.tight_layout(w_pad=1.6, h_pad=1.2)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'Figure4_iron_dependence.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'Figure4_iron_dependence.pdf')

    # transition points in real units, for the legend and the text
    print('\nTransitions (iron uptake capacity, mmol gDW-1 h-1):')
    for pair in pairs:
        g = df[df.pair == pair].sort_values('axis_value_lumen')
        m, v = g['mode'].tolist(), (-g.axis_value_lumen).tolist()
        for i in range(1, len(m)):
            if m[i] != m[i - 1]:
                print(f'  {pair}: {m[i-1]} -> {m[i]} at {v[i]:.5f}')


if __name__ == '__main__':
    main()
