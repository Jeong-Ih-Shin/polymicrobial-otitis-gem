"""Figure 3 - iron availability determines the interaction mode.

Panel order follows the coauthor's request, so that the text can walk
through a-d in sequence:

    a  C. parapsilosis + P. aeruginosa   boost factors
    b  C. parapsilosis + P. aeruginosa   absolute growth rates
    c  C. parapsilosis + S. aureus       boost factors
    d  C. parapsilosis + S. aureus       absolute growth rates

The absolute-rate panels are not decoration. Along the oxygen axis the
boost factor rose while both joint growth rates stayed exactly constant,
so the apparent coupling there was a denominator effect; showing the
absolute rates is what separates a real response from that artefact.

The x axis is the iron uptake bound assigned in the medium definition,
not a measured concentration, and the two medium values are labelled on
the axis so that a reader outside the field can see what they are.

Input  : results/axis_scan_iron.csv
Output : figures/Figure3.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
}
LOW, HIGH = 0.90, 1.10
MEDIA = [(0.005, 'Iron level of CSOM', 45, 0),
         (0.010, 'Iron level of OE', 40, 0)]
PAIRS = [('CP+PA', 'CP + PA'), ('CP+SA', 'CP + SA')]
XLABEL = ('Iron uptake bound assigned in the medium\n'
          '(mmol gDW$^{-1}$ h$^{-1}$)')


def shade(ax, x, modes):
    start = 0
    for i in range(1, len(modes) + 1):
        if i == len(modes) or modes[i] != modes[start]:
            ax.axvspan(x[start], x[min(i, len(x) - 1)],
                       color=MODE_COLORS.get(modes[start], '#EEEEEE'),
                       alpha=0.18, lw=0, zorder=0)
            start = i


def mark_media(ax, label=True):
    """Mark the single iron value assigned to each medium.

    The label is offset and joined to the line by a leader, so that it
    reads as pointing at one position on the axis rather than spanning a
    region of it.
    """
    for val, name, dx, _ in MEDIA:
        ax.axvline(val, color='#1A1A1A', lw=1.1, zorder=4)
        if label:
            # label sits on the same baseline as the panel letter, with a
            # leader down to the line so it still reads as one position
            ax.annotate(name,
                        xy=(val, 1.0), xycoords=('data', 'axes fraction'),
                        xytext=(dx, 1.06), textcoords=('offset points',
                                                       'axes fraction'),
                        fontsize=6.4, fontweight='bold',
                        ha='center', va='bottom', linespacing=1.15,
                        annotation_clip=False,
                        arrowprops=dict(arrowstyle='-', lw=0.7,
                                        color='#1A1A1A',
                                        shrinkA=1, shrinkB=1),
                        zorder=6)


def main():
    df = pd.read_csv(RESULTS / 'axis_scan_iron.csv')
    df = df[df.base_env == BASE_ENV]

    fig, axes = plt.subplots(2, 2, figsize=(7.09, 5.6), sharex=True)
    panel = iter('abcd')

    for r, (pair, title) in enumerate(PAIRS):
        g = df[df.pair == pair].sort_values('axis_value_lumen')
        x = -g.axis_value_lumen.values
        modes = g['mode'].tolist()

        # ---------------------------------------------- boost factors
        ax = axes[r, 0]
        shade(ax, x, modes)
        ax.axhline(1.0, color='#555555', lw=0.7)
        for h in (LOW, HIGH):
            ax.axhline(h, color='#999999', lw=0.6, ls=':')
        ax.plot(x, g.cp_boost, color=CP_C, lw=1.9, zorder=3, label='CP')
        ax.plot(x, g.bac_boost, color=BAC_C, lw=1.9, zorder=3,
                label=pair.split('+')[1])
        mark_media(ax, label=True)
        ax.set_title(title, fontsize=8.5, fontweight='bold', pad=16,
                     loc='center')
        ax.set_ylabel('Boost factor\n(joint / solo growth)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(-0.12, 1.06, next(panel), transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

        # ---------------------------------------------- absolute rates
        ax = axes[r, 1]
        ax.set_title(title, fontsize=8.5, fontweight='bold', pad=16, loc='center')
        shade(ax, x, modes)
        ax.plot(x, g.cp_free, color=CP_C, lw=1.9, zorder=3)
        ax.plot(x, g.cp_solo, color=CP_C, lw=1.4, ls='--', alpha=0.85, zorder=3)
        ax.plot(x, g.bac_free, color=BAC_C, lw=1.9, zorder=3)
        ax.plot(x, g.bac_solo, color=BAC_C, lw=1.4, ls='--', alpha=0.85, zorder=3)
        mark_media(ax, label=True)
        ax.set_ylabel('Growth rate (h$^{-1}$)', fontsize=7.5)
        ax.tick_params(labelsize=6.5)
        ax.spines[['top', 'right']].set_visible(False)
        ax.text(-0.12, 1.06, next(panel), transform=ax.transAxes,
                fontsize=11, fontweight='bold', va='bottom')

        # one row title spanning both panels of the pair
        mid = (axes[r, 0].get_position().x0 + axes[r, 1].get_position().x1) / 2
        fig.text(mid, 0, '', fontsize=1)     # placeholder, set after layout
        axes[r, 0].set_title(title, fontsize=8.5, fontweight='bold',
                             pad=16, loc='center')

    for c in (0, 1):
        axes[1, c].set_xlabel(XLABEL, fontsize=7)

    # a single legend outside the axes, so no panel is obscured
    present = [m for m in MODE_COLORS if m in set(df['mode'])]
    handles = [
        Line2D([], [], color=CP_C, lw=1.9, label='CP'),
        Line2D([], [], color=BAC_C, lw=1.9, label='Bacterial partner (PA or SA)'),
        Line2D([], [], color='#666666', lw=2.4, label='Joint culture (solid lines)'),
        Line2D([], [], color='#666666', lw=1.4, ls=(0, (4, 2)), label='Solo culture (dashed lines)'),
        Line2D([], [], color='#999999', lw=0.6, ls=':',
               label='Boost factor thresholds (0.90, 1.10)'),
    ] + [Patch(facecolor=MODE_COLORS[m], alpha=0.5, label=m) for m in present]

    fig.tight_layout(w_pad=2.0, h_pad=2.0, rect=(0, 0.10, 1, 1))
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False,
               fontsize=6.6, bbox_to_anchor=(0.5, 0.005),
               handlelength=1.8, columnspacing=1.6, labelspacing=0.5)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'Figure3.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'Figure3.pdf')

    print('\nTransitions (iron uptake bound, mmol gDW-1 h-1):')
    for pair, _ in PAIRS:
        g = df[df.pair == pair].sort_values('axis_value_lumen')
        m, v = g['mode'].tolist(), (-g.axis_value_lumen).tolist()
        found = False
        for i in range(1, len(m)):
            if m[i] != m[i - 1]:
                print(f'  {pair}: {m[i-1]} -> {m[i]} at {v[i]:.5f}')
                found = True
        if not found:
            print(f'  {pair}: no transition ({m[0]} throughout)')


if __name__ == '__main__':
    main()
