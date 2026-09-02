"""Figure 5 - predicted bacterial response to antifungal intervention.

Panel a : dose-response, bacterial growth vs degree of fungal suppression
Panel b : bacterial growth under fungistatic vs fungicidal removal

Input  : results/antifungal_simulation.csv
Output : figures/Figure5_antifungal.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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

COND_COLOR = {
    ('CSOM', 'CP+PA'): '#3B7EA1',
    ('CSOM', 'CP+SA'): '#4C9F70',
    ('OE',   'CP+PA'): '#B3472F',
    ('OE',   'CP+SA'): '#B0763F',
}
LABEL = {'CP+PA': '$\\it{P.\\ aeruginosa}$', 'CP+SA': '$\\it{S.\\ aureus}$'}


def main():
    df = pd.read_csv(RESULTS / 'antifungal_simulation.csv')
    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.0),
                             gridspec_kw={'width_ratios': [1.25, 1]})

    # ------------------------------------------------------- panel a
    ax = axes[0]
    dose = df[df.scenario == 'dose_response']
    for (env, pair), g in dose.groupby(['env', 'pair']):
        g = g.sort_values('fungal_fraction')
        # x axis: degree of suppression, 0 = untreated, 1 = fully suppressed
        ax.plot(1 - g.fungal_fraction, 100 * g.bac_rel_to_coinfection,
                color=COND_COLOR.get((env, pair), '#666666'), lw=1.8,
                label=f'{env}, {LABEL[pair]}')
    ax.axhline(100, color='#888888', lw=0.7, ls=':')
    ax.set_xlabel('Degree of fungal suppression\n'
                  '(0 = untreated co-infection, 1 = growth fully arrested)',
                  fontsize=7)
    ax.set_ylabel('Predicted bacterial growth\n(% of co-infection)',
                  fontsize=7.5)
    ax.set_xlim(0, 1)
    ax.tick_params(labelsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6.0, frameon=False, loc='best')
    ax.text(0.01, 1.04, 'a', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom')

    # ------------------------------------------------------- panel b
    ax = axes[1]
    sub = df[df.scenario.isin(['fungistatic', 'fungicidal'])].copy()
    sub['cond'] = sub.env + '\n' + sub.pair
    conds = sorted(sub.cond.unique())
    x = np.arange(len(conds))
    w = 0.36
    for k, (scen, hatch) in enumerate([('fungistatic', ''),
                                       ('fungicidal', '///')]):
        vals = [100 * sub[(sub.cond == c) &
                          (sub.scenario == scen)].bac_rel_to_coinfection.mean()
                for c in conds]
        cols = [COND_COLOR.get(tuple(c.split('\n')), '#666666') for c in conds]
        ax.bar(x + (k - 0.5) * w, vals, w, color=cols,
               alpha=1.0 if k == 0 else 0.55, hatch=hatch,
               edgecolor='white', linewidth=0.6,
               label=scen)
    ax.axhline(100, color='#888888', lw=0.7, ls=':')
    ax.set_xticks(x)
    ax.set_xticklabels([c.replace('CP+', '') for c in conds], fontsize=6.2)
    ax.set_ylabel('Predicted bacterial growth\n(% of co-infection)',
                  fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6.2, frameon=False)
    ax.text(0.02, 1.04, 'b', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom')

    fig.tight_layout(w_pad=2.0)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'Figure5_antifungal.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'Figure5_antifungal.pdf')


if __name__ == '__main__':
    main()
