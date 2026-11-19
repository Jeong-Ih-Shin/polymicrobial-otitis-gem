"""Figure 4 - predicted bacterial growth after removal of C. parapsilosis.

Replaces the earlier version, which was drawn from 08 and included a
fungistatic scenario now known to be a structural artefact of pinning the
fungal biomass at zero while leaving fungal metabolism intact.

Panel a : bacterial growth as a function of remaining fungal metabolic
          activity. The response is flat over almost the whole range and
          falls only once fungal activity is nearly abolished, so the
          x axis is drawn on a scale that resolves the region near zero.
Panel b : bacterial growth at complete elimination, by condition.

Input  : results/fungal_activity_budget.csv   (from 14)
Output : figures/Figure4.{pdf,png}
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

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

COND_COLOR = {
    ('CSOM', 'CP+PA'): '#3B7EA1',
    ('CSOM', 'CP+SA'): '#4C9F70',
    ('OE',   'CP+PA'): '#B3472F',
    ('OE',   'CP+SA'): '#B0763F',
}
SHORT = {'CP+PA': 'PA', 'CP+SA': 'SA'}

# 세균별로 묶어 환경 간 비교가 보이도록 한 순서
ORDER = [('CSOM', 'CP+PA'), ('OE', 'CP+PA'),
         ('CSOM', 'CP+SA'), ('OE', 'CP+SA')]


def main():
    src = RESULTS / 'fungal_activity_budget.csv'
    if not src.exists():
        raise SystemExit(f'{src} not found — run 14_fungal_activity_budget.py')
    df = pd.read_csv(src)
    if 'fungal_activity' not in df.columns:
        raise SystemExit('unexpected columns; expected fungal_activity')

    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.0),
                             gridspec_kw={'width_ratios': [1.35, 1]})

    # -------------------------------------------------------- panel a
    ax = axes[0]
    for env, pair in ORDER:
        g = df[(df.env == env) & (df.pair == pair)]
        if not len(g):
            continue
        g = g.sort_values('fungal_activity', ascending=False)
        ax.plot(100 * g.fungal_activity, 100 * g.bac_rel_to_coinfection,
                color=COND_COLOR.get((env, pair), '#666666'), lw=1.8,
                marker='o', ms=2.2, label=f'{env}, {SHORT.get(pair, pair)}')
    ax.axhline(100, color='#888888', lw=0.7, ls=':')
    ax.set_xlabel('Remaining fungal metabolic activity\n'
                  '(% of unsuppressed joint-model value)', fontsize=7)
    ax.set_ylabel('Predicted bacterial growth\n(% of unsuppressed joint-model value)',
                  fontsize=7.5)
    ax.set_xlim(102, -2)          # untreated on the left, removal on right
    ax.tick_params(labelsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6.0, frameon=False, loc='lower left')
    ax.text(-0.14, 1.04, 'a', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom')

    # inset resolving the region near complete removal
    near = df[df.fungal_activity <= 0.20]
    if len(near):
        ins = ax.inset_axes([0.44, 0.18, 0.34, 0.46])
        for env, pair in ORDER:
            g = near[(near.env == env) & (near.pair == pair)]
            if not len(g):
                continue
            g = g.sort_values('fungal_activity', ascending=False)
            ins.plot(100 * g.fungal_activity, 100 * g.bac_rel_to_coinfection,
                     color=COND_COLOR.get((env, pair), '#666666'), lw=1.4,
                     marker='o', ms=2.0)
        ins.set_xlim(20, -1)
        ins.tick_params(labelsize=5.2)
        ins.set_title('near complete removal', fontsize=5.6, pad=2)
        for s in ('top', 'right'):
            ins.spines[s].set_visible(False)

    # -------------------------------------------------------- panel b
    ax = axes[1]
    recs = []
    for (env, pair), g in df.groupby(['env', 'pair']):
        g = g.sort_values('fungal_activity', ascending=False)
        hit = g[g.bac_rel_to_coinfection < 0.999]
        if not len(hit):
            continue
        recs.append(dict(env=env, pair=pair,
                         cond=env + '\n' + pair.replace('CP+', ''),
                         thr=100 * hit.iloc[0].fungal_activity))
    thr_df = (pd.DataFrame(recs)
              .set_index(['env', 'pair'])
              .reindex(pd.MultiIndex.from_tuples(ORDER, names=['env', 'pair']))
              .dropna(subset=['thr'])
              .reset_index())
    x = np.arange(len(thr_df))
    cols = [COND_COLOR.get((r.env, r.pair), '#666666')
            for _, r in thr_df.iterrows()]
    vals = thr_df.thr.values
    ax.bar(x, vals, 0.62, color=cols, edgecolor='white', linewidth=0.6)
    for xi, v in zip(x, vals):
        ax.text(xi, v + 0.6, f'{v:.1f}%', ha='center', fontsize=6.4,
                fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(thr_df.cond, fontsize=6.2)
    ax.set_ylim(0, max(vals) * 1.25)
    ax.set_ylabel('Fungal activity at which bacterial growth\n'
                  'first declines (% of unsuppressed value)',
                  fontsize=7.5)
    ax.tick_params(labelsize=6.5)
    ax.spines[['top', 'right']].set_visible(False)
    ax.text(-0.16, 1.04, 'b', transform=ax.transAxes,
            fontsize=11, fontweight='bold', va='bottom')

    fig.tight_layout(w_pad=2.0)
    for ext in ('pdf', 'png'):
        fig.savefig(FIGS / f'Figure4.{ext}',
                    bbox_inches='tight', dpi=300)
    print('Saved:', FIGS / 'Figure4.pdf')

    # where the response first departs from the untreated value
    print('\nFirst departure from 100% (grid resolution '
          f'{df.fungal_activity.diff().abs().max():.3f}):')
    for (env, pair), g in df.groupby(['env', 'pair']):
        g = g.sort_values('fungal_activity', ascending=False)
        hit = g[g.bac_rel_to_coinfection < 0.999]
        first = f's={hit.iloc[0].fungal_activity:.3f}' if len(hit) else 'never'
        final = 100 * g.iloc[-1].bac_rel_to_coinfection
        print(f'  {env:5s} {pair:6s}  {first:<10} '
              f'final {final:.1f}%')


if __name__ == '__main__':
    main()
