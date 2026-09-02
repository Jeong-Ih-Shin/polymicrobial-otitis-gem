"""Build Table 5 and Supplementary Table S2 from the result files.

Table 5   component attribution for the CSOM-to-OE interaction switch
Table S2  absolute solo and joint growth rates behind every boost factor

Each table is written twice: a CSV for the record and a markdown file that
can be pasted into Word and converted with Insert > Table > Convert Text
to Table.

Inputs:
    results/component_attribution.csv   (from 10)
    results/pair_results.csv            (from 01)
    results/axis_scan_iron.csv          (from 12, optional)
    results/fungal_activity_budget.csv  (from 14, optional)

Outputs:
    results/Table5_component_attribution.{csv,md}
    results/TableS2_absolute_growth.{csv,md}
"""

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / 'results'

GROUP_LABEL = {
    'oxygen'      : 'Oxygen',
    'mucin'       : 'Mucin-associated substrates',
    'ammonium'    : 'Ammonium',
    'fermentation': 'Fermentation products',
    'iron'        : 'Inorganic iron',
    'amino_acids' : 'Amino acids',
    'sugars'      : 'Sugars',
}
PAIR_LABEL = {'CP+PA': 'C. parapsilosis + P. aeruginosa',
              'CP+SA': 'C. parapsilosis + S. aureus'}


def to_markdown(df, path, title, note=''):
    lines = [f'**{title}**', '']
    lines.append('| ' + ' | '.join(df.columns) + ' |')
    lines.append('| ' + ' | '.join('---' for _ in df.columns) + ' |')
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join(str(v) for v in r.values) + ' |')
    if note:
        lines += ['', note]
    path.write_text('\n'.join(lines), encoding='utf-8')


# ------------------------------------------------------------- Table 5
def table5():
    src = RESULTS / 'component_attribution.csv'
    if not src.exists():
        print(f'  missing {src.name} — run 10_component_attribution.py')
        return
    df = pd.read_csv(src)

    rows = []
    for pair in ['CP+PA', 'CP+SA']:
        g = df[df.pair == pair]
        ref = g[g.analysis == 'reference'].set_index('component')['mode']
        m_csom = ref.get('CSOM_reference', '-')
        m_oe = ref.get('OE_reference', '-')

        suf = g[g.analysis == 'sufficiency'].set_index('component')['mode']
        nec = g[g.analysis == 'necessity'].set_index('component')['mode']

        for grp in sorted(set(suf.index) | set(nec.index)):
            s, n = suf.get(grp, '-'), nec.get(grp, '-')
            if m_csom == m_oe:
                verdict = 'not applicable'
            else:
                is_s = (s == m_oe)
                is_n = (n == m_csom)
                verdict = ('Necessary and sufficient' if is_s and is_n else
                           'Necessary' if is_n else
                           'Sufficient' if is_s else
                           'Partial effect' if s not in (m_csom, '-') else
                           'No effect')
            rows.append({
                'Organism pair': PAIR_LABEL[pair],
                'Component group': GROUP_LABEL.get(grp, grp),
                'CSOM reference': m_csom,
                'OE reference': m_oe,
                'Sufficiency (CSOM + group at OE value)': s,
                'Necessity (OE with group reverted)': n,
                'Interpretation': verdict,
            })

    t = pd.DataFrame(rows)
    t.to_csv(RESULTS / 'Table5_component_attribution.csv', index=False)
    to_markdown(
        t, RESULTS / 'Table5_component_attribution.md',
        'Table 5. Attribution of the CSOM-to-OE interaction switch to '
        'individual medium components.',
        'Sufficiency was tested by setting a single component group to its '
        'OE value in the CSOM medium; necessity by reverting a single group '
        'to its CSOM value in the OE medium. A group is interpreted as '
        'necessary when its reversion restores the CSOM interaction class.')
    print(f'  Table 5: {len(t)} rows')
    for _, r in t[t.Interpretation.str.contains('Necessary')].iterrows():
        print(f'    {r["Organism pair"]}: {r["Component group"]} '
              f'-> {r["Interpretation"]}')


# ------------------------------------------------------------ Table S2
def table_s2():
    src = RESULTS / 'pair_results.csv'
    if not src.exists():
        print(f'  missing {src.name} — run 01_simulate.py')
        return
    df = pd.read_csv(src)

    rows = []
    for _, r in df.iterrows():
        rows.append({
            'Environment': r['env'],
            'Organism pair': PAIR_LABEL.get(r['pair'], r['pair']),
            'C. parapsilosis solo (h-1)': f"{r['cp_solo']:.4f}",
            'C. parapsilosis joint (h-1)': f"{r['cp_free']:.4f}",
            'C. parapsilosis boost': f"{r['cp_boost']:.2f}x",
            'Bacterium solo (h-1)': f"{r['bac_solo']:.4f}",
            'Bacterium joint (h-1)': f"{r['bac_free']:.4f}",
            'Bacterium boost': f"{r['bac_boost']:.2f}x",
        })
    t = pd.DataFrame(rows)

    note = ('Boost factors are ratios of joint to solo growth. Absolute '
            'rates are reported so that a change in the ratio can be '
            'distinguished from a change in the solo reference alone.')

    iron = RESULTS / 'axis_scan_iron.csv'
    if iron.exists():
        g = pd.read_csv(iron)
        g = g[g.base_env == 'OE']
        extra = []
        for pair, gg in g.groupby('pair'):
            gg = gg.sort_values('axis_value_lumen')
            lo, hi = gg.iloc[0], gg.iloc[-1]
            extra.append(
                f'{PAIR_LABEL.get(pair, pair)}: across the iron scan, '
                f'joint growth ranged {gg.cp_free.min():.3f}-'
                f'{gg.cp_free.max():.3f} (fungus) and '
                f'{gg.bac_free.min():.3f}-{gg.bac_free.max():.3f} '
                f'(bacterium), while bacterial solo growth ranged '
                f'{gg.bac_solo.min():.3f}-{gg.bac_solo.max():.3f}.')
        note += ' ' + ' '.join(extra)

    t.to_csv(RESULTS / 'TableS2_absolute_growth.csv', index=False)
    to_markdown(
        t, RESULTS / 'TableS2_absolute_growth.md',
        'Supplementary Table S2. Absolute solo and joint growth rates '
        'underlying every reported boost factor.', note)
    print(f'  Table S2: {len(t)} rows')


def main():
    print('Building manuscript tables...')
    table5()
    table_s2()
    print(f'\nWritten to {RESULTS}')
    for f in sorted(RESULTS.glob('Table*')):
        print(f'  {f.name}')


if __name__ == '__main__':
    main()
