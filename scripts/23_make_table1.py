"""Build Table 1 from the medium definitions themselves.

The previous Table 1 listed 14 components taken from a console-output
list in build_oe_media.py, which was written to spot-check a few values
rather than to define a table. That selection turned out to be
inconsistent: acetate differs between the two media but was absent,
ferric dicitrate was absent although iron is central to the results, and
formate appeared although it was not in the original list.

This script instead reads the two media and emits every exchange bound
that differs, so the table cannot drift away from the code. Component
groups match those used in the attribution analysis (Table 4), so the
group sizes quoted in the Table 4 footnote are derived from the same
source.

Output: results/Table1_media_components.{csv,md}
"""

import pickle
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)

# display names for the metabolites that differ; anything missing falls
# back to its exchange id so nothing is silently dropped
NAMES = {
    'glc_D': 'D-glucose', 'gal': 'D-galactose', 'man': 'D-mannose',
    'gam': 'D-glucosamine',
    'o2': 'Oxygen', 'nh4': 'Ammonium',
    'lac_L': 'L-lactate', 'ac': 'Acetate', 'for': 'Formate',
    'succ': 'Succinate', 'ppa': 'Propionate',
    'acgam': 'N-acetylglucosamine', 'acnam': 'N-acetylneuraminate',
    'fuc_L': 'L-fucose',
    'fe2': 'Ferrous iron', 'fe3': 'Ferric iron',
    'fe3dcit': 'Ferric dicitrate',
    'ala_L': 'L-alanine', 'arg_L': 'L-arginine', 'asn_L': 'L-asparagine',
    'asp_L': 'L-aspartate', 'cys_L': 'L-cysteine', 'gln_L': 'L-glutamine',
    'glu_L': 'L-glutamate', 'gly': 'Glycine', 'his_L': 'L-histidine',
    'ile_L': 'L-isoleucine', 'leu_L': 'L-leucine', 'lys_L': 'L-lysine',
    'met_L': 'L-methionine', 'phe_L': 'L-phenylalanine',
    'pro_L': 'L-proline', 'ser_L': 'L-serine', 'thr_L': 'L-threonine',
    'trp_L': 'L-tryptophan', 'tyr_L': 'L-tyrosine', 'val_L': 'L-valine',
}

GROUPS = {
    'Oxygen'                     : ['o2'],
    'Nitrogen source'            : ['nh4'],
    'Iron source'                : ['fe2', 'fe3', 'fe3dcit'],
    'Mucin-associated substrate' : ['acgam', 'acnam', 'fuc_L'],
    'Fermentation product'       : ['lac_L', 'ac', 'for', 'succ', 'ppa'],
    'Sugar'                      : ['glc_D', 'gal', 'man', 'gam'],
    'Host-derived amino acid'    : ['ala_L', 'arg_L', 'asn_L', 'asp_L',
                                    'cys_L', 'gln_L', 'glu_L', 'gly',
                                    'his_L', 'ile_L', 'leu_L', 'lys_L',
                                    'met_L', 'phe_L', 'pro_L', 'ser_L',
                                    'thr_L', 'trp_L', 'tyr_L', 'val_L'],
}
ORDER = ['Oxygen', 'Sugar', 'Fermentation product',
         'Mucin-associated substrate', 'Host-derived amino acid',
         'Nitrogen source', 'Iron source']


def load(name):
    with open(DATA / name, 'rb') as f:
        return pickle.load(f)


def fmt(v):
    if v == 0:
        return '0'
    s = f'{v:.4f}'.rstrip('0').rstrip('.')
    return s


def main():
    csom, oe = load('csom_lumen.pckl'), load('oe_lumen.pckl')

    met_group = {m: g for g, mets in GROUPS.items() for m in mets}

    rows, ungrouped = [], []
    for key in set(csom) | set(oe):
        if not key.startswith('EX_u_'):
            continue
        met = key[5:]
        c, o = float(csom.get(key, 0.0)), float(oe.get(key, 0.0))
        if abs(c - o) < 1e-12:
            continue
        g = met_group.get(met)
        if g is None:
            ungrouped.append(met)
            g = 'Other'
        rows.append({'Category': g, 'Metabolite': NAMES.get(met, met),
                     'Exchange ID': met, 'CSOM': fmt(c), 'OE': fmt(o),
                     '_sort': ORDER.index(g) if g in ORDER else 99})

    df = (pd.DataFrame(rows)
          .sort_values(['_sort', 'Metabolite'])
          .drop(columns='_sort')
          .reset_index(drop=True))

    # blank repeated category labels, as in the published table
    shown, cat = [], None
    for c in df.Category:
        shown.append('' if c == cat else c)
        cat = c
    df['Category'] = shown

    df.to_csv(RESULTS / 'Table1_media_components.csv', index=False)

    counts = (pd.DataFrame(rows).Category.value_counts()
              .reindex(ORDER).dropna().astype(int))
    note = ('Values are lower bounds for metabolite uptake '
            '(mmol gDW-1 h-1); negative values indicate permitted uptake '
            'and 0 indicates that uptake was not permitted. All components '
            f'whose bounds differed between the two media are listed '
            f'({len(df)} in total); components with identical bounds are '
            f'not shown. The media comprised {len(csom)} (CSOM) and '
            f'{len(oe)} (OE) exchange reactions in total. CSOM, chronic '
            'suppurative otitis media; OE, otitis externa.')

    lines = ['**Table 1. Components that differed between the two simulated '
             'environmental media**', '']
    lines.append('| ' + ' | '.join(df.columns) + ' |')
    lines.append('|' + '|'.join(' --- ' for _ in df.columns) + '|')
    for _, r in df.iterrows():
        lines.append('| ' + ' | '.join(str(v) for v in r.values) + ' |')
    lines += ['', note]
    (RESULTS / 'Table1_media_components.md').write_text('\n'.join(lines),
                                                        encoding='utf-8')

    print(f'Rows written: {len(df)}')
    print('\nGroup sizes (for the Table 4 footnote):')
    for g, n in counts.items():
        print(f'  {g:<28} {n}')
    if ungrouped:
        print('\n! not assigned to a group:', sorted(set(ungrouped)))
    print(f'\nSaved: {RESULTS}/Table1_media_components.csv and .md')


if __name__ == '__main__':
    main()
