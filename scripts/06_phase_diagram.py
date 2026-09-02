"""2D phase diagram: oxygen x mucin availability.

Maps the ecological interaction class over a plane defined by oxygen uptake
and mucin-associated substrate availability, so that the two media reported
in the manuscript appear as two points within a continuous landscape rather
than as arbitrary choices.

Mucin axis
----------
Mucin availability is treated as a single ecological axis: the uptake bounds
of all mucin-associated substrates are scaled together by one factor, from 0
(fully non-mucosal, as in the OE medium) to 1 (CSOM level).

Important limitation to state in the manuscript
-----------------------------------------------
Only these two axes vary. All other components are held at the values of the
base medium, so the OE and CSOM reference points marked on the plane are
projections, not exact reproductions of the two published media.

Output: results/phase_diagram.csv
"""

import sys
import copy
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)

# ----------------------------------------------------------------- config
N_O2    = 21             # oxygen grid points
N_MUC   = 21             # mucin grid points
BASE    = 'CSOM'         # scaffold medium; mucin present here
LOW, HIGH = 0.90, 1.10

O2_KEY_LUMEN = None      # override if auto-detection is wrong
O2_KEY_AGORA = None

# Substrate name fragments treated as the mucin axis.
MUCIN_PATTERNS = ['acnam', 'neuac', 'neuraminate', 'sialic',
                  'fuc_l', 'fuc-l', 'fucose',
                  'acgam', 'glcnac', 'acetylglucosamine']

PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]
DIETS = {
    'CSOM': ('csom_lumen.pckl', 'csom_agora.pckl'),
    'OE':   ('oe_lumen.pckl',   'oe_agora.pckl'),
}


def load_diet(f):
    with open(DATA / f, 'rb') as fh:
        return pickle.load(fh)


def find_o2_key(diet, override=None):
    if override is not None:
        return override
    hits = sorted([k for k in diet
                   if 'o2' in k.lower() or 'oxygen' in k.lower()], key=len)
    if not hits:
        raise KeyError('no oxygen-like key; set O2_KEY_* manually')
    return hits[0]


def find_mucin_keys(diet):
    keys = [k for k in diet
            if any(p in k.lower() for p in MUCIN_PATTERNS)]
    return sorted(keys)


def classify(f, b):
    sgn = lambda v: '+' if v > HIGH else ('-' if v < LOW else 'o')
    sf, sb = sgn(f), sgn(b)
    table = {('+', '+'): 'Mutualism',    ('-', '-'): 'Competition',
             ('o', 'o'): 'Neutralism',   ('+', '-'): 'Parasitism',
             ('-', '+'): 'Parasitism',   ('+', 'o'): 'Commensalism',
             ('o', '+'): 'Commensalism', ('-', 'o'): 'Amensalism',
             ('o', '-'): 'Amensalism'}
    return table[(sf, sb)]


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def main():
    print('Loading CP iDC1003...')
    cp_raw = cobra.io.read_sbml_model(str(MODELS / 'cp' / 'iDC1003_processed.xml'))
    mapping = pd.read_csv(DATA / 'id_mapping_cp.csv', index_col=0)
    cp = cp_join_lumen_rxns(cp_raw.copy(), mapping)
    for leak in ['EX_CO2__extr', 'EX_Oxygen__extr']:
        try:
            cp.reactions.get_by_id(leak).bounds = (0, 0)
        except KeyError:
            pass
    cp.objective = 'e_Biomass__cyto'
    cp.solver = 'glpk'

    diets = {e: (load_diet(l), load_diet(a)) for e, (l, a) in DIETS.items()}
    base_lum, base_ago = diets[BASE]

    k_lum = find_o2_key(base_lum, O2_KEY_LUMEN)
    k_ago = find_o2_key(base_ago, O2_KEY_AGORA)
    muc_lum = find_mucin_keys(base_lum)
    muc_ago = find_mucin_keys(base_ago)
    print(f'O2 keys: lumen={k_lum!r} agora={k_ago!r}')
    print(f'Mucin keys (lumen, n={len(muc_lum)}): {muc_lum}')
    print(f'Mucin keys (agora, n={len(muc_ago)}): {muc_ago}')
    if not muc_lum:
        raise SystemExit('No mucin substrates found — extend MUCIN_PATTERNS.')

    o2_c, o2_o = base_lum[k_lum], diets['OE'][0][k_lum]
    lo, hi = sorted((o2_c, o2_o))
    o2_grid_l = np.linspace(lo, hi, N_O2)
    ratio = ((diets['OE'][1][k_ago] - base_ago[k_ago]) /
             (o2_o - o2_c)) if o2_o != o2_c else 1.0
    o2_grid_a = base_ago[k_ago] + (o2_grid_l - o2_c) * ratio
    muc_grid = np.linspace(0.0, 1.0, N_MUC)

    muc_base_l = {k: base_lum[k] for k in muc_lum}
    muc_base_a = {k: base_ago[k] for k in muc_ago}

    rows = []
    for p in PAIRS:
        print(f'\n[{p["bac_label"]}] building joint model once...')
        bac = cobra.io.read_sbml_model(
            str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
        bac_j = change_ids(bac.copy())
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
            bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                    base_lum, copy, 'agora')
        fu_b, bac_b = biomass_rxns(joint)
        cp_solo_model, bac_solo_model = cp.copy(), bac.copy()

        total = N_O2 * N_MUC
        for i, (ol, oa) in enumerate(zip(o2_grid_l, o2_grid_a)):
            for j, mf in enumerate(muc_grid):
                d_lum = dict(base_lum)
                d_ago = dict(base_ago)
                d_lum[k_lum] = float(ol)
                d_ago[k_ago] = float(oa)
                for k, v in muc_base_l.items():
                    d_lum[k] = float(v * mf)
                for k, v in muc_base_a.items():
                    d_ago[k] = float(v * mf)

                cp_solo = constrain_model(
                    cp_solo_model, d_lum, 'EX_u').optimize().objective_value or 0.0
                bac_solo = constrain_model(
                    bac_solo_model, d_ago, 'EX_').optimize().objective_value or 0.0

                constrain_model(joint, d_lum, 'EX_u')
                fu_b.bounds = (0.001, 1000)
                bac_b.bounds = (0.001, 1000)
                joint.objective = {fu_b: 1.0, bac_b: 1.0}
                sol = joint.optimize()

                if (sol.status != 'optimal' or cp_solo <= 1e-9
                        or bac_solo <= 1e-9):
                    fb = bb = np.nan
                    mode = 'infeasible'
                else:
                    fb = sol.fluxes[fu_b.id] / cp_solo
                    bb = sol.fluxes[bac_b.id] / bac_solo
                    mode = classify(fb, bb)

                rows.append({
                    'pair'      : f'CP+{p["bac_label"]}',
                    'o2_lumen'  : ol,
                    'mucin_frac': mf,
                    'cp_boost'  : fb,
                    'bac_boost' : bb,
                    'mode'      : mode,
                })
            done = (i + 1) * N_MUC
            print(f'  {done}/{total} points', end='\r')
        print()

    df = pd.DataFrame(rows)
    out = RESULTS / 'phase_diagram.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    # reference points for plotting
    ref = pd.DataFrame([
        {'env': 'CSOM', 'o2_lumen': o2_c, 'mucin_frac': 1.0},
        {'env': 'OE',   'o2_lumen': o2_o, 'mucin_frac': 0.0},
    ])
    ref.to_csv(RESULTS / 'phase_reference_points.csv', index=False)

    print('\n=== Mode composition ===')
    for pair, g in df.groupby('pair'):
        print(f'\n{pair}:')
        for m, n in g['mode'].value_counts().items():
            print(f'  {m:<14} {n:4d}  ({100*n/len(g):.1f}%)')


if __name__ == '__main__':
    main()
