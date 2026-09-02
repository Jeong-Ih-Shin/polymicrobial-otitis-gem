"""1D oxygen-gradient scan of bacterial-fungal interaction mode.

Sweeps the oxygen uptake bound across a continuum spanning the CSOM and OE
values, recomputing solo and joint growth at every point, and records the
boost factors and ecological interaction class along the gradient.

Rationale: the manuscript reports interaction mode at two discrete
environments. If mode switching is driven by oxygen availability, the
transition should be locatable on a continuous axis rather than asserted
from two points.

Design notes
------------
* CP solo is constrained with the *lumen* diet, the bacterium solo with the
  *agora* diet (as in 01_simulate.py). Both must be perturbed together or
  the boost factor becomes internally inconsistent. This is enforced below.
* Joint models are expensive to assemble, so each pair is built once and
  only exchange bounds are re-applied per grid point.

Output: results/oxygen_scan.csv
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

try:
    from functions import interaction_type
    HAVE_ITYPE = True
except ImportError:                                   # optional cross-check
    HAVE_ITYPE = False

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)

# ----------------------------------------------------------------- config
N_POINTS   = 41          # grid resolution along the oxygen axis
PAD        = 0.15        # extend range 15% beyond the CSOM-OE span
LOW, HIGH  = 0.90, 1.10  # manuscript thresholds on the raw boost factor

# Set these manually if auto-detection picks the wrong key.
O2_KEY_LUMEN = None
O2_KEY_AGORA = None

PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]

DIETS = {
    'CSOM': ('csom_lumen.pckl', 'csom_agora.pckl'),
    'OE':   ('oe_lumen.pckl',   'oe_agora.pckl'),
}


# ------------------------------------------------------------- utilities
def load_diet(fname):
    with open(DATA / fname, 'rb') as f:
        return pickle.load(f)


def find_o2_key(diet, override=None):
    """Locate the oxygen exchange key in a diet dictionary."""
    if override is not None:
        if override not in diet:
            raise KeyError(f'{override!r} not present in diet')
        return override
    hits = [k for k in diet
            if 'o2' in k.lower() or 'oxygen' in k.lower()]
    # prefer the shortest match: 'EX_o2(e)' over 'EX_o2s(e)' etc.
    hits.sort(key=len)
    if not hits:
        raise KeyError('no oxygen-like key found; set O2_KEY_* manually')
    if len(hits) > 1:
        print(f'  ! multiple oxygen-like keys {hits}; using {hits[0]!r}')
    return hits[0]


def classify(f, b):
    """Ecological class from raw boost factors, per the manuscript."""
    sgn = lambda v: '+' if v > HIGH else ('-' if v < LOW else 'o')
    sf, sb = sgn(f), sgn(b)
    table = {('+', '+'): 'Mutualism',    ('-', '-'): 'Competition',
             ('o', 'o'): 'Neutralism',   ('+', '-'): 'Parasitism',
             ('-', '+'): 'Parasitism',   ('+', 'o'): 'Commensalism',
             ('o', '+'): 'Commensalism', ('-', 'o'): 'Amensalism',
             ('o', '-'): 'Amensalism'}
    return table[(sf, sb)], sf, sb


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def solve_joint(joint, fu_b, bac_b):
    fu_b.lower_bound, fu_b.upper_bound = 0.001, 1000
    bac_b.lower_bound, bac_b.upper_bound = 0.001, 1000
    joint.objective = {fu_b: 1.0, bac_b: 1.0}
    sol = joint.optimize()
    if sol.status != 'optimal':
        return None, None
    return sol.fluxes[fu_b.id], sol.fluxes[bac_b.id]


# ------------------------------------------------------------------ main
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

    diets = {env: (load_diet(l), load_diet(a)) for env, (l, a) in DIETS.items()}

    k_lum = find_o2_key(diets['CSOM'][0], O2_KEY_LUMEN)
    k_ago = find_o2_key(diets['CSOM'][1], O2_KEY_AGORA)
    print(f'Oxygen key (lumen): {k_lum!r}\nOxygen key (agora): {k_ago!r}')

    o2_csom_l, o2_oe_l = diets['CSOM'][0][k_lum], diets['OE'][0][k_lum]
    o2_csom_a, o2_oe_a = diets['CSOM'][1][k_ago], diets['OE'][1][k_ago]
    print(f'  lumen bounds  CSOM={o2_csom_l}  OE={o2_oe_l}')
    print(f'  agora bounds  CSOM={o2_csom_a}  OE={o2_oe_a}')

    # Grid spans CSOM->OE with padding at both ends. Bounds are uptake
    # lower bounds (negative), so ordering is handled by min/max.
    lo, hi = sorted((o2_csom_l, o2_oe_l))
    span = hi - lo
    grid_l = np.linspace(lo - PAD * span, hi + PAD * span, N_POINTS)
    # scale the agora bound proportionally so the two diets stay coupled
    ratio = (o2_oe_a - o2_csom_a) / (o2_oe_l - o2_csom_l) if o2_oe_l != o2_csom_l else 1.0
    grid_a = o2_csom_a + (grid_l - o2_csom_l) * ratio

    rows = []
    for base_env in ['CSOM', 'OE']:
        diet_lum_base, diet_ago_base = diets[base_env]
        print(f'\n=== Base medium: {base_env} (all components except O2 fixed) ===')

        for p in PAIRS:
            print(f'  [{p["bac_label"]}] building joint model once...')
            bac = cobra.io.read_sbml_model(
                str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
            bac_j = change_ids(bac.copy())
            bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
                bac_j, cp, Metabolite, Reaction, 'agora')
            joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                        diet_lum_base, copy, 'agora')
            fu_b, bac_b = biomass_rxns(joint)

            cp_solo_model = cp.copy()
            bac_solo_model = bac.copy()

            for i, (ol, oa) in enumerate(zip(grid_l, grid_a)):
                d_lum = dict(diet_lum_base); d_lum[k_lum] = float(ol)
                d_ago = dict(diet_ago_base); d_ago[k_ago] = float(oa)

                cp_solo = constrain_model(
                    cp_solo_model, d_lum, 'EX_u').optimize().objective_value or 0.0
                bac_solo = constrain_model(
                    bac_solo_model, d_ago, 'EX_').optimize().objective_value or 0.0

                constrain_model(joint, d_lum, 'EX_u')
                fu_free, bac_free = solve_joint(joint, fu_b, bac_b)

                if fu_free is None or cp_solo <= 1e-9 or bac_solo <= 1e-9:
                    mode, sf, sb = 'infeasible', '.', '.'
                    fb = bb = np.nan
                else:
                    fb = fu_free / cp_solo
                    bb = bac_free / bac_solo
                    mode, sf, sb = classify(fb, bb)

                row = {
                    'base_env'  : base_env,
                    'pair'      : f'CP+{p["bac_label"]}',
                    'o2_lumen'  : ol,
                    'o2_agora'  : oa,
                    'cp_solo'   : cp_solo,
                    'bac_solo'  : bac_solo,
                    'cp_free'   : fu_free,
                    'bac_free'  : bac_free,
                    'cp_boost'  : fb,
                    'bac_boost' : bb,
                    'cp_sign'   : sf,
                    'bac_sign'  : sb,
                    'mode'      : mode,
                }

                # cross-check against the project's own classifier
                if HAVE_ITYPE and np.isfinite(fb) and fb > 0 and bb > 0:
                    import math
                    try:
                        alt, _ = interaction_type(math.log2(fb), math.log2(bb), math)
                        row['mode_interaction_type'] = alt
                        row['mode_agrees'] = alt.split(' [')[0] == mode
                    except Exception:
                        pass

                rows.append(row)
                if i % 10 == 0:
                    print(f'    o2={ol:8.3f}  CP={fb:6.3f}x  '
                          f'bac={bb:6.3f}x  {mode}')

    df = pd.DataFrame(rows)
    out = RESULTS / 'oxygen_scan.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    if 'mode_agrees' in df.columns:
        bad = df[df.mode_agrees == False]
        if len(bad):
            print(f'\n!! {len(bad)} grid points where the manuscript thresholds '
                  f'and interaction_type() disagree — check functions.py')
        else:
            print('\nClassification consistent with interaction_type() '
                  'at all grid points.')

    # report transition points
    print('\n=== Mode transitions ===')
    for (env, pair), g in df.groupby(['base_env', 'pair']):
        g = g.sort_values('o2_lumen')
        prev = None
        for _, r in g.iterrows():
            if prev is not None and r['mode'] != prev:
                print(f'  {env} / {pair}: {prev} -> {r["mode"]} '
                      f'near o2={r["o2_lumen"]:.3f}')
            prev = r['mode']
        if prev is not None and g['mode'].nunique() == 1:
            print(f'  {env} / {pair}: no transition ({prev} throughout)')


if __name__ == '__main__':
    main()
