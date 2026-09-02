"""Which cross-feeding links actually carry the bacterial dependency?

Three independent ways of suppressing the fungus as a whole (08, 11, 14) all
gave the same answer: bacterial growth is unaffected until fungal metabolism
is almost entirely gone. That points to a dependency resting on a small
number of exchanged metabolites rather than on the overall scale of fungal
activity.

This script tests that directly. For each fungus-to-bacterium metabolite
identified in 02_cross_feeding.py, secretion of that single metabolite into
the shared lumen is blocked and the bacterium is re-maximised.

    single    one link blocked at a time
    greedy    links added cumulatively, each step choosing the one that
              reduces bacterial growth most, until the fungicidal level is
              reached or MAX_GREEDY links are blocked

Blocking sets the upper bound of the fungal lumen reaction to zero, which
removes secretion while leaving uptake intact. That is the in silico
analogue of inhibiting export or synthesis of one compound, not of killing
the fungus.

Output: results/link_knockout.csv
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
from cp_helpers import (cp_join_lumen_rxns, make_joint_model_cp,
                        categorize_joint_rxn)

try:
    from metabolite_names import label
except ImportError:
    def label(m):
        return m

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'

MAX_GREEDY = 8
ENVS = {
    'CSOM': ('csom_lumen.pckl', 'csom_agora.pckl'),
    'OE':   ('oe_lumen.pckl',   'oe_agora.pckl'),
}
PAIRS = {'CP+PA': 'Pseudomonas_aeruginosa_NCGM2_S1',
         'CP+SA': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'}


def load(f):
    with open(DATA / f, 'rb') as fh:
        return pickle.load(fh)


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def fungal_lumen_map(joint):
    """agora metabolite id -> fungal-side lumen reaction."""
    out = {}
    for r in joint.reactions:
        try:
            if categorize_joint_rxn(r) != 'fungal':
                continue
        except Exception:
            continue
        for m in r.metabolites:
            if m.compartment == 'u':
                mid = m.id[:-2] if m.id.endswith('_u') else m.id
                out[mid] = r
                break
    return out


def max_bac(joint, fu_b, bac_b, cp_level):
    fu_b.bounds = (cp_level, cp_level)
    bac_b.bounds = (0.0, 1000.0)
    joint.objective = bac_b
    sol = joint.optimize()
    return sol.fluxes[bac_b.id] if sol.status == 'optimal' else np.nan


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

    cf = pd.read_csv(RESULTS / 'cross_feeding.csv')
    pair_results = pd.read_csv(RESULTS / 'pair_results.csv')
    rows = []

    for env, (lf, af) in ENVS.items():
        diet_lum = load(lf)
        print(f'\n=== {env} ===')

        for label_pair, model_name in PAIRS.items():
            base = pair_results[(pair_results.env == env) &
                                (pair_results.pair == label_pair)]
            if base.empty:
                continue
            cp_free = base.iloc[0]['cp_free']

            print(f'  [{label_pair}] assembling joint model...')
            bac = cobra.io.read_sbml_model(
                str(MODELS / 'agora' / f'{model_name}.xml'))
            bac_j = change_ids(bac.copy())
            bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
                bac_j, cp, Metabolite, Reaction, 'agora')
            joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                        diet_lum, copy, 'agora')
            joint = constrain_model(joint, diet_lum, 'EX_u')
            fu_b, bac_b = biomass_rxns(joint)
            fmap = fungal_lumen_map(joint)

            ref = max_bac(joint, fu_b, bac_b, cp_free)

            # full elimination, for scale
            saved_all = {r.id: r.bounds for r in fmap.values()}
            for r in fmap.values():
                r.bounds = (0.0, 0.0)
            cidal = max_bac(joint, fu_b, bac_b, 0.0)
            for r in fmap.values():
                r.bounds = saved_all[r.id]

            print(f'    reference {ref:.4f} | full elimination {cidal:.4f} '
                  f'({100*cidal/ref:.1f}%)')

            mets = (cf[(cf.env == env) & (cf.pair == label_pair) &
                       (cf.direction == 'fungus_to_bac')]
                    .sort_values('mag', ascending=False)['metabolite'].tolist())
            mets = [m for m in mets if m in fmap]
            print(f'    fungus-to-bacterium links tested: {len(mets)}')

            # ------------------------------------------------ single
            singles = {}
            for m in mets:
                r = fmap[m]
                old = r.bounds
                r.bounds = (old[0], 0.0)          # block secretion only
                v = max_bac(joint, fu_b, bac_b, cp_free)
                r.bounds = old
                singles[m] = v
                rows.append({
                    'env': env, 'pair': label_pair, 'analysis': 'single',
                    'n_blocked': 1, 'metabolite': label(m),
                    'blocked_set': label(m),
                    'bac_growth': v,
                    'bac_rel_to_coinfection': v / ref if ref else np.nan,
                })

            top = sorted(singles.items(), key=lambda kv: kv[1])[:5]
            for m, v in top:
                print(f'      {label(m):<24} {v:.4f} ({100*v/ref:.1f}%)')

            # ------------------------------------------------ greedy
            blocked, saved = [], {}
            remaining = list(mets)
            for step in range(1, MAX_GREEDY + 1):
                best, best_v = None, np.inf
                for m in remaining:
                    r = fmap[m]
                    old = r.bounds
                    r.bounds = (old[0], 0.0)
                    v = max_bac(joint, fu_b, bac_b, cp_free)
                    r.bounds = old
                    if v < best_v:
                        best, best_v = m, v
                if best is None:
                    break
                r = fmap[best]
                saved[best] = r.bounds
                r.bounds = (r.bounds[0], 0.0)
                blocked.append(best)
                remaining.remove(best)
                rows.append({
                    'env': env, 'pair': label_pair, 'analysis': 'greedy',
                    'n_blocked': step, 'metabolite': label(best),
                    'blocked_set': ' + '.join(label(x) for x in blocked),
                    'bac_growth': best_v,
                    'bac_rel_to_coinfection': best_v / ref if ref else np.nan,
                })
                print(f'      greedy {step}: +{label(best):<20} '
                      f'{best_v:.4f} ({100*best_v/ref:.1f}%)')
                if best_v <= cidal * 1.01:
                    print('      reached full-elimination level')
                    break
            for m, b in saved.items():
                fmap[m].bounds = b

            rows.append({
                'env': env, 'pair': label_pair, 'analysis': 'full_elimination',
                'n_blocked': len(fmap), 'metabolite': 'all fungal lumen',
                'blocked_set': 'all', 'bac_growth': cidal,
                'bac_rel_to_coinfection': cidal / ref if ref else np.nan,
            })

    df = pd.DataFrame(rows)
    out = RESULTS / 'link_knockout.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    print('\n=== Minimal link sets ===')
    for (env, pair), g in df.groupby(['env', 'pair']):
        full = g[g.analysis == 'full_elimination'].bac_rel_to_coinfection
        tgt = full.iloc[0] if len(full) else np.nan
        gr = g[g.analysis == 'greedy'].sort_values('n_blocked')
        hit = gr[gr.bac_rel_to_coinfection <= tgt * 1.05] if len(gr) else gr
        if len(hit):
            h = hit.iloc[0]
            print(f'  {env:5s} {pair:6s}  {int(h.n_blocked)} link(s) reach '
                  f'{100*h.bac_rel_to_coinfection:.1f}% '
                  f'(elimination {100*tgt:.1f}%): {h.blocked_set}')
        else:
            print(f'  {env:5s} {pair:6s}  no link set within '
                  f'{MAX_GREEDY} reached the elimination level')


if __name__ == '__main__':
    main()
