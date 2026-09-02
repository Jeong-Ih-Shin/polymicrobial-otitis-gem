"""Why some link knockouts were infeasible, and whether the result survives.

In 15 the fungal biomass was pinned at its free-optimum value while a link
was blocked. Blocking Cys-Gly secretion then returned 'infeasible' for two
conditions, which does not mean the intervention fails; it means the fungus
cannot sustain that exact growth rate without secreting the compound. The
pinned constraint is also questionable on its own terms: an intervention
that removes a metabolic route would be expected to change fungal growth,
not to hold it fixed.

Each single-link knockout is therefore repeated under three regimes:

    fixed     fungal biomass = cp_free                (as in 15)
    capped    fungal biomass in [0, cp_free]          fungus may slow down
    free      fungal biomass in [0, 1000]             fungus re-optimises

A result that holds across all three is not an artefact of the constraint.
The script also reports, for every link, the highest fungal growth still
attainable once that link is blocked, which is what makes the infeasibility
interpretable.

Output: results/knockout_regimes.csv
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

RESULTS = ROOT / 'results'
MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'

ENVS = {'CSOM': 'csom_lumen.pckl', 'OE': 'oe_lumen.pckl'}
PAIRS = {'CP+PA': 'Pseudomonas_aeruginosa_NCGM2_S1',
         'CP+SA': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'}


def clean(name):
    """Strip LaTeX markup that metabolite_names.label() adds for figures."""
    return (str(name).replace('$_', '').replace('$', '')
            .replace('^', '').replace('{', '').replace('}', ''))


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def fungal_lumen_map(joint):
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


def solve(joint, obj_rxn):
    joint.objective = obj_rxn
    sol = joint.optimize()
    return (sol.fluxes[obj_rxn.id], sol.status) if sol.status == 'optimal' \
        else (np.nan, sol.status)


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

    for env, lf in ENVS.items():
        with open(DATA / lf, 'rb') as fh:
            diet_lum = pickle.load(fh)
        print(f'\n=== {env} ===')

        for pair, model_name in PAIRS.items():
            base = pair_results[(pair_results.env == env) &
                                (pair_results.pair == pair)]
            if base.empty:
                continue
            cp_free = float(base.iloc[0]['cp_free'])

            print(f'  [{pair}] assembling joint model...')
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

            # untreated reference under each regime
            refs = {}
            for regime, fb in [('fixed', (cp_free, cp_free)),
                               ('capped', (0.0, cp_free)),
                               ('free', (0.0, 1000.0))]:
                fu_b.bounds = fb
                bac_b.bounds = (0.0, 1000.0)
                refs[regime] = solve(joint, bac_b)[0]
            print(f'    reference  fixed={refs["fixed"]:.4f}  '
                  f'capped={refs["capped"]:.4f}  free={refs["free"]:.4f}')

            mets = (cf[(cf.env == env) & (cf.pair == pair) &
                       (cf.direction == 'fungus_to_bac')]
                    .sort_values('mag', ascending=False)['metabolite'].tolist())
            mets = [m for m in mets if m in fmap]

            for m in mets:
                r = fmap[m]
                old = r.bounds
                r.bounds = (old[0], 0.0)

                # how much can the fungus still grow without this outlet?
                fu_b.bounds = (0.0, 1000.0)
                cp_max, _ = solve(joint, fu_b)

                rec = {'env': env, 'pair': pair, 'metabolite': clean(m),
                       'raw_id': m, 'cp_free_reference': cp_free,
                       'cp_max_after_block': cp_max,
                       'blocks_fungal_optimum': (np.isnan(cp_max)
                                                 or cp_max < cp_free - 1e-6)}

                for regime, fb in [('fixed', (cp_free, cp_free)),
                                   ('capped', (0.0, cp_free)),
                                   ('free', (0.0, 1000.0))]:
                    fu_b.bounds = fb
                    bac_b.bounds = (0.0, 1000.0)
                    v, st = solve(joint, bac_b)
                    rec[f'bac_{regime}'] = v
                    rec[f'status_{regime}'] = st
                    ref = refs[regime]
                    rec[f'rel_{regime}'] = v / ref if ref else np.nan

                r.bounds = old
                rows.append(rec)

            g = pd.DataFrame([x for x in rows
                              if x['env'] == env and x['pair'] == pair])
            g = g.sort_values('rel_free')
            print('    top links by effect under the free regime:')
            for _, x in g.head(4).iterrows():
                flag = ' (blocks fungal optimum)' if x.blocks_fungal_optimum else ''
                print(f'      {x.metabolite:<22} '
                      f'fixed={100*x.rel_fixed:6.1f}%  '
                      f'capped={100*x.rel_capped:6.1f}%  '
                      f'free={100*x.rel_free:6.1f}%{flag}')

    df = pd.DataFrame(rows)
    out = RESULTS / 'knockout_regimes.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    print('\n=== Links that were infeasible under the pinned regime ===')
    bad = df[df.status_fixed != 'optimal']
    if bad.empty:
        print('  none')
    for _, x in bad.iterrows():
        print(f'  {x.env:5s} {x.pair:6s} {x.metabolite:<22} '
              f'fungal growth attainable after block: {x.cp_max_after_block:.4f} '
              f'(reference {x.cp_free_reference:.4f}) -> '
              f'capped={100*x.rel_capped:.1f}%  free={100*x.rel_free:.1f}%')

    print('\n=== Robustness of the leading link across regimes ===')
    for (env, pair), g in df.groupby(['env', 'pair']):
        g = g.sort_values('rel_free')
        if g.empty:
            continue
        x = g.iloc[0]
        agree = all(np.isfinite([x.rel_capped, x.rel_free])) and \
            abs(x.rel_capped - x.rel_free) < 0.05
        print(f'  {env:5s} {pair:6s}  best={x.metabolite:<20} '
              f'capped={100*x.rel_capped:5.1f}%  free={100*x.rel_free:5.1f}%  '
              f'{"consistent" if agree else "regime-dependent"}')


if __name__ == '__main__':
    main()
