"""Does the cross-feeding result depend on pinning both biomass reactions?

02_cross_feeding.py fixes both biomass fluxes at their free-optimum values
before running pFBA. Diagnostics in 16 showed that pinning fungal biomass
manufactures dependencies that vanish the moment the fungus is allowed to
grow even slightly slower, so the same constraint may also be shaping the
exchange table that the manuscript reports.

Three constraint regimes are compared, all at the same medium and with the
same joint objective used to define the free optimum:

    pinned    both biomass fluxes fixed at the free optimum   (as in 02)
    capped    both allowed anywhere in [0, free optimum]
    joint     no pinning; the 1:1 joint objective is optimised and pFBA is
              run on that solution

For each regime the fungus-to-bacterium and bacterium-to-fungus metabolite
sets are recorded, then compared. Metabolites appearing under all three are
robust; those appearing only under pinning are constraint artefacts and
should not be presented as predicted cross-feeding.

Output: results/crossfeeding_regimes.csv
        results/crossfeeding_regime_comparison.csv
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

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'

THRESHOLD = 0.01
ENV_DIET = {'CSOM': 'csom_lumen.pckl', 'OE': 'oe_lumen.pckl'}


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def classify_exchange(joint, sol):
    """Reproduce the direction assignment used in 02_cross_feeding.py."""
    fu_side, bac_side = {}, {}
    for r in joint.reactions:
        try:
            cat = categorize_joint_rxn(r)
        except Exception:
            continue
        if cat in ('fungal', 'BAC'):
            for m in r.metabolites:
                if m.compartment == 'u':
                    mid = m.id[:-2] if m.id.endswith('_u') else m.id
                    (fu_side if cat == 'fungal' else bac_side)[mid] = \
                        sol.fluxes[r.id]
                    break
    out = {}
    for met in set(fu_side) | set(bac_side):
        ff, bf = fu_side.get(met, 0.0), bac_side.get(met, 0.0)
        if abs(ff) < THRESHOLD and abs(bf) < THRESHOLD:
            continue
        if ff > THRESHOLD and bf < -THRESHOLD:
            d = 'fungus_to_bac'
        elif ff < -THRESHOLD and bf > THRESHOLD:
            d = 'bac_to_fungus'
        elif ff > THRESHOLD and bf > THRESHOLD:
            d = 'both_produce'
        elif ff < -THRESHOLD and bf < -THRESHOLD:
            d = 'both_consume'
        elif abs(bf) < THRESHOLD:
            d = 'fungus_only'
        elif abs(ff) < THRESHOLD:
            d = 'bac_only'
        else:
            d = 'unclassified'
        out[met] = (d, ff, bf, (abs(ff) + abs(bf)) / 2)
    return out


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

    pair_results = pd.read_csv(RESULTS / 'pair_results.csv')
    rows = []

    for _, base in pair_results.iterrows():
        env, pair = base['env'], base['pair']
        cp_free, bac_free = float(base['cp_free']), float(base['bac_free'])
        print(f'\n=== {env} / {pair} ===')

        with open(DATA / ENV_DIET[env], 'rb') as fh:
            diet_lum = pickle.load(fh)

        bac = cobra.io.read_sbml_model(
            str(MODELS / 'agora' / f'{base["bac_model"]}.xml'))
        bac_j = change_ids(bac.copy())
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
            bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                    diet_lum, copy, 'agora')
        joint = constrain_model(joint, diet_lum, 'EX_u')
        fu_b, bac_b = biomass_rxns(joint)

        for regime in ('pinned', 'capped', 'joint'):
            if regime == 'pinned':
                fu_b.bounds = (cp_free, cp_free)
                bac_b.bounds = (bac_free, bac_free)
                joint.objective = fu_b
            elif regime == 'capped':
                fu_b.bounds = (0.0, cp_free)
                bac_b.bounds = (0.0, bac_free)
                joint.objective = fu_b
            else:
                fu_b.bounds = (0.001, 1000.0)
                bac_b.bounds = (0.001, 1000.0)
                joint.objective = {fu_b: 1.0, bac_b: 1.0}

            try:
                sol = cobra.flux_analysis.pfba(joint)
            except Exception as e:
                print(f'  {regime:<7} pFBA failed: {e}')
                continue

            ex = classify_exchange(joint, sol)
            n_f2b = sum(1 for v in ex.values() if v[0] == 'fungus_to_bac')
            n_b2f = sum(1 for v in ex.values() if v[0] == 'bac_to_fungus')
            print(f'  {regime:<7} CP={sol.fluxes[fu_b.id]:.4f} '
                  f'bac={sol.fluxes[bac_b.id]:.4f}  '
                  f'f2b={n_f2b} b2f={n_b2f}')

            for met, (d, ff, bf, mag) in ex.items():
                rows.append({'env': env, 'pair': pair, 'regime': regime,
                             'metabolite': met, 'direction': d,
                             'fungus_flux': round(ff, 4),
                             'bac_flux': round(bf, 4),
                             'mag': round(mag, 4),
                             'cp_growth': sol.fluxes[fu_b.id],
                             'bac_growth': sol.fluxes[bac_b.id]})

    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / 'crossfeeding_regimes.csv', index=False)
    print(f'\nSaved: {RESULTS}/crossfeeding_regimes.csv  ({len(df)} rows)')

    # ------------------------------------------------ comparison
    comp = []
    directional = ('fungus_to_bac', 'bac_to_fungus')
    for (env, pair), g in df.groupby(['env', 'pair']):
        sets = {r: set(g[(g.regime == r) & (g.direction.isin(directional))]
                       .metabolite) for r in ('pinned', 'capped', 'joint')}
        allm = set().union(*sets.values())
        for m in sorted(allm):
            comp.append({
                'env': env, 'pair': pair, 'metabolite': m,
                'in_pinned': m in sets['pinned'],
                'in_capped': m in sets['capped'],
                'in_joint': m in sets['joint'],
                'robust': all(m in s for s in sets.values()),
                'pinned_only': (m in sets['pinned']
                                and m not in sets['capped']
                                and m not in sets['joint']),
            })
    cdf = pd.DataFrame(comp)
    cdf.to_csv(RESULTS / 'crossfeeding_regime_comparison.csv', index=False)

    print('\n=== Directional exchange sets ===')
    for (env, pair), g in cdf.groupby(['env', 'pair']):
        n_p = int(g.in_pinned.sum())
        n_c = int(g.in_capped.sum())
        n_j = int(g.in_joint.sum())
        n_r = int(g.robust.sum())
        n_o = int(g.pinned_only.sum())
        print(f'  {env:5s} {pair:6s}  pinned={n_p:3d} capped={n_c:3d} '
              f'joint={n_j:3d} | robust={n_r:3d} pinned-only={n_o:3d}')
        if n_o:
            print('      pinned-only: '
                  + ', '.join(g[g.pinned_only].metabolite.tolist()))

    print('\n=== Verdict ===')
    tot_o = int(cdf.pinned_only.sum())
    tot_p = int(cdf.in_pinned.sum())
    if tot_p:
        print(f'  {tot_o}/{tot_p} directional exchanges reported under pinning '
              f'({100*tot_o/tot_p:.0f}%) do not survive relaxation.')
    if tot_o == 0:
        print('  The published exchange table is not an artefact of pinning.')
    else:
        print('  Table 4 should be rebuilt from the robust set.')


if __name__ == '__main__':
    main()
