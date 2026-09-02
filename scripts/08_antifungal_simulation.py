"""In silico antifungal intervention: bacterial response to fungal removal.

Asks what the joint model predicts for the bacterial partner when the fungus
is suppressed or eliminated, i.e. the metabolic consequence of adding
antifungal therapy to a bacterial-fungal co-infection.

Three scenarios, all evaluated inside the *joint* model so that lumen
topology and medium are identical across conditions. Comparing against the
bacterium's solo model instead would confound fungal removal with the
difference between the agora and lumen diets.

  1. dose_response  fungal biomass constrained to a fraction f of its
                    free-optimum value, f from 1.0 down to 0.0
  2. fungistatic    fungal biomass fixed at 0; fungal metabolism and lumen
                    exchange still permitted (growth arrest without lysis)
  3. fungicidal     fungal biomass fixed at 0 and every fungal-side lumen
                    exchange blocked (metabolic elimination)

In each case the bacterial biomass reaction is maximised alone.

Output: results/antifungal_simulation.csv
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
RESULTS.mkdir(exist_ok=True)

N_DOSE = 21                      # points on the suppression axis
ENVS = {
    'CSOM': ('csom_lumen.pckl', 'csom_agora.pckl'),
    'OE':   ('oe_lumen.pckl',   'oe_agora.pckl'),
}
PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]


def load(f):
    with open(DATA / f, 'rb') as fh:
        return pickle.load(fh)


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def fungal_lumen_rxns(joint):
    """Fungal-side lumen exchange reactions, via the project's classifier."""
    out = []
    for r in joint.reactions:
        try:
            if categorize_joint_rxn(r) == 'fungal':
                out.append(r)
        except Exception:
            continue
    return out


def maximise_bacterium(joint, fu_b, bac_b, fungal_level):
    """Fix fungal biomass, maximise the bacterium, return its flux."""
    fu_b.bounds = (fungal_level, fungal_level)
    bac_b.bounds = (0.0, 1000.0)
    joint.objective = bac_b
    sol = joint.optimize()
    if sol.status != 'optimal':
        return np.nan
    return sol.fluxes[bac_b.id]


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

    for env, (lum_f, ago_f) in ENVS.items():
        diet_lum, diet_ago = load(lum_f), load(ago_f)
        print(f'\n=== {env} ===')

        for p in PAIRS:
            label = f'CP+{p["bac_label"]}'
            base = pair_results[(pair_results.env == env) &
                                (pair_results.pair == label)]
            if base.empty:
                print(f'  {label}: absent from pair_results.csv — skip')
                continue
            base = base.iloc[0]

            print(f'  [{label}] assembling joint model...')
            bac = cobra.io.read_sbml_model(
                str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
            bac_j = change_ids(bac.copy())
            bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
                bac_j, cp, Metabolite, Reaction, 'agora')
            joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                        diet_lum, copy, 'agora')
            joint = constrain_model(joint, diet_lum, 'EX_u')
            fu_b, bac_b = biomass_rxns(joint)

            cp_free, bac_free = base['cp_free'], base['bac_free']

            # reference: bacterium maximised with the fungus at free optimum
            bac_ref = maximise_bacterium(joint, fu_b, bac_b, cp_free)

            # ---------------------------------------- 1. dose response
            for frac in np.linspace(1.0, 0.0, N_DOSE):
                val = maximise_bacterium(joint, fu_b, bac_b, cp_free * frac)
                rows.append({
                    'env': env, 'pair': label, 'scenario': 'dose_response',
                    'fungal_fraction': frac,
                    'bac_growth': val,
                    'bac_rel_to_coinfection': val / bac_ref if bac_ref else np.nan,
                    'bac_rel_to_joint_free': val / bac_free if bac_free else np.nan,
                })

            # ---------------------------------------- 2. fungistatic
            static = maximise_bacterium(joint, fu_b, bac_b, 0.0)

            # ---------------------------------------- 3. fungicidal
            blocked = fungal_lumen_rxns(joint)
            saved = {r.id: r.bounds for r in blocked}
            for r in blocked:
                r.bounds = (0.0, 0.0)
            cidal = maximise_bacterium(joint, fu_b, bac_b, 0.0)
            for r in blocked:
                r.bounds = saved[r.id]

            for name, val in [('fungistatic', static), ('fungicidal', cidal)]:
                rows.append({
                    'env': env, 'pair': label, 'scenario': name,
                    'fungal_fraction': 0.0,
                    'bac_growth': val,
                    'bac_rel_to_coinfection': val / bac_ref if bac_ref else np.nan,
                    'bac_rel_to_joint_free': val / bac_free if bac_free else np.nan,
                })

            print(f'    co-infection ref : {bac_ref:.4f}')
            print(f'    fungistatic      : {static:.4f} '
                  f'({100*static/bac_ref:.1f}% of co-infection)'
                  if bac_ref else '')
            print(f'    fungicidal       : {cidal:.4f} '
                  f'({100*cidal/bac_ref:.1f}% of co-infection)'
                  if bac_ref else '')
            print(f'    n fungal lumen reactions blocked: {len(blocked)}')

    df = pd.DataFrame(rows)
    out = RESULTS / 'antifungal_simulation.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    print('\n=== Predicted bacterial growth after fungal removal ===')
    sub = df[df.scenario.isin(['fungistatic', 'fungicidal'])]
    for (env, pair), g in sub.groupby(['env', 'pair']):
        line = '  '.join(f'{r.scenario}={100*r.bac_rel_to_coinfection:.1f}%'
                         for _, r in g.iterrows())
        print(f'  {env:5s} {pair:6s}  {line}')


if __name__ == '__main__':
    main()
