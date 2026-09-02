"""Graded antifungal effect via proportional scaling of fungal capacity.

Replaces the fungistatic scenario in 08, which was not a valid model of
growth arrest: fixing fungal biomass at zero leaves the fungal network fully
active and in fact frees the carbon otherwise spent on biomass, so the
bacterium loses nothing. That scenario returned exactly 100% in every
condition, which is a structural artefact rather than a result.

Here the flux capacity of every fungal reaction is scaled by a factor s,
from 1.0 (untreated co-infection) down to 0.0 (complete metabolic
elimination, equivalent to the fungicidal scenario). Intermediate values
represent partial suppression of fungal metabolic activity, which is the
closer analogue of a fungistatic agent.

Reaction partition inside the joint model
-----------------------------------------
    ids starting 'model2_'  bacterial
    ids starting 'EX_u_'    system boundary (the medium)
    everything else         fungal

Output: results/fungistatic_scaling.csv
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

N_STEPS = 21
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


def fungal_reactions(joint):
    return [r for r in joint.reactions
            if not r.id.startswith('model2_') and not r.id.startswith('EX_u_')]


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

    rows = []
    for env, (lf, af) in ENVS.items():
        diet_lum, diet_ago = load(lf), load(af)
        print(f'\n=== {env} ===')

        for p in PAIRS:
            label = f'CP+{p["bac_label"]}'
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

            fung = fungal_reactions(joint)
            base_bounds = {r.id: r.bounds for r in fung}
            print(f'    fungal reactions scaled: {len(fung)}')

            ref = None
            for s in np.linspace(1.0, 0.0, N_STEPS):
                for r in fung:
                    lb, ub = base_bounds[r.id]
                    r.bounds = (lb * s, ub * s)
                bac_b.bounds = (0.0, 1000.0)
                joint.objective = bac_b
                sol = joint.optimize()

                if sol.status != 'optimal':
                    bacg, fung_g = np.nan, np.nan
                else:
                    bacg = sol.fluxes[bac_b.id]
                    fung_g = sol.fluxes[fu_b.id]
                if ref is None:
                    ref = bacg

                rows.append({
                    'env': env, 'pair': label,
                    'scenario': 'metabolic_scaling',
                    'fungal_fraction': s,
                    'fungal_growth': fung_g,
                    'bac_growth': bacg,
                    'bac_rel_to_coinfection': bacg / ref if ref else np.nan,
                })

            for r in fung:
                r.bounds = base_bounds[r.id]

            g = pd.DataFrame([r for r in rows
                              if r['env'] == env and r['pair'] == label])
            half = g[g.fungal_fraction.between(0.45, 0.55)]
            print(f'    s=1.00 -> {ref:.4f}')
            if not half.empty:
                v = half.iloc[0]
                print(f'    s=0.50 -> {v.bac_growth:.4f} '
                      f'({100*v.bac_rel_to_coinfection:.1f}%)')
            last = g.iloc[-1]
            print(f'    s=0.00 -> {last.bac_growth:.4f} '
                  f'({100*last.bac_rel_to_coinfection:.1f}%)')

    df = pd.DataFrame(rows)
    out = RESULTS / 'fungistatic_scaling.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    # where does suppression start to matter?
    print('\n=== Suppression needed for a 25% / 50% bacterial reduction ===')
    for (env, pair), g in df.groupby(['env', 'pair']):
        g = g.sort_values('fungal_fraction', ascending=False)
        line = []
        for target in (0.75, 0.50):
            hit = g[g.bac_rel_to_coinfection <= target]
            line.append(f'{int((1-target)*100)}%: '
                        + (f's={hit.iloc[0].fungal_fraction:.2f}'
                           if not hit.empty else 'not reached'))
        print(f'  {env:5s} {pair:6s}  ' + '   '.join(line))


if __name__ == '__main__':
    main()
