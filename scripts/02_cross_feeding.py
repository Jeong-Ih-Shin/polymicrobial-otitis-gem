"""Cross-feeding direction analysis at free-optimum biomass (pFBA).
Output: results/cross_feeding.csv  (env × pair × metabolite)"""

import sys, pickle, copy
from pathlib import Path
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp, categorize_joint_rxn

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'
THRESHOLD = 0.01

ENV_DIET = {
    'CSOM': 'csom_lumen.pckl',
    'OE'  : 'oe_lumen.pckl',
}

pair_results = pd.read_csv(RESULTS / 'pair_results.csv')

# Load CP once
cp_raw = cobra.io.read_sbml_model(str(MODELS / 'cp' / 'iDC1003_processed.xml'))
mapping = pd.read_csv(DATA / 'id_mapping_cp.csv', index_col=0)
cp = cp_join_lumen_rxns(cp_raw.copy(), mapping)
for leak in ['EX_CO2__extr', 'EX_Oxygen__extr']:
    try: cp.reactions.get_by_id(leak).bounds = (0, 0)
    except KeyError: pass
cp.objective = 'e_Biomass__cyto'
cp.solver = 'glpk'

rows = []
for _, row in pair_results.iterrows():
    print(f'\n[{row["env"]} / {row["pair"]}] pFBA at free optimum...')
    with open(DATA / ENV_DIET[row['env']], 'rb') as f:
        diet_lumen = pickle.load(f)

    bac = cobra.io.read_sbml_model(str(MODELS / 'agora' / f'{row["bac_model"]}.xml'))
    bac_j = change_ids(bac.copy())
    bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
    joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets, diet_lumen, copy, 'agora')
    joint = constrain_model(joint, diet_lumen, 'EX_u')

    fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
    bac_b = next(r for r in joint.reactions
                 if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
    fu_b.lower_bound = fu_b.upper_bound = row['cp_free']
    bac_b.lower_bound = bac_b.upper_bound = row['bac_free']
    joint.objective = fu_b

    try:
        sol = cobra.flux_analysis.pfba(joint)
    except Exception as e:
        print(f'  pFBA failed: {e}'); continue

    fu_side, bac_side, sys_side = {}, {}, {}
    for r in joint.reactions:
        cat = categorize_joint_rxn(r)
        if cat == 'system':
            sys_side[r.id[5:]] = sol.fluxes[r.id]
        elif cat in ('fungal', 'BAC'):
            for m in r.metabolites:
                if m.compartment == 'u':
                    agora = m.id[:-2] if m.id.endswith('_u') else m.id
                    (fu_side if cat == 'fungal' else bac_side)[agora] = sol.fluxes[r.id]
                    break

    for met in sorted(set(fu_side) | set(bac_side)):
        ff = fu_side.get(met, 0); bf = bac_side.get(met, 0); sf = sys_side.get(met, 0)
        if abs(ff) < THRESHOLD and abs(bf) < THRESHOLD: continue
        if ff > THRESHOLD and bf < -THRESHOLD: d = 'fungus_to_bac'
        elif ff < -THRESHOLD and bf > THRESHOLD: d = 'bac_to_fungus'
        elif ff > THRESHOLD and bf > THRESHOLD: d = 'both_produce'
        elif ff < -THRESHOLD and bf < -THRESHOLD: d = 'both_consume'
        elif abs(bf) < THRESHOLD: d = 'fungus_only'
        elif abs(ff) < THRESHOLD: d = 'bac_only'
        else: d = 'unclassified'
        rows.append({
            'env'         : row['env'],
            'pair'        : row['pair'],
            'metabolite'  : met,
            'fungus_flux' : round(ff, 3),
            'bac_flux'    : round(bf, 3),
            'sys_flux'    : round(sf, 3),
            'direction'   : d,
            'mag'         : round((abs(ff) + abs(bf)) / 2, 3),
        })

df = pd.DataFrame(rows)
df.to_csv(RESULTS / 'cross_feeding.csv', index=False)
print(f'\nSaved {len(df)} entries: {RESULTS}/cross_feeding.csv')

print('\n=== Direction counts (env × pair) ===')
for (env, pn), grp in df.groupby(['env', 'pair']):
    print(f'\n{env} / {pn}:')
    for d, n in grp.direction.value_counts().items():
        print(f'  {d:<20} n={n}')
