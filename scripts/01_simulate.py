"""Solo and joint biomass for CP + bacterium pairs across environments.
Environments: CSOM (chronic suppurative OM), OE (otitis externa)
Pairs: CP+PA, CP+SA
Output: results/pair_results.csv  (env × pair rows)"""

import sys, pickle, copy
from pathlib import Path
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

ENVS = [
    {'name': 'CSOM', 'lumen': 'csom_lumen.pckl', 'agora': 'csom_agora.pckl'},
    {'name': 'OE',   'lumen': 'oe_lumen.pckl',   'agora': 'oe_agora.pckl'},
]
PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]

# Load CP once
print('Loading CP iDC1003...')
cp_raw = cobra.io.read_sbml_model(str(MODELS / 'cp' / 'iDC1003_processed.xml'))
mapping = pd.read_csv(DATA / 'id_mapping_cp.csv', index_col=0)
cp = cp_join_lumen_rxns(cp_raw.copy(), mapping)
for leak in ['EX_CO2__extr', 'EX_Oxygen__extr']:
    try: cp.reactions.get_by_id(leak).bounds = (0, 0)
    except KeyError: pass
cp.objective = 'e_Biomass__cyto'
cp.solver = 'glpk'

results = []
for env in ENVS:
    print(f'\n=== Environment: {env["name"]} ===')
    with open(DATA / env['lumen'], 'rb') as f: diet_lumen = pickle.load(f)
    with open(DATA / env['agora'], 'rb') as f: diet_agora = pickle.load(f)

    cp_solo = constrain_model(cp.copy(), diet_lumen, 'EX_u').optimize().objective_value or 0
    print(f'CP solo: {cp_solo:.4f}')

    for p in PAIRS:
        print(f'  [{p["bac_label"]}] {p["bac_model"]}')
        bac = cobra.io.read_sbml_model(str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
        bac_solo = constrain_model(bac.copy(), diet_agora, 'EX_').optimize().objective_value or 0
        print(f'    bac solo: {bac_solo:.4f}')
        if bac_solo < 1e-4:
            print('    bacterium does not grow alone — skip')
            continue

        bac_j = change_ids(bac.copy())
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets, diet_lumen, copy, 'agora')
        joint = constrain_model(joint, diet_lumen, 'EX_u')

        fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
        bac_b = next(r for r in joint.reactions
                     if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
        fu_b.lower_bound = 0.001; fu_b.upper_bound = 1000
        bac_b.lower_bound = 0.001; bac_b.upper_bound = 1000
        joint.objective = {fu_b: 1.0, bac_b: 1.0}
        sol = joint.optimize()

        fu_free = sol.fluxes[fu_b.id]
        bac_free = sol.fluxes[bac_b.id]
        fu_boost = fu_free / cp_solo if cp_solo > 0 else 0
        bac_boost = bac_free / bac_solo if bac_solo > 0 else 0
        print(f'    Free optimum: CP={fu_free:.3f} ({fu_boost:.2f}x), bac={bac_free:.3f} ({bac_boost:.2f}x)')

        results.append({
            'env'       : env['name'],
            'pair'      : f'CP+{p["bac_label"]}',
            'bac_label' : p['bac_label'],
            'bac_model' : p['bac_model'],
            'cp_solo'   : cp_solo,
            'bac_solo'  : bac_solo,
            'cp_free'   : fu_free,
            'bac_free'  : bac_free,
            'cp_boost'  : fu_boost,
            'bac_boost' : bac_boost,
        })

pd.DataFrame(results).to_csv(RESULTS / 'pair_results.csv', index=False)
print(f'\nSaved: {RESULTS}/pair_results.csv')
print('\n=== Summary ===')
print(pd.DataFrame(results)[['env', 'pair', 'cp_solo', 'bac_solo', 'cp_free', 'bac_free',
                             'cp_boost', 'bac_boost']].to_string(index=False))
