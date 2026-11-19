"""전면의 꺾임점을 사전식 최적화로 정확히 구한다.
1단계: 세균 biomass 최대화. 2단계: 그 값에 고정하고 진균 biomass 최대화."""
import sys, pickle, copy
from pathlib import Path
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp

MODELS, DATA, RESULTS = ROOT/'models', ROOT/'data', ROOT/'results'

ENVS = [{'name':'CSOM','lumen':'csom_lumen.pckl','agora':'csom_agora.pckl'},
        {'name':'OE',  'lumen':'oe_lumen.pckl',  'agora':'oe_agora.pckl'}]
PAIRS = [{'bac_label':'PA','bac_model':'Pseudomonas_aeruginosa_NCGM2_S1'},
         {'bac_label':'SA','bac_model':'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'}]

cp_raw = cobra.io.read_sbml_model(str(MODELS/'cp'/'iDC1003_processed.xml'))
mapping = pd.read_csv(DATA/'id_mapping_cp.csv', index_col=0)
cp = cp_join_lumen_rxns(cp_raw.copy(), mapping)
for leak in ['EX_CO2__extr','EX_Oxygen__extr']:
    try: cp.reactions.get_by_id(leak).bounds = (0,0)
    except KeyError: pass
cp.objective = 'e_Biomass__cyto'
cp.solver = 'glpk'

rows = []
for env in ENVS:
    with open(DATA/env['lumen'],'rb') as f: diet_lumen = pickle.load(f)
    with open(DATA/env['agora'],'rb') as f: diet_agora = pickle.load(f)
    cp_solo = constrain_model(cp.copy(), diet_lumen, 'EX_u').optimize().objective_value or 0

    for p in PAIRS:
        bac = cobra.io.read_sbml_model(str(MODELS/'agora'/f'{p["bac_model"]}.xml'))
        bac_solo = constrain_model(bac.copy(), diet_agora, 'EX_').optimize().objective_value or 0
        bac_j = change_ids(bac.copy())
        bac_j, lr, lm = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lr, lm, diet_lumen, copy, 'agora')
        joint = constrain_model(joint, diet_lumen, 'EX_u')

        fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
        bac_b = next(r for r in joint.reactions
                     if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
        fu_b.bounds = (0.0, 1000); bac_b.bounds = (0.0, 1000)

        joint.objective = bac_b
        bac_max = joint.optimize().objective_value or 0
        bac_b.bounds = (bac_max, bac_max)
        joint.objective = fu_b
        fu_at_kink = joint.optimize().objective_value or 0

        rows.append(dict(env=env['name'], pair=f'CP+{p["bac_label"]}',
                         cp_growth=fu_at_kink, bac_growth=bac_max,
                         cp_boost=fu_at_kink/cp_solo, bac_boost=bac_max/bac_solo))
        print(f'{env["name"]:5s} CP+{p["bac_label"]}  '
              f'CP {fu_at_kink/cp_solo:.4f} / bac {bac_max/bac_solo:.4f}')

pd.DataFrame(rows).to_csv(RESULTS/'pareto_kink.csv', index=False)
print(f'\nSaved: {RESULTS}/pareto_kink.csv')
