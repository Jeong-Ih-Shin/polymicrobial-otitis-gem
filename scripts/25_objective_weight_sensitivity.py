"""biomass 목적함수 가중치를 바꿔가며 상호작용 분류가 유지되는지 본다."""
import sys, pickle, copy, math
from pathlib import Path
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns, interaction_type
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp

MODELS, DATA, RESULTS = ROOT/'models', ROOT/'data', ROOT/'results'
RESULTS.mkdir(exist_ok=True)

ENVS = [{'name':'CSOM','lumen':'csom_lumen.pckl','agora':'csom_agora.pckl'},
        {'name':'OE',  'lumen':'oe_lumen.pckl',  'agora':'oe_agora.pckl'}]
PAIRS = [{'bac_label':'PA','bac_model':'Pseudomonas_aeruginosa_NCGM2_S1'},
         {'bac_label':'SA','bac_model':'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'}]
WEIGHTS = [(1,1),(1,2),(2,1),(1,3),(3,1),(1,5),(5,1)]

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
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets, diet_lumen, copy, 'agora')
        joint = constrain_model(joint, diet_lumen, 'EX_u')

        fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
        bac_b = next(r for r in joint.reactions
                     if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
        fu_b.bounds = (0.001, 1000); bac_b.bounds = (0.001, 1000)

        for wf, wb in WEIGHTS:
            joint.objective = {fu_b: float(wf), bac_b: float(wb)}
            sol = joint.optimize()
            fu_free, bac_free = sol.fluxes[fu_b.id], sol.fluxes[bac_b.id]
            fu_boost = fu_free/cp_solo if cp_solo > 0 else 0
            bac_boost = bac_free/bac_solo if bac_solo > 0 else 0
            ccf = math.log2(fu_boost) if fu_boost > 0 else float('-inf')
            bcf = math.log2(bac_boost) if bac_boost > 0 else float('-inf')
            mode, _ = interaction_type(ccf, bcf, math)
            rows.append(dict(env=env['name'], pair=f'CP+{p["bac_label"]}',
                             w_fungus=wf, w_bacterium=wb,
                             cp_boost=round(fu_boost,4), bac_boost=round(bac_boost,4),
                             mode=mode))
            print(f'{env["name"]:5s} CP+{p["bac_label"]} w={wf}:{wb}  '
                  f'CP={fu_boost:.3f} bac={bac_boost:.3f}  {mode}')

df = pd.DataFrame(rows)
df.to_csv(RESULTS/'weight_sensitivity.csv', index=False)
print('\n=== 가중치에 따라 분류가 달라진 조합 ===')
chk = df.groupby(['env','pair'])['mode'].nunique()
print(chk[chk > 1] if (chk > 1).any() else '없음 — 모든 가중치에서 분류 동일')
