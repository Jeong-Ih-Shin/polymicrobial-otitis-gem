"""두 biomass의 트레이드오프 전면(Pareto front)을 계산한다.
진균 biomass를 격자로 고정하고 각 지점에서 세균 biomass를 최대화한다.
목적함수 가중치를 쓰지 않으므로 가중치 선택에 의존하지 않는다."""
import sys, pickle, copy
from pathlib import Path
import numpy as np
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp

MODELS, DATA, RESULTS = ROOT/'models', ROOT/'data', ROOT/'results'
RESULTS.mkdir(exist_ok=True)
N_POINTS = 41

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

rows, anchors = [], []
for env in ENVS:
    with open(DATA/env['lumen'],'rb') as f: diet_lumen = pickle.load(f)
    with open(DATA/env['agora'],'rb') as f: diet_agora = pickle.load(f)
    cp_solo = constrain_model(cp.copy(), diet_lumen, 'EX_u').optimize().objective_value or 0

    for p in PAIRS:
        tag = f'{env["name"]} / CP+{p["bac_label"]}'
        bac = cobra.io.read_sbml_model(str(MODELS/'agora'/f'{p["bac_model"]}.xml'))
        bac_solo = constrain_model(bac.copy(), diet_agora, 'EX_').optimize().objective_value or 0
        bac_j = change_ids(bac.copy())
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets, diet_lumen, copy, 'agora')
        joint = constrain_model(joint, diet_lumen, 'EX_u')

        fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
        bac_b = next(r for r in joint.reactions
                     if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
        fu_b.bounds = (0.0, 1000); bac_b.bounds = (0.0, 1000)

        # 두 축의 최대치 (상대 조직은 제약하지 않음)
        joint.objective = fu_b
        fu_max = joint.optimize().objective_value or 0
        joint.objective = bac_b
        bac_max = joint.optimize().objective_value or 0

        # 1:1 동등가중 자유최적 (기존 논문의 지점)
        fu_b.bounds = (0.001, 1000); bac_b.bounds = (0.001, 1000)
        joint.objective = {fu_b: 1.0, bac_b: 1.0}
        s = joint.optimize()
        anchors.append(dict(env=env['name'], pair=f'CP+{p["bac_label"]}',
                            point='equal_weight_optimum',
                            cp_growth=s.fluxes[fu_b.id], bac_growth=s.fluxes[bac_b.id],
                            cp_solo=cp_solo, bac_solo=bac_solo))
        fu_b.bounds = (0.0, 1000); bac_b.bounds = (0.0, 1000)

        # 전면 계산
        joint.objective = bac_b
        for fu_fix in np.linspace(0.0, fu_max, N_POINTS):
            fu_b.bounds = (fu_fix, fu_fix)
            sol = joint.optimize()
            bac_val = sol.objective_value if sol.status == 'optimal' else np.nan
            rows.append(dict(env=env['name'], pair=f'CP+{p["bac_label"]}',
                             cp_growth=fu_fix, bac_growth=bac_val,
                             cp_boost=fu_fix/cp_solo if cp_solo else np.nan,
                             bac_boost=(bac_val/bac_solo) if bac_solo else np.nan))
        fu_b.bounds = (0.0, 1000)
        print(f'{tag}: CP max={fu_max:.3f} (solo {cp_solo:.3f}), '
              f'bac max={bac_max:.3f} (solo {bac_solo:.3f})')

pd.DataFrame(rows).to_csv(RESULTS/'pareto_front.csv', index=False)
pd.DataFrame(anchors).to_csv(RESULTS/'pareto_anchors.csv', index=False)
print(f'\nSaved: {RESULTS}/pareto_front.csv, pareto_anchors.csv')
