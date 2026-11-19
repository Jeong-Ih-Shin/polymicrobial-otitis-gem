"""진균 '제거'의 두 정의를 구분해 비교한다.
 A) 1:1 자유최적 (기준)
 B) CP biomass = 0 (Pareto 전면의 왼쪽 끝) — 진균 교환은 자유
 C) CP의 lumen 교환 전면 차단 (Figure 4의 완전 제거)
B에서 진균 교환 플럭스가 살아있는지 수치로 확인한다."""
import sys, pickle, copy
from pathlib import Path
import pandas as pd
import cobra
from cobra import Metabolite, Reaction

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp, categorize_joint_rxn

MODELS, DATA, RESULTS = ROOT/'models', ROOT/'data', ROOT/'results'
THRESHOLD = 0.01

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

def fungal_lumen_rxns(model):
    out = []
    for r in model.reactions:
        if categorize_joint_rxn(r) == 'fungal':
            if any(m.compartment == 'u' for m in r.metabolites):
                out.append(r)
    return out

rows = []
for env in ENVS:
    with open(DATA/env['lumen'],'rb') as f: diet_lumen = pickle.load(f)
    with open(DATA/env['agora'],'rb') as f: diet_agora = pickle.load(f)

    for p in PAIRS:
        tag = f'{env["name"]} / CP+{p["bac_label"]}'
        bac = cobra.io.read_sbml_model(str(MODELS/'agora'/f'{p["bac_model"]}.xml'))
        bac_j = change_ids(bac.copy())
        bac_j, lr, lm = tune_bac_lumen_rxns(bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lr, lm, diet_lumen, copy, 'agora')
        joint = constrain_model(joint, diet_lumen, 'EX_u')

        fu_b = joint.reactions.get_by_id('e_Biomass__cyto')
        bac_b = next(r for r in joint.reactions
                     if r.id.startswith('model2_biomass') and not r.id.startswith('model2_EX_'))
        fx = fungal_lumen_rxns(joint)

        # A) 1:1 자유최적
        fu_b.bounds = (0.001, 1000); bac_b.bounds = (0.001, 1000)
        joint.objective = {fu_b: 1.0, bac_b: 1.0}
        sA = joint.optimize()
        bacA = sA.fluxes[bac_b.id]
        actA = sum(1 for r in fx if abs(sA.fluxes[r.id]) > THRESHOLD)

        # B) CP biomass = 0, 교환은 자유
        fu_b.bounds = (0.0, 0.0); bac_b.bounds = (0.0, 1000)
        joint.objective = bac_b
        sB = joint.optimize()
        bacB = sB.objective_value or 0
        actB = sum(1 for r in fx if abs(sB.fluxes[r.id]) > THRESHOLD)
        magB = sum(abs(sB.fluxes[r.id]) for r in fx)

        # C) CP lumen 교환 전면 차단
        saved = {r.id: r.bounds for r in fx}
        for r in fx: r.bounds = (0.0, 0.0)
        fu_b.bounds = (0.0, 1000)
        joint.objective = bac_b
        sC = joint.optimize()
        bacC = sC.objective_value or 0
        for r in fx: r.bounds = saved[r.id]

        rows.append(dict(env=env['name'], pair=f'CP+{p["bac_label"]}',
                         n_fungal_lumen_rxns=len(fx),
                         bac_A=round(bacA,4), bac_B=round(bacB,4), bac_C=round(bacC,4),
                         pct_B=round(100*bacB/bacA,1), pct_C=round(100*bacC/bacA,1),
                         active_exch_A=actA, active_exch_B=actB,
                         total_exch_flux_B=round(magB,2)))
        print(f'{tag}')
        print(f'  진균 lumen 교환 반응 {len(fx)}개')
        print(f'  A 1:1 자유최적      : bac {bacA:.4f} (100.0%), 활성 교환 {actA}개')
        print(f'  B CP biomass=0      : bac {bacB:.4f} ({100*bacB/bacA:.1f}%), '
              f'활성 교환 {actB}개, 총 |flux| {magB:.2f}')
        print(f'  C lumen 교환 차단   : bac {bacC:.4f} ({100*bacC/bacA:.1f}%)\n')

df = pd.DataFrame(rows)
df.to_csv(RESULTS/'removal_definitions.csv', index=False)
print(df.to_string(index=False))
