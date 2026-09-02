"""Graded antifungal effect as a constraint on total fungal metabolic activity.

Two earlier attempts failed, both for the same underlying reason: partial
fungal suppression was expressed as a limit on individual reactions, which
says nothing about how much flux those reactions actually carry.

    08  fungal biomass fixed at zero. The network stayed fully active and
        the carbon previously spent on biomass was freed, so the bacterium
        lost nothing and every condition returned exactly 100%.
    11  every fungal reaction bound multiplied by s. Nominal bounds are
        +/-1000 while real fluxes are of order 1-10, so nothing binds until
        s falls below roughly 0.01 and the response is a step at zero.

Here a single constraint caps the summed flux carried by the fungal network
at a fraction s of its untreated value. The fungus may redistribute flux
within that budget, so s reads as the fraction of fungal metabolic activity
remaining, which is what a partially effective antifungal would leave.

    sum over fungal reactions of (forward + reverse) <= s * reference_total

Sanity checks printed at the end:
    s = 1  should reproduce the untreated co-infection growth
    s = 0  should reproduce the fungicidal result from 08

Output: results/fungal_activity_budget.csv
"""

import sys
import copy
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import cobra
from cobra import Metabolite, Reaction
from optlang.symbolics import Zero

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from functions import constrain_model, change_ids, tune_bac_lumen_rxns
from cp_helpers import cp_join_lumen_rxns, make_joint_model_cp

MODELS  = ROOT / 'models'
DATA    = ROOT / 'data'
RESULTS = ROOT / 'results'
RESULTS.mkdir(exist_ok=True)

N_STEPS = 41
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

    pair_results = pd.read_csv(RESULTS / 'pair_results.csv')
    rows, checks = [], []

    for env, (lf, af) in ENVS.items():
        diet_lum, diet_ago = load(lf), load(af)
        print(f'\n=== {env} ===')

        for p in PAIRS:
            label = f'CP+{p["bac_label"]}'
            base = pair_results[(pair_results.env == env) &
                                (pair_results.pair == label)]
            if base.empty:
                print(f'  {label}: absent from pair_results.csv — skip')
                continue
            cp_free = base.iloc[0]['cp_free']

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

            # ---- untreated reference: fungus at its joint free optimum,
            #      bacterium maximised, fungal activity measured by pFBA
            fu_b.bounds = (cp_free, cp_free)
            bac_b.bounds = (0.0, 1000.0)
            joint.objective = bac_b
            bac_ref = joint.optimize().fluxes[bac_b.id]

            pf = cobra.flux_analysis.pfba(joint)
            ref_total = float(sum(abs(pf.fluxes[r.id]) for r in fung))
            print(f'    reference bacterial growth : {bac_ref:.4f}')
            print(f'    reference fungal flux total: {ref_total:.2f}')

            # ---- add the activity budget constraint
            expr = Zero
            for r in fung:
                expr += r.forward_variable + r.reverse_variable
            budget = joint.problem.Constraint(expr, lb=0, ub=ref_total,
                                              name='fungal_activity_budget')
            joint.add_cons_vars(budget)

            # the fungus is now free to grow or not within its budget
            fu_b.bounds = (0.0, 1000.0)

            for s in np.concatenate([np.linspace(1.0, 0.1, 19), np.linspace(0.1, 0.0, 21)]):
                budget.ub = s * ref_total
                bac_b.bounds = (0.0, 1000.0)
                joint.objective = bac_b
                sol = joint.optimize()
                if sol.status != 'optimal':
                    bg = fg = np.nan
                else:
                    bg = sol.fluxes[bac_b.id]
                    fg = sol.fluxes[fu_b.id]
                rows.append({
                    'env': env, 'pair': label,
                    'scenario': 'activity_budget',
                    'fungal_activity': s,
                    'fungal_growth': fg,
                    'bac_growth': bg,
                    'bac_rel_to_coinfection': bg / bac_ref if bac_ref else np.nan,
                })

            joint.remove_cons_vars([budget])

            g = pd.DataFrame([r for r in rows
                              if r['env'] == env and r['pair'] == label])
            at1 = g.iloc[0]
            at0 = g.iloc[-1]
            checks.append({
                'env': env, 'pair': label,
                's1_pct': 100 * at1.bac_rel_to_coinfection,
                's0_pct': 100 * at0.bac_rel_to_coinfection,
            })
            for probe in (0.75, 0.50, 0.25):
                near = g.iloc[(g.fungal_activity - probe).abs().argmin()]
                print(f'    s={probe:.2f} -> {near.bac_growth:.4f} '
                      f'({100*near.bac_rel_to_coinfection:.1f}%)')
            print(f'    s=0.00 -> {at0.bac_growth:.4f} '
                  f'({100*at0.bac_rel_to_coinfection:.1f}%)')

    df = pd.DataFrame(rows)
    out = RESULTS / 'fungal_activity_budget.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    print('\n=== Sanity checks ===')
    print('  s=1 should be 100%; s=0 should match 08 fungicidal')
    for c in checks:
        print(f'  {c["env"]:5s} {c["pair"]:6s}  '
              f's=1: {c["s1_pct"]:6.1f}%   s=0: {c["s0_pct"]:6.1f}%')

    print('\n=== Fungal activity remaining at a given bacterial reduction ===')
    for (env, pair), g in df.groupby(['env', 'pair']):
        g = g.sort_values('fungal_activity', ascending=False)
        parts = []
        for target in (0.75, 0.50):
            hit = g[g.bac_rel_to_coinfection <= target]
            parts.append(f'{int((1-target)*100)}% drop at s='
                         + (f'{hit.iloc[0].fungal_activity:.3f}'
                            if not hit.empty else 'never'))
        print(f'  {env:5s} {pair:6s}  ' + '   '.join(parts))


if __name__ == '__main__':
    main()
