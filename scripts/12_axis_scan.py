"""Scan a multi-key medium axis and separate real coupling from ratio artefacts.

The oxygen scan (05) produced a rising boost factor that turned out to be an
artefact: joint growth was identical at every grid point and only the solo
denominator moved. Any axis scan therefore has to report absolute growth
rates alongside the ratio, and say plainly whether the joint solution
responded at all.

Component attribution (10) identified iron as the only group whose reversal
abolishes the CSOM -> OE switch, so iron is the default axis here. Set
AXIS to any group in GROUPS to scan a different one.

Output: results/axis_scan_<axis>.csv
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

# ----------------------------------------------------------------- config
AXIS      = 'iron'
N_POINTS  = 81
PAD       = 0.25          # extend 25% beyond the CSOM-OE span
LOW, HIGH = 0.90, 1.10
FLAT_TOL  = 1e-6          # joint fluxes varying less than this count as flat

GROUPS = {
    'oxygen'      : ['_o2', 'o2('],
    'mucin'       : ['acgam', 'acnam', 'fuc_l'],
    'ammonium'    : ['nh4'],
    'fermentation': ['lac_', '_for', '_ac', 'succ', 'ppa'],
    'iron'        : ['fe2', 'fe3'],
    'amino_acids' : ['ala_l', 'arg_l', 'asn_l', 'asp_l', 'cys_l', 'gln_l',
                     'glu_l', '_gly', 'gly(', 'his_l', 'ile_l', 'leu_l',
                     'lys_l', 'met_l', 'phe_l', 'pro_l', 'ser_l', 'thr_l',
                     'trp_l', 'tyr_l', 'val_l'],
    'sugars'      : ['glc_d', 'gal', 'man', 'gam'],
}

PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]


def load(f):
    with open(DATA / f, 'rb') as fh:
        return pickle.load(fh)


def axis_keys(diet, axis):
    pats = GROUPS[axis]
    return sorted(k for k in diet if any(p in k.lower() for p in pats))


def classify(f, b):
    sgn = lambda v: '+' if v > HIGH else ('-' if v < LOW else 'o')
    sf, sb = sgn(f), sgn(b)
    return {('+', '+'): 'Mutualism',    ('-', '-'): 'Competition',
            ('o', 'o'): 'Neutralism',   ('+', '-'): 'Parasitism',
            ('-', '+'): 'Parasitism',   ('+', 'o'): 'Commensalism',
            ('o', '+'): 'Commensalism', ('-', 'o'): 'Amensalism',
            ('o', '-'): 'Amensalism'}[(sf, sb)]


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


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

    diets = {
        'CSOM': (load('csom_lumen.pckl'), load('csom_agora.pckl')),
        'OE':   (load('oe_lumen.pckl'),   load('oe_agora.pckl')),
    }

    kl = axis_keys(diets['CSOM'][0], AXIS)
    ka = axis_keys(diets['CSOM'][1], AXIS)
    print(f"\nAxis '{AXIS}'")
    print(f'  lumen keys: {kl}')
    print(f'  agora keys: {ka}')
    for k in kl:
        print(f'    {k:<18} CSOM={diets["CSOM"][0].get(k, 0):>8}  '
              f'OE={diets["OE"][0].get(k, 0):>8}')

    # alpha 0 = CSOM values, 1 = OE values, padded on both sides
    alphas = np.linspace(-PAD, 1 + PAD, N_POINTS)

    rows = []
    for base_env in ['CSOM', 'OE']:
        base_l, base_a = diets[base_env]
        print(f'\n=== Base medium: {base_env} '
              f'(all components except {AXIS} fixed) ===')

        for p in PAIRS:
            label = f'CP+{p["bac_label"]}'
            print(f'  [{label}] assembling joint model...')
            bac = cobra.io.read_sbml_model(
                str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
            bac_j = change_ids(bac.copy())
            bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
                bac_j, cp, Metabolite, Reaction, 'agora')
            joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                        base_l, copy, 'agora')
            fu_b, bac_b = biomass_rxns(joint)
            cp_m, bac_m = cp.copy(), bac.copy()

            for a in alphas:
                d_l, d_a = dict(base_l), dict(base_a)
                for k in kl:
                    c, o = diets['CSOM'][0].get(k, 0.0), diets['OE'][0].get(k, 0.0)
                    d_l[k] = float(c + a * (o - c))
                for k in ka:
                    c, o = diets['CSOM'][1].get(k, 0.0), diets['OE'][1].get(k, 0.0)
                    d_a[k] = float(c + a * (o - c))

                cp_solo = constrain_model(
                    cp_m, d_l, 'EX_u').optimize().objective_value or 0.0
                bac_solo = constrain_model(
                    bac_m, d_a, 'EX_').optimize().objective_value or 0.0

                constrain_model(joint, d_l, 'EX_u')
                fu_b.bounds = (0.001, 1000)
                bac_b.bounds = (0.001, 1000)
                joint.objective = {fu_b: 1.0, bac_b: 1.0}
                sol = joint.optimize()

                if sol.status != 'optimal' or cp_solo <= 1e-9 or bac_solo <= 1e-9:
                    rows.append({'base_env': base_env, 'pair': label,
                                 'alpha': a, 'cp_solo': cp_solo,
                                 'bac_solo': bac_solo, 'cp_free': np.nan,
                                 'bac_free': np.nan, 'cp_boost': np.nan,
                                 'bac_boost': np.nan, 'mode': 'infeasible'})
                    continue

                cf, bf = sol.fluxes[fu_b.id], sol.fluxes[bac_b.id]
                cb, bb = cf / cp_solo, bf / bac_solo
                rows.append({
                    'base_env': base_env, 'pair': label, 'alpha': a,
                    'axis_value_lumen': d_l[kl[0]] if kl else np.nan,
                    'cp_solo': cp_solo, 'bac_solo': bac_solo,
                    'cp_free': cf, 'bac_free': bf,
                    'cp_boost': cb, 'bac_boost': bb,
                    'mode': classify(cb, bb),
                })

            g = pd.DataFrame([r for r in rows
                              if r['base_env'] == base_env
                              and r['pair'] == label])
            print(f"    modes: {g['mode'].value_counts().to_dict()}")
            print(f"    cp_free  range {g.cp_free.min():.4f} - {g.cp_free.max():.4f}")
            print(f"    bac_free range {g.bac_free.min():.4f} - {g.bac_free.max():.4f}")
            print(f"    bac_solo range {g.bac_solo.min():.4f} - {g.bac_solo.max():.4f}")

    df = pd.DataFrame(rows)
    out = RESULTS / f'axis_scan_{AXIS}.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    # ---------------------------------------------- artefact diagnostics
    print('\n=== Is the joint solution actually responding? ===')
    for (env, pair), g in df.groupby(['base_env', 'pair']):
        sp_c = g.cp_free.max() - g.cp_free.min()
        sp_b = g.bac_free.max() - g.bac_free.min()
        flat = sp_c < FLAT_TOL and sp_b < FLAT_TOL
        verdict = ('FLAT - boost changes are a denominator artefact'
                   if flat else 'joint solution responds to the axis')
        print(f'  {env:5s} {pair:6s}  d(cp_free)={sp_c:.2e}  '
              f'd(bac_free)={sp_b:.2e}  -> {verdict}')

    print('\n=== Mode transitions ===')
    for (env, pair), g in df.groupby(['base_env', 'pair']):
        g = g.sort_values('alpha')
        modes, al = g['mode'].tolist(), g['alpha'].tolist()
        found = False
        for i in range(1, len(modes)):
            if modes[i] != modes[i - 1]:
                print(f'  {env} / {pair}: {modes[i-1]} -> {modes[i]} '
                      f'at alpha={al[i]:.3f}')
                found = True
        if not found:
            print(f'  {env} / {pair}: no transition ({modes[0]} throughout)')


if __name__ == '__main__':
    main()
