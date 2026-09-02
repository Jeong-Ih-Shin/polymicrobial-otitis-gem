"""Which medium component drives the CSOM -> OE interaction switch?

The oxygen scan (05) and the oxygen x mucin phase diagram (06) both failed
to reproduce the switch, so the responsible component lies elsewhere in the
medium definition. This script attributes the switch directly.

Three analyses
--------------
1. morph            All components that differ between the two media are
                    interpolated together, CSOM -> OE, and the mode is
                    recorded along the way.
2. sufficiency      Start from CSOM; swap ONE component group to its OE
                    value. Does the OE mode appear? A group that alone
                    reproduces the switch is sufficient.
3. necessity        Start from OE; revert ONE component group to its CSOM
                    value. Does the CSOM mode return? A group whose removal
                    abolishes the switch is necessary.

A component that is both sufficient and necessary is the driver.

Output: results/component_attribution.csv
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

N_MORPH = 21
LOW, HIGH = 0.90, 1.10

# Component groups. Any differing key not matched here becomes its own
# single-key group, so attribution always covers the full difference set.
GROUPS = {
    'oxygen'      : ['_o2', 'o2('],
    'mucin'       : ['acgam', 'acnam', 'fuc_l'],
    'ammonium'    : ['nh4'],
    'fermentation': ['lac_', 'lac-', 'for(', '_for', 'ac(', '_ac', 'succ',
                     'etoh', 'ppa', 'but', 'pyr'],
    'iron'        : ['fe2', 'fe3', 'pheme', 'sheme'],
    'amino_acids' : ['ala_l', 'arg_l', 'asn_l', 'asp_l', 'cys_l', 'gln_l',
                     'glu_l', 'gly', 'his_l', 'ile_l', 'leu_l', 'lys_l',
                     'met_l', 'phe_l', 'pro_l', 'ser_l', 'thr_l', 'trp_l',
                     'tyr_l', 'val_l'],
    'sugars'      : ['glc_d', 'gal', 'man', 'gam'],
}

PAIRS = [
    {'bac_label': 'PA', 'bac_model': 'Pseudomonas_aeruginosa_NCGM2_S1'},
    {'bac_label': 'SA', 'bac_model': 'Staphylococcus_aureus_subsp_aureus_USA300_FPR3757'},
]


def load(f):
    with open(DATA / f, 'rb') as fh:
        return pickle.load(fh)


def classify(f, b):
    sgn = lambda v: '+' if v > HIGH else ('-' if v < LOW else 'o')
    sf, sb = sgn(f), sgn(b)
    return {('+', '+'): 'Mutualism',    ('-', '-'): 'Competition',
            ('o', 'o'): 'Neutralism',   ('+', '-'): 'Parasitism',
            ('-', '+'): 'Parasitism',   ('+', 'o'): 'Commensalism',
            ('o', '+'): 'Commensalism', ('-', 'o'): 'Amensalism',
            ('o', '-'): 'Amensalism'}[(sf, sb)]


def group_of(key):
    k = key.lower()
    for g, pats in GROUPS.items():
        if any(p in k for p in pats):
            return g
    return f'single:{key}'


def diff_keys(a, b):
    """Keys whose uptake bound differs between two diets (absent = 0)."""
    out = {}
    for k in set(a) | set(b):
        va, vb = float(a.get(k, 0.0)), float(b.get(k, 0.0))
        if abs(va - vb) > 1e-12:
            out[k] = (va, vb)
    return out


def biomass_rxns(joint):
    fu = joint.reactions.get_by_id('e_Biomass__cyto')
    bac = next(r for r in joint.reactions
               if r.id.startswith('model2_biomass')
               and not r.id.startswith('model2_EX_'))
    return fu, bac


def evaluate(joint, fu_b, bac_b, cp_model, bac_model, d_lum, d_ago):
    """Return (cp_boost, bac_boost, mode) for one medium definition."""
    cp_solo = constrain_model(cp_model, d_lum, 'EX_u').optimize().objective_value or 0.0
    bac_solo = constrain_model(bac_model, d_ago, 'EX_').optimize().objective_value or 0.0
    constrain_model(joint, d_lum, 'EX_u')
    fu_b.bounds = (0.001, 1000)
    bac_b.bounds = (0.001, 1000)
    joint.objective = {fu_b: 1.0, bac_b: 1.0}
    sol = joint.optimize()
    if sol.status != 'optimal' or cp_solo <= 1e-9 or bac_solo <= 1e-9:
        return np.nan, np.nan, 'infeasible'
    fb = sol.fluxes[fu_b.id] / cp_solo
    bb = sol.fluxes[bac_b.id] / bac_solo
    return fb, bb, classify(fb, bb)


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

    csom_l, csom_a = load('csom_lumen.pckl'), load('csom_agora.pckl')
    oe_l, oe_a = load('oe_lumen.pckl'), load('oe_agora.pckl')

    dl = diff_keys(csom_l, oe_l)
    da = diff_keys(csom_a, oe_a)
    print(f'\nDiffering keys: lumen n={len(dl)}, agora n={len(da)}')

    groups = sorted({group_of(k) for k in dl} | {group_of(k) for k in da})
    print(f'Component groups ({len(groups)}):')
    for g in groups:
        ks = [k for k in dl if group_of(k) == g]
        print(f'  {g:<28} lumen keys: {ks}')

    rows = []
    for p in PAIRS:
        label = f'CP+{p["bac_label"]}'
        print(f'\n[{label}] assembling joint model...')
        bac = cobra.io.read_sbml_model(
            str(MODELS / 'agora' / f'{p["bac_model"]}.xml'))
        bac_j = change_ids(bac.copy())
        bac_j, lumen_rxns, lumen_mets = tune_bac_lumen_rxns(
            bac_j, cp, Metabolite, Reaction, 'agora')
        joint = make_joint_model_cp(cp, bac_j, lumen_rxns, lumen_mets,
                                    csom_l, copy, 'agora')
        fu_b, bac_b = biomass_rxns(joint)
        cp_m, bac_m = cp.copy(), bac.copy()

        def ev(d_lum, d_ago):
            return evaluate(joint, fu_b, bac_b, cp_m, bac_m, d_lum, d_ago)

        # reference endpoints
        for name, (dlum, dago) in [('CSOM_reference', (csom_l, csom_a)),
                                   ('OE_reference', (oe_l, oe_a))]:
            fb, bb, mode = ev(dict(dlum), dict(dago))
            rows.append({'pair': label, 'analysis': 'reference',
                         'component': name, 'alpha': np.nan,
                         'cp_boost': fb, 'bac_boost': bb, 'mode': mode})
            print(f'  {name:<16} CP={fb:.3f}x bac={bb:.3f}x  {mode}')

        # ------------------------------------------------ 1. morph
        print('  morphing CSOM -> OE ...')
        for alpha in np.linspace(0, 1, N_MORPH):
            d_lum, d_ago = dict(csom_l), dict(csom_a)
            for k, (v_c, v_o) in dl.items():
                d_lum[k] = v_c + alpha * (v_o - v_c)
            for k, (v_c, v_o) in da.items():
                d_ago[k] = v_c + alpha * (v_o - v_c)
            fb, bb, mode = ev(d_lum, d_ago)
            rows.append({'pair': label, 'analysis': 'morph',
                         'component': 'all', 'alpha': alpha,
                         'cp_boost': fb, 'bac_boost': bb, 'mode': mode})

        # ------------------------------------------------ 2. sufficiency
        print('  sufficiency (CSOM + one group at OE value) ...')
        for g in groups:
            d_lum, d_ago = dict(csom_l), dict(csom_a)
            for k, (_, v_o) in dl.items():
                if group_of(k) == g:
                    d_lum[k] = v_o
            for k, (_, v_o) in da.items():
                if group_of(k) == g:
                    d_ago[k] = v_o
            fb, bb, mode = ev(d_lum, d_ago)
            rows.append({'pair': label, 'analysis': 'sufficiency',
                         'component': g, 'alpha': np.nan,
                         'cp_boost': fb, 'bac_boost': bb, 'mode': mode})
            print(f'    +{g:<28} {mode}')

        # ------------------------------------------------ 3. necessity
        print('  necessity (OE with one group reverted to CSOM) ...')
        for g in groups:
            d_lum, d_ago = dict(oe_l), dict(oe_a)
            for k, (v_c, _) in dl.items():
                if group_of(k) == g:
                    d_lum[k] = v_c
            for k, (v_c, _) in da.items():
                if group_of(k) == g:
                    d_ago[k] = v_c
            fb, bb, mode = ev(d_lum, d_ago)
            rows.append({'pair': label, 'analysis': 'necessity',
                         'component': g, 'alpha': np.nan,
                         'cp_boost': fb, 'bac_boost': bb, 'mode': mode})
            print(f'    -{g:<28} {mode}')

    df = pd.DataFrame(rows)
    out = RESULTS / 'component_attribution.csv'
    df.to_csv(out, index=False)
    print(f'\nSaved: {out}  ({len(df)} rows)')

    # ------------------------------------------------ verdict
    print('\n=== Attribution ===')
    for pair, g in df.groupby('pair'):
        ref = g[g.analysis == 'reference'].set_index('component')['mode']
        m_csom, m_oe = ref.get('CSOM_reference'), ref.get('OE_reference')
        print(f'\n{pair}:  CSOM={m_csom}  OE={m_oe}')
        if m_csom == m_oe:
            print('  no switch between references — nothing to attribute')
            continue
        suf = g[g.analysis == 'sufficiency']
        nec = g[g.analysis == 'necessity']
        s_hits = set(suf[suf['mode'] == m_oe]['component'])
        n_hits = set(nec[nec['mode'] == m_csom]['component'])
        print(f'  sufficient alone : {sorted(s_hits) or "none"}')
        print(f'  necessary        : {sorted(n_hits) or "none"}')
        both = s_hits & n_hits
        if both:
            print(f'  DRIVER           : {sorted(both)}')
        elif not s_hits and not n_hits:
            print('  no single group accounts for the switch — '
                  'likely a combination effect; inspect the morph curve')

        mo = g[g.analysis == 'morph'].sort_values('alpha')
        modes = mo['mode'].tolist()
        alphas = mo['alpha'].tolist()
        for i in range(1, len(modes)):
            if modes[i] != modes[i - 1]:
                print(f'  morph switch     : {modes[i-1]} -> {modes[i]} '
                      f'at alpha={alphas[i]:.2f}')
                break
        else:
            print('  morph switch     : none along the interpolation')


if __name__ == '__main__':
    main()
