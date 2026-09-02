"""Build OM media for 3 clinical scenarios (acute / ome / csom)."""
import pickle
import copy

# Scaffold: Western diet (has full cofactors, minerals, AAs)
with open('dat/calb_western_diet.pckl', 'rb') as f:
    western = pickle.load(f)

# Remove irrelevant plant-derived carbons
REMOVE = ['arab_L', 'cellb', 'drib', 'fru', 'glcn', 'lcts', 'malt',
          'sucr', 'xyl_D', 'pect', 'strch1', 'strch2', 'starch1200',
          'rmn', 'tre']

SCENARIOS = {
    'acute': {
        'carbons': {'glc_D': -1.0, 'acgam': -0.1, 'acnam': -0.05, 'fuc_L': -0.05,
                    'gal': -0.1, 'man': -0.1, 'gam': -0.1,
                    'lac_L': -0.1, 'ac': -0.1, 'for': -0.1},
        'iron': {'fe2': -0.05, 'fe3': -0.05, 'fe3dcit': -0.05},
        'o2': -10.0, 'nh4': -1.0, 'aa_bound': -0.3,
    },
    'ome': {
        'carbons': {'glc_D': -0.3, 'acgam': -0.5, 'acnam': -0.3, 'fuc_L': -0.3,
                    'gal': -0.3, 'man': -0.2, 'gam': -0.3,
                    'lac_L': -0.5, 'ac': -0.3, 'for': -0.1},
        'iron': {'fe2': -0.02, 'fe3': -0.02, 'fe3dcit': -0.02},
        'o2': -3.0, 'nh4': -1.0, 'aa_bound': -0.2,
    },
    'csom': {
        'carbons': {'glc_D': -0.05, 'acgam': -0.5, 'acnam': -0.5, 'fuc_L': -0.5,
                    'gal': -0.3, 'man': -0.3, 'gam': -0.3,
                    'lac_L': -2.0, 'ac': -0.5, 'for': -0.3,
                    'succ': -0.2, 'ppa': -0.2},
        'iron': {'fe2': -0.005, 'fe3': -0.005, 'fe3dcit': -0.005},
        'o2': -0.5, 'nh4': -1.0, 'aa_bound': -0.2,
    },
}

AA_LIST = ['ala_L', 'arg_L', 'asn_L', 'asp_L', 'cys_L', 'gln_L', 'glu_L', 'gly',
           'his_L', 'ile_L', 'leu_L', 'lys_L', 'met_L', 'phe_L', 'pro_L',
           'ser_L', 'thr_L', 'trp_L', 'tyr_L', 'val_L']

def build(name):
    s = SCENARIOS[name]
    m = copy.deepcopy(western)
    for k in REMOVE:
        m.pop(f'EX_u_{k}', None)
    for met, b in s['carbons'].items():
        m[f'EX_u_{met}'] = b
    for met, b in s['iron'].items():
        m[f'EX_u_{met}'] = b
    m['EX_u_o2'] = s['o2']
    m['EX_u_nh4'] = s['nh4']
    for aa in AA_LIST:
        if f'EX_u_{aa}' in m:
            m[f'EX_u_{aa}'] = s['aa_bound']
    return m

def to_agora(d):
    return {k.replace('EX_u_', 'EX_') + '(e)': v for k, v in d.items()}

print('Building OM media...')
for sc in ['acute', 'ome', 'csom']:
    calb = build(sc)
    agora = to_agora(calb)
    with open(f'dat/calb_om_{sc}.pckl', 'wb') as f:
        pickle.dump(calb, f)
    with open(f'dat/agora_om_{sc}.pckl', 'wb') as f:
        pickle.dump(agora, f)
    print(f'  [{sc}] {len(calb)} entries')
    print(f'    glc_D={calb.get("EX_u_glc_D")}, lac_L={calb.get("EX_u_lac_L")}, '
          f'acgam={calb.get("EX_u_acgam")}, o2={calb.get("EX_u_o2")}, '
          f'fe2={calb.get("EX_u_fe2")}')
print('\nDone!')
