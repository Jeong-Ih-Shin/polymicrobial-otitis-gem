"""Build OE (otitis externa) media — derived from CSOM scaffold,
modified for external ear canal (skin surface) environment.

OE characteristics (vs CSOM):
  + Aerobic (external surface exposed to air): o2 -0.5 → -8.0
  + Lower glucose (skin < mucosa): glc_D -0.05 → -0.1
  - No mucin-derived sugars (acgam ↓, acnam removed, fuc_L removed)
  + Sweat metabolites: lactate (lower), glycine, alanine, higher NH4
  + Cerumen iron-restricted (lactoferrin) but less than CSOM
  - No fermentation byproducts (succ, ppa removed)

References:
  - Cerumen composition: Burkhart CN, Burkhart CG. Otolaryngol Head Neck Surg. 2000
  - Sweat metabolites: Harker M, et al. Br J Dermatol. 2014
  - OE microbiology: Schaefer P, Baugh RF. Am Fam Physician. 2012

Limitation: cerumen lipids (squalene, wax esters, fatty acids) are not
GEM-tractable in AGORA exchange namespace and are therefore omitted.
"""

import pickle
import copy
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'

# Load CSOM as scaffold
with open(DATA / 'csom_lumen.pckl', 'rb') as f:
    oe = copy.deepcopy(pickle.load(f))

# 1. Remove mucin-derived carbons (skin has no mucin)
for met in ['acnam', 'fuc_L']:
    oe.pop(f'EX_u_{met}', None)

# 2. Close fermentation byproducts (CSOM-specific anaerobic environment)
for met in ['succ', 'ppa', 'for']:
    key = f'EX_u_{met}'
    if key in oe and oe[key] < 0:
        oe[key] = 0

# 3. Adjust carbons for skin environment
CARBONS_OE = {
    'glc_D': -0.1,        # low skin glucose
    'lac_L': -0.5,        # sweat lactate
    'ac':    -0.2,        # skin metabolism
    'acgam': -0.05,       # minimal residual
    'gam':   -0.05,
    'man':   -0.05,
    'gal':   -0.05,
    'gly':   -0.1,        # sweat amino acid
    'ala_L': -0.1,        # sweat amino acid
}
for met, b in CARBONS_OE.items():
    oe[f'EX_u_{met}'] = b

# 4. Iron — cerumen lactoferrin → moderate restriction
for met in ['fe2', 'fe3', 'fe3dcit']:
    oe[f'EX_u_{met}'] = -0.01

# 5. Oxygen — aerobic
oe['EX_u_o2'] = -8.0

# 6. NH4 — sweat urea/ammonia
oe['EX_u_nh4'] = -1.5

# 7. AA bound — skin < mucosa (don't override gly, ala_L set above)
AA_LIST = ['ala_L', 'arg_L', 'asn_L', 'asp_L', 'cys_L', 'gln_L', 'glu_L', 'gly',
           'his_L', 'ile_L', 'leu_L', 'lys_L', 'met_L', 'phe_L', 'pro_L',
           'ser_L', 'thr_L', 'trp_L', 'tyr_L', 'val_L']
for aa in AA_LIST:
    key = f'EX_u_{aa}'
    if key in oe and aa not in ('gly', 'ala_L'):
        oe[key] = -0.1

# 8. Convert to AGORA namespace (EX_<met>(e))
def to_agora(d):
    return {k.replace('EX_u_', 'EX_') + '(e)': v for k, v in d.items()}

# Save
with open(DATA / 'oe_lumen.pckl', 'wb') as f:
    pickle.dump(oe, f)
with open(DATA / 'oe_agora.pckl', 'wb') as f:
    pickle.dump(to_agora(oe), f)

print(f'Saved OE media:')
print(f'  data/oe_lumen.pckl  ({len(oe)} entries)')
print(f'  data/oe_agora.pckl  ({len(oe)} entries)')
print()
print('Key value comparison (CSOM vs OE):')
with open(DATA / 'csom_lumen.pckl', 'rb') as f:
    csom = pickle.load(f)
keys = ['glc_D', 'lac_L', 'ac', 'acgam', 'acnam', 'fuc_L', 'succ', 'ppa',
        'gly', 'ala_L', 'fe2', 'fe3', 'o2', 'nh4']
print(f'  {"metabolite":<12s}  {"CSOM":>8s}  {"OE":>8s}')
print(f'  {"-"*12}  {"-"*8}  {"-"*8}')
for k in keys:
    key = f'EX_u_{k}'
    c = csom.get(key, '—')
    o = oe.get(key, '—')
    c_str = f'{c:>8.3f}' if isinstance(c, (int, float)) else f'{c:>8s}'
    o_str = f'{o:>8.3f}' if isinstance(o, (int, float)) else f'{o:>8s}'
    marker = '←' if c != o else ' '
    print(f'  {k:<12s}  {c_str}  {o_str}  {marker}')
