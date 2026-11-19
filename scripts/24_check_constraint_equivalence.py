"""Solo(agora diet)와 joint(lumen diet)의 영양 제약이 동등한지 검증한다.
agora 네임스페이스 'EX_<met>(e)'와 lumen 네임스페이스 'EX_u_<met>'를
대사물 이름으로 정규화한 뒤 상한/하한 값을 비교한다."""
import pickle, re
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA, RESULTS = ROOT/'data', ROOT/'results'

PAIRS = {'CSOM': ('csom_agora.pckl', 'csom_lumen.pckl'),
         'OE'  : ('oe_agora.pckl',   'oe_lumen.pckl')}

def norm_agora(k):   # EX_fuc_L(e) -> fuc_L
    return re.sub(r'^EX_', '', re.sub(r'\(e\)$', '', k))

def norm_lumen(k):   # EX_u_fuc_L -> fuc_L
    return re.sub(r'^EX_u_', '', k)

rows = []
for env, (fa, fl) in PAIRS.items():
    with open(DATA/fa, 'rb') as f: a = pickle.load(f)
    with open(DATA/fl, 'rb') as f: l = pickle.load(f)
    an = {norm_agora(k): v for k, v in a.items()}
    ln = {norm_lumen(k): v for k, v in l.items()}
    if len(an) != len(a) or len(ln) != len(l):
        print(f'[경고] {env}: 정규화 후 키가 줄었다 (충돌). '
              f'agora {len(a)}->{len(an)}, lumen {len(l)}->{len(ln)}')
    for k in sorted(set(an) | set(ln)):
        va, vl = an.get(k), ln.get(k)
        if va != vl:
            rows.append(dict(env=env, metabolite=k, agora=va, lumen=vl,
                             reason='agora_only' if vl is None
                                    else 'lumen_only' if va is None
                                    else 'value_differs'))

df = pd.DataFrame(rows)
RESULTS.mkdir(exist_ok=True)
df.to_csv(RESULTS/'constraint_equivalence.csv', index=False)
print(f'불일치 {len(df)}건')
if len(df):
    print(df.groupby(['env','reason']).size())
    print(df.to_string(index=False))
else:
    print('두 네임스페이스의 섭취 제약이 모든 대사물에서 동일하다.')
