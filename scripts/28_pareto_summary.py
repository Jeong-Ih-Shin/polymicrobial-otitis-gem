"""Pareto front 요약과 Figure S3 생성."""
import sys
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT/'src'))
import style  # noqa: F401
RESULTS, FIGS = ROOT/'results', ROOT/'figures'
FIGS.mkdir(exist_ok=True)

LINE   = '#264653'
MARKER = '#E76F51'
MUT    = '#2A9D8F'

df  = pd.read_csv(RESULTS/'pareto_front.csv')
anc = pd.read_csv(RESULTS/'pareto_anchors.csv')
kdf = pd.read_csv(RESULTS/'pareto_kink.csv')

def kink(env, pair):
    """사전식 최적화로 구한 정확한 꺾임점 (30_exact_kink.py)."""
    return kdf[(kdf.env==env) & (kdf.pair==pair)].iloc[0]

print('=== 전면 요약 ===')
for (env, pair), g in df.groupby(['env','pair'], sort=False):
    g = g.dropna(subset=['bac_boost'])
    k = kink(env, pair)
    a = anc[(anc.env==env) & (anc.pair==pair)].iloc[0]
    eq_cp, eq_bac = a.cp_growth/a.cp_solo, a.bac_growth/a.bac_solo
    print(f'\n{env} / {pair}')
    print(f'  전면의 꺾임점   : CP {k.cp_boost:.3f} / bac {k.bac_boost:.3f}')
    print(f'  1:1 자유최적    : CP {eq_cp:.3f} / bac {eq_bac:.3f}')
    print(f'  두 지점 일치 여부: {"일치" if abs(k.cp_boost-eq_cp) < 0.05 else "불일치 — 확인 필요"}')
    print(f'  세균 비용 없이 CP가 도달 가능한 최대: solo의 {k.cp_boost:.2f}배')

ENVS, PAIRS = ['CSOM','OE'], ['CP+PA','CP+SA']
fig, axes = plt.subplots(2, 2, figsize=(9.5, 8.5))

for i, env in enumerate(ENVS):
    for j, pair in enumerate(PAIRS):
        ax = axes[i, j]
        g = df[(df.env==env) & (df.pair==pair)].dropna(subset=['bac_boost'])
        k = kink(env, pair)
        a = anc[(anc.env==env) & (anc.pair==pair)].iloc[0]

        xmax = max(g.cp_boost.max()*1.08, 1.15)
        ymax = g.bac_boost.max()*1.18

        # 둘 다 solo를 넘는 영역
        ax.add_patch(plt.Rectangle((1, 1), xmax-1, ymax-1,
                                   fc=MUT, alpha=0.10, ec='none', zorder=0))
        ax.axvline(1.0, color='grey', lw=0.9, ls=':', zorder=1)
        ax.axhline(1.0, color='grey', lw=0.9, ls=':', zorder=1)

        ax.plot(g.cp_boost, g.bac_boost, '-', color=LINE, lw=2.6, zorder=3)
        ax.plot(k.cp_boost, k.bac_boost, 'o', ms=10, mfc=MARKER,
                mec='white', mew=1.6, zorder=5)
        ax.annotate(f'{k.cp_boost:.2f}', xy=(k.cp_boost, k.bac_boost),
                    xytext=(-38, 9), textcoords='offset points',
                    fontsize=10, fontweight='bold', color=MARKER)

        ax.text(-0.14, 1.06, 'abcd'[i*2+j], transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='left')
        ax.set_title(f'{env}: {pair}', fontsize=12, fontweight='bold')
        ax.set_xlabel('C. parapsilosis boost factor')
        ax.set_ylabel('Bacterial boost factor')
        ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)

handles = [
    Line2D([], [], color=LINE, lw=2.6, label='Trade-off front'),
    Line2D([], [], color=MARKER, marker='o', ls='none', ms=9, mec='white',
           mew=1.4, label='Greatest fungal growth at no cost to the bacterium'),
    Line2D([], [], color='grey', lw=0.9, ls=':', label='Solo growth (boost factor 1.0)'),
    Patch(fc=MUT, alpha=0.10, ec='none', label='Both partners above solo growth'),
]
fig.legend(handles=handles, loc='lower center', ncol=2, frameon=False,
           fontsize=10, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout(rect=[0, 0.07, 1, 1])
plt.savefig(FIGS/'FigureS3_pareto_front.png', dpi=300, bbox_inches='tight')
plt.savefig(FIGS/'FigureS3_pareto_front.pdf', bbox_inches='tight')
print(f'\nSaved: {FIGS}/FigureS3_pareto_front.png')
