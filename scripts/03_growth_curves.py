"""Generate growth curve figure (CSOM/OE × CP+PA/CP+SA, 2×2 grid).

Reads `results/pair_results.csv` (produced by 01_simulate.py) and renders
logistic growth curves for solo and co-culture conditions side-by-side.
Areas where co-culture exceeds solo growth are shaded.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
import style  # noqa: F401 — applies global matplotlib style

RESULTS = ROOT / 'results'
FIGS    = ROOT / 'figures'
FIGS.mkdir(exist_ok=True)

# Colors
CP_COLOR  = '#2A9D8F'
BAC_COLOR = '#264653'

# Grid layout
ENVS  = ['CSOM', 'OE']
PAIRS = ['CP+PA', 'CP+SA']

# Logistic growth parameters
K_CARRYING = 2.0   # carrying capacity (proxy OD)
B0_INITIAL = 0.01  # initial biomass
T_MAX_HR   = 20    # time horizon for simulation


def logistic(t, mu, K=K_CARRYING, B0=B0_INITIAL):
    """Logistic growth: B(t) = K / (1 + (K/B0 - 1) * exp(-μt))."""
    return K / (1 + (K / B0 - 1) * np.exp(-mu * t))


def main():
    df_pair = pd.read_csv(RESULTS / 'pair_results.csv')
    t = np.linspace(0, T_MAX_HR, 400)

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    for i, pair_name in enumerate(PAIRS):
        for j, env in enumerate(ENVS):
            ax = axes[i, j]
            row = df_pair[(df_pair.env == env) & (df_pair.pair == pair_name)]
            if len(row) == 0:
                ax.axis('off')
                continue
            row = row.iloc[0]

            # Panel label (A, B, C, D)
            ax.text(-0.12, 1.08, 'abcd'[i * 2 + j],
                    transform=ax.transAxes,
                    fontsize=18, fontweight='bold', va='top', ha='left')

            bac_label = row['bac_label']

            # Compute trajectories
            cp_solo   = logistic(t, row['cp_solo'])
            cp_joint  = logistic(t, row['cp_free'])
            bac_solo  = logistic(t, row['bac_solo'])
            bac_joint = logistic(t, row['bac_free'])

            # Plot solo (dashed) and joint (solid)
            ax.plot(t, cp_solo,   '--', color=CP_COLOR,  alpha=0.5, lw=2,
                    label='CP (solo)')
            ax.plot(t, cp_joint,  '-',  color=CP_COLOR,  alpha=1.0, lw=2.8,
                    label='CP (co-culture)')
            ax.plot(t, bac_solo,  '--', color=BAC_COLOR, alpha=0.5, lw=2,
                    label=f'{bac_label} (solo)')
            ax.plot(t, bac_joint, '-',  color=BAC_COLOR, alpha=1.0, lw=2.8,
                    label=f'{bac_label} (co-culture)')

            # Highlight co-culture boost regions
            if row['cp_boost'] > 1.05:
                ax.fill_between(t, cp_solo, cp_joint,
                                where=(cp_joint > cp_solo),
                                color=CP_COLOR, alpha=0.1)
            if row['bac_boost'] > 1.05:
                ax.fill_between(t, bac_solo, bac_joint,
                                where=(bac_joint > bac_solo),
                                color=BAC_COLOR, alpha=0.1)

            ax.set_title(f'{env}: {pair_name}', fontsize=13, fontweight='bold')
            ax.set_xlabel('Time (h)')
            ax.set_ylabel('Biomass (proxy OD)')
            ax.set_xlim(0, T_MAX_HR)
            ax.set_ylim(0, 2.2)
            ax.grid(alpha=0.3)

            leg = ax.legend(loc='lower right', framealpha=0.9,
                            prop={'size': 9, 'weight': 'bold'})
            for txt_obj in leg.get_texts():
                txt_obj.set_color('#1E293B')

            txt = (f"μ boost\n"
                   f"Fungus: {row['cp_boost']:.2f}×\n"
                   f"Bacterium: {row['bac_boost']:.2f}×")
            ax.text(0.5, 2.05, txt, fontsize=9, color='#1E293B',
                    ha='left', va='top', fontweight='bold',
                    bbox=dict(boxstyle='round,pad=0.3', fc='white',
                              ec='#CBD5E1', alpha=0.9))

    plt.tight_layout()
    out_png = FIGS / 'FigureS1.png'
    out_pdf = FIGS / 'FigureS1.pdf'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.savefig(out_pdf, bbox_inches='tight')
    plt.close()
    print(f'Saved: {out_png}')
    print(f'Saved: {out_pdf}')


if __name__ == '__main__':
    main()
