"""Generate `results/mechanism_summary.csv`.

Integrates cross-feeding directionality with interaction-mode classification
(based on biomass coupling factors) into a per-condition table. Each row
contains the interaction mode, key metabolic axis, mechanistic interpretation,
and a clinical relevance statement.

Inputs:
    results/pair_results.csv   (from 01_simulate.py)
    results/cross_feeding.csv  (from 02_cross_feeding.py)

Output:
    results/mechanism_summary.csv
"""
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'src'))
from metabolite_names import label
from functions import interaction_type

RESULTS = ROOT / 'results'

# ----------------------------------------------------------
# Per-condition mechanistic interpretation
# (manually curated from cross-feeding directionality + biology)
# ----------------------------------------------------------
INTERPRETATION = {
    ('CSOM', 'CP+PA'): {
        'Key_axis': 'Fermentation-respiration coupling',
        'Mechanism': (
            "PA's formate/succinate to CP carbon; "
            "CP's O2/serine/fumarate to PA respiration"
        ),
        'Clinical_message': (
            'In the CSOM anaerobic environment, CP serves as an O2 and amino acid '
            'donor while PA recycles fermentation byproducts -- a true syntrophy. '
            'Single-target therapy may be insufficient.'
        ),
        'Key_metabolites': {'for', 'succ', 'o2', 'fum'},
    },
    ('CSOM', 'CP+SA'): {
        'Key_axis': 'Nitrogen cycling',
        'Mechanism': "SA's NH4 to CP amino acid synthesis; CP's serine to SA",
        'Clinical_message': (
            'SA-CP co-infection is sustained by nitrogen cycling -- SA NH4 '
            'secretion fuels CP amino acid synthesis. Worth considering fungal '
            'identification even when only SA is cultured.'
        ),
        'Key_metabolites': {'nh4', 'ser_L'},
    },
    ('OE', 'CP+PA'): {
        'Key_axis': 'Aerobic carbon monopoly',
        'Mechanism': (
            "PA exploits OE's O2; CP forced to donate fumarate/serine/O2 "
            "at biosynthetic cost"
        ),
        'Clinical_message': (
            'PA thrives in the OE aerobic environment even solo, while CP is '
            'suppressed by donating carbon, oxygen, and amino acids. Explains '
            'the clinical PA-dominance in OE through metabolic parasitism '
            '(PA exploits CP-derived substrates while suppressing CP growth). '
            'PA-targeted therapy is priority.'
        ),
        'Key_metabolites': {'o2', 'fum', 'ser_L'},
    },
    ('OE', 'CP+SA'): {
        'Key_axis': 'Nucleoside salvage',
        'Mechanism': (
            'dcyt <-> deoxyuridine exchange -- DNA precursor swap reduces '
            'de novo synthesis cost'
        ),
        'Clinical_message': (
            'In the nutrient-limited OE environment, SA-CP exploit nucleoside '
            'salvage (DNA precursor exchange) as a cooperative axis. Saves ATP '
            'cost of de novo synthesis, enabling rapid growth of both -- '
            'a chronicity risk.'
        ),
        'Key_metabolites': {'dcyt', 'duri'},
    },
}


def main():
    cf = pd.read_csv(RESULTS / 'cross_feeding.csv')
    cf = cf[cf.direction.isin(['fungus_to_bac', 'bac_to_fungus'])].copy()
    pair_df = pd.read_csv(RESULTS / 'pair_results.csv')

    summary_rows = []
    for env in ['CSOM', 'OE']:
        for pair in ['CP+PA', 'CP+SA']:
            sub = cf[(cf.env == env) & (cf.pair == pair)]
            info = pair_df[(pair_df.env == env) & (pair_df.pair == pair)].iloc[0]
            top_b2f = sub[sub.direction == 'bac_to_fungus'].nlargest(3, 'mag')
            top_f2b = sub[sub.direction == 'fungus_to_bac'].nlargest(3, 'mag')

            # Coupling factors (log2 of solo→joint biomass ratio)
            ccf = (math.log2(info.cp_boost)  if info.cp_boost  > 0
                   else float('-inf'))
            bcf = (math.log2(info.bac_boost) if info.bac_boost > 0
                   else float('-inf'))
            mode, _ = interaction_type(ccf, bcf, math)

            interp = INTERPRETATION[(env, pair)]

            summary_rows.append({
                'Environment'      : env,
                'Pair'             : pair,
                'CP_boost'         : f"{info.cp_boost:.2f}x",
                'Bac_boost'        : f"{info.bac_boost:.2f}x",
                'Mode'             : mode,
                'Key_axis'         : interp['Key_axis'],
                'Mechanism'        : interp['Mechanism'],
                'Clinical_message' : interp['Clinical_message'],
                'Key_metabolites'  : ', '.join(
                    sorted(label(m) for m in interp['Key_metabolites'])),
                'Top_Bac_to_Fungus': '; '.join(
                    [f"{label(r.metabolite)} ({r.mag:.1f})"
                     for _, r in top_b2f.iterrows()]),
                'Top_Fungus_to_Bac': '; '.join(
                    [f"{label(r.metabolite)} ({r.mag:.1f})"
                     for _, r in top_f2b.iterrows()]),
            })

    summary_df = pd.DataFrame(summary_rows)
    out_path = RESULTS / 'mechanism_summary.csv'
    summary_df.to_csv(out_path, index=False)
    print(f'Saved: {out_path}')


if __name__ == '__main__':
    main()
