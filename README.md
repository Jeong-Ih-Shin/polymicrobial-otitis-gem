# polymicrobial-otitis-gem

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20469574.svg)](https://doi.org/10.5281/zenodo.20469574)

In silico modeling of *Candida parapsilosis* interactions with *Pseudomonas aeruginosa*
and *Staphylococcus aureus* in chronic otitis environments (CSOM and OE), using
genome-scale metabolic models.

## Overview

This repository contains the genome-scale metabolic modeling pipeline for predicting
microbial cross-feeding patterns and resulting interaction outcomes
(mutualism / parasitism / commensalism / amensalism) in two clinically distinct
otitis environments:

- **CSOM** (Chronic Suppurative Otitis Media) — anaerobic middle ear environment
- **OE** (Otitis Externa) — aerobic external auditory canal environment

For each environment, two clinically observed polymicrobial pairings are simulated:

- CP + PA (*C. parapsilosis* + *P. aeruginosa*)
- CP + SA (*C. parapsilosis* + *S. aureus*)

## Outputs

| File | Description |
|------|-------------|
| `figures/fig_growth_curves.{png,pdf}` | Predicted growth curves for solo vs co-culture (CSOM/OE × CP+PA/CP+SA) |
| `results/mechanism_summary.csv` | Integrated cross-feeding analysis with interaction mode, key axis metabolites, and clinical interpretation per condition |

## Reproduce

### Requirements

- Python ≥ 3.9
- Linear programming solver (GLPK installed by default with cobra; gurobi/cplex optional for speed)

See `requirements.txt` for Python dependencies.

### Setup

```bash
git clone <repository-url>
cd polymicrobial-otitis-gem
pip install -r requirements.txt
```

### Run pipeline

The pipeline consists of four sequential steps:

```bash
# 1. Solo and joint biomass simulation
python scripts/01_simulate.py

# 2. Cross-feeding directionality analysis
python scripts/02_cross_feeding.py

# 3. Generate growth curve figure
python scripts/03_growth_curves.py

# 4. Generate mechanism summary table
python scripts/04_mechanism_summary.py
```

Or run all at once:

```bash
bash scripts/run_all.sh
```

### Pipeline data flow

```
scripts/01_simulate.py          →  results/pair_results.csv     (intermediate)
scripts/02_cross_feeding.py     →  results/cross_feeding.csv    (intermediate)
scripts/03_growth_curves.py     →  figures/fig_growth_curves.{png,pdf}
scripts/04_mechanism_summary.py →  results/mechanism_summary.csv
```

Intermediate CSVs (`pair_results.csv`, `cross_feeding.csv`) are excluded from version
control but are regenerated locally when the pipeline runs.

## Repository structure

```
polymicrobial-otitis-gem/
├── data/                  Input pickles (lumen exchange profiles for CSOM/OE diets)
├── models/                Genome-scale metabolic models
│   ├── agora/             Bacterial models (AGORA universe)
│   └── cp/                Fungal model (iDC1003-derived processed)
├── src/                   Utility modules
│   ├── cp_helpers.py        — C. parapsilosis model preparation
│   ├── functions.py         — Joint model construction, interaction classification
│   ├── metabolite_names.py  — Display name mapping for metabolites
│   └── style.py             — Matplotlib styling
├── scripts/               Sequential analysis pipeline
├── figures/               Output figures
└── results/               Output tables
```

## Modeling approach

1. **Solo and joint simulation** (`01_simulate.py`)
   Loads AGORA bacterial models (PA, SA) and the processed *C. parapsilosis* model.
   Constrains each model with the appropriate lumen diet (CSOM anaerobic / OE aerobic)
   and computes biomass for both monoculture and joint co-culture via FBA.
   The biomass ratio (joint vs solo) provides the growth-coupling factor used for
   interaction-mode classification.

2. **Cross-feeding analysis** (`02_cross_feeding.py`)
   Decomposes the joint-model lumen exchange flux into directional pairs (bacterium →
   fungus, fungus → bacterium), identifying which metabolites are produced by one
   organism and consumed by the other.

3. **Growth curve visualization** (`03_growth_curves.py`)
   Translates each computed growth rate into a logistic growth curve
   (μ, carrying capacity K = 2.0, initial biomass B₀ = 0.01) and displays solo
   vs co-culture trajectories side-by-side. Co-culture areas exceeding solo are
   shaded.

4. **Mechanism summary** (`04_mechanism_summary.py`)
   Combines the cross-feeding directionality with the interaction mode classification
   (Mutualism, Parasitism, Commensalism, Amensalism, Neutralism) into a per-condition
   summary table. Each row includes the key metabolic axis, mechanistic interpretation,
   and a clinical relevance statement.

## License

See `LICENSE`.
