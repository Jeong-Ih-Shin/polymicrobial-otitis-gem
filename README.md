# polymicrobial-otitis-gem

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20469574.svg)](https://doi.org/10.5281/zenodo.20469574)

In silico modeling of *Candida parapsilosis* interactions with *Pseudomonas aeruginosa*
and *Staphylococcus aureus* in chronic otitis environments (CSOM and OE), using
genome-scale metabolic models.

## Overview

This repository contains the genome-scale metabolic modeling pipeline used to predict
cross-feeding patterns and interaction outcomes (mutualism / commensalism / parasitism /
amensalism / neutralism) in two clinically distinct otologic environments:

- **CSOM** (chronic suppurative otitis media) — chronically infected middle ear
- **OE** (otitis externa) — chronically inflamed external auditory canal

For each environment, two clinically observed pairings are simulated:

- CP + PA (*C. parapsilosis* + *P. aeruginosa*)
- CP + SA (*C. parapsilosis* + *S. aureus*)

### Main finding

Both pairs are mutualistic under CSOM-like conditions, whereas the CP + PA pair shifts
to parasitism under OE-like conditions. Systematic attribution across all 37 exchange
bounds that differ between the two media identified **inorganic iron availability** as
the component required for this shift. Oxygen and mucin-associated substrates, which
distinguish the two anatomical environments conspicuously, did not alter the interaction
class in either direction.

Two points of methodological caution are built into the pipeline and are worth noting
before reusing it:

- **Boost factors are ratios.** Along the oxygen axis the boost factor rose steadily
  while joint growth rates stayed constant to fourteen decimal places; the apparent
  coupling was entirely a change in the solo denominator. Absolute solo and joint
  growth rates are therefore recorded alongside every ratio.
- **Pinning a biomass flux can manufacture dependencies.** Constraints that hold an
  organism at a fixed growth rate while an intervention is applied produce apparent
  requirements that disappear once the organism is allowed to slow down (see
  `16_knockout_regimes.py`).

## Reproduce

### Requirements

- Python ≥ 3.9
- Linear programming solver (GLPK ships with cobrapy; gurobi/cplex optional for speed)

See `requirements.txt` for Python dependencies.

### Setup

```bash
git clone <repository-url>
cd polymicrobial-otitis-gem
pip install -r requirements.txt
bash scripts/setup.sh          # checks for external model files
```

### Pipeline

Scripts are numbered in dependency order. Steps 00 build the model and media inputs;
01–04 produce the core simulation; 05–17 are the analyses added for the revised
manuscript; 18–23 render the figures and tables.

```bash
# 00 — inputs (run once)
python scripts/00_preprocess_cp.py      # iDC1003 + 4 transport reactions, CP id mapping
python scripts/00_build_om_media.py     # otitis media scenarios, incl. CSOM
python scripts/00_build_oe_media.py     # OE medium, derived from the CSOM scaffold

# 01-04 — core simulation
python scripts/01_simulate.py           # solo and joint biomass
python scripts/02_cross_feeding.py      # exchange direction by pFBA
python scripts/03_growth_curves.py      # Figure S1
python scripts/04_mechanism_summary.py  # per-condition summary

# 05-06 — is oxygen responsible? (negative result)
python scripts/05_oxygen_scan.py
python scripts/06_phase_diagram.py

# 10, 12 — what is responsible? (iron)
python scripts/10_component_attribution.py
python scripts/12_axis_scan.py

# 08, 11, 14 — what happens if the fungus is removed?
python scripts/08_antifungal_simulation.py
python scripts/11_fungistatic_scaling.py
python scripts/14_fungal_activity_budget.py

# 15-17 — robustness of the knockout and cross-feeding results
python scripts/15_link_knockout.py
python scripts/16_knockout_regimes.py
python scripts/17_crossfeeding_regimes.py

# 18-23 — manuscript figures and tables
python scripts/21_plot_figure1.py       # Figure 1  framework schematic
python scripts/22_plot_figure2.py       # Figure 2  interaction heatmap
python scripts/13b_plot_figure3_iron.py # Figure 3  iron dependence
python scripts/19_plot_figure4_removal.py  # Figure 4  fungal removal
python scripts/20_plot_figureS2.py      # Figure S2 oxygen and mucin
python scripts/18_make_tables.py        # Tables 4 and S2
python scripts/23_make_table1.py        # Table 1
```

`scripts/superseded/` holds earlier plotting scripts that have been replaced; they are
kept for provenance and are not part of the current pipeline.

### Three implementations of partial fungal suppression

`08`, `11` and `14` each suppress the fungus differently, and the first two do not
constrain the network in practice — fixing biomass at zero leaves fungal metabolism
active, and the nominal reaction bounds (±1000) exceed the fluxes actually carried
(order 1–10). All three converge at complete removal, which is what validates the
implementation. Reported values come from `14`.

## Outputs

| File | Description |
|------|-------------|
| `results/pair_results.csv` | Solo and joint biomass, boost factors, interaction class |
| `results/cross_feeding.csv` | Exchange direction per metabolite (pFBA) |
| `results/component_attribution.csv` | Sufficiency and necessity per component group |
| `results/axis_scan_iron.csv` | Iron gradient scan with absolute growth rates |
| `results/fungal_activity_budget.csv` | Bacterial growth vs remaining fungal activity |
| `results/Table1_media_components.{csv,md}` | Components differing between the media |
| `results/Table5_component_attribution.{csv,md}` | Attribution table |
| `results/TableS2_absolute_growth.{csv,md}` | Absolute rates behind every boost factor |
| `figures/Figure1–4.{png,pdf}` | Main manuscript figures |
| `figures/FigureS1–S2.{png,pdf}` | Supplementary figures |

## Repository structure

```
polymicrobial-otitis-gem/
├── data/                  Medium definitions (pickled dictionaries) and CP id mapping
├── models/                Genome-scale metabolic models
│   ├── agora/             Bacterial models (AGORA 1.03)
│   └── cp/                iDC1003, raw and processed
├── src/                   Utility modules
│   ├── cp_helpers.py        — C. parapsilosis model preparation, joint model assembly
│   ├── functions.py         — Medium constraints, interaction classification
│   ├── metabolite_names.py  — Display names
│   └── style.py             — Matplotlib styling (Arial, Type 42 embedding)
├── scripts/               Analysis pipeline, numbered in dependency order
│   └── superseded/        Replaced scripts, kept for provenance
├── figures/               Output figures
└── results/               Output tables
```

## Media definitions

Each medium is a dictionary mapping system-boundary exchange reactions to lower-bound
uptake rates (mmol gDW⁻¹ h⁻¹); negative values permit uptake and 0 does not. Two
namespaces are kept in parallel: `EX_u_*` for the shared lumen and `EX_*(e)` for AGORA
solo simulations.

The CSOM medium was adapted from a host-microbiome modeling framework and is one of a
set representing the otitis media disease spectrum; the acute and effusion conditions
in `00_build_om_media.py` were not analysed here. The OE medium was derived from the
CSOM scaffold. Iron was assigned in graded fashion across the set, most restricted in
CSOM; these are relative rather than measured values, and the predicted role of iron is
conditional on them.

## License

See `LICENSE`.
