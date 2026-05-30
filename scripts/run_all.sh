#!/usr/bin/env bash
# Full pipeline: simulate -> cross-feeding -> growth curves -> mechanism summary
#
# Run from anywhere; the script changes to the repository root.

set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

echo "=== Step 1/4: solo + joint biomass simulation ==="
python scripts/01_simulate.py

echo
echo "=== Step 2/4: cross-feeding directionality ==="
python scripts/02_cross_feeding.py

echo
echo "=== Step 3/4: growth curve figure ==="
python scripts/03_growth_curves.py

echo
echo "=== Step 4/4: mechanism summary table ==="
python scripts/04_mechanism_summary.py

echo
echo "=== Pipeline complete ==="
echo "Outputs:"
echo "  figures/fig_growth_curves.{png,pdf}"
echo "  results/mechanism_summary.csv"
echo
echo "Intermediates (regenerated each run, not tracked):"
echo "  results/pair_results.csv"
echo "  results/cross_feeding.csv"
