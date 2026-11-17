#!/usr/bin/env bash
# Regenerate every paper figure produced by this repository (Figs 4-11).
set -euo pipefail
cd "$(dirname "$0")"
for script in fig4/fig4_*.py fig5/fig5_*.py fig6/fig6_*.py fig7/fig7_*.py \
              fig8/fig8_*.py fig9/fig9_*.py fig10/fig10_*.py fig11/fig11_*.py; do
    echo "== ${script}"
    python "${script}"
done
