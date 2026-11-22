#!/usr/bin/env bash
# Figure 2: force-directed layout (Graphviz neato) of the iGEA network of the three GO collections.
# Input: src.dot (written by ../figure_1_and_dot_sources/make_dots.py; Fig2_iGEA_GO_network_q0.10.dot).
# Requires: graphviz (neato), rsvg-convert, Liberation Serif font, python3.
# Layout chosen by score.py (fewest edges crossing term boxes, compact area): L=0.25, font 8.6, seed 25, mode KK.
# Usage: bash make_figure_2.sh [<output folder>]   (default: current folder; intermediate files go to a temporary folder)
set -euo pipefail
OUT="$(cd "${1:-.}" && pwd)"; HERE="$(cd "$(dirname "$0")" && pwd)"; W="$(mktemp -d)"
cp "$HERE"/build2.py "$HERE"/final.py "$HERE"/score.py "$HERE"/src.dot "$W"/ && cd "$W"
python3 final.py 0.25 8.6 25 KK          # writes f3.svg (boxes with 2 pt corners, bold serif labels, colour key)
rsvg-convert -w 4488 f3.svg -b white -o "$OUT/Figure 2 - iGEA GO network.png"                 # 600 dpi at 7.48 in
rsvg-convert -f pdf -w 718.1 f3.svg -o "$OUT/Figure 2 - iGEA GO network.pdf"                 # 7.48 in wide
rm -rf "$W"
