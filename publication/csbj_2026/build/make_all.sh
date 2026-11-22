#!/usr/bin/env bash
# Rebuild every supplementary item and figure of the CSBJ 2026 paper from the archived analysis outputs,
# then compare the rebuilt text with the files in ../supplement (the versions submitted with the paper).
# Usage: bash make_all.sh [<output folder>]          (default: ../rebuilt, ignored by git)
# Needs: Python >= 3.12 with the repository requirements + markdown, leidenalg, adjustText (see ../environment.txt);
#        Graphviz; for Figure 2 rsvg-convert; for PDFs Node.js with Playwright. Steps whose tools are missing are skipped.
set -uo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; PUB="$(dirname "$HERE")"; A="$PUB/analyses"; S="$HERE/supplement"; F="$HERE/figures"
O="$(mkdir -p "${1:-$PUB/rebuilt}" && cd "${1:-$PUB/rebuilt}" && pwd)"; mkdir -p "$O/supplement" "$O/figures" "$O/figure_1_dots"
PY="${PYTHON:-python3}"; export PYTHONHASHSEED=0
run() { echo "== $*"; "$@" || echo "!! FAILED: $*"; }

# Supplementary tables and note (order matters: Table 2 reads Tables 4 and 5)
run $PY "$S/build_supp_table_01.py" "$O/supplement"
run $PY "$S/build_supp_tables_04_05_07_txt.py" "$A/step09_end_to_end/ora_regular" "$A/step09_end_to_end" "$F/figure_1_and_dot_sources" "$O/supplement"
run $PY "$S/build_supp_table_06_txt.py" "$A/step09_end_to_end/scaling_all.tsv" "$O/supplement"
run $PY "$S/build_supp_tables_04_07_html.py" "$O/supplement" "$S/st16.css"
$PY "$S/collection_table_data.py" "$O/supplement" > "$O/collection_table_data.tsv" 2>/dev/null || echo "!! FAILED: collection_table_data.py"
run $PY "$S/build_supp_table_02.py" "$O/collection_table_data.tsv" "$O/supplement"
run $PY "$S/build_supp_table_03.py" "$O/supplement"
run $PY "$S/build_supp_table_08_robustness.py" "$A" "$O/supplement"
GENERATED_DATE=2026-11-20 LLM_DIR=llm_gpt-6.2-sol run $PY "$S/build_llm_supp.py" "$A/step09_end_to_end" "$O/supplement"
run $PY "$S/build_supp_table_14_three_screens.py" "$A" "$O/supplement"
run $PY "$S/build_supp_note_01.py" "$A" "$O/supplement"
# written by hand (no builder): Supplementary Table 1 text version, figure legends, index; Table 13 is the DAVID report of the first submission
for f in "Supp Table 01 - iGEA pseudocode.txt" "Supp Figures - legends.txt" "Supplementary Material - index.txt"; do cp "$PUB/supplement/$f" "$O/supplement/"; done

# Figures
run $PY "$F/figure_1_and_dot_sources/make_dots.py" "$O/figure_1_dots"     # Figure 1 was drawn by the author from Fig1A/Fig1B DOT files
if command -v rsvg-convert >/dev/null && command -v neato >/dev/null; then run bash "$F/figure_2/make_figure_2.sh" "$O/figures"; else echo "-- skipped Figure 2 (needs rsvg-convert and graphviz)"; fi
( cd "$O/figures" && run $PY "$F/figure_3/plot_figure_3.py" && run $PY "$F/figure_4/plot_figure_4.py" )
for s in supp_figure_01 supp_figures_02_03 supp_figure_04 supp_figure_05 supp_figure_06; do run $PY "$F/supp_figures/$s.py" "$O/figures"; done

# PDFs of the supplementary tables
if node -e "require('playwright')" 2>/dev/null; then run node "$HERE/render_pdfs.js" "$O/supplement"; else echo "-- skipped PDFs (needs Node.js with Playwright)"; fi

# Compare with the submitted versions (text content)
echo "== comparison with ../supplement"; $PY "$HERE/compare_outputs.py" "$O/supplement" "$PUB/supplement"
