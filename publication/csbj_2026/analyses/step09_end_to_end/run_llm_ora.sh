#!/bin/bash
# ORA baseline (R2-4): same prompt on the FDR-filtered ORA GO network (ora_go_fdr0.05/, built by build_ora_llm_input.py).
# Usage (Terminal):  export OPENAI_API_KEY="sk-..."   then   bash run_llm_ora.sh [model] [folder: ora_go_fdr0.05 | ora_all_fdr0.05]
MODEL="${1:-gpt-6.2-sol}"
DIR="${2:-ora_go_fdr0.05}"
C="$HOME/Library/Mobile Documents/com~apple~CloudDocs/AION Labs/15 Coding Projects/Projects/iterative_enrichment/code/llm_client.py"
cd "$(dirname "$0")" || exit 1
/usr/bin/python3 "$C" "$DIR" --model "$MODEL" --repeats 3 --out-name "llm_$MODEL" && echo "All done."
