#!/bin/bash
# Two-turn LLM interpretation (prompt v5) of the iGEA networks of the König and Zhou screens, same model and protocol
# as the Brass network (step09_end_to_end/HIV_dependency_factors.symbols_*/llm_gpt-6.2-sol/): 3 independent first turns
# without experimental information, then the second turn revealing the experiment. The context sentence is identical
# for the three screens so that the interpretations can be compared.
# Usage (Terminal):  export OPENAI_API_KEY="sk-..."   then   bash run_llm_screens.sh [model]
MODEL="${1:-gpt-6.2-sol}"
OUT="llm_$MODEL"
C="$HOME/Library/Mobile Documents/com~apple~CloudDocs/AION Labs/15 Coding Projects/Projects/iterative_enrichment/code/llm_client.py"
cd "$(dirname "$0")" || exit 1
HIV="an siRNA screen in human cells for host factors required for HIV infection (HIV dependency factors)"
for d in Konig_2008_HIV_dependency_factors.symbols_* Zhou_2008_HIV_dependency_factors.symbols_*; do
  [ -f "$d/llm_prompt.txt" ] || continue
  echo "== first answer: $d (3 x $MODEL)"
  /usr/bin/python3 "$C" "$d" --model "$MODEL" --repeats 3 --out-name "$OUT" || exit 1
done
for d in Konig_2008_HIV_dependency_factors.symbols_* Zhou_2008_HIV_dependency_factors.symbols_*; do
  [ -f "$d/llm_prompt.txt" ] || continue
  echo "== second turn: $d"
  /usr/bin/python3 "$C" "$d" --model "$MODEL" --context "$HIV" --out-name "$OUT" || exit 1
done
echo "All done."
