#!/bin/bash
# Final LLM run (prompt v5): first answer for every run folder, then a second turn revealing the experiment
# (not for the random controls). Responses go to <run folder>/llm/. Finished calls are skipped, so the
# script can be rerun after an interruption.
# Usage (Terminal):  export OPENAI_API_KEY="sk-..."   then   bash run_llm_test_hiv.sh [model] [answer folder]
# The answers of each model go to their own subfolder of every run folder (default: llm_<model>; the
# gpt-6-astra answers of 20 Nov are in llm/).
MODEL="${1:-gpt-6-astra}"
OUT="${2:-llm_$MODEL}"
C="$HOME/Library/Mobile Documents/com~apple~CloudDocs/AION Labs/15 Coding Projects/Projects/iterative_enrichment/code/llm_client.py"
cd "$(dirname "$0")" || exit 1

HIV="an siRNA screen in human cells for host factors required for HIV infection (HIV dependency factors)"
HYPOXIA="a meta analysis of genes that are consistently regulated by mRNA expression when cells are exposed to hypoxia"
VACCINATION="genes whose expression in blood 3 days after seasonal influenza vaccination correlated positively with the antibody response (hemagglutination-inhibition titer) at day 28, from a human systems-vaccinology study"

first() {  # folder pattern, repeats
  for d in $1; do
    [ -f "$d/llm_prompt.txt" ] || continue
    echo "== first answer: $d ($2 x $MODEL)"
    /usr/bin/python3 "$C" "$d" --model "$MODEL" --repeats "$2" --out-name "$OUT" || exit 1
  done
}
second() {  # folder pattern, context
  for d in $1; do
    [ -f "$d/llm_prompt.txt" ] || continue
    echo "== second turn: $d"
    /usr/bin/python3 "$C" "$d" --model "$MODEL" --context "$2" --out-name "$OUT" || exit 1
  done
}

# Test of one model on the HIV list only (3 first answers + 3 second turns), before running everything.
first "HIV_dependency_factors.symbols_2*" 3
second "HIV_dependency_factors.symbols_2*" "$HIV"
echo "HIV test done: answers in HIV_dependency_factors.symbols_*/$OUT"
