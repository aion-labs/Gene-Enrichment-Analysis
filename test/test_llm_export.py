import os, sys, subprocess, shutil
import pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))
from network_modules import analyze_network_modules
from llm_export import export_for_llm, network_stats
from test_network_modules import _toy


def _stats(r=0, p=1e-4):
    table = pd.DataFrame([dict(Metric=m, Observed=10, **{"Null Mean": 2.0, "Random >= Observed": r, "Empirical p": p})
                          for m in ["Cluster Size", "Number of Genes", "Number of Terms", "Number of Edges",
                                    "Average Edges per Gene"]])
    return network_stats(table, dict(n_lists=1000, q_threshold=0.1, expected_chance_terms=3.2),
                         dict(n_terms_total=8, n_terms_kept=6), 50)


def _terms():
    return [dict(t, iteration=1, p=2.5e-5, q=0.01) for t in _toy()]


def test_dot_structure_and_prompt():
    terms = _terms()
    mods = analyze_network_modules(terms, n_seeds=5)
    dot, prompt = export_for_llm(terms, mods, _stats())
    assert dot.count("subgraph cluster_C") == 2
    assert '"X" [type=gene, degree=2, community=' in dot and "bridge=true" in dot
    assert 'p="2.5e-05"' in dot and 'network_p="<0.001"' in dot
    # bridge gene X is declared outside every cluster
    tail = dot[dot.rindex("  }"):]
    assert '"X" [type=gene' in tail
    assert "STATISTICAL CONTEXT" in prompt and "COMMUNITIES" in prompt and dot in prompt
    assert "no random list reached this value" in prompt
    if shutil.which("dot"):
        assert subprocess.run(["dot", "-Tsvg", "-Kfdp"], input=dot.encode(), capture_output=True).returncode == 0


def test_without_statistics_or_communities():
    dot, prompt = export_for_llm(_terms())
    assert "subgraph cluster" not in dot
    assert "No statistics were computed" in prompt and "COMMUNITIES" not in prompt


def test_not_significant_instruction():
    terms = _terms()
    _, prompt = export_for_llm(terms, analyze_network_modules(terms, n_seeds=3), _stats(r=400, p=0.4))
    assert "p = 0.4" in prompt and "not distinguishable from what random gene lists produce" in prompt and "The biological story" in prompt
