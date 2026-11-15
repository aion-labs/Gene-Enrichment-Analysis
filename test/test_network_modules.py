import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "code"))
from network_modules import analyze_network_modules


def _toy():
    # two groups of three terms (three collections each) sharing genes; gene X bridges the groups
    a = ["A1", "A2", "A3", "A4"]
    b = ["B1", "B2", "B3", "B4"]
    return [dict(library=f"L{i}", term=f"TA{i}", genes=a + (["X"] if i == 1 else [])) for i in range(1, 4)] + \
           [dict(library=f"L{i}", term=f"TB{i}", genes=b + (["X"] if i == 2 else [])) for i in range(1, 4)]


def test_two_communities_and_bridge():
    m = analyze_network_modules(_toy(), min_terms=6, n_seeds=10)
    comm = dict(zip(m.terms.Term, m.terms.Community))
    assert len({comm[f"TA{i}"] for i in range(1, 4)}) == 1
    assert len({comm[f"TB{i}"] for i in range(1, 4)}) == 1
    assert comm["TA1"] != comm["TB1"]
    g = m.genes.set_index("Gene")
    assert g.loc["X", "Bridge"] and g.loc["X", "Degree"] == 2 and abs(g.loc["X", "Participation"] - 0.5) < 1e-9
    assert g.loc["A1", "Degree"] == 3 and not g.loc["A1", "Bridge"] and g.loc["A1", "Participation"] == 0
    assert m.summary["n_communities"] == 2 and m.summary["partition_stability_mean_ARI"] == 1.0


def test_deterministic():
    a = analyze_network_modules(_toy(), n_seeds=5)
    b = analyze_network_modules(_toy(), n_seeds=5)
    assert a.terms.equals(b.terms) and a.genes.equals(b.genes)


def test_small_network_uses_components():
    m = analyze_network_modules([_toy()[i] for i in (1, 2, 3, 5)], min_terms=6)
    assert m.summary["method"].startswith("connected components")
    assert m.summary["n_communities"] == 2
