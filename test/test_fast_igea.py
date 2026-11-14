"""Equivalence and reproducibility tests for the fast engine and the per-run null.

Run from the repository root:  python -m pytest test/test_fast_igea.py -q
"""
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))

from background_gene_set import BackgroundGeneSet  # noqa: E402
from fast_igea import FISHER, HYPERGEOMETRIC, FastCollection  # noqa: E402
from gene_set import GeneSet  # noqa: E402
from gene_set_library import GeneSetLibrary  # noqa: E402
from iter_enrichment import IterativeEnrichment  # noqa: E402
from per_run_null import compare_to_null, generate_null  # noqa: E402

LIBS = ["h.all", "c2.cp.biocarta", "c2.cp.kegg_medicus", "c2.cp.pid", "c2.cp.reactome",
        "c2.cp.wikipathways", "c2.cp.kegg_legacy", "c5.go.bp", "c5.go.cc", "c5.go.mf"]


@pytest.fixture(scope="module")
def data():
    bg = BackgroundGeneSet(str(ROOT / "data/backgrounds/all_genes.txt"), skip_validation=True)
    libs = [GeneSetLibrary(str(ROOT / f"data/libraries/{f}.v2025.1.Hs.symbols.gmt"), name=f) for f in LIBS]
    hiv = [l.strip() for l in open(ROOT / "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()]
    return bg, libs, hiv


def _reference(genes, lib, bg, method):
    gs = GeneSet(genes, bg.genes, "t", hgcn=False, format=False)
    res = IterativeEnrichment(gs, lib, bg, 10, 600, method, p_threshold=0.01, min_overlap=3,
                              max_iterations=30, use_multiprocessing=False, save_iteration_results=False).results
    return [(r["Iteration"], r["Term_raw"], r["iteration p-value"], sorted(r["Genes removed for next iteration"])) for r in res]


def _fast(genes, lib, bg, method):
    col = FastCollection.from_library(lib, bg.genes)
    return [(i, t, p, sorted(g)) for i, t, p, g in col.iterate(genes, 0.01, 3, 30, method)]


@pytest.mark.parametrize("method", [FISHER, HYPERGEOMETRIC])
def test_hiv_identical(data, method):
    bg, libs, hiv = data
    for lib in libs:
        assert _fast(hiv, lib, bg, method) == _reference(hiv, lib, bg, method), lib.name


def test_random_lists_identical(data):
    bg, libs, _ = data
    pool = sorted(bg.genes)
    rnd = random.Random(11)
    for size in (300, 600):
        genes = rnd.sample(pool, size)
        for lib in libs:
            assert _fast(genes, lib, bg, FISHER) == _reference(genes, lib, bg, FISHER), (size, lib.name)


def test_null_reproducible_and_parallel_safe(data):
    bg, libs, _ = data
    a = generate_null(libs, bg.genes, 272, n_lists=40, seed=7, n_jobs=1)
    b = generate_null(libs, bg.genes, 272, n_lists=40, seed=7, n_jobs=2)
    c = generate_null(libs, bg.genes, 272, n_lists=40, seed=8, n_jobs=1)
    cols = ["cluster_size", "n_genes", "n_terms", "n_edges", "avg_edges_per_gene"]
    assert a[cols].equals(b[cols])
    assert not a[cols].equals(c[cols])
    t = compare_to_null({"cluster_size": 241, "n_genes": 179, "n_terms": 62, "n_edges": 412,
                         "avg_edges_per_gene": 412 / 179}, a)
    assert ((t["Empirical p"] >= 1 / 41) & (t["Empirical p"] <= 1)).all()


def test_two_layer_statistics(data):
    from per_run_null import run_statistics, signal_adjusted_q, ChanceReference
    bg, libs, hiv = data
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False)
    obs = [(c, t, it, p, g) for c, lib in enumerate(libs)
           for it, t, p, g in FastCollection.from_library(lib, bg.genes).iterate(gs.genes, 0.01, 3, None)]
    a = run_statistics(obs, libs, bg.genes, gs.size, n_lists=200, grid_lists=100, seed=3, n_jobs=1)
    b = run_statistics(obs, libs, bg.genes, gs.size, n_lists=200, grid_lists=100, seed=3, n_jobs=2)
    cols = ["cluster_size", "n_genes", "n_terms", "n_edges", "n_terms_kept"]
    assert a[1][cols].equals(b[1][cols])                       # reproducible, parallel-safe
    assert a[2]["q-value"].tolist() == b[2]["q-value"].tolist()
    terms, observed = a[2], a[3]
    assert observed["n_terms_kept"] == int(terms["Kept"].sum()) <= len(obs)
    assert ((terms["q-value"] >= 0) & (terms["q-value"] <= 1)).all()
    # q is monotone in p within a collection
    for _, g in terms.groupby("Library"):
        g = g.sort_values("p-value")
        assert (g["q-value"].diff().dropna() >= -1e-12).all()
    # random networks are filtered: kept terms per random list well below selected terms
    assert a[1]["n_terms_kept"].mean() < a[1]["n_terms_total"].mean()
    # without the filter all observed terms are kept
    c = run_statistics(obs, libs, bg.genes, gs.size, n_lists=50, seed=3, n_jobs=1, q_threshold=None)
    assert c[3]["n_terms_kept"] == len(obs)
