import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, time, random, logging
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fast_igea import Collection
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from iter_enrichment import IterativeEnrichment
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt", skip_validation=True)
F = ["h.all","c2.cp.biocarta","c2.cp.kegg_medicus","c2.cp.pid","c2.cp.reactome","c2.cp.wikipathways","c2.cp.kegg_legacy","c5.go.bp","c5.go.cc","c5.go.mf"]
t = time.time(); fast = {f: Collection(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", bg.genes) for f in F}; print("build", round(time.time() - t, 1), "s")
libs = {f: GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt") for f in F}
hiv = [l.strip() for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()]
bgl = sorted(bg.genes); random.seed(7)
lists = [("HIV", hiv)] + [(f"rand{s}_{i}", random.sample(bgl, s)) for s in (300, 1000) for i in range(int(sys.argv[1]))]
tf = to = 0; mism = []; maxrel = 0
for name, g in lists:
    gs = GeneSet(g, bg.genes, name, hgcn=False, format=False)
    for f in F:
        t0 = time.time(); a = fast[f].igea(g); tf += time.time() - t0
        t0 = time.time(); b = IterativeEnrichment(gs, libs[f], bg, 10, 600, "Fisher's Exact Test", p_threshold=0.01, min_overlap=3, max_iterations=30, use_multiprocessing=False, save_iteration_results=False).results; to += time.time() - t0
        sa = [(r[1], sorted(r[3])) for r in a]; sb = [(r["Term_raw"], sorted(r["Genes removed for next iteration"])) for r in b]
        if sa != sb: mism.append((name, f, len(sa), len(sb)))
        for r, q in zip(a, b): maxrel = max(maxrel, abs(r[2] - q["iteration p-value"]) / q["iteration p-value"])
print(f"lists: {len(lists)}  collections x lists: {len(lists)*len(F)}  mismatched paths: {len(mism)} {mism[:5]}")
print(f"max relative p difference: {maxrel:.2e}")
print(f"time current code: {to:.1f}s   fast engine: {tf:.2f}s   speed-up: {to/tf:.0f}x")
