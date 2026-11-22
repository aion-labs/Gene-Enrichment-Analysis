"""Are terms that appear only in subsamples substitutes for full-list terms? For every term kept in a subsample but
not in the full-list network, the highest gene-set Jaccard with any full-list kept term (same 100 subsamples, seed 2026)."""
import sys, os, logging
from multiprocessing import get_context
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s06_stability as s6
R = s6.R
if __name__ == "__main__":
    OUT = sys.argv[1]
    from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
    from gene_converter import GeneConverter; from fast_igea import FastCollection; import per_run_null as prn
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter(); B = sorted(bg.genes)
    hiv = sorted(GeneSet([m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m], bg.genes, "HIV", hgcn=False, format=False).genes)
    libs = [GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n) for n, f in s6.ALL]; names = [n for n, _ in s6.ALL]
    cols = [FastCollection.from_library(L, B) for L in libs]; n = int(round(0.9 * len(hiv)))
    sets = {t["name"]: set(t["genes"]) for L in libs for t in L.library}
    grid = prn.grid_sizes(n); NREF = 5000; GL = 500
    tasks = [(n, i) for i in range(NREF)] + [(s, i) for s in grid if s != n for i in range(GL)]
    ctx = get_context("spawn")
    with ctx.Pool(os.cpu_count(), initializer=prn._init_gen, initargs=(cols, B, prn.DEFAULT_SEED, s6.PARAMS)) as p:
        gen = p.map(prn._gen, tasks, chunksize=50)
    ref = prn.ChanceReference(len(cols)); ref.add(n, gen[:NREF]); pos = NREF
    for s in grid:
        if s != n: ref.add(s, gen[pos:pos + GL]); pos += GL
    rng = np.random.default_rng(2026); subs = [sorted(rng.choice(hiv, n, replace=False)) for _ in range(100)]
    with ctx.Pool(os.cpu_count(), initializer=s6.init_sub, initargs=(cols, names, ref, n)) as p:
        res = p.map(s6.run_list, subs, chunksize=5)
    full = pd.read_csv(os.path.join(os.path.dirname(OUT), "step05_two_layer/hiv_10k_q0.10_terms.tsv"), sep="\t")
    fk = list(full[full.Kept.astype(bool)].Term)
    best = []
    for r in res:
        for t in r["kept"] - set(fk):
            j = max(len(sets[t] & sets[f]) / len(sets[t] | sets[f]) for f in fk); best.append((t, j))
    b = pd.DataFrame(best, columns=["term", "max_jaccard_with_full_term"]); b.to_csv(f"{OUT}/s06_new_term_similarity.tsv", sep="\t", index=False)
    print("new-term occurrences:", len(b), "| distinct:", b.term.nunique())
    print("max Jaccard with a full-list term: median %.2f; >=0.3: %.0f%%; >=0.5: %.0f%%" % (b.max_jaccard_with_full_term.median(), 100*(b.max_jaccard_with_full_term>=0.3).mean(), 100*(b.max_jaccard_with_full_term>=0.5).mean()))
    print(b.groupby("term").size().sort_values(ascending=False).head(10).to_string())
