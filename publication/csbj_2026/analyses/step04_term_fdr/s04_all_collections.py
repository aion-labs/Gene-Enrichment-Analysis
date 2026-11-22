"""HIV, all 10 collections: ORA (BH FDR < 0.05) vs iGEA (empirical q < 0.05), gain/loss of terms and genes.
Random lists = the per-run null lists (default_rng([2026, 272, i]), sorted validated background), 10,000.
iGEA q per collection: mean # random iGEA terms with p' <= p / # real iGEA terms with p' <= p, monotone.
Usage: python s04_all_collections.py <out_dir> [N]"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging, itertools
from multiprocessing import get_context
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
P = [("Hallmark","h.all"),("Reactome","c2.cp.reactome"),("KEGG Legacy","c2.cp.kegg_legacy"),("KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("BioCarta","c2.cp.biocarta"),("PID","c2.cp.pid"),("WikiPathways","c2.cp.wikipathways"),("GO BP","c5.go.bp"),("GO CC","c5.go.cc"),("GO MF","c5.go.mf")]
W = {}
def init(cols, B, n): W.update(cols=cols, B=B, n=n)
def one(i):
    g = [W["B"][j] for j in np.random.default_rng([2026, W["n"], i]).choice(len(W["B"]), W["n"], replace=False)]
    return [(c, p) for c, col in enumerate(W["cols"]) for _, _, p, _ in col.iterate(g, 0.01, 3, None)]
def mean_jaccard(sets):
    if len(sets) < 2: return np.nan
    return float(np.mean([len(a & b) / len(a | b) for a, b in itertools.combinations(sets, 2)]))
if __name__ == "__main__":
    OUT = sys.argv[1]; N = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False); B = sorted(set(bg.genes)); H = set(gs.genes)
    cols = [FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=n), B) for n, f in P]
    with get_context("spawn").Pool(os.cpu_count(), initializer=init, initargs=(cols, B, gs.size)) as pool:
        res = pool.map(one, range(N), chunksize=50)
    rnd = pd.DataFrame([x for r in res for x in r], columns=["c", "p"])
    terms, summary = [], []
    for c, ((name, _), col) in enumerate(zip(P, cols)):
        hs = lambda t: {col.genes[j] for j in col.members[t] if col.genes[j] in H}
        # ORA (exact p, BH over all tested sets, as in iGEA's Enrichment)
        x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]
        n = x.sum(); k = col.A @ x
        p_all = np.array([exact_pvalue("Fisher's Exact Test", int(k[t]), int(col.K[t]), int(n), col.M) for t in range(len(col.term_names))])
        fdr = multipletests(p_all, method="fdr_bh")[1]
        ora_fdr = [t for t in range(len(p_all)) if p_all[t] < 0.01 and k[t] >= 3 and fdr[t] < 0.05]
        # iGEA + empirical q
        ig = col.iterate(gs.genes, 0.01, 3, None); tid = {nm: i for i, nm in enumerate(col.term_names)}
        rp = np.sort(rnd.p[rnd.c == c].values); op = np.array([p for _, _, p, _ in ig])
        qraw = [np.searchsorted(rp, p, side="right") / N / (op <= p).sum() for p in op]
        order = np.argsort(op); qs = np.empty(len(op))
        if len(op): qs[order] = np.minimum.accumulate(np.clip(np.array(qraw)[order], 0, 1)[::-1])[::-1]
        ig_keep = [(it, tid[t], p, set(g), q) for (it, t, p, g), q in zip(ig, qs) if q < 0.05]
        ora_set = {t for t in ora_fdr}; ig_set = {t for _, t, _, _, _ in ig_keep}
        g_ora = set().union(*[hs(t) for t in ora_fdr]) if ora_fdr else set()
        g_ig = set().union(*[g for *_, g, _ in ig_keep]) if ig_keep else set()
        g_ig_all = set().union(*[set(g) for *_, g in ig]) if ig else set()
        lost = [col.term_names[t] for t in ora_fdr if len(hs(t) & g_ig) < 0.5 * len(hs(t))]
        summary.append({"Collection": name, "ORA FDR<0.05": len(ora_set), "iGEA terms": len(ig), "iGEA q<0.05": len(ig_set),
            "Both": len(ora_set & ig_set), "iGEA q<0.05 only": len(ig_set - ora_set), "ORA only": len(ora_set - ig_set),
            "Genes ORA FDR": len(g_ora), "Genes iGEA q<0.05": len(g_ig), "Genes only ORA": len(g_ora - g_ig), "Genes only iGEA": len(g_ig - g_ora),
            "ORA terms mostly uncovered (<50% genes)": len(lost),
            "Mean Jaccard ORA FDR": round(mean_jaccard([hs(t) for t in ora_fdr]), 3),
            "Mean Jaccard iGEA q<0.05": round(mean_jaccard([g for *_, g, _ in ig_keep]), 3),
            "Random iGEA terms per list": round((rnd.c == c).sum() / N, 2)})
        for it, t, p, g, q in [(a, tid[b], cc, set(d), qq) for (a, b, cc, d), qq in zip(ig, qs)]:
            terms.append({"Collection": name, "Source": "iGEA", "Term": col.term_names[t], "Iteration": it, "p": p, "q / FDR": q,
                          "Kept": q < 0.05, "In ORA FDR<0.05": t in ora_set, "Genes": len(g), "Gene list": ", ".join(sorted(g))})
        for t in ora_fdr:
            terms.append({"Collection": name, "Source": "ORA", "Term": col.term_names[t], "Iteration": None, "p": p_all[t], "q / FDR": fdr[t],
                          "Kept": True, "In iGEA q<0.05": t in ig_set, "Genes": len(hs(t)), "Mostly uncovered by iGEA": col.term_names[t] in lost,
                          "Gene list": ", ".join(sorted(hs(t)))})
        globals().setdefault("G", {"ora": set(), "ig": set(), "ig_all": set()})
        G["ora"] |= g_ora; G["ig"] |= g_ig; G["ig_all"] |= g_ig_all
    s = pd.DataFrame(summary); tot = s.drop(columns=["Collection", "Mean Jaccard ORA FDR", "Mean Jaccard iGEA q<0.05", "Random iGEA terms per list"]).sum()
    s.to_csv(f"{OUT}/hiv_all_ora_vs_igea_summary.tsv", sep="\t", index=False); pd.DataFrame(terms).to_csv(f"{OUT}/hiv_all_ora_vs_igea_terms.tsv", sep="\t", index=False)
    pd.set_option("display.width", 300)
    print(s.to_string(index=False)); print("\nColumn sums:", tot.to_dict())
    print(f"Distinct genes (union over collections): ORA FDR {len(G['ora'])}, iGEA q<0.05 {len(G['ig'])}, all iGEA {len(G['ig_all'])}, "
          f"only ORA {len(G['ora'] - G['ig'])}, only iGEA q<0.05 {len(G['ig'] - G['ora'])}, both {len(G['ora'] & G['ig'])}")
