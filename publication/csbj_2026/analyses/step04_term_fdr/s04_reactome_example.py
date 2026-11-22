"""Worked example of the per-term empirical q-value, HIV x Reactome.
Random lists are the same as in the per-run null (list i = default_rng([2026, 272, i]) on the sorted background)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, logging
import numpy as np, pandas as pd
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection
OUT = sys.argv[1]; N = 10000; SEED = 2026
bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
gs = GeneSet(hiv, bg.genes, "HIV", hgcn=False, format=False)
B = sorted(set(bg.genes))
col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/c2.cp.reactome.v2025.1.Hs.symbols.gmt", name="Reactome"), B)
obs = col.iterate(gs.genes, 0.01, 3, None)
rnd = []
for i in range(N):
    g = [B[j] for j in np.random.default_rng([SEED, gs.size, i]).choice(len(B), gs.size, replace=False)]
    for it, t, p, rem in col.iterate(g, 0.01, 3, None):
        rnd.append((i, it, t, p))
rnd = pd.DataFrame(rnd, columns=["list", "iteration", "term", "p"])
rnd.to_csv(f"{OUT}/hiv_reactome_random_terms.tsv", sep="\t", index=False)
rp = np.sort(rnd.p.values)
rows = []
op = np.array([p for _, _, p, _ in obs])
for it, t, p, rem in obs:
    exp_false = np.searchsorted(rp, p, side="right") / N      # random terms with p' <= p, per list
    n_obs = int((op <= p).sum())                               # real terms with p' <= p
    rows.append({"Iter": it, "Term": t.replace("REACTOME_", "").replace("_", " ").lower(), "Genes": len(rem), "p": p,
                 "Real terms with p' <= p": n_obs, "Random terms with p' <= p (total over 10,000 lists)": int(exp_false * N),
                 "Expected false (per list)": exp_false, "q (raw)": exp_false / n_obs})
d = pd.DataFrame(rows).sort_values("p")
d["q"] = np.minimum.accumulate(d["q (raw)"][::-1].clip(upper=1))[::-1]   # monotone, like BH
d.sort_values("Iter").to_csv(f"{OUT}/hiv_reactome_qvalues.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(d.sort_values("Iter").to_string(index=False, formatters={"p": "{:.2e}".format, "Expected false (per list)": "{:.4f}".format, "q (raw)": "{:.3f}".format, "q": "{:.3f}".format}))
print(f"\nRandom lists: {N}; Reactome terms selected per random list: mean {len(rnd)/N:.2f}; "
      f"lists with >=1 Reactome term: {rnd.list.nunique()/N:.1%}")
print("Random-term p-value distribution (quantiles 1%, 5%, 25%, 50%):", np.quantile(rnd.p, [.01, .05, .25, .5]).round(5))
print("Random terms by iteration:", rnd.iteration.value_counts().sort_index().head(8).to_dict())
