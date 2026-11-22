"""Hub genes vs annotation: gene degree in the network vs number of gene sets (10 collections,
size 10-600) containing the gene. iGEA (q<0.10) vs ORA (BH FDR<0.05), HIV."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, re, collections, pandas as pd
from scipy.stats import spearmanr
R = (str(_REPO) + "/")
src = open("../step05_two_layer/run_two_layer.py").read()
P = re.findall(r'\("([^"]+)","(c\d[^"]+|h\.[^"]+)"\)', src)
sys.path.insert(0, R + "code"); import logging; logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet
_here = os.getcwd(); os.chdir(os.path.expanduser("~/s01wd"))  # needs data/recent_release/Homo_sapiens.gene_info (validated background)
bg = set(BackgroundGeneSet(R + "data/backgrounds/all_genes.txt").genes); os.chdir(_here)
print("background", len(bg))
ann = collections.Counter()
for name, f in P:
    for l in open(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt"):
        c = l.rstrip("\n").split("\t")
        if 10 <= len(c) - 2 <= 600:
            for g in set(c[2:]) & bg: ann[g] += 1
print(len(P), "collections")
def deg(rows):
    d = collections.Counter()
    for gl in rows:
        for g in gl.split(","): d[g.strip()] += 1
    return d
import glob
igea = pd.read_csv(glob.glob("../step09_end_to_end/HIV_dependency_factors.symbols_*/HIV_dependency_factors.symbols_term_qvalues.tsv")[0], sep="\t"); igea = igea[igea.Kept]
ora = pd.read_csv("../step04_term_fdr/hiv_final_ora_vs_igea_terms.tsv", sep="\t"); ora = ora[ora.Source == "ORA"]
out = []
for name, d in (("iGEA q<0.10", deg(igea.Genes)), ("ORA FDR<0.05", deg(ora["Gene list"]))):
    g = sorted(d); x = [d[k] for k in g]; y = [ann[k] for k in g]
    rho, p = spearmanr(x, y)
    top = sorted(g, key=lambda k: (-d[k], k))[:10]
    out.append(dict(network=name, terms=len(igea) if name.startswith("i") else len(ora), genes=len(g),
                    max_degree=max(x), median_degree=pd.Series(x).median(), spearman_rho=round(rho, 3), p=f"{p:.1e}",
                    top10_hubs=", ".join(f"{k} ({d[k]}; {ann[k]} sets)" for k in top)))
    pd.DataFrame({"Gene": g, "Degree": x, "Annotation sets": y}).to_csv(f"hiv_degree_vs_annotation_{name.split()[0]}.tsv", sep="\t", index=False)
o = pd.DataFrame(out); o.to_csv("hiv_hub_annotation_bias.tsv", sep="\t", index=False)
print(o.drop(columns="top10_hubs").to_string(index=False))
for r in out: print(r["network"], ":", r["top10_hubs"])
