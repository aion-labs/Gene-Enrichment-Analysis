"""HIV, 10 collections: ORA (p < 0.01, overlap >= 3, BH FDR < 0.05 per collection) vs the final iGEA network
(signal-adjusted q < 0.10, CLI end-to-end run in step09). One universe: input and gene sets restricted to the background.
Usage: python s04_final_ora_vs_igea.py <out_dir>"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import sys, os, glob, itertools, logging
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
HERE = os.path.dirname(os.path.abspath(__file__))
P = [("H: Hallmark Gene Sets","h.all"),("C2: Reactome Pathways","c2.cp.reactome"),("C2: KEGG Legacy","c2.cp.kegg_legacy"),("C2: KEGG MEDICUS","c2.cp.kegg_medicus"),
     ("C2: BioCarta","c2.cp.biocarta"),("C2: Pathway Interaction Database","c2.cp.pid"),("C2: WikiPathways","c2.cp.wikipathways"),
     ("C5: Gene Ontology: Biological Process","c5.go.bp"),("C5: Gene Ontology: Cellular Component","c5.go.cc"),("C5: Gene Ontology: Molecular Function","c5.go.mf")]
def mj(sets):
    return round(float(np.mean([len(a & b) / len(a | b) for a, b in itertools.combinations(sets, 2)])), 3) if len(sets) > 1 else np.nan
if __name__ == "__main__":
    OUT = sys.argv[1]
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    H = set(GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False).genes); B = sorted(bg.genes)
    ig = pd.read_csv(glob.glob(f"{HERE}/../step09_end_to_end/HIV_dependency_factors.symbols_*/HIV_dependency_factors.symbols_term_qvalues.tsv")[0], sep="\t")
    ig["Kept"] = ig["Kept"].astype(bool)
    rows, terms, G = [], [], {"ora": set(), "ig": set()}
    for name, f in P:
        col = FastCollection.from_library(GeneSetLibrary(R + f"data/libraries/{f}.v2025.1.Hs.symbols.gmt", name=name), B)
        x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]; n = x.sum(); k = col.A @ x
        p = np.array([exact_pvalue("Fisher's Exact Test", int(k[t]), int(col.K[t]), int(n), col.M) for t in range(len(col.term_names))])
        fdr = multipletests(p, method="fdr_bh")[1]
        sel = [t for t in range(len(p)) if p[t] < 0.01 and k[t] >= 3 and fdr[t] < 0.05]
        hs = {col.term_names[t]: {col.genes[j] for j in col.members[t] if col.genes[j] in H} for t in sel}
        igc = ig[(ig.Library == name) & ig.Kept]; igs = {r.Term: set(str(r.Genes).split(", ")) for r in igc.itertuples()}
        g_ora = set().union(*hs.values()) if hs else set(); g_ig = set().union(*igs.values()) if igs else set()
        lost = [t for t, g in hs.items() if len(g & g_ig) < 0.5 * len(g)]
        rows.append({"Collection": name, "ORA FDR<0.05": len(hs), "iGEA q<0.10": len(igs), "Both": len(set(hs) & set(igs)),
                     "Genes ORA": len(g_ora), "Genes iGEA": len(g_ig), "Genes only ORA": len(g_ora - g_ig), "Genes only iGEA": len(g_ig - g_ora),
                     "ORA terms mostly uncovered by iGEA (<50% genes)": len(lost), "Mean Jaccard ORA": mj(list(hs.values())), "Mean Jaccard iGEA": mj(list(igs.values()))})
        terms += [{"Collection": name, "Source": "ORA", "Term": t, "p": p[col.term_names.index(t)], "FDR": fdr[col.term_names.index(t)],
                   "In iGEA network": t in igs, "Mostly uncovered by iGEA": t in lost, "Genes": len(g), "Gene list": ", ".join(sorted(g))} for t, g in hs.items()]
        terms += [{"Collection": name, "Source": "iGEA", "Term": t, "In ORA FDR<0.05": t in hs, "Genes": len(g), "Gene list": ", ".join(sorted(g))} for t, g in igs.items()]
        G["ora"] |= g_ora; G["ig"] |= g_ig
    s = pd.DataFrame(rows); s.to_csv(f"{OUT}/hiv_final_ora_vs_igea_summary.tsv", sep="\t", index=False)
    pd.DataFrame(terms).to_csv(f"{OUT}/hiv_final_ora_vs_igea_terms.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250); print(s.to_string(index=False))
    print("sums:", s.drop(columns=["Collection", "Mean Jaccard ORA", "Mean Jaccard iGEA"]).sum().to_dict())
    print(f"distinct genes: ORA {len(G['ora'])}, iGEA {len(G['ig'])}, both {len(G['ora'] & G['ig'])}, only ORA {len(G['ora'] - G['ig'])}, only iGEA {len(G['ig'] - G['ora'])}")
    lost = pd.DataFrame(terms); lost = lost[(lost.Source == "ORA") & (lost["Mostly uncovered by iGEA"] == True)]
    print("ORA terms mostly uncovered:", "; ".join(lost.Term))
    print("ORA-only genes:", ", ".join(sorted(G["ora"] - G["ig"])))
