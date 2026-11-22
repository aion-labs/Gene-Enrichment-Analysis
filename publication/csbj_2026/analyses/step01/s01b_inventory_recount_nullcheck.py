"""Step 1b — collection inventory (10 collections), 5/4/2 recount, stored-null consistency check.
Run from a working dir with data/recent_release: python s01b_... <repo> <out_dir>"""
import sys, os, json, logging
from pathlib import Path
import numpy as np, pandas as pd
REPO = Path(sys.argv[1]).resolve(); OUT = Path(sys.argv[2]).resolve()
sys.path.insert(0, str(REPO / "code")); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet
from gene_set_library import GeneSetLibrary
from network_connectivity_benchmark import NetworkConnectivityAnalyzer
LIBS = [("H: Hallmark Gene Sets","h.all","H: Hallmark Gene Sets","MSigDB H"),
 ("C2: Reactome Pathways","c2.cp.reactome","C2: CP: Reactome Pathways","MSigDB C2:CP:REACTOME"),
 ("C2: KEGG Legacy","c2.cp.kegg_legacy","C2: CP: KEGG Legacy","MSigDB C2:CP:KEGG_LEGACY"),
 ("C2: KEGG MEDICUS","c2.cp.kegg_medicus","C2: CP: KEGG MEDICUS","MSigDB C2:CP:KEGG_MEDICUS"),
 ("C2: BioCarta","c2.cp.biocarta","C2: CP: BioCarta","MSigDB C2:CP:BIOCARTA"),
 ("C2: Pathway Interaction Database","c2.cp.pid","C2: CP: Pathway Interaction Database","MSigDB C2:CP:PID"),
 ("C2: WikiPathways","c2.cp.wikipathways","C2: CP: WikiPathways","MSigDB C2:CP:WIKIPATHWAYS"),
 ("C5: Gene Ontology: Biological Process","c5.go.bp","C5: GO: Biological Process","MSigDB C5:GO:BP"),
 ("C5: Gene Ontology: Cellular Component","c5.go.cc","C5: GO: Cellular Component","MSigDB C5:GO:CC"),
 ("C5: Gene Ontology: Molecular Function","c5.go.mf","C5: GO: Molecular Function","MSigDB C5:GO:MF")]
bg = BackgroundGeneSet(str(REPO / "data/backgrounds/all_genes.txt"))
hiv = {l.strip() for l in open(REPO / "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()}
ora = pd.read_csv(OUT / "s01_ora_all_tested.tsv", sep="\t"); ig = pd.read_csv(OUT / "s01_igea_terms.tsv", sep="\t")
inv, rec = [], []
for disp, f, pname, src in LIBS:
    lib = GeneSetLibrary(str(REPO / f"data/libraries/{f}.v2025.1.Hs.symbols.gmt"), name=disp)
    kept = [t for t in lib.library if 10 <= t["size"] <= 600]
    kgenes = set().union(*[set(t["genes"]) for t in kept])
    M = len(bg.genes & kgenes); n_in = len(hiv & kgenes)
    maxposs = sum(len(hiv & set(t["genes"])) >= 3 for t in kept)
    o = ora[ora.Library == disp]; osig = o[(o.p < 0.01) & (o.overlap >= 3)]; ofdr = set(osig[osig.fdr < 0.05].Term)
    iset = set(ig[ig.Library == disp].Term.str.replace(": ", "_", n=1).str.replace(" ", "_"))
    # iGEA Term is formatted ("GOBP: X Y"); map back to raw names via ORA table
    raw = {t.replace("_", " ").replace(" ", "_") for t in o.Term}
    fmt2raw = {}
    for t in o.Term:
        s = t.replace("_", " "); w = s.split(); fmt2raw[f"{w[0]}: {' '.join(w[1:])}" if len(w) > 1 else s] = t
    iset = {fmt2raw.get(t, t) for t in ig[ig.Library == disp].Term}
    if iset == ofdr: cat = "identical"
    elif ofdr <= iset: cat = "iGEA adds terms (all FDR terms captured)"
    else: cat = "limited overlap"
    inv.append({"Collection": disp, "Source (release)": f"{src} (v2025.1.Hs)", "GMT file": f"{f}.v2025.1.Hs.symbols.gmt",
        "Gene sets (total)": lib.num_terms, "Gene sets 10-600": len(kept), "Distinct genes (all sets)": lib.size,
        "Distinct genes (sets 10-600)": len(kgenes), "Collection background M (all_genes ∩ collection)": M,
        "HIV genes in collection (n)": n_in, "Max possible terms (sets 10-600, overlap>=3)": maxposs,
        "ORA p<0.01 (overlap>=3)": len(osig), "ORA FDR<0.05": len(ofdr), "iGEA terms": len(iset),
        "In ORA comparison": "yes", "In iGEA network": "yes", "In random-list null": "yes"})
    rec.append({"Collection": disp, "ORA FDR<0.05": len(ofdr), "iGEA": len(iset), "shared": len(iset & ofdr),
                "FDR terms missed by iGEA": len(ofdr - iset), "iGEA-only": len(iset - ofdr), "category": cat})
inv = pd.DataFrame(inv); rec = pd.DataFrame(rec)
inv.to_csv(OUT / "s01_collection_inventory.csv", index=False); rec.to_csv(OUT / "s01_recount_ora_vs_igea.csv", index=False)
print(inv.drop(columns=["GMT file","Source (release)"]).to_string()); print(rec.to_string()); print(rec.category.value_counts())
# ---- stored null: do any runs contain iterations AFTER a p>=0.01 iteration was dropped (gaps)?
df = pd.read_csv(REPO / "permutations/merged_permutation_results.tsv", sep="\t", usecols=["Iteration","Library","gene_list_size","filename","iteration p-value"])
g = df.groupby(["gene_list_size","filename","Library"]).Iteration.agg(["max","count"])
g["gap"] = g["max"] > g["count"]
gap_by_size = g.groupby(level=0).gap.mean().round(4)
runs_with_gap = g.groupby(level=[0,1]).gap.any().groupby(level=0).mean().round(4)
res = {"rows": len(df), "max_p_stored": float(df["iteration p-value"].max()), "max_iteration": int(df.Iteration.max()),
       "libraries": sorted(df.Library.unique()), "frac_library_runs_with_gap_by_size": gap_by_size.to_dict(),
       "frac_random_lists_with_any_gap_by_size": runs_with_gap.to_dict()}
# impact at size 300 (10 network collections): null with vs without post-gap iterations
d = df[(df.gene_list_size == 300) & df.Library.isin([l[2] for l in LIBS])].copy()
full = pd.read_csv(REPO / "permutations/merged_permutation_results.tsv", sep="\t")
full = full[(full.gene_list_size == 300) & full.Library.isin([l[2] for l in LIBS])]
def prefix(x):
    x = x.sort_values("Iteration"); ok = (x.Iteration.values == np.arange(1, len(x) + 1)); 
    k = len(x) if ok.all() else int(np.argmin(ok)); return x.iloc[:k]
trunc = full.groupby(["filename","Library"], group_keys=False).apply(prefix)
def stats(t):
    out = []
    for fn, x in t.groupby("filename"):
        a = NetworkConnectivityAnalyzer()
        a.add_igea_results([{"Term": r.Term, "Iteration": r.Iteration, "Library": r.Library,
            "Genes removed for next iteration": [s.strip() for s in str(r["Genes removed for next iteration"]).split(",") if s.strip()]} for _, r in x.iterrows()])
        a.compute_metrics(include_library_diversity=True); c = a.get_clusters()
        out.append({"filename": fn, "size": c[0]["size"] if c else 0, "genes": c[0]["n_genes"] if c else 0, "terms": c[0]["n_terms"] if c else 0, "edges": c[0]["n_edges"] if c else 0})
    o = pd.DataFrame(out).set_index("filename").reindex(sorted(full.filename.unique())).fillna(0)
    return o
a, b = stats(full), stats(trunc)
a.to_csv(OUT / "s01_null300_as_stored.csv"); b.to_csv(OUT / "s01_null300_truncated.csv")
res["size300_null_as_used"] = a.agg(["mean","std"]).round(2).to_dict(); res["size300_null_stopping_rule_respected"] = b.agg(["mean","std"]).round(2).to_dict()
json.dump(res, open(OUT / "s01_null_check.json", "w"), indent=2, default=str); print(json.dumps(res, indent=2, default=str))
