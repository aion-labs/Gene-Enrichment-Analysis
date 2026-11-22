"""(1) Parameter sensitivity (R5-M10, R1-1): one parameter changed at a time, for HIV and the random control.
The random lists are always generated with the changed parameters (per-run statistics; 2,000 lists, 200 per
reduced size). Usage: s10a_sensitivity.py <dataset> <config index>...   rows appended to s10a_parameter_sensitivity.tsv"""
import sys, os, json
import pandas as pd
from common import *
CONFIGS = [("default", {}), ("p < 0.005", dict(p_threshold=0.005)), ("p < 0.05", dict(p_threshold=0.05)),
           ("min overlap 2", dict(min_overlap=2)), ("min overlap 5", dict(min_overlap=5)),
           ("set size 15-500", dict(min_term_size=15, max_term_size=500)), ("set size 5-1000", dict(min_term_size=5, max_term_size=1000)),
           ("max 10 iterations", dict(max_iterations=10))]
OUT = "s10a_parameter_sensitivity.tsv"

def main(ds, idx):
    g = genes_of(ds); tag = ds.replace(" ", "_")
    for i in idx:
        name, par = CONFIGS[i]
        recs, table, null, terms, obs, attrs = stats(g, n_lists=2000, **par)
        kept = terms[terms.Kept == True]
        ks = sorted(f"{a}|{b}" for a, b in zip(kept.Library, kept.Term))
        cg = sorted(set(", ".join(kept.Genes).split(", "))) if len(kept) else []
        if name == "default":
            json.dump(dict(terms=ks, genes=cg), open(f"s10a_default_{tag}.json", "w"))
        d = json.load(open(f"s10a_default_{tag}.json"))
        J = lambda a, b: round(len(set(a) & set(b)) / max(1, len(set(a) | set(b))), 3)
        t = table.set_index("Metric")
        row = dict(dataset=ds, setting=name, selected=len(terms), kept=len(kept), chance_terms_per_list=attrs["expected_chance_terms"],
                   cluster_genes=int(obs["n_genes"]), cluster_terms=int(obs["n_terms"]),
                   p_cluster_size=float(t.loc["Cluster Size", "Empirical p"]), p_edges_per_gene=float(t.loc["Average Edges per Gene", "Empirical p"]),
                   jaccard_kept_terms_vs_default=J(ks, d["terms"]), jaccard_genes_vs_default=J(cg, d["genes"]))
        pd.DataFrame([row]).to_csv(OUT, sep="\t", index=False, mode="a", header=not os.path.exists(OUT))
        print(row, flush=True)

if __name__ == "__main__":
    main(sys.argv[1], [int(x) for x in sys.argv[2:]])
