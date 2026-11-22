"""Graphviz DOT files for the regenerated HIV figures (paper colouring), plus PNG/SVG previews.
Inputs: q-filtered iGEA terms from the two-layer statistics (q < 0.10, 10,000 random lists):
  full network  -> analyses/step09_end_to_end/HIV_*/..._term_qvalues.tsv (CLI, 10 collections)
  GO network    -> hiv_go_terms.tsv   (random lists analysed with the 3 GO collections)
  GO BP         -> hiv_gobp_terms.tsv (random lists analysed with GO BP only)
ORA (Figure 1A): GO BP sets 10-600, hypergeometric test (upper tail), overlap >= 3, p < 0.01, BH FDR over all tested sets.
Figure 1B colours: iGEA terms shared with ORA at FDR < 0.05 vs iGEA-only (corrected 21 Nov; the first version compared with FDR < 0.10).
Usage: python make_dots.py [<output folder>]   (default: this folder)"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import os, sys, re, textwrap, itertools, logging
import numpy as np, pandas as pd
from statsmodels.stats.multitest import multipletests
R = (str(_REPO) + "/"); sys.path.insert(0, R + "code"); logging.disable(logging.WARNING)
from background_gene_set import BackgroundGeneSet; from gene_set import GeneSet; from gene_set_library import GeneSetLibrary
from gene_converter import GeneConverter; from fast_igea import FastCollection, exact_pvalue
from enrichment import format_term_name
F = os.path.dirname(os.path.abspath(__file__)); V7 = os.path.dirname(os.path.dirname(F))
import glob
FULL = glob.glob(str(_PUB / "analyses" / "step09_end_to_end" / "HIV_dependency_factors.symbols_2*" / "HIV_dependency_factors.symbols_term_qvalues.tsv"))[0]  # CLI end-to-end run
OUT = sys.argv[1] if len(sys.argv) > 1 else F   # output folder for DOT files and Fig1_counts.tsv (inputs are read from this folder)
# paper colours: saturated = alias.json (Supp Fig 2); pastel = Figure 2 / Supp Fig 3
SAT = {"H: Hallmark Gene Sets": "#FE456B", "C2: BioCarta": "#A35E9A", "C2: KEGG MEDICUS": "#29996D", "C2: Pathway Interaction Database": "#116B1A",
       "C2: Reactome Pathways": "#5A3F91", "C2: WikiPathways": "#0005FD", "C2: KEGG Legacy": "#B17500",
       "C5: Gene Ontology: Biological Process": "#C50A00", "C5: Gene Ontology: Cellular Component": "#009AD8", "C5: Gene Ontology: Molecular Function": "#9A00C3"}
PASTEL = {"C5: Gene Ontology: Biological Process": ("#F1948A", "#EC7063"), "C5: Gene Ontology: Cellular Component": ("#85C1E9", "#5DADE2"),
          "C5: Gene Ontology: Molecular Function": ("#D7BDE2", "#BB8FCE")}            # (fill, edge)
FIG1 = {"shared": ("#F5C518", "#1F1F1F"), "ora": ("#3E78B2", "white"), "igea": ("#EF3E4A", "white")}

def q(s): return '"' + s.replace('"', "'") + '"'
def label(raw, width=24): return "\\n".join(textwrap.wrap(format_term_name(raw), width))
def header(name, layout="neato", extra=""):
    return [f"graph {name} {{", f'  graph [layout={layout}, overlap=false, splines=false, outputorder=edgesfirst, bgcolor="white", pad=0.3, sep="+6"{extra}];',
            '  node [fontname="Helvetica-Bold"];', '  edge [penwidth=0.8];']
def term_node(tid, raw, fill, font, border="#1F1F1F", width=24, fs=14):
    return f'  {q(tid)} [label="{label(raw, width)}", shape=box, style="rounded,filled", fillcolor="{fill}", fontcolor="{font}", color="{border}", penwidth=1.2, fontsize={fs}, margin="0.12,0.06"];'
def gene_node(g, mode):
    if mode == "point": return f'  {q("g_" + g)} [label="", shape=point, width=0.11, color="black", fillcolor="black", tooltip="{g}"];'
    return f'  {q("g_" + g)} [label="{g}", shape=plaintext, fontname="Helvetica", fontsize=10, fontcolor="black", margin=0];'

def network_dot(terms, name, palette, genes="label", legend=False):
    t = terms[terms.Kept].copy()
    L = header(name); allg = set()
    for _, r in t.iterrows():
        tid = "t_" + r.Term
        if palette == "pastel":
            fill, _ = PASTEL[r.Library]; L.append(term_node(tid, r.Term, fill, "#1F1F1F"))
        else:
            L.append(term_node(tid, r.Term, SAT[r.Library], "white", border=SAT[r.Library], width=34, fs=12))
        allg |= set(str(r.Genes).split(", "))
    L += [gene_node(g, genes) for g in sorted(allg)]
    for _, r in t.iterrows():
        col = PASTEL[r.Library][1] if palette == "pastel" else SAT[r.Library]
        L += [f'  {q("g_" + g)} -- {q("t_" + r.Term)} [color="{col}"];' for g in str(r.Genes).split(", ")]
    if legend:
        rows = "".join(f'<TR><TD BGCOLOR="{SAT[l]}"><FONT COLOR="white">{l}</FONT></TD></TR>' for l in SAT if l in set(t.Library))
        L.append(f'  legend [shape=plaintext, fontname="Helvetica", fontsize=11, label=<<TABLE BORDER="1" CELLBORDER="0" CELLSPACING="3" CELLPADDING="3">{rows}</TABLE>>];')
    return "\n".join(L + ["}"]) + "\n", len(t), len(allg)

def fig1_dot(sets, colours, name, clusters):
    """sets: {raw term: genes}; colours: {raw term: key in FIG1}; clusters: draw connected components in grey ellipses (fdp)."""
    L = header(name, layout="fdp" if clusters else "neato")
    comps = []
    if clusters:   # connected components of the term-gene graph
        terms = list(sets); parent = {t: t for t in terms}
        def find(x):
            while parent[x] != x: x = parent[x]
            return x
        for a, b in itertools.combinations(terms, 2):
            if sets[a] & sets[b]: parent[find(a)] = find(b)
        groups = {}
        for t in terms: groups.setdefault(find(t), []).append(t)
        comps = list(groups.values())
    else:
        comps = [[t] for t in sets]
    placed = set()
    for k, comp in enumerate(comps):
        genes = sorted(set().union(*[sets[t] for t in comp]) - placed); placed |= set(genes)
        if clusters:
            L.append(f'  subgraph cluster_{k} {{ style="filled,rounded"; fillcolor="#F2F2F2"; color="#F2F2F2"; label="";')
        for t in comp:
            fill, font = FIG1[colours[t]]; L.append(term_node("t_" + t, t, fill, font, width=22))
        L += [gene_node(g, "label") for g in genes]
        if clusters: L.append("  }")
    for t, gs in sets.items():
        L += [f'  {q("g_" + g)} -- {q("t_" + t)} [color="#BDBDBD"];' for g in sorted(gs)]
    return "\n".join(L + ["}"]) + "\n"

if __name__ == "__main__":
    out = []
    full = pd.read_csv(FULL, sep="\t"); go = pd.read_csv(f"{F}/hiv_go_terms.tsv", sep="\t"); gobp = pd.read_csv(f"{F}/hiv_gobp_terms.tsv", sep="\t")
    for df in (full, go, gobp): df["Kept"] = df["Kept"].astype(bool)
    d, nt, ng = network_dot(full, "SuppFig2_iGEA_HIV_full_network", "sat", "label", legend=True)
    open(f"{OUT}/SuppFig2_iGEA_HIV_full_network_q0.10.dot", "w").write(d); out.append(("Supp Fig 2 (10 collections)", nt, ng))
    d, nt, ng = network_dot(go, "Fig2_iGEA_GO_network", "pastel", "point")
    open(f"{OUT}/Fig2_iGEA_GO_network_q0.10.dot", "w").write(d); out.append(("Figure 2 (GO, dots)", nt, ng))
    d, nt, ng = network_dot(go, "SuppFig3_iGEA_GO_network_genes", "pastel", "label")
    open(f"{OUT}/SuppFig3_iGEA_GO_network_genes_q0.10.dot", "w").write(d); out.append(("Supp Fig 3 (GO, gene names)", nt, ng))
    # Figure 1: ORA (GO BP, BH FDR) vs iGEA (GO BP, q < 0.10, GO BP-only random lists)
    bg = BackgroundGeneSet(R + "data/backgrounds/all_genes.txt"); conv = GeneConverter()
    hiv = [m for m in (conv.validate_and_map_symbol(l.strip()) for l in open(R + "data/gene_lists/HIV_dependency_factors.symbols.txt") if l.strip()) if m]
    H = set(GeneSet(hiv, bg.genes, "HIV", hgcn=True, format=False).genes)  # genes outside the background removed (one universe)
    col = FastCollection.from_library(GeneSetLibrary(R + "data/libraries/c5.go.bp.v2025.1.Hs.symbols.gmt", name="GO BP"), sorted(bg.genes))
    x = np.zeros(len(col.genes)); [x.__setitem__(col.index[g], 1) for g in H if g in col.index]; n = x.sum(); k = col.A @ x
    p = np.array([exact_pvalue("Fisher's Exact Test", int(k[i]), int(col.K[i]), int(n), col.M) for i in range(len(col.term_names))])
    fdr = multipletests(p, method="fdr_bh")[1]
    ig = {r.Term: set(str(r.Genes).split(", ")) for _, r in gobp[gobp.Kept].iterrows()}
    rows = []
    for thr in (0.05, 0.10):
        sel = [i for i in range(len(p)) if p[i] < 0.01 and k[i] >= 3 and fdr[i] < thr]
        ora = {col.term_names[i]: {col.genes[j] for j in col.members[i] if col.genes[j] in H} for i in sel}
        shared = set(ora) & set(ig)
        tag = f"FDR{thr:.2f}"
        open(f"{OUT}/Fig1A_ORA_GOBP_{tag}.dot", "w").write(fig1_dot(ora, {t: "shared" if t in shared else "ora" for t in ora}, f"Fig1A_ORA_GOBP_{tag.replace('.', '')}", False))
        rows.append({"ORA threshold": f"BH FDR < {thr}", "ORA terms": len(ora), "ORA genes": len(set().union(*ora.values())),
                     "iGEA terms (q < 0.10)": len(ig), "iGEA genes": len(set().union(*ig.values())), "shared terms": len(shared),
                     "shared": "; ".join(sorted(format_term_name(t) for t in shared)), "ORA only": "; ".join(sorted(format_term_name(t) for t in set(ora) - shared)),
                     "iGEA only": "; ".join(sorted(format_term_name(t) for t in set(ig) - shared))})
        if thr == 0.05: shared05 = shared
    open(f"{OUT}/Fig1B_iGEA_GOBP_q0.10.dot", "w").write(fig1_dot(ig, {t: "shared" if t in shared05 else "igea" for t in ig}, "Fig1B_iGEA_GOBP", False))
    pd.DataFrame(rows).to_csv(f"{OUT}/Fig1_counts.tsv", sep="\t", index=False)
    for r in out: print("%-30s terms %3d  genes %3d" % r)
    print(pd.DataFrame(rows).drop(columns=["shared", "ORA only", "iGEA only"]).to_string(index=False))
    for r in rows: print(f"\n{r['ORA threshold']}\n  shared: {r['shared']}\n  ORA only: {r['ORA only']}\n  iGEA only: {r['iGEA only']}")
