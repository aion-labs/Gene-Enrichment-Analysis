"""Compare the iGEA networks of the three HIV siRNA screens (Brass, König, Zhou): genes, kept terms, communities.
Community similarity = Jaccard index of their annotation footprints: the union of the full gene sets (background-restricted,
as analysed) of the terms in each community. This compares biology, not the screen hits themselves."""
import glob, os, sys, itertools
import numpy as np, pandas as pd
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "step10_extras"))
from common import LIBS, BACKGROUND, FastCollection
H = os.path.dirname(os.path.abspath(__file__))
RUNS = {"Brass": glob.glob(f"{H}/../step09_end_to_end/HIV_dependency_factors.symbols_2026*/")[0],
        "König": glob.glob(f"{H}/Konig_2008_*/")[0], "Zhou": glob.glob(f"{H}/Zhou_2008_*/")[0]}
cols = {l.name: FastCollection.from_library(l, BACKGROUND, 10, 600) for l in LIBS}
members = {}
for n, c in cols.items():
    for t, name in enumerate(c.term_names):
        members[(n, name)] = {c.genes[j] for j in c.members[t]}
def key(lib, term):   # combined tables show display names ("GOBP: X Y"); collection files use GOBP_X_Y
    raw = term.replace(": ", "_").replace(" ", "_")
    for k in [(lib, term), (lib, raw)]:
        if k in members: return k
    cand = [k for k in members if k[0] == lib and k[1].replace("_", " ").upper() == term.replace(":", "").replace("_", " ").upper().replace("  ", " ")]
    return cand[0] if cand else None
data = {}
for s, d in RUNS.items():
    ct = pd.read_csv(d + "combined_iterative_results.tsv", sep="\t")
    kept = ct[ct["In network"] == True]
    genes_in = set(open(glob.glob(d + "*statistical_report.txt")[0]).read() and [])  # placeholder
    mt = pd.read_csv(d + "network_modules_terms.tsv", sep="\t"); mc = pd.read_csv(d + "network_modules_communities.tsv", sep="\t")
    terms = {(r.Library, r.Term) for r in kept.itertuples()}
    tkeys = {t: key(*t) for t in terms}
    missing = [t for t, k in tkeys.items() if k is None]
    assert not missing, missing[:3]
    comm = {}
    for r in mt.itertuples():
        comm.setdefault(r.Community, set()).add(tkeys[(r.Library, r.Term)])
    data[s] = dict(terms=set(tkeys.values()), comm=comm, lead=dict(zip(mc.Community, mc["Lead term"])),
                   hubs=dict(zip(mc.Community, mc["Top hub genes (degree)"])), nterms=dict(zip(mc.Community, mc.Terms)))
# input genes (validated) from the runs' term tables are not the inputs; read the validated input lists
from common import genes_of, DATA
import common
inputs = {}
IN = os.path.join(H, "..", "inputs")
for s, f in [("Brass", "HIV_dependency_factors.symbols.txt"), ("König", "Konig_2008_HIV_dependency_factors.symbols.txt"), ("Zhou", "Zhou_2008_HIV_dependency_factors.symbols.txt")]:
    common.DATA[s] = f; inputs[s] = set(genes_of(s))
rows = []
for a, b in itertools.combinations(RUNS, 2):
    ga, gb = inputs[a], inputs[b]; ta, tb = data[a]["terms"], data[b]["terms"]
    rows.append(dict(pair=f"{a}–{b}", genes_a=len(ga), genes_b=len(gb), shared_genes=len(ga & gb), gene_jaccard=round(len(ga & gb) / len(ga | gb), 3),
                     terms_a=len(ta), terms_b=len(tb), shared_terms=len(ta & tb), term_jaccard=round(len(ta & tb) / len(ta | tb), 3)))
pairs = pd.DataFrame(rows); pairs.to_csv(f"{H}/s11b_pairwise_overlap.tsv", sep="\t", index=False); print(pairs.to_string(index=False))
allg = inputs["Brass"] & inputs["König"] & inputs["Zhou"]; allt = data["Brass"]["terms"] & data["König"]["terms"] & data["Zhou"]["terms"]
print("genes in all three:", len(allg), sorted(allg)); print("terms in all three:", len(allt), sorted(t[1] for t in allt))
# community footprints and cross-screen similarity
fp = {(s, c): set().union(*[members[t] for t in ts]) for s in data for c, ts in data[s]["comm"].items()}
cm = []
for (s1, c1), (s2, c2) in itertools.combinations(fp, 2):
    if s1 == s2: continue
    A, B = fp[(s1, c1)], fp[(s2, c2)]
    sh = data[s1]["comm"][c1] & data[s2]["comm"][c2]
    cm.append(dict(screen_1=s1, community_1=c1, screen_2=s2, community_2=c2, footprint_jaccard=round(len(A & B) / len(A | B), 3),
                   shared_terms=len(sh), shared_term_names="; ".join(sorted(t[1] for t in sh))[:200]))
cm = pd.DataFrame(cm).sort_values("footprint_jaccard", ascending=False); cm.to_csv(f"{H}/s11b_community_similarity.tsv", sep="\t", index=False)
summ = []
for s in data:
    for c in sorted(data[s]["comm"], key=lambda x: int(x[1:])):
        m = cm[((cm.screen_1 == s) & (cm.community_1 == c)) | ((cm.screen_2 == s) & (cm.community_2 == c))]
        best = {}
        for o in data:
            if o == s: continue
            mo = m[(m.screen_1 == o) | (m.screen_2 == o)].iloc[0]
            oc = mo.community_1 if mo.screen_1 == o else mo.community_2
            best[o] = (oc, mo.footprint_jaccard)
        summ.append(dict(screen=s, community=c, terms=data[s]["nterms"][c], lead_term=data[s]["lead"][c], hubs=data[s]["hubs"][c],
                         **{f"best_{o}": f"{v[0]} ({v[1]})" for o, v in best.items()}))
pd.DataFrame(summ).to_csv(f"{H}/s11b_community_best_matches.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 55)
print(pd.DataFrame(summ).drop(columns=["hubs"]).to_string(index=False))
print(cm.head(25).drop(columns=["shared_term_names"]).to_string(index=False))
