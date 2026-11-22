"""Do LLM responses group terms the way Leiden communities do?
For each response, terms cited on the same line (bullet, table row or paragraph) are 'co-mentioned'.
Compare co-mention rates for term pairs within the same Leiden community vs different communities."""
import re, glob, os, itertools, pandas as pd
L = "../llm_package/responses"
JOB = {"hiv_go_igea": "hiv_go", "hiv_all_igea": "hiv", "hypoxia_all_igea": "hypoxia", "vaccination_all_igea": "vaccination"}
norm = lambda s: re.sub(r"[^A-Z0-9]", "", s.upper())
rows = []
for job, net in JOB.items():
    t = pd.read_csv(f"{net}_terms.tsv", sep="\t")
    comm = dict(zip(t.Term.map(norm), t.Community))
    keys = sorted(comm, key=len, reverse=True)
    for f in sorted(glob.glob(f"{L}/{job}__r*.md")):
        lines = open(f).read().split("\n")
        co = set(); seen = set()
        for ln in lines:
            s = norm(ln); found = []
            for k in keys:
                if k in s: found.append(k); s = s.replace(k, "|")
            seen.update(found)
            co.update(frozenset(p) for p in itertools.combinations(sorted(set(found)), 2))
        w = [frozenset(p) for p in itertools.combinations(sorted(seen), 2) if comm[p[0]] == comm[p[1]]]
        b = [frozenset(p) for p in itertools.combinations(sorted(seen), 2) if comm[p[0]] != comm[p[1]]]
        rw = sum(p in co for p in w) / max(1, len(w)); rb = sum(p in co for p in b) / max(1, len(b))
        rows.append(dict(response=os.path.basename(f)[:-3], terms_cited=len(seen), of_terms=len(comm),
                         within_pairs=len(w), within_comention=round(rw, 3), between_pairs=len(b),
                         between_comention=round(rb, 3), ratio=round(rw / rb, 1) if rb else None))
o = pd.DataFrame(rows); o.to_csv("llm_vs_communities.tsv", sep="\t", index=False); print(o.to_string(index=False))
