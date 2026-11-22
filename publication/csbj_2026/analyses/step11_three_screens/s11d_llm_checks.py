"""Checks and section extracts of the two-turn LLM answers for the three screens (Brass, König, Zhou)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, gzip, os, re, json
import pandas as pd
H = os.path.dirname(os.path.abspath(__file__)); R = (str(_REPO) + "/")
RUNS = {"Brass": glob.glob(f"{H}/../step09_end_to_end/HIV_dependency_factors.symbols_2026*/")[0],
        "König": glob.glob(f"{H}/Konig_2008_*/")[0], "Zhou": glob.glob(f"{H}/Zhou_2008_*/")[0]}
SYM = set()
with gzip.open(R + "data/recent_release/Homo_sapiens.gene_info.gz", "rt") as f:
    next(f)
    for l in f: SYM.add(l.split("\t")[2])
STOP = set("GO ER RNA DNA ATP NOT SET CAT MET CAN AND FOR ARE HAS WAS ITS THE ALL ANY ONE TWO PC CD IL TNF MHC HLA GPI ECM HIF NF ERK MAPK PI3K TOR AKT NFKB KEGG PID WP RNAI CRISPR SEQ FDR NA ID PH CO DE HIV SIRNA MRNA COPII COPI TGN NPC ESCRT LTR TAT REV NEF VPR GAG ENV".split())
def section(t, name):
    m = re.search(rf"#+\s*\d*\.?\s*\**{name}[^\n]*\n(.*?)(?=\n#+\s|\Z)", t, re.S | re.I)
    return m.group(1).strip() if m else ""
rows, ext = [], {}
for s, d in RUNS.items():
    dot = open(d + "network_llm.dot").read()
    genes = set(re.findall(r'^\s*"([^"]+)" \[type=gene', dot, re.M))
    comms = list(pd.read_csv(d + "network_modules_communities.tsv", sep="\t").Community)
    words = set(w for t in re.findall(r'label="([^"]+)"', dot) for w in re.findall(r"[A-Z0-9]+", t))
    for f in sorted(glob.glob(d + "llm_gpt-6.2-sol/*_r*.md")):
        t = open(f).read(); turn = "2" if "followup" in f else "1"; rep = re.search(r"_r(\d)", f).group(1)
        toks = set(re.findall(r"(?<![A-Za-z0-9-])[A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]+)?(?![a-z])", t))
        ment = {x for x in toks if x in SYM and x not in STOP}
        rows.append(dict(screen=s, turn=turn, rep=rep, words=len(t.split()), genes_named=len(ment & genes),
                         outside=" ".join(sorted(ment - genes - words)), communities=str(sum(1 for c in comms if re.search(r"\b" + c + r"\b", t))) + "/" + str(len(comms)),
                         names_HIV_or_virus=bool(re.search(r"\bHIV\b|virus|viral|infect", t, re.I))))
        key = (s, turn, rep)
        if turn == "1":
            ext[key] = dict(summary=section(t, "Executive summary"), context=section(t, "Estimated experimental context"),
                            sections=re.findall(r"^#+\s*(C\d+[^\n]*)", t, re.M))
        else:
            ext[key] = dict(summary=section(t, "Revised summary"), new=section(t, "What may be new"), priorities=section(t, "Priorities"),
                            fits=section(t, "What fits"), notfit=section(t, "What does not fit"))
pd.DataFrame(rows).to_csv(f"{H}/s11d_llm_checks.tsv", sep="\t", index=False)
json.dump({"|".join(k): v for k, v in ext.items()}, open(f"{H}/s11d_llm_sections.json", "w"), indent=1, ensure_ascii=False)
pd.set_option("display.width", 200); print(pd.DataFrame(rows).to_string(index=False))
