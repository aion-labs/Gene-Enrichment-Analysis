"""Checks of the v5 LLM answers: invented genes, communities covered, structure, word counts,
significance statement, and the experimental-context sections (for reading)."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, gzip, os, re, csv
import pandas as pd
R = (str(_REPO) + "/")
LLM_DIR = os.environ.get("LLM_DIR", "llm")                       # answer folder (one per model)
SUFFIX = "" if LLM_DIR == "llm" else "_" + LLM_DIR.replace("llm_", "")
SYM = set()
with gzip.open(R + "data/recent_release/Homo_sapiens.gene_info.gz", "rt") as f:
    next(f)
    for l in f: SYM.add(l.split("\t")[2])
STOP = set("""GO ER RNA DNA ATP NOT SET CAT MET CAN AND FOR ARE HAS WAS ITS THE ALL ANY ONE TWO MAX MIN PC CD IL TNF MHC HLA
GPI ECM ROS RNS HIF NF ERK MAPK PI3K TOR AKT NFKB KEGG PID WP RNAI CRISPR SEQ FDR NA ID PH CO DE LOW IMPACT SHE BAD
CLOCK C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 H HIV HAI PBMC SIRNA MRNA COPII COPI TGN NPC ESCRT""".split())
RUNS = {"HIV": "HIV_dependency_factors.symbols_*", "HIV GO": "hiv_go_only/HIV*", "Random q<0.10": "random_control_q0.10/*",
        "Random q=1": "random_control_q1/*", "Hypoxia": "hypoxia-genes.symbols_*", "Vaccination": "vaccination_*"}
rows, ctx = [], []
for name, pat in RUNS.items():
    d = [x for x in glob.glob(pat) if os.path.isfile(x + "/network_llm.dot")][0]
    dot = open(f"{d}/network_llm.dot").read()
    genes = set(re.findall(r'^\s*"([^"]+)" \[type=gene', dot, re.M))
    comms = list(pd.read_csv(f"{d}/network_modules_communities.tsv", sep="\t").Community)
    termwords = set(w for t in re.findall(r'label="([^"]+)"', dot) for w in re.findall(r"[A-Z0-9]+", t))
    for f in sorted(glob.glob(f"{d}/{LLM_DIR}/*_r*.md")):
        t = open(f).read(); kind = "second" if "followup" in f else "first"
        toks = set(re.findall(r"(?<![A-Za-z0-9-])[A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]+)?(?![a-z])", t))
        ment = {x for x in toks if x in SYM and x not in STOP}
        out = sorted(x for x in ment - genes if x not in termwords)
        heads = re.findall(r"^#+\s*(.+)$", t, re.M)
        comm_cov = [c for c in comms if re.search(rf"\b{c}\b", t)]
        first_sent = re.split(r"(?<=[.!?])\s", re.sub(r"^#.*\n", "", t.strip(), count=1).strip())[0]
        rows.append(dict(run=name, turn=kind, rep=re.search(r"_r(\d)", f).group(1), words=len(t.split()),
                         genes_named=len(ment & genes), not_in_network=" ".join(out),
                         communities=f"{len(comm_cov)}/{len(comms)}", sections=" | ".join(h.strip("* ") for h in heads if not h.startswith("C"))[:160],
                         opening=first_sent[:160]))
        if kind == "first":
            m = re.search(r"#+\s*Estimated experimental context\s*\n(.*?)(?=\n#+\s|\Z)", t, re.S)
            ctx.append(dict(run=name, rep=rows[-1]["rep"], context=(m.group(1).strip()[:700] if m else "(missing)")))
o = pd.DataFrame(rows); o.to_csv(f"llm_v5_checks{SUFFIX}.tsv", sep="\t", index=False)
pd.DataFrame(ctx).to_csv(f"llm_v5_experimental_context{SUFFIX}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(o[["run","turn","rep","words","genes_named","not_in_network","communities"]].to_string(index=False))
