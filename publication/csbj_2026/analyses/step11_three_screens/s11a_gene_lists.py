"""Gene lists of the three genome-wide siRNA screens for HIV dependency factors (as published), mapped to current
NCBI symbols through Entrez Gene IDs (gene_info; discontinued IDs through gene_history).
Brass 2008: the list used throughout the paper (inputs/HIV_dependency_factors.symbols.txt).
König 2008 (Cell 135:49): Table S3 (NIHMS73737-supplement-04.xls), 295 genes (unique Gene_IDs).
Zhou 2008 (Cell Host Microbe 4:495): Table S1 (mmc2.xls), the 390 primary hits minus those not expressed in activated
T cells or macrophages ("Not expressed", "No expression data", "Expression unknown"), i.e. the published 311 host factors."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import gzip, os, re
import pandas as pd
R = (str(_REPO) + "/data/recent_release/")
H = os.path.dirname(os.path.abspath(__file__)); KZ = os.path.join(H, "..", "..", "Konig and Zhou"); OUT = os.path.join(H, "..", "inputs")
sym = {}
with gzip.open(R + "Homo_sapiens.gene_info.gz", "rt") as f:
    next(f)
    for l in f:
        p = l.split("\t"); sym[p[1]] = p[2]
hist = {}
with gzip.open(R + "gene_history_9606.gz", "rt") as f:
    next(f)
    for l in f:
        p = l.split("\t")
        if p[0] == "9606": hist[p[2]] = p[1]
def to_symbol(gid):
    g = str(gid).strip().split(".")[0]
    seen = set()
    while g not in sym and g in hist and g not in seen and hist[g] != "-":
        seen.add(g); g = hist[g]
    return sym.get(g)
rep = []
k = pd.read_excel(os.path.join(KZ, "Konig", "NIHMS73737-supplement-04.xls"), header=0)
k = k[pd.to_numeric(k.Gene_ID, errors="coerce").notna()]
kid = sorted(set(k.Gene_ID.astype(int).astype(str)))
z = pd.read_excel(os.path.join(KZ, "Zhou", "mmc2.xls"), header=0)
note = z.iloc[:, 6].astype(str)
z = z[~note.str.contains("Not expressed|No expression data|Expression unknown")]
zid_raw = z.iloc[:, 2].astype(str).tolist()
for name, ids, fallback in [("Konig", kid, None), ("Zhou", zid_raw, z.iloc[:, 1].astype(str).tolist())]:
    out, miss = [], []
    for i, gid in enumerate(ids):
        s = to_symbol(gid) if re.fullmatch(r"\d+(\.0)?", gid) else None
        if s: out.append(s)
        else: miss.append(gid if fallback is None else f"{fallback[i]} ({gid})")
    out = sorted(set(out))
    open(os.path.join(OUT, f"{name}_2008_HIV_dependency_factors.symbols.txt"), "w").write("\n".join(out) + "\n")
    rep.append(dict(screen=name, published=len(ids), mapped_symbols=len(out), unmapped=len(miss), unmapped_ids="; ".join(miss)))
    print(name, len(ids), "->", len(out), "symbols; unmapped:", miss)
pd.DataFrame(rep).to_csv(os.path.join(H, "s11a_gene_list_mapping.tsv"), sep="\t", index=False)
