"""Predefined theme scoring for the ORA vs iGEA LLM comparison (defined 20 Nov, before the 10-collection ORA run).
Themes = the seven communities of the HIV 10-collection iGEA network; marker genes = their hub genes.
For each network: which themes it contains (>= 1 marker gene in the network). For each answer (first turn):
which themes it discusses (>= 1 marker gene named), genes named, words, words in the community sections."""
import glob, re, os, sys
import pandas as pd
THEMES = {
    "Golgi vesicle transport": ["GOSR2", "RAB2A", "TMED2", "ARF1", "COG2", "COG3", "COG4", "RAB1B"],
    "EGFR/AKT/NF-kB signalling, endocytosis": ["EGFR", "EGF", "AKT1", "RELA", "IKBKG", "AP2M1"],
    "Mediator transcription": ["MED14", "MED4", "MED6", "MED7", "MED28"],
    "Nuclear pore and export": ["NUP107", "NUP133", "NUP155", "NUP160", "NUP85", "NUP153", "RANBP2"],
    "RNA exosome, RNA processing": ["EXOSC3", "EXOSC5", "MPHOSPH6"],
    "ER glycosylation, GPI anchor": ["PIGY", "PIGH", "PIGK", "DPM1", "MGAT1", "STT3A", "DDOST"],
    "Cullin-RING ubiquitin ligases": ["FBXW11", "DCAF1", "KLHDC2", "KBTBD7", "ANKRD9"],
}
RUNS = {"ORA GO": "ora_go_fdr0.05", "iGEA GO": "hiv_go_only/HIV_dependency_factors.symbols_*",
        "ORA 10 collections": "ora_all_fdr0.05", "iGEA 10 collections": "HIV_dependency_factors.symbols_*"}
LLM = os.environ.get("LLM_DIR", "llm_gpt-6.2-sol")
rows = []
for name, pat in RUNS.items():
    d = [x for x in glob.glob(pat) if os.path.isfile(x + "/network_llm.dot")][0]
    dot = open(d + "/network_llm.dot").read()
    genes = set(re.findall(r'^\s*"([^"]+)" \[type=gene', dot, re.M))
    present = [t for t, m in THEMES.items() if genes & set(m)]
    rows.append(dict(network=name, answer="(network)", themes=len(present), of=len(THEMES),
                     detail="; ".join(present), genes_named=len(genes), words=None))
    for f in sorted(glob.glob(f"{d}/{LLM}/response_r*.md")):
        t = open(f).read()
        named = {g for g in genes if re.search(rf"(?<![A-Za-z0-9-]){re.escape(g)}(?![A-Za-z0-9])", t)}
        disc = [th for th in present if named & set(THEMES[th])]
        rows.append(dict(network=name, answer=os.path.basename(f)[:-3], themes=len(disc), of=len(present),
                         detail="; ".join(th for th in present if th not in disc) or "-",
                         genes_named=len(named), words=len(t.split())))
o = pd.DataFrame(rows); o.to_csv(f"theme_scores_{LLM}.tsv", sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 90)
print(o.to_string(index=False))
