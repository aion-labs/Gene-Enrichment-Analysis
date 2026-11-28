"""Supplementary Table 14 (HTML): iGEA networks of three independent genome-wide siRNA screens for HIV dependency factors.
Inputs (read from publication/csbj_2026/analyses): the Brass network (step09_end_to_end, HIV run), the König and Zhou networks
and the comparison tables (step11_three_screens: s11a mapping, s11b overlap and community similarity, s11d LLM checks, s11e ORA,
s11h annotation-matched null with Benjamini-Hochberg correction, K = 10/25/50). R2 (Dec 2026): evidence = matched null + BH.
Usage: python build_supp_table_14_three_screens.py [<analyses dir>] [<out dir>]"""
import sys, os, re, json, html
import pandas as pd
import glob
_A = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "analyses")
_O = sys.argv[2] if len(sys.argv) > 2 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "supplement")
OUT = os.path.join(_O, "Supp Table 14 - Three independent HIV screens.html")
B = os.path.join(_A, "step11_three_screens")
RUN = {"Brass": glob.glob(os.path.join(_A, "step09_end_to_end", "HIV_dependency_factors.symbols_2*"))[0],
       "König": glob.glob(os.path.join(B, "Konig_2008_HIV_dependency_factors.symbols_2*"))[0],
       "Zhou": glob.glob(os.path.join(B, "Zhou_2008_HIV_dependency_factors.symbols_2*"))[0]}
S = ["Brass", "König", "Zhou"]; F = {"Brass": "Brass", "König": "Konig", "Zhou": "Zhou"}
REF = {"Brass": "Brass et al. 2008, Science", "König": "König et al. 2008, Cell", "Zhou": "Zhou et al. 2008, Cell Host Microbe"}
NAME = {("Brass", "C1"): "Golgi transport", ("Brass", "C2"): "EGFR/AKT/NF-κB signaling, endocytosis", ("Brass", "C3"): "Mediator",
        ("Brass", "C4"): "Nuclear pore, export", ("Brass", "C5"): "RNA exosome, RNA processing", ("Brass", "C6"): "ER glycosylation, GPI anchor",
        ("Brass", "C7"): "Cullin-RING ligases", ("König", "C1"): "mRNA splicing", ("König", "C2"): "Proteasome",
        ("König", "C3"): "Mediator, transcription coregulators", ("König", "C4"): "Nuclear pore, import", ("König", "C5"): "RNA polymerase II, DNA repair",
        ("König", "C6"): "Nucleases", ("König", "C7"): "Tyrosine phosphatases", ("König", "C8"): "Perinuclear ER",
        ("Zhou", "C1"): "NF-κB/AKT signaling, autophagy", ("Zhou", "C2"): "RNA polymerase II, entry receptors", ("Zhou", "C3"): "Mediator",
        ("Zhou", "C4"): "Mitochondrial complex I", ("Zhou", "C5"): "K48 ubiquitination"}
THEMES = [("Mediator / transcription coregulators", "#2a78d6", {"Brass": ["C3"], "König": ["C3"], "Zhou": ["C3"]}),
          ("Ubiquitin–proteasome", "#1baf7a", {"Brass": ["C7"], "König": ["C2"], "Zhou": ["C5"]}),
          ("RNA processing and RNA polymerase", "#6250d6", {"Brass": ["C5"], "König": ["C1", "C5"], "Zhou": ["C2"]}),
          ("Nuclear pore and transport", "#eb6834", {"Brass": ["C4"], "König": ["C4"], "Zhou": []}),
          ("NF-κB / AKT signaling", "#c2416b", {"Brass": ["C2"], "König": [], "Zhou": ["C1"]})]
SPECIFIC = {"Brass": ["C1", "C6"], "König": ["C6", "C7", "C8"], "Zhou": ["C4"]}
comm = {s: pd.read_csv(f"{RUN[s]}/network_modules_communities.tsv", sep="\t").set_index("Community") for s in S}
rep = {s: open(glob.glob(f"{RUN[s]}/*_statistical_report.txt")[0]).read() for s in S}
summ = {s: dict(l.rstrip("\n").split("\t") for l in open(f"{RUN[s]}/network_modules_summary.tsv")) for s in S}
mp = pd.read_csv(f"{B}/s11a_gene_list_mapping.tsv", sep="\t").set_index("screen")
ov = pd.read_csv(f"{B}/s11b_pairwise_overlap.tsv", sep="\t")
cs = pd.read_csv(f"{B}/s11b_community_similarity.tsv", sep="\t")
ora = pd.read_csv(f"{B}/s11e_ora_community_similarity.tsv", sep="\t")
# annotation-matched null (s11h): K = 25 main, K = 10 and 50 sensitivity; BH within each method
KEYS = ["screen_1", "community_1", "screen_2", "community_2"]
mt = pd.read_csv(f"{B}/s11h_igea_pairs_matched.tsv", sep="\t")
for k in ("K10", "K50"):
    mt = mt.merge(pd.read_csv(f"{B}/s11h_igea_pairs_matched_{k}.tsv", sep="\t")[KEYS + ["matched_BH_FDR"]].rename(columns={"matched_BH_FDR": f"FDR_{k}"}), on=KEYS)
mt = mt.merge(cs[KEYS + ["shared_terms"]], on=KEYS, how="left")
if mt.shared_terms.isna().any():   # pairs stored in the other orientation in s11b
    sw = cs.rename(columns={"screen_1": "screen_2", "community_1": "community_2", "screen_2": "screen_1", "community_2": "community_1"})
    mt = mt.drop(columns="shared_terms").merge(pd.concat([cs, sw])[KEYS + ["shared_terms"]].drop_duplicates(KEYS), on=KEYS, how="left")
assert len(mt) == 131 and mt.shared_terms.notna().all()
mo = pd.read_csv(f"{B}/s11h_ora_pairs_matched.tsv", sep="\t"); assert len(mo) == 16
N05, N10 = int((mt.matched_BH_FDR < 0.05).sum()), int((mt.matched_BH_FDR < 0.10).sum())
ORA10 = int((mo.matched_BH_FDR < 0.10).sum()); NZERO = int((mt.footprint_jaccard == 0).sum())
def fdr_of(s1, c1, s2, c2):
    m = mt[((mt.screen_1 == s1) & (mt.community_1 == c1) & (mt.screen_2 == s2) & (mt.community_2 == c2)) |
           ((mt.screen_1 == s2) & (mt.community_1 == c2) & (mt.screen_2 == s1) & (mt.community_2 == c1))]
    return m.iloc[0]
def qtxt(q): return "&lt; 0.001" if q < 0.001 else f"{q:.3f}" if q < 0.1 else f"{q:.2f}"
def qcls(q): return "sig" if q < 0.05 else "bord" if q < 0.10 else "ns"
llm = json.load(open(f"{B}/s11d_llm_sections.json")); chk = pd.read_csv(f"{B}/s11d_llm_checks.tsv", sep="\t")
e = html.escape
def num(pat, s, cast=float):
    m = re.search(pat, rep[s]); return cast(m.group(1)) if m else None
def hubs(s, c): return re.sub(r" \(\d+\)", "", comm[s].loc[c, "Top hub genes (degree)"])
def match(s1, c1, s2, c2):
    m = cs[((cs.screen_1 == s1) & (cs.community_1 == c1) & (cs.screen_2 == s2) & (cs.community_2 == c2)) |
           ((cs.screen_1 == s2) & (cs.community_1 == c2) & (cs.screen_2 == s1) & (cs.community_2 == c1))]
    return m.iloc[0] if len(m) else None
def ptxt(p): return "&lt; 0.001" if p < 0.001 else f"{p:.3f}"
# ---------- screen table
srows = []
for s in S:
    pub = {"Brass": "275 (list of the paper)", "König": "295", "Zhou": "311 (expressed hits)"}[s]
    srows.append(f"<tr><td><b>{s}</b><div class='sub'>{REF[s]}</div></td><td class='n'>{pub}</td><td class='n'>{num(r'Validated input genes: (\d+)', s, int)}</td>"
                 f"<td class='n'>{re.search(r'Terms kept \(q < 0.1\): (\d+ of \d+)', rep[s]).group(1)}</td>"
                 f"<td class='n'>{num(r'Largest cluster: \d+ nodes, (\d+) genes', s, int)} / {num(r'genes, (\d+) terms', s, int)}</td>"
                 f"<td><span class='pb sig'>&lt; 0.001</span></td><td class='n'>{summ[s]['n_communities']}</td><td class='n'>{float(summ[s]['modularity']):.2f}</td></tr>")
# ---------- theme grid
def cell(s, cs_, col):
    if not cs_: return "<td class='none'>—</td>"
    out = "".join(f"<div class='cc' style='border-left-color:{col}'><b>{c}</b> {e(NAME[(s, c)])}<div class='sub'>{e(hubs(s, c))}</div></div>" for c in cs_)
    return f"<td>{out}</td>"
grid = []
for t, col, m in THEMES:
    ev = []
    for a, b in [("Brass", "König"), ("Brass", "Zhou"), ("König", "Zhou")]:
        best = None
        for c1 in m[a]:
            for c2 in m[b]:
                r = match(a, c1, b, c2)
                if r is not None and (best is None or r.footprint_jaccard > best.footprint_jaccard): best = r
        if best is not None:
            q = fdr_of(a, best.community_1 if best.screen_1 == a else best.community_2, b, best.community_2 if best.screen_2 == b else best.community_1).matched_BH_FDR
            ev.append(f"{a[0]}–{b[0]} {best.footprint_jaccard:.2f} <span class='{'' if q < 0.05 else 'muted'}'>(FDR {qtxt(q)})</span>")
    n = sum(1 for s in S if m[s])
    grid.append(f"<tr><td class='th' style='border-left:5px solid {col}'><b>{e(t)}</b><div class='sub'>{n} of 3 screens</div></td>"
                + "".join(cell(s, m[s], col) for s in S) + f"<td class='ev'>{'<br>'.join(ev)}</td></tr>")
grid.append("<tr><td class='th' style='border-left:5px solid #9aa1ab'><b>Screen-specific</b><div class='sub'>no pair with FDR &lt; 0.10</div></td>"
            + "".join(cell(s, SPECIFIC[s], "#9aa1ab") for s in S) + "<td class='ev muted'>—</td></tr>")
# ---------- overlap
ovr = "".join(f"<tr><td>{r.pair}</td><td class='n'>{r.shared_genes}</td><td class='n'>{r.gene_jaccard:.3f}</td><td class='n'>{r.shared_terms}</td><td class='n'>{r.term_jaccard:.3f}</td></tr>" for r in ov.itertuples())
# ---------- matches
def prow(r):
    return (f"<tr><td>{r.screen_1} {r.community_1}<div class='sub'>{e(NAME[(r.screen_1, r.community_1)])}</div></td><td>{r.screen_2} {r.community_2}<div class='sub'>{e(NAME[(r.screen_2, r.community_2)])}</div></td>"
            f"<td class='n'>{r.footprint_jaccard:.3f}</td><td class='n'>{r.matched_median_jaccard:.3f}</td><td class='n'>{int(r.shared_terms)}</td><td class='n'>{ptxt(r.matched_p)}</td>"
            f"<td><span class='pb {qcls(r.matched_BH_FDR)}'>{qtxt(r.matched_BH_FDR)}</span></td></tr>")
mts = mt.sort_values(["matched_p", "footprint_jaccard"], ascending=[True, False])
mrows = "".join(prow(r) for r in mts[mts.matched_BH_FDR < 0.25].itertuples())

morows = "".join(f"<tr><td>{r.screen_1} {r.community_1}<div class='sub'>{int(r.terms_1)} terms</div></td><td>{r.screen_2} {r.community_2}<div class='sub'>{int(r.terms_2)} terms</div></td>"
                 f"<td class='n'>{r.footprint_jaccard:.3f}</td><td class='n'>{r.matched_median_jaccard:.3f}</td><td class='n'>{ptxt(r.matched_p)}</td><td><span class='pb {qcls(r.matched_BH_FDR)}'>{qtxt(r.matched_BH_FDR)}</span></td></tr>"
                 for r in mo.sort_values("matched_p").itertuples())
nfkb = fdr_of("Brass", "C2", "Zhou", "C1"); golgi = fdr_of("Brass", "C1", "König", "C4")
# ---------- ORA
oraT = {"Brass": (87, 6.37, 2), "König": (246, 17.32, 2), "Zhou": (38, 4.99, 3)}
igT = {s: (int(re.search(r"Terms kept \(q < 0.1\): (\d+)", rep[s]).group(1)), round(num(r"([\d.]+) edges per gene", s), 2), int(summ[s]["n_communities"])) for s in S}
orar = "".join(f"<tr><td><b>{s}</b></td><td class='n'>{oraT[s][0]}</td><td class='n'>{igT[s][0]}</td><td class='n'>{oraT[s][1]:.1f}</td><td class='n'>{igT[s][1]:.1f}</td>"
               f"<td class='n'>{oraT[s][2]}</td><td class='n'>{igT[s][2]}</td></tr>" for s in S)

# ---------- LLM
def first_sent(t, n=1):
    t = re.sub(r"\*\*", "", t.strip()); parts = re.split(r"(?<=[.!?])\s+", t); return " ".join(parts[:n])
def bold_heads(t): return [re.sub(r"\.$", "", h) for h in re.findall(r"\*\*([^*]{8,140})\*\*", t) if len(h) >= 25][:2]
lrows = []
for s in S:
    guesses = "".join(f"<li>{e(first_sent(llm[f'{s}|1|{r}']['context']))}</li>" for r in "123")
    rev = e(first_sent(llm[f"{s}|2|1"]["summary"], 2))
    new = "".join(f"<li>{e(h)}</li>" for h in bold_heads(llm[f"{s}|2|1"]["new"]) + bold_heads(llm[f"{s}|2|2"]["new"]))
    c = chk[chk.screen == s]
    lrows.append(f"<div class='llm'><h3>{s}</h3><div class='lg'><div><div class='lab'>Turn 1 · blind guess of the experiment (3 repeats)</div><ul>{guesses}</ul></div>"
                 f"<div><div class='lab'>Turn 2 · revised summary (repeat 1)</div><p>{rev}</p><div class='lab'>Turn 2 · what may be new (repeats 1–2)</div><ul>{new or '<li class=muted>see Supplementary Table 14B/C</li>'}</ul></div></div>"
                 f"<p class='sub'>Checks: every answer covers every community ({c.communities.iloc[0]}); genes named per answer {c.genes_named.min()}–{c.genes_named.max()}; "
                 f"genes not in the network: {', '.join(sorted(set(' '.join(c.outside.fillna('')).split()))) or 'none'}.</p></div>")
CSS = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "st16.css")).read()
doc = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Table 14 — Three independent HIV screens</title><style>{CSS}</style></head><body><main>
<p class="kicker">Supplementary Table 14</p><h1>iGEA networks of three independent genome-wide siRNA screens for HIV dependency factors</h1>
<p class="lead">The hit lists of three genome-wide siRNA screens published in 2008 overlap little at the gene level [Bushman et al. 2009]. Each list was
analyzed separately with iGEA (ten collections, default settings, 10,000 random gene lists) and the resulting networks were compared at the
level of genes, terms and communities. Each cross-screen pair of communities was tested against an annotation-matched null and the
p-values were corrected for multiple testing (Benjamini–Hochberg over all pairs; section 4). The same comparison was made with ORA. The similarity of all communities is shown in Supplementary Figure 6.</p>
<div class="cards"><div class="card c"><b>Genes: 2–3% shared</b><span>Pairwise Jaccard 0.017–0.029; three genes in all screens (MED6, MED7, RELA).</span></div>
<div class="card a"><b>Terms: 4–8% shared</b><span>Pairwise Jaccard 0.035–0.080; two terms in all screens, both Mediator.</span></div>
<div class="card b"><b>Communities: two themes in all three screens</b><span>{N05} of 131 community pairs match at FDR &lt; 0.05 (Mediator; ubiquitin–proteasome), {N10 - N05} more at FDR &lt; 0.10. ORA: {ORA10} of 16 pairs.</span></div></div>
<h2>1. Networks of the three screens</h2><div class="tbl"><table><thead><tr><th>Screen</th><th class="n">Published hits</th><th class="n">Input genes<br><span>validated, in background</span></th>
<th class="n">Terms kept<br><span>q &lt; 0.10</span></th><th class="n">Largest cluster<br><span>genes / terms</span></th><th>Network p<br><span>cluster size</span></th><th class="n">Communities</th><th class="n">Modularity</th></tr></thead><tbody>{''.join(srows)}</tbody></table></div>
<h2>2. Shared and screen-specific communities</h2>
<p class="note">Each cell lists the community of that screen's network (ID, short name, top hub genes). Themes group communities with similar footprints (Supplementary Figure 6). Evidence: footprint Jaccard of the best pair for each pair of screens and its FDR (section 4); FDR &lt; 0.05 in black.</p>
<div class="tbl"><table class="grid"><thead><tr><th>Theme</th><th>Brass</th><th>König</th><th>Zhou</th><th>Evidence (footprint Jaccard)</th></tr></thead><tbody>{''.join(grid)}</tbody></table></div>
<h2>3. Overlap of genes and terms</h2><div class="tbl"><table class="narrow"><thead><tr><th>Screens</th><th class="n">Shared input genes</th><th class="n">Gene Jaccard</th><th class="n">Shared kept terms</th><th class="n">Term Jaccard</th></tr></thead><tbody>{ovr}</tbody></table></div>
<h2>4. Community pairs between screens</h2>
<p class="note">The footprint of a community is the union of the full gene sets of its terms (all background genes annotated to them), so the comparison is between the biology the communities describe, not between the screen hits; similarity is the Jaccard index of two footprints. To assess whether two communities overlap more than expected by chance, their Jaccard was compared with that of 50,000 pairs of random communities with the same number of terms, each term replaced by a gene set of the same collection and similar size (one of the 25 closest in size). All 131 cross-screen pairs (7 × 8 + 7 × 5 + 8 × 5) were tested and the p-values corrected with the Benjamini–Hochberg procedure. {N05} pairs have FDR &lt; 0.05 and {N10} FDR &lt; 0.10; results were similar with size windows of 10 and 50 gene sets. Pairs are not independent (a community takes part in several pairs), and hits that recur across siRNA screens in general may contribute to shared themes. All pairs and the code: publication/csbj_2026/analyses/step11_three_screens (s11h).</p>
<div class="tbl"><table><thead><tr><th>Community</th><th>Community</th><th class="n">Footprint Jaccard</th><th class="n">Chance<br><span>median of random pairs</span></th><th class="n">Identical terms</th><th class="n">p</th><th>FDR</th></tr></thead><tbody>{mrows}</tbody></table></div>
<p class="note">Pairs with FDR &lt; 0.25, ordered by p-value.</p>
<h2>5. The same comparison with ORA</h2>
<p class="note">ORA with the same gene sets, background and minimum overlap (Benjamini-Hochberg FDR &lt; 0.05 within each collection), communities with the same procedure.</p>
<div class="tbl"><table class="narrow"><thead><tr><th>Screen</th><th class="n">Terms<br><span>ORA</span></th><th class="n">Terms<br><span>iGEA</span></th><th class="n">Edges per gene<br><span>ORA</span></th><th class="n">Edges per gene<br><span>iGEA</span></th><th class="n">Communities<br><span>ORA</span></th><th class="n">Communities<br><span>iGEA</span></th></tr></thead><tbody>{orar}</tbody></table></div>
<p class="note">ORA terms overlap more between screens (Jaccard 0.05–0.14; six terms in all three), because the same broad and nested terms recur, but each ORA network collapses into two or three large communities, giving 16 cross-screen pairs. Tested in the same way, none of the 16 pairs matched (FDR &lt; 0.10): these large communities overlap as much as random communities of the same size (for example Brass C2, 54 terms, with König C1, 124 terms: Jaccard 0.305, chance 0.307). Because the two methods yield different numbers of pairs, the counts of matching pairs are not compared.</p>
<h2>6. LLM interpretation of the three networks</h2>
<p class="note">Same model (gpt-6.2-sol), prompt (Supplementary Table 3) and protocol as for the other networks; in turn 2 the same sentence described the experiment for all three screens: “an siRNA screen in human cells for host factors required for HIV infection (HIV dependency factors)”. Full answers: Supplementary Table 10 (Brass) and 14B–C (König, Zhou).</p>
{''.join(lrows)}
<h2>Notes</h2><ul class="notes">
<li><b>Gene lists.</b> Brass: the list used throughout the paper. König: Table S3 of König et al. 2008 (295 Entrez Gene IDs; 292 current symbols). Zhou: Table S1 of Zhou et al. 2008, the 390 primary hits minus those not expressed in activated T cells or macrophages (the published 311 host factors; 311 current symbols). Identifiers mapped to current NCBI symbols through Entrez Gene IDs (gene_info, gene_history); input genes outside the background are removed by iGEA.</li>
<li><b>iGEA.</b> Hypergeometric test (upper tail), p &lt; 0.01, overlap ≥ 3, gene sets of 10–600 genes, ≤ 30 iterations, ten collections, 10,000 random gene lists per screen (seed 2026), terms kept at empirical q &lt; 0.10; Leiden communities on the term network (best of 100 seeds).</li>
<li><b>Themes.</b> Themes are clusters of the footprint similarity matrix that contain communities from at least two screens (Supplementary Figure 6); statistical support is given per pair (section 4). Mediator and ubiquitin–proteasome match in all three screens (FDR &lt; 0.05 for every pair); nuclear pore (Brass–König) and RNA polymerase II (König–Zhou) match at FDR 0.052; RNA processing (Brass C5 with König C1 and Zhou C2) only nominally (FDR 0.09–0.13); NF-κB/AKT (Brass C2–Zhou C1) shares hub genes (AKT1, RELA) but not more annotation than expected by chance (FDR {nfkb.matched_BH_FDR:.2f}). The Brass Golgi community resembles the König nuclear pore community (Jaccard 0.19) through small-GTPase terms, but not beyond chance (FDR {golgi.matched_BH_FDR:.2f}).</li>
<li><b>Comparison with Bushman et al. 2009.</b> The meta-analysis found shared functional categories by annotating the pooled screens; iGEA recovers the same shared themes (Mediator, proteasome, RNA polymerase and splicing, nuclear pore) from each screen separately, as communities, and keeps screen-specific biology visible. ER/Golgi transport, shared in the meta-analysis, is a community only in the Brass network.</li></ul>
<footer>Code and data: iGEA repository, publication/csbj_2026/analyses/step11_three_screens (s11a gene lists, s11b comparison, s11h annotation-matched null and multiple testing, s11i similarity heatmap (Supplementary Figure 6), s11d LLM checks, s11e ORA); this page: publication/csbj_2026/build/supplement/build_supp_table_14_three_screens.py.</footer>
</main></body></html>"""
open(OUT, "w").write(doc); print("written", OUT)
