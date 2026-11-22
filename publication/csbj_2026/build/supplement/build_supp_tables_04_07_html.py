"""Supplementary Tables 4-7 as HTML (from the tab-separated / text versions in the same folder).
Usage: python build_supp_tables_04_07_html.py <folder with the .txt tables> <css file>"""
import sys, os, re, html
import pandas as pd
D, CSSF = sys.argv[1], sys.argv[2]
CSS = open(CSSF).read() + """
table.data td,table.data th{font-size:12px;padding:5px 7px;vertical-align:top}td.num{text-align:right;white-space:nowrap;font-variant-numeric:tabular-nums}
td.genes{font-size:11px;color:#3a3f47;max-width:420px}tr.sig td{background:#f2f8f5}tr.notkept td{color:#8a909a}
h2 .count{font-weight:400;color:var(--muted);font-size:14px}a{color:#1f5fae;text-decoration:none}
td.sec{background:#f6f7f9;font-weight:600}@media print{.tbl{overflow:visible}table{break-inside:auto}tr{break-inside:avoid}h2,h3{break-after:avoid}}table.metrics td{font-size:13px}
"""
e = lambda s: html.escape(str(s), quote=False)
def fmt_p(x):
    try: x = float(x)
    except: return e(x)
    return "0" if x == 0 else (f"{x:.2e}" if x < 1e-3 else f"{x:.4f}".rstrip("0").rstrip("."))
def page(num, title, lead, body):
    return (f"<!doctype html><html lang='en'><head><meta charset='utf-8'><title>Supplementary Table {num} — {e(title)}</title>"
            f"<style>{CSS}</style></head><body><main><p class='kicker'>Supplementary Table {num}</p><h1>{e(title)}</h1>"
            f"<p class='lead'>{lead}</p>{body}</main></body></html>")
def term_link(t, url):
    return f"<a href='{e(url)}'>{e(t)}</a>" if isinstance(url, str) and url.startswith("http") else e(t)
ORDER = ["H: Hallmark Gene Sets", "C2: BioCarta", "C2: KEGG Legacy", "C2: KEGG MEDICUS", "C2: Pathway Interaction Database",
         "C2: Reactome Pathways", "C2: WikiPathways", "C5: Gene Ontology: Biological Process", "C5: Gene Ontology: Cellular Component",
         "C5: Gene Ontology: Molecular Function"]
outs = {}
# ---------------- Table 4: ORA
o = pd.read_csv(f"{D}/Supp Table 04 - ORA results.txt", sep="\t")
body = (f"<div class='cards'><div class='card a'><b>{len(o)} terms with p &lt; 0.01</b><span>hypergeometric test (upper tail), at least 3 input genes</span></div>"
        f"<div class='card b'><b>{(o.FDR < 0.05).sum()} terms with FDR &lt; 0.05</b><span>Benjamini–Hochberg within each collection (shaded rows)</span></div></div>")
for c in ORDER:
    s = o[o.Collection == c]
    if s.empty: continue
    body += f"<h2>{e(c)} <span class='count'>{len(s)} terms, {(s.FDR < 0.05).sum()} with FDR &lt; 0.05</span></h2><div class='tbl'><table class='data'>"
    body += "<tr><th>Rank</th><th>Term</th><th>Overlap</th><th>p-value</th><th>FDR</th><th>Genes</th></tr>"
    for r in s.itertuples():
        body += (f"<tr class='{'sig' if r.FDR < 0.05 else ''}'><td class='num'>{r.Rank}</td><td>{term_link(r.Term, r.Description)}</td>"
                 f"<td class='num'>{e(r._5)}</td><td class='num'>{fmt_p(r._6)}</td><td class='num'>{fmt_p(r.FDR)}</td><td class='genes'>{e(r.Genes)}</td></tr>")
    body += "</table></div>"
outs["Supp Table 04 - ORA results.html"] = page(4, "ORA results for the HIV dependency factors (ten collections)",
    "Over-representation analysis of the 265 HIV dependency factors in the analysis background (23,391 genes), with the same gene sets "
    "(10–600 genes after restriction to the background), background and minimum overlap (3 genes) as iGEA. All terms with p &lt; 0.01 are "
    "listed; shaded rows passed FDR &lt; 0.05 (Benjamini–Hochberg within each collection, over all gene sets in the size range). "
    "Overlap = input genes in the term / genes in the term. Term names link to MSigDB.", body)
# ---------------- Table 5: iGEA
g = pd.read_csv(f"{D}/Supp Table 05 - iGEA results.txt", sep="\t")
kept = g["In network"].astype(str) == "True"
body = (f"<div class='cards'><div class='card a'><b>{len(g)} terms selected</b><span>iterative selection, p &lt; 0.01, at most 30 iterations per collection</span></div>"
        f"<div class='card b'><b>{kept.sum()} terms kept</b><span>empirical q &lt; 0.10 (10,000 random gene lists); they form the network</span></div></div>")
for c in ORDER:
    s = g[g.Collection == c]
    if s.empty: continue
    k = (s["In network"].astype(str) == "True").sum()
    body += f"<h2>{e(c)} <span class='count'>{len(s)} selected, {k} kept</span></h2><div class='tbl'><table class='data'>"
    body += ("<tr><th>Iter.</th><th>Term</th><th>Overlap</th><th>p-value</th><th>q-value</th><th>Kept</th><th>Comm.</th>"
             "<th>Genes assigned (removed for the next iteration)</th><th>Full-list overlap</th><th>Full-list p</th><th>Full-list FDR</th></tr>")
    for r in s.itertuples(index=False):
        rr = r._asdict() if hasattr(r, "_asdict") else None
    for _, r in s.iterrows():
        kk = str(r["In network"]) == "True"
        body += (f"<tr class='{'' if kk else 'notkept'}'><td class='num'>{r['Iteration']}</td><td>{term_link(r['Term'], r['Description'])}</td>"
                 f"<td class='num'>{e(r['Iteration overlap'])}</td><td class='num'>{fmt_p(r['Iteration p-value'])}</td>"
                 f"<td class='num'>{fmt_p(r['q-value (empirical FDR)'])}</td><td>{'yes' if kk else 'no'}</td>"
                 f"<td>{e(r['Community']) if isinstance(r['Community'], str) else '–'}</td><td class='genes'>{e(r['Genes removed for next iteration'])}</td>"
                 f"<td class='num'>{e(r['Full list overlap'])}</td><td class='num'>{fmt_p(r['Full list p-value'])}</td><td class='num'>{fmt_p(r['Full list FDR (BH)'])}</td></tr>")
    body += "</table></div>"
outs["Supp Table 05 - iGEA results.html"] = page(5, "iGEA results for the HIV dependency factors (ten collections)",
    "Terms selected by iGEA for the 265 HIV dependency factors, in order of selection within each collection. Overlap and p-value refer to the "
    "iteration in which the term was selected (genes of earlier terms of the same collection removed); the full-list columns give the same term "
    "tested on the complete list, as in ORA. q-value: signal-adjusted empirical FDR from 10,000 random gene lists of the same size. Terms with "
    "q &lt; 0.10 are kept and form the network (gray rows: not kept). Comm.: community in the network (Leiden).", body)
# ---------------- Table 6: ORA vs iGEA as collections are added
L = open(f"{D}/Supp Table 06 - ORA vs iGEA as collections are added.txt").read().splitlines()
start = next(i for i, l in enumerate(L) if l.startswith("Dataset\t"))
notes = " ".join(l.strip() for l in L[1:start] if l.strip())
t6 = pd.read_csv(f"{D}/Supp Table 06 - ORA vs iGEA as collections are added.txt", sep="\t", skiprows=start)
cols = [c for c in t6.columns if c not in ("Dataset", "Input genes")]
body = "<div class='tbl'><table class='data'><tr>" + "".join(f"<th>{e(c)}</th>" for c in cols) + "</tr>"
prev = None
for _, r in t6.iterrows():
    if r["Dataset"] != prev:
        body += f"<tr><td class='sec' colspan='{len(cols)}'>{e(r['Dataset'])} ({r['Input genes']} input genes)</td></tr>"; prev = r["Dataset"]
    cls = "sig" if r["Method"] == "iGEA" else ""
    body += f"<tr class='{cls}'>" + "".join(f"<td class='{'num' if isinstance(v, (int, float)) and not isinstance(v, bool) else ''}'>{e(v)}</td>" for v in r[cols].values) + "</tr>"
body += "</table></div><p class='note'>Shaded rows: iGEA.</p>"
outs["Supp Table 06 - ORA vs iGEA as collections are added.html"] = page(6, "ORA and iGEA networks as gene-set collections are added", e(notes), body)
# ---------------- Table 7: network statistical report
T = open(f"{D}/Supp Table 07 - Network statistical report.txt").read()
parts = re.split(r"#{20,}\n(PART [A-Z]\..*?)\n#{20,}\n", T)
body = ""
for k in range(1, len(parts), 2):
    title, txt = parts[k], parts[k + 1]
    m = lambda pat: (re.search(pat, txt).group(1) if re.search(pat, txt) else "")
    body += f"<h2>{e(title)}</h2>"
    body += (f"<p class='note'><b>Input:</b> {e(m(r'Validated input genes: (\d+)'))} genes. <b>Largest cluster:</b> {e(m(r'Largest cluster: (.*)'))}.<br>"
             f"<b>Random lists:</b> {e(m(r'Random lists: (.*?),'))}, {e(m(r'each of (\d+) genes'))} genes each, from the analysis background "
             f"({e(m(r'background \((\d+) genes\)'))} genes); seed {e(m(r'seed (\d+)'))}. Parameters: {e(m(r'  (p_value_method=.*)'))}.<br>"
             f"<b>Term level:</b> a random list of this size selects on average {e(m(r'selects on average ([\d.]+) terms'))} terms by chance; "
             f"terms kept (q &lt; 0.1): {e(m(r'Terms kept \(q < 0.1\): (\d+ of \d+)'))}.</p>")
    mt = re.search(r"Metric .*?\n\s*-+\n(.*?)\n\n", txt, re.S).group(1).splitlines()
    body += "<h3>Network level (largest connected cluster)</h3><div class='tbl'><table class='data metrics narrow'><tr><th>Metric</th><th>Observed</th><th>Random mean</th><th>Random SD</th><th>Random max</th><th>z</th><th>Random ≥ observed</th><th>Empirical p</th></tr>"
    for l in mt:
        f = re.match(r"\s*(.+?)\s{2,}(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$", l).groups()
        p = float(f[7]); ps = "&lt; 0.0001" if p < 1.0001e-4 else fmt_p(p)
        body += f"<tr><td>{e(f[0])}</td>" + "".join(f"<td class='num'>{e(v)}</td>" for v in f[1:7]) + f"<td class='num'>{ps}</td></tr>"
    body += "</table></div><p class='note'>z = (observed − random mean)/SD; empirical p = (r + 1)/(n + 1), r = random lists with a value ≥ observed. The metrics describe the same cluster and are not independent tests.</p>"
    tl = txt.split("\nTERMS\n", 1)[1].splitlines()[1:]
    body += "<h3>Terms</h3><div class='tbl'><table class='data'><tr><th>Collection</th><th>Iter.</th><th>p-value</th><th>q-value</th><th>Kept</th><th>Term</th></tr>"
    for l in tl:
        mm = re.match(r"\s+(.+?)\s{2,}(\d+)\s+(\S+)\s+(\S+)\s+(yes|no)\s+(\S+)", l)
        if not mm: continue
        c, it, p, q, kp, term = mm.groups()
        body += f"<tr class='{'' if kp == 'yes' else 'notkept'}'><td>{e(c)}</td><td class='num'>{it}</td><td class='num'>{e(p)}</td><td class='num'>{e(q)}</td><td>{kp}</td><td>{e(term)}</td></tr>"
    body += "</table></div>"
outs["Supp Table 07 - Network statistical report.html"] = page(7, "Statistical report of the iGEA network (per-run random gene lists)",
    "The report produced by iGEA for the HIV dependency factors: (A) with the ten collections (the network of Supplementary Figure 2) and (B) with "
    "the three Gene Ontology collections, with random lists analyzed with the same three collections (the network of Figure 2). Each random list "
    "was analyzed with the same collections and parameters as the input and filtered with the same term-level q-value before its network was built. "
    "The iGEA software names the test \"Fisher's Exact Test\" in its report; it is the one-sided test, identical to the hypergeometric test (upper tail) described in the Methods.", body)
for n, h in outs.items():
    open(os.path.join(D, n), "w").write(h); print("written", n, len(h))
