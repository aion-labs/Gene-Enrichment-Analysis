"""Reader-friendly supplementary documents for the LLM interpretations (one HTML per dataset, later printed to PDF).

Each document: key facts of the network, the community table, the first answer (experiment unknown to the model)
and the second turn (experiment revealed) for repeat 1, then a consistency table across repeats and repeats 2-3
in an appendix. Gene symbols are highlighted: coloured by community when they are in the network, dashed grey when
the model names a gene that is not in the network (outside knowledge).

Usage: [LLM_DIR=llm_gpt-6.2-sol] python build_llm_supp.py <step09_end_to_end folder> <output folder>
"""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import csv
import glob
import gzip
import html
import os
import re
import sys
from datetime import date

import markdown

REPO = (str(_REPO) + "/")
LLM_DIR = os.environ.get("LLM_DIR", "llm")   # answer folder in each run folder (one per model)

NUM = {'LLM_HIV_GO_network': ('9', 'Supp Table 09 - LLM answers, HIV iGEA network, GO collections'), 'LLM_HIV_10_collections': ('10', 'Supp Table 10 - LLM answers, HIV iGEA network, ten collections'), 'LLM_HIV_GO_ORA_baseline': ('11A', 'Supp Table 11A - LLM answers, HIV ORA network, GO collections'), 'LLM_HIV_ORA_10_collections': ('11B', 'Supp Table 11B - LLM answers, HIV ORA network, ten collections'), 'LLM_Random_control': ('12', 'Supp Table 12 - LLM answers, random control list'), 'LLM_Hypoxia': ('15', 'Supp Table 15 - LLM answers, hypoxia signature'), 'LLM_Vaccination': ('16', 'Supp Table 16 - LLM answers, vaccination signature'), 'LLM_Konig_screen': ('14B', 'Supp Table 14B - LLM answers, König screen'), 'LLM_Zhou_screen': ('14C', 'Supp Table 14C - LLM answers, Zhou screen')}  # stem -> (Supplementary Table label, output file name)

DATASETS = [  # (file stem, title, subtitle, run folder patterns)
    ("LLM_HIV_GO_network", "HIV dependency factors: GO network (Figure 2)",
     "iGEA network built from the three Gene Ontology collections", ["hiv_go_only/HIV_dependency_factors.symbols_*"]),
    ("LLM_HIV_GO_ORA_baseline", "HIV dependency factors: GO network from ORA (baseline)",
     "ORA network of the three Gene Ontology collections (BH FDR < 0.05), same prompt as for iGEA", ["ora_go_fdr0.05"]),
    ("LLM_HIV_ORA_10_collections", "HIV dependency factors: ORA network, 10 collections (baseline)",
     "ORA network of all ten gene-set collections (BH FDR < 0.05), same prompt as for iGEA", ["ora_all_fdr0.05"]),
    ("LLM_Konig_screen", "HIV dependency factors, König et al. 2008: 10 collections",
     "iGEA network of the König screen (Supplementary Table 14), all ten gene-set collections", ["../step11_three_screens/Konig_2008_HIV_dependency_factors.symbols_*"]),
    ("LLM_Zhou_screen", "HIV dependency factors, Zhou et al. 2008: 10 collections",
     "iGEA network of the Zhou screen (Supplementary Table 14), all ten gene-set collections", ["../step11_three_screens/Zhou_2008_HIV_dependency_factors.symbols_*"]),
    ("LLM_HIV_10_collections", "HIV dependency factors: 10 collections",
     "iGEA network built from all ten gene-set collections", ["HIV_dependency_factors.symbols_*"]),
    ("LLM_Hypoxia", "Hypoxia signature: 10 collections",
     "iGEA network built from all ten gene-set collections", ["hypoxia-genes.symbols_*"]),
    ("LLM_Vaccination", "Influenza vaccination signature: 10 collections",
     "iGEA network built from all ten gene-set collections", ["vaccination_NAKAYA_HAI_28DY_3DY_POSITIVE_*"]),
    ("LLM_Random_control", "Random gene list (negative control)",
     "265 genes drawn at random from the background (264 validated); networks with and without the term-level filter",
     ["random_control_q0.10/random_control_*", "random_control_q1/random_control_*"]),
]

PALETTE = [("#E3EEFB", "#1F5FA8"), ("#FBE4EC", "#A3245A"), ("#E4F3E6", "#2C7A3A"), ("#FFF0DE", "#B5560B"),
           ("#EFE6F6", "#6A3796"), ("#DDF3F4", "#16777E"), ("#FBF4D6", "#866A00"), ("#EFE9E6", "#6B4B3E"),
           ("#E8ECEF", "#3E4C57"), ("#F6E3E1", "#9B3B2E")]

STOP = set("""GO ER RNA DNA ATP NOT SET CAT MET CAN AND FOR ARE HAS WAS ITS THE ALL ANY ONE TWO MAX MIN PC CD IL TNF MHC
HLA GPI ECM ROS RNS HIF NF ERK MAPK PI3K TOR AKT NFKB KEGG PID WP RNAI CRISPR SEQ FDR NA ID PH CO DE LOW IMPACT SHE BAD
CLOCK HIV HAI PBMC SIRNA MRNA COPII COPI TGN NPC ESCRT TCR BCR NK HDF HDFS DC PCR""".split())

METRIC_LABEL = {"Cluster Size": "Nodes in largest cluster", "Number of Genes": "Genes in largest cluster",
                "Number of Terms": "Terms in largest cluster", "Number of Edges": "Edges in largest cluster",
                "Average Edges per Gene": "Average edges per gene"}


def gene_symbols():
    s = set()
    with gzip.open(REPO + "data/recent_release/Homo_sapiens.gene_info.gz", "rt") as f:
        next(f)
        for line in f:
            s.add(line.split("\t")[2])
    return s


SYM = gene_symbols()


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def fmt_p(r, p, n):
    return f"&lt; {1 / n:.0g}".replace("e-0", "e-") if int(r) == 0 else f"{float(p):.2g}"


def load_run(d):
    """Everything needed from one end-to-end run folder."""
    run = {"dir": d}
    run["communities"] = read_tsv(f"{d}/network_modules_communities.tsv")
    run["genes"] = {r["Gene"]: r for r in read_tsv(f"{d}/network_modules_genes.tsv")}
    run["summary"] = dict(r for r in csv.reader(open(f"{d}/network_modules_summary.tsv"), delimiter="\t"))
    run["ora"] = os.path.exists(f"{d}/ora_terms.tsv")       # ORA baseline: no random-list statistics
    if run["ora"]:
        terms = read_tsv(f"{d}/ora_terms.tsv")
        run["stats"], run["n_selected"], run["n_kept"] = None, len(terms), len(terms)
        run["n_genes"] = len({g.strip() for t in terms for g in t["Genes"].split(",")})
        run["input_genes"], run["n_lists"], run["q"] = 265, 0, "FDR < 0.05"
    else:
        stats = glob.glob(f"{d}/*_statistics.tsv")[0]
        run["stats"] = read_tsv(stats)
        terms = read_tsv(glob.glob(f"{d}/*_term_qvalues.tsv")[0])
        run["n_selected"] = len(terms)
        run["n_kept"] = sum(t["Kept"] == "True" for t in terms)
        rep = open(glob.glob(f"{d}/*_statistical_report.txt")[0]).read()
        run["input_genes"] = int(re.search(r"Validated input genes: (\d+)", rep).group(1))
        run["n_lists"] = int(re.search(r"Random lists: (\d+)", rep).group(1))
        m = re.search(r"Terms kept \(q < ([\d.]+)\)", rep)
        run["q"] = f"{float(m.group(1)):.2f}" if m else "1 (no filter)"
    run["manifest"] = {r["repeat"]: r for r in read_tsv(f"{d}/{LLM_DIR}/manifest.tsv")}
    fm = f"{d}/{LLM_DIR}/followup_manifest.tsv"
    run["followup"] = {r["repeat"]: r for r in read_tsv(fm)} if os.path.exists(fm) else {}
    run["dot_terms"] = set(re.findall(r'label="([^"]+)"', open(f"{d}/network_llm.dot").read()))
    return run


def term_words(run):
    return {w for t in run["dot_terms"] for w in re.findall(r"[A-Z0-9]+", t)}


def comm_colour(run, c):
    ids = [x["Community"] for x in run["communities"]]
    return PALETTE[ids.index(c) % len(PALETTE)] if c in ids else ("#EEE", "#555")


def chip(run, c, text=None):
    bg, fg = comm_colour(run, c)
    return f'<span class="chip" style="background:{bg};color:{fg};border-color:{fg}33">{html.escape(text or c)}</span>'


def highlight(run, html_text):
    """Mark gene symbols and community IDs inside the text (never inside tags)."""
    words = term_words(run)
    comms = {x["Community"] for x in run["communities"]}
    out = []
    for piece in re.split(r"(<[^>]+>)", html_text):
        if piece.startswith("<"):
            out.append(piece)
            continue

        def gene(m):
            g = m.group(1)
            if g in run["genes"]:
                bg, fg = comm_colour(run, run["genes"][g]["Community"])
                br = " bridge" if run["genes"][g]["Bridge"] == "True" else ""
                return f'<span class="g{br}" style="background:{bg};color:{fg}">{g}</span>'
            if g in SYM and g not in STOP and g not in words and not re.fullmatch(r"C\d+", g):
                return f'<span class="g out" title="not in the network">{g}</span>'
            return g

        piece = re.sub(r"(?<![A-Za-z0-9-])([A-Z][A-Z0-9]{1,9}(?:-[A-Z0-9]+)?)(?![A-Za-z0-9]|-[A-Z0-9])", gene, piece)
        piece = re.sub(r"(?<![A-Za-z0-9-])(C\d{1,2})(?![A-Za-z0-9:])",
                       lambda m: chip(run, m.group(1)) if m.group(1) in comms else m.group(1), piece)
        out.append(piece)
    return "".join(out)


def md(run, text, shift=1):
    """Markdown -> HTML with headings shifted down (## -> h3) and genes highlighted."""
    h = markdown.markdown(text, extensions=["tables", "sane_lists"])
    for lvl in range(5, 0, -1):
        h = re.sub(rf"<(/?)h{lvl}>", rf"<\1h{min(lvl + shift, 6)}>", h)
    return highlight(run, h)


def read(path):
    return open(path).read().strip() if os.path.exists(path) else None


def facts(run, label=None):
    if run["ora"]:
        head = f'<div class="facts-label">{html.escape(label)}</div>' if label else ""
        cards = [("Input genes", f'{run["input_genes"]}', "validated, in the background"),
                 ("Terms in network", f'{run["n_kept"]}', "ORA, BH FDR &lt; 0.05"),
                 ("Network genes", f'{run["n_genes"]}', "each term with all its genes"),
                 ("Network p", "not tested", "ORA has no network test"),
                 ("Communities", run["summary"].get("n_communities", "–"), f'modularity {run["summary"].get("modularity", "–")}')]
        return head + '<div class="facts">' + "".join(
            f'<div class="card"><div class="k">{k}</div><div class="v">{v}</div><div class="s">{s_}</div></div>' for k, v, s_ in cards) + "</div>"
    s = {r["Metric"]: r for r in run["stats"]}
    cs = s["Cluster Size"]
    p = fmt_p(cs["Random >= Observed"], cs["Empirical p"], run["n_lists"])
    cards = [
        ("Input genes", f'{run["input_genes"]}', "validated, in the background"),
        ("Terms in network", f'{run["n_kept"]} / {run["n_selected"]}', f'kept at q &lt; {run["q"]}' if "no filter" not in run["q"] else "no term filter (q = 1)"),
        ("Largest cluster", f'{int(float(s["Number of Genes"]["Observed"]))} genes', f'{int(float(s["Number of Terms"]["Observed"]))} terms, {int(float(s["Number of Edges"]["Observed"]))} edges'),
        ("Network p", p, f'vs {run["n_lists"]:,} random lists'),
        ("Communities", run["summary"].get("n_communities", "–"), f'modularity {run["summary"].get("modularity", "–")}' if run["summary"].get("modularity") not in (None, "", "None") else "connected components"),
    ]
    head = f'<div class="facts-label">{html.escape(label)}</div>' if label else ""
    return head + '<div class="facts">' + "".join(
        f'<div class="card"><div class="k">{k}</div><div class="v">{v}</div><div class="s">{s_}</div></div>' for k, v, s_ in cards) + "</div>"


def num(x):
    x = float(x)
    return f"{int(x)}" if x == int(x) else f"{x:.2f}"


def stats_table(run):
    if run["ora"]:
        return ('<p class="small">Over-representation analysis of the full gene list (one-sided Fisher\'s exact test, '
                'overlap of at least 3 genes, Benjamini-Hochberg FDR within each collection). The prompt is the iGEA prompt '
                '(Supplementary Table 3) with only the method sentences changed: ORA introduction, FDR instead of the '
                'random-list statistics, and gene degree = number of terms. Compare with the iGEA network of the same '
                'collections (' + ('Supplementary Table 10' if 'ora_all' in run['dir'] else 'Supplementary Table 9') + ').</p>')
    rows = "".join(
        f'<tr><td>{METRIC_LABEL.get(r["Metric"], r["Metric"])}</td><td class="n">{num(r["Observed"])}</td>'
        f'<td class="n">{float(r["Null Mean"]):.2f}</td><td class="n">{fmt_p(r["Random >= Observed"], r["Empirical p"], run["n_lists"])}</td></tr>'
        for r in run["stats"])
    return ('<table class="stats"><thead><tr><th>Network metric</th><th>This network</th><th>Random mean</th><th>Empirical p</th></tr></thead>'
            f"<tbody>{rows}</tbody></table>")


def community_table(run):
    rows = []
    for c in run["communities"]:
        hubs = re.sub(r"([A-Z0-9-]+) \((\d+)\)", lambda m: f'{m.group(1)}<sub>{m.group(2)}</sub>', c["Top hub genes (degree)"])
        rows.append(f'<tr><td>{chip(run, c["Community"])}</td><td class="lead">{html.escape(c["Lead term"])}</td>'
                    f'<td class="n">{c["Terms"]}</td><td class="n">{c["Genes"]}</td><td>{highlight(run, hubs)}</td><td class="n">{c["Bridge genes"]}</td></tr>')
    return ('<table class="comm"><thead><tr><th>Community</th><th>Lead term (most genes)</th><th>Terms</th><th>Genes</th>'
            '<th>Top hub genes (degree)</th><th>Bridge genes</th></tr></thead><tbody>' + "".join(rows) + "</tbody></table>")


def meta_line(m, turn):
    when = m["started_utc"][:16].replace("T", " ")
    toks = f'{int(m["input_tokens"]):,} in / {int(m["output_tokens"]):,} out tokens' if m.get("input_tokens") else ""
    return (f'<div class="meta">{turn} · model <b>{m["returned_model"]}</b> (requested {m["requested_model"]}) · {when} UTC · '
            f'{float(m["seconds"]):.0f} s · {toks} · response id <code>{m["response_id"]}</code></div>')


def answer_block(run, rep, first_title="Turn 1 · interpretation without experimental information"):
    d = run["dir"]
    parts = []
    t1 = read(f"{d}/{LLM_DIR}/response_r{rep}.md")
    if t1:
        parts.append(f'<section class="turn"><h2 class="turn-h"><span class="tag t1">Turn 1</span>{first_title[9:] if first_title.startswith("Turn 1 · ") else first_title}</h2>'
                     + meta_line(run["manifest"][rep], f"Repeat {rep}") + f'<div class="answer">{md(run, t1)}</div></section>')
    elif os.path.exists(f"{d}/{LLM_DIR}/response_r{rep}.refused.json"):
        parts.append('<p class="refused">Turn 1 was declined by the provider\'s content policy; the request is archived.</p>')
    t2 = read(f"{d}/{LLM_DIR}/followup_r{rep}.md")
    if t2 and rep in run["followup"]:
        ctx = run["followup"][rep]["context"]
        parts.append(f'<section class="turn"><h2 class="turn-h"><span class="tag t2">Turn 2</span>after the experiment was revealed</h2>'
                     f'<div class="reveal"><div class="rl">Message sent to the model</div><p>Here is how the gene list was generated: '
                     f'<b>{html.escape(ctx)}</b></p><p class="small">The model was asked to rebuild its interpretation: what fits, '
                     'what does not, what may be new, and priorities (prompt in Supplementary Table 3).</p></div>'
                     + meta_line(run["followup"][rep], f"Repeat {rep}") + f'<div class="answer">{md(run, t2)}</div></section>')
    elif os.path.exists(f"{d}/{LLM_DIR}/followup_r{rep}.refused.json"):
        parts.append('<p class="refused">Turn 2 was declined by the provider\'s content policy; the request is archived.</p>')
    return "".join(parts)


def section_titles(text):
    """Community names chosen by the model: '### C1 — Name'."""
    out = {}
    for m in re.finditer(r"^#{2,4}\s*\**(C\d{1,2})\**\s*[—–:-]+\s*(.+?)\**\s*$", text or "", re.M):
        out.setdefault(m.group(1), m.group(2).strip("* "))
    return out


def first_sentences(text, heading, n=1):
    m = re.search(rf"^#+\s*\**{heading}.*?\n(.*?)(?=\n#+\s|\Z)", text or "", re.S | re.I | re.M)
    if not m:
        return None
    body = re.sub(r"\s+", " ", re.sub(r"[*_#>]", "", m.group(1))).strip()
    sents = re.split(r"(?<=[.!?])\s+(?=[A-Z])", body)
    return " ".join(sents[:n])


def consistency(run, reps):
    d = run["dir"]
    texts = {r: read(f"{d}/{LLM_DIR}/response_r{r}.md") for r in reps}
    names = {r: section_titles(t) for r, t in texts.items()}
    head = "".join(f"<th>Repeat {r}</th>" for r in reps)
    rows = []
    for c in run["communities"]:
        cid = c["Community"]
        rows.append(f"<tr><td>{chip(run, cid)}</td>" + "".join(f"<td>{html.escape(names[r].get(cid, '–'))}</td>" for r in reps) + "</tr>")
    words = "".join(f'<td class="n">{len((texts[r] or "").split()):,}</td>' for r in reps)
    rows.append(f'<tr class="sub"><td>Words</td>{words}</tr>')
    ctx = [first_sentences(texts[r], "Estimated experimental context") for r in reps]
    if any(ctx):
        rows.append('<tr class="sub"><td>Experimental context (first sentence)</td>' + "".join(
            f'<td class="small">{highlight(run, html.escape(x or "– (not proposed)"))}</td>' for x in ctx) + "</tr>")
    return (f'<table class="cons"><thead><tr><th>Community</th>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table>'
            '<p class="small">Community names are the model\'s own section headings in each independent repeat (same prompt, '
            'no memory between calls).</p>')


CSS = """
@page { size: A4; margin: 16mm 15mm 18mm 15mm; }
:root { --ink:#1d2329; --muted:#5f6b76; --rule:#d9dee3; --accent:#1F5FA8; }
* { box-sizing: border-box; }
body { font-family: Georgia, 'Source Serif Pro', 'DejaVu Serif', serif; color: var(--ink); font-size: 10.4pt; line-height: 1.5;
       max-width: 182mm; margin: 0 auto; padding: 18px 0 40px; background: white; }
h1, h2, h3, h4, h5, .card .k, .facts-label, .tag, th, .meta, .chip, .g, .rl, .small, .lead { font-family: 'Helvetica Neue', Helvetica, Arial, 'DejaVu Sans', sans-serif; }
.kicker { text-transform: uppercase; letter-spacing: .08em; font: 600 8.5pt 'Helvetica Neue', Arial, sans-serif; color: var(--accent); margin: 0 0 4px; }
h1 { font-size: 19pt; line-height: 1.2; margin: 0 0 4px; }
.subtitle { color: var(--muted); font-size: 10.5pt; margin: 0 0 14px; }
h2 { font-size: 13pt; margin: 22px 0 8px; padding-bottom: 4px; border-bottom: 1px solid var(--rule); }
.answer h3 { font-size: 11.5pt; margin: 16px 0 6px; color: #23303b; }
.answer h4 { font-size: 10.6pt; margin: 12px 0 4px; }
.answer p { margin: 0 0 8px; text-align: left; }
.answer ul, .answer ol { margin: 0 0 8px 18px; padding: 0; }
.answer li { margin: 0 0 3px; }
.facts { display: grid; grid-template-columns: repeat(5, 1fr); gap: 8px; margin: 6px 0 10px; }
.facts-label { font-size: 9pt; font-weight: 600; color: var(--muted); margin-top: 8px; }
.card { border: 1px solid var(--rule); border-radius: 8px; padding: 8px 10px; background: #fafbfc; }
.card .k { font-size: 7.6pt; text-transform: uppercase; letter-spacing: .05em; color: var(--muted); }
.card .v { font: 600 14pt 'Helvetica Neue', Arial, sans-serif; margin: 2px 0; }
.card .s { font-size: 7.8pt; color: var(--muted); line-height: 1.3; }
table { border-collapse: collapse; width: 100%; font-size: 8.8pt; margin: 6px 0 10px; }
th { text-align: left; font-size: 7.8pt; text-transform: uppercase; letter-spacing: .04em; color: var(--muted); border-bottom: 1.5px solid var(--ink); padding: 4px 6px; }
td { border-bottom: 1px solid var(--rule); padding: 4px 6px; vertical-align: top; }
td.n { text-align: right; font-variant-numeric: tabular-nums; white-space: nowrap; }
td.lead { font-size: 8pt; color: #33414d; }
tr { break-inside: avoid; }
tr.sub td { background: #fafbfc; }
.chip { display: inline-block; font-size: 7.8pt; font-weight: 700; padding: 0 6px; border-radius: 9px; border: 1px solid; line-height: 1.55; white-space: nowrap; }
.g { font-size: 8.4pt; font-weight: 600; padding: 0 3px; border-radius: 3px; white-space: nowrap; }
.g.bridge { box-shadow: inset 0 -2px 0 currentColor; }
.g.out { background: transparent; color: #6b7680; border: 1px dashed #9aa5ae; }
.legend { display: flex; flex-wrap: wrap; gap: 6px 14px; font-size: 8.2pt; color: var(--muted); margin: 4px 0 10px; font-family: 'Helvetica Neue', Arial, sans-serif; }
.turn { margin-top: 10px; }
.turn-h { display: flex; align-items: center; gap: 8px; break-after: avoid; }
.tag { font-size: 8pt; font-weight: 700; color: white; padding: 2px 8px; border-radius: 4px; }
.tag.t1 { background: #1F5FA8; } .tag.t2 { background: #A3245A; }
.meta { font-size: 7.6pt; color: var(--muted); margin: -2px 0 10px; }
.meta code { font-size: 7pt; }
.reveal { border-left: 4px solid #A3245A; background: #fbf3f6; padding: 8px 12px; margin: 4px 0 10px; border-radius: 0 6px 6px 0; }
.reveal p { margin: 2px 0; } .rl { font-size: 7.6pt; text-transform: uppercase; letter-spacing: .05em; color: #A3245A; font-weight: 700; }
.small { font-size: 8.4pt; color: var(--muted); }
.refused { font-style: italic; color: #9B3B2E; }
.howto { border: 1px solid var(--rule); border-radius: 8px; padding: 8px 12px; background: #fcfcfd; font-size: 8.8pt; margin: 8px 0 12px; }
.howto p { margin: 2px 0; }
.page-break { break-before: page; }
.appendix-h { color: var(--muted); }
footer { margin-top: 24px; padding-top: 8px; border-top: 1px solid var(--rule); font-size: 7.6pt; color: var(--muted); font-family: 'Helvetica Neue', Arial, sans-serif; }
"""


def legend(run):
    chips = " ".join(chip(run, c["Community"]) for c in run["communities"])
    return (f'<div class="legend"><span>Gene colours follow the communities: {chips}</span>'
            '<span><span class="g bridge" style="background:#eee;color:#333">GENE</span> underlined = bridge gene</span>'
            '<span><span class="g out">GENE</span> named by the model, not in the network</span></div>')


def build(step09, stem, title, subtitle, patterns, out_dir):
    runs = []
    for pat in patterns:
        d = [x for x in glob.glob(os.path.join(step09, pat)) if os.path.isfile(os.path.join(x, "network_llm.dot"))]
        if d:
            runs.append(load_run(sorted(d)[-1]))
    if not runs:
        return None
    body = [f'<p class="kicker">Supplementary Table {NUM[stem][0]} · LLM interpretation of {"an ORA" if runs[0]["ora"] else "an iGEA"} network</p><h1>{html.escape(title)}</h1>'
            f'<p class="subtitle">{html.escape(subtitle)}</p>']
    second = any(r["followup"] for r in runs)
    turn2 = ('In <b>Turn 2</b> the experiment was revealed and the model rebuilt its interpretation. Each request was '
             'independent and stateless (earlier turns were resent; nothing was stored by the provider).' if second else
             ('No second turn was sent for this baseline. ' if runs[0]["ora"] else 'No second turn was sent for this negative control. ')
             + 'Each request was independent and stateless (nothing was stored by the provider).')
    body.append('<div class="howto"><p><b>What is shown.</b> The network (statistics, communities and the network itself in DOT '
                'format) was sent to a large language model with the prompt in Supplementary Table 3. In <b>Turn 1</b> the '
                f'model did not know how the gene list was generated and was asked to estimate it. {turn2} Repeat 1 is '
                'shown in full; repeats 2 and 3 follow in the appendix.</p></div>')
    for i, run in enumerate(runs):
        label = None
        if len(runs) > 1:
            label = ("A. Default analysis: terms filtered at q < 0.10" if i == 0 else
                     "B. Without the term-level filter (q = 1): all selected terms in the network")
            cls = ' class="page-break"' if i else ""
            body.append(f"<h2{cls}>{html.escape(label)}</h2>")
        body.append(facts(run))
        body.append(stats_table(run))
        body.append("<h2>Communities</h2>" if len(runs) == 1 else "<h3>Communities</h3>")
        body.append(legend(run))
        body.append(community_table(run))
        body.append(answer_block(run, "1"))
    reps = sorted(runs[0]["manifest"], key=int)
    body.append('<h2 class="page-break appendix-h">Appendix · consistency across independent repeats</h2>')
    for i, run in enumerate(runs):
        if len(runs) > 1:
            body.append(f'<h3>{"A. q &lt; 0.10" if i == 0 else "B. q = 1"}</h3>')
        body.append(consistency(run, sorted(run["manifest"], key=int)))
    for rep in reps[1:]:
        for i, run in enumerate(runs):
            suffix = "" if len(runs) == 1 else (" · A. q &lt; 0.10" if i == 0 else " · B. q = 1")
            body.append(f'<h2 class="page-break appendix-h">Appendix · repeat {rep}{suffix}</h2>' + answer_block(run, rep))
    m = runs[0]["manifest"]["1"]
    body.append(f'<footer>iGEA (iterative gene set enrichment analysis). Model {m["returned_model"]} via the OpenAI Responses API, '
                f'store = false. ' + ('Network: ORA terms at BH FDR &lt; 0.05 with all their genes. ' if runs[0]["ora"] else
                f'Network: largest-cluster empirical p from {runs[0]["n_lists"]:,} random gene lists of the same size, '
                'background, collections and parameters. ') + 'Communities: Leiden (modularity) on terms linked by shared genes. '
                f'Raw requests and responses (JSON) are provided with the code. Generated {os.environ.get("GENERATED_DATE", date.today().isoformat())}.</footer>')
    doc = (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>Supplementary Table {NUM[stem][0]} — {html.escape(title)}</title>'
           f'<style>{CSS}</style></head><body>{"".join(body)}</body></html>')
    path = os.path.join(out_dir, NUM[stem][1] + ".html")
    open(path, "w").write(doc)
    return path


if __name__ == "__main__":
    step09, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    only = sys.argv[3:]
    for stem, title, sub, pats in DATASETS:
        if only and stem not in only:
            continue
        p = build(step09, stem, title, sub, pats, out)
        print(stem, "->", p)
