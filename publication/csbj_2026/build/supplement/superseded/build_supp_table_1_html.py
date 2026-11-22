import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import html, re

KW = r"\b(OUTPUT|FOR|WHILE|END|IF|BREAK|REPEAT|STOP|RETURN|FUNCTION|APPEND)\b"

def hl(code):
    out, band = [], None
    for line in code.strip("\n").split("\n"):
        if line.startswith("@@"):
            band, label = line[2:].split("|")
            if label:
                out.append(f'<span class="ln {band} bandh">{html.escape(label)}</span>')
            continue
        cm = ""
        if "//" in line:
            i = line.index("//")
            line, cm = line[:i], line[i:]
        s = html.escape(line)
        s = re.sub(r"&quot;([^&]*)&quot;", r'<span class="str">"\1"</span>', s)
        s = re.sub(KW, r'<span class="kw">\1</span>', s)
        s = re.sub(r"(←|∩|∪|&lt;=|&gt;=|&lt;|&gt;|≤|≥)", r'<span class="op">\1</span>', s)
        s = re.sub(r"(?<![\w-])(\d[\d,.]*)(?![\w])", r'<span class="num">\1</span>', s)
        if cm:
            s += f'<span class="cm">{html.escape(cm)}</span>'
        s = s.replace("\t", "    ")
        out.append(f'<span class="ln {band}">{s or " "}</span>' if band else s)
    return "\n".join(out) if not any(o.startswith('<span class="ln') for o in out) else "".join(out)

def block(title, tag, cls, code, note=""):
    n = f'<p class="note">{note}</p>' if note else ""
    return f'''<section class="step {cls}"><div class="steph"><span class="tag">{tag}</span><h3>{title}</h3></div>{n}<pre><code>{hl(code)}</code></pre></section>'''

P1 = """
@@b1|Prepare the collection
G ← G ∩ B                                   // input genes outside the background are removed (one gene universe)
results ← []
FOR Collection in C:
	terms ← {T ∩ B : T in Collection}           // every gene set restricted to the background
	terms ← {T in terms : min_gene_set_size <= |T| <= max_gene_set_size}
	U ← union of all genes in terms              // collection-specific universe
	M ← |U|                                      // collection-specific background size
	remaining_genes ← G ∩ U                      // each collection starts from the full input
	iteration ← 1
@@b2|Select the most significant term
	WHILE (|remaining_genes| > 0 and iteration <= max_iterations):
		N ← |remaining_genes|
		FOR T in terms:
			k ← |T ∩ remaining_genes|              // overlap
			IF k < min_overlap: skip T
			table ← [[k, |T| - k], [N - k, M - |T| - N + k]]
			p(T) ← Fisher's exact test, one-sided (over-representation), on table
			       // = hypergeom.sf(k - 1, M, |T|, N) = P(overlap >= k); under-represented sets get large p
		IF no term was tested: BREAK
		top_term ← term with the lowest p(T)    // ties: order of the gene sets in the collection file
		top_pvalue ← p(top_term)
		IF top_pvalue >= p_threshold: BREAK
@@b3|Record the term and remove its genes
		top_overlap_genes ← top_term ∩ remaining_genes
		record ← {
			Collection: Collection.name,
			Iteration: iteration,
			Term: top_term.name,
			"iteration p-value": top_pvalue,
			"iteration overlapping genes": |top_overlap_genes|,
			"Genes removed for next iteration": top_overlap_genes
		}
		results.APPEND(record)
		remaining_genes ← remaining_genes - top_overlap_genes   // each gene: at most one term per collection
		iteration ← iteration + 1
@@b0|
	END WHILE
END FOR
OUTPUT: results (selected terms of all collections)
"""
P2_RAND = """
@@b1|Random lists of the input size (network test and first q-value round)
n ← |G|                                             // HIV: 265
FOR i in 1..n_lists:                                // default 10,000
	R_i ← n genes drawn from B without replacement      // generator seeded with (seed, n, i)
	random_results_i ← Part 1 applied to R_i            // all collections, all iterations, same parameters
@@b3|Reduced-size random lists (later q-value rounds only)
reduced_sizes ← n - s, n - 2s, ..., not below n/4      // s = n/20, rounded up
                                                    // HIV: 251, 237, ..., 97, ..., 69 (14 sizes)
FOR size in reduced_sizes:
	FOR j in 1..max(200, n_lists/10):                   // 10,000 lists: 1,000 per size
		R ← size genes drawn from B without replacement  // seeded with (seed, size, j)
		                                                  // a fresh draw, not a subset of the lists above
		Part 1 applied to R                               // all collections, all iterations, same parameters
@@b0|
// All random lists are drawn and analysed here, once, before any q-value is computed.
"""
P2_Q = """
FUNCTION q_values(selected_terms, list_size):
	// The rounds below are not iterations of Part 1: they start after the terms of
	// all collections have been selected, and they change only the q-values.
	signal_genes ← {}                                     // round 1: every gene is treated as noise
	REPEAT (at most 10 rounds):
		noise_size ← list_size - |signal_genes|            // genes not explained by accepted terms
		ref_size ← smallest of n and reduced_sizes that is >= noise_size
		           // rounded UP: can only over-estimate chance terms
		FOR each collection c, FOR each term t selected in c with p-value p:
			expected ← mean number of terms with p-value <= p selected in c from the random lists of ref_size
			           // when a random list of size n is scored, that list itself is left out
			observed ← number of terms with p-value <= p selected in c from the scored list
			q_raw(t) ← min(1, expected / observed)
		q ← q_raw made monotone in p within each collection   // as Benjamini-Hochberg
		new_signal ← all genes assigned to terms with q < accept_level, pooled across collections
		             // accept_level = q_threshold (0.10 when q_threshold = 1)
		IF new_signal = signal_genes: STOP                 // the accepted terms no longer change
		signal_genes ← new_signal
	RETURN q

q ← q_values(results, n)                               // the input list
kept ← terms in results with q < q_threshold          // q_threshold = 1: all terms kept
// Every random list of size n is scored by the same function in step 2d.
"""
P2_NET = """
network ← gene-term graph: an edge between each gene and the term it was assigned to (no edge weights)
cluster ← largest connected component of network
metrics ← cluster size (genes + terms), genes, terms, edges, edges per gene
          // number of collections in the cluster: descriptive
"""
P2_TEST = """
FOR i in 1..n_lists:
	q_i ← q_values(random_results_i, n)
	kept_i ← terms with q_i < q_threshold                // q_threshold = 1: all terms kept
	metrics_i ← metrics of the largest connected component built from kept_i
FOR each metric:
	r ← number of random lists with metrics_i >= metrics
	empirical p ← (r + 1) / (n_lists + 1)
	z-score ← (metrics - mean(metrics_i)) / SD(metrics_i)      // descriptive only
"""
P2_COMM = """
communities ← Leiden (modularity) on the term-term projection of network
              // terms linked by shared genes, weighted by their number; best of 100 seeds
              // networks with fewer than 6 terms: connected components
OUTPUT: q-value of every selected term, kept terms, network with communities, empirical p-values
"""

PARAMS1 = [("B", "Background gene set", "23,391 validated human gene symbols"),
           ("G", "Input genes", "—"),
           ("C", "Gene-set collections", "10 MSigDB collections (v2025.1)"),
           ("p_threshold", "p-value cutoff for selecting a term; also the stopping criterion", "0.01"),
           ("min_overlap", "Minimum overlap between the remaining input genes and a gene set", "3"),
           ("min/max_gene_set_size", "Gene-set size range, counted after restriction to B", "10–600"),
           ("max_iterations", "Maximum iterations per collection", "30"),
           ("test", "Fisher's exact test, one-sided (over-representation); the hypergeometric test gives the same p-values", "Fisher, one-sided")]
PARAMS2 = [("n_lists", "Number of random gene lists", "10,000 (web app: 1,000)"),
           ("q_threshold", "Term-level q-value threshold; 1 keeps all terms (q-values still reported)", "0.10"),
           ("seed", "Random seed", "2026")]

def ptable(rows):
    r = "".join(f"<tr><td><code>{html.escape(a)}</code></td><td>{html.escape(b)}</td><td class='def'>{html.escape(c)}</td></tr>" for a, b, c in rows)
    return f"<table><thead><tr><th>Parameter</th><th>Meaning</th><th>Default</th></tr></thead><tbody>{r}</tbody></table>"

CSS = """
:root{--ink:#22252b;--muted:#5d6470;--line:#e3e6ea;--bg:#fbfbfa;
--a:#2a78d6;--a-bg:#eef5fd;--b:#1baf7a;--b-bg:#ecf8f3;--c:#eb6834;--c-bg:#fdf1ec;--d:#6250d6;--d-bg:#f1effc}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif}
main{max-width:1000px;margin:0 auto;padding:40px 28px 60px}
.kicker{font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0}
h1{font-size:26px;margin:4px 0 12px}
h2{font-size:19px;margin:38px 0 10px;padding-bottom:6px;border-bottom:2px solid var(--line)}
h3{font-size:15px;margin:0}
p{margin:8px 0}
.lead{color:#3a3f47;max-width:860px}
.flow{display:flex;gap:10px;align-items:stretch;margin:22px 0 6px;flex-wrap:wrap}
.box{flex:1 1 180px;border-radius:10px;padding:12px 14px;border:1px solid var(--line);background:#fff}
.box b{display:block;font-size:13.5px;margin-bottom:3px}
.box span{font-size:12.5px;color:var(--muted)}
.box.a{border-top:4px solid var(--a)} .box.b{border-top:4px solid var(--b)}
.box.c{border-top:4px solid var(--c)} .box.d{border-top:4px solid var(--d)}
.arrow{align-self:center;color:#9aa1ab;font-size:20px}
table{border-collapse:collapse;width:100%;background:#fff;border:1px solid var(--line);border-radius:8px;overflow:hidden;font-size:14px;margin:10px 0}
th{background:#f3f4f6;text-align:left;font-weight:600;padding:8px 12px;border-bottom:1px solid var(--line)}
td{padding:7px 12px;border-bottom:1px solid var(--line);vertical-align:top}
tr:last-child td{border-bottom:none}
td code{font-size:13px;color:#7a3e9d}
td.def{white-space:nowrap;font-weight:600}
.step{border-radius:10px;margin:14px 0;border:1px solid var(--line);overflow:hidden;background:#fff}
.steph{display:flex;align-items:center;gap:10px;padding:10px 14px}
.tag{font-size:11.5px;font-weight:700;color:#fff;border-radius:999px;padding:2px 9px;white-space:nowrap}
.a .steph{background:var(--a-bg)} .a .tag{background:var(--a)} .a{border-left:5px solid var(--a)}
.b .steph{background:var(--b-bg)} .b .tag{background:var(--b)} .b{border-left:5px solid var(--b)}
.c .steph{background:var(--c-bg)} .c .tag{background:var(--c)} .c{border-left:5px solid var(--c)}
.d .steph{background:var(--d-bg)} .d .tag{background:var(--d)} .d{border-left:5px solid var(--d)}
.note{margin:10px 14px 0;color:#3a3f47;font-size:14px}
pre{margin:0;padding:12px 16px 14px;overflow-x:auto;background:#fff}
code{font-family:"SF Mono",Menlo,Consolas,"Liberation Mono",monospace}
pre code{font-size:13px;line-height:1.6;color:#2b2f36}
.kw{color:#b4237a;font-weight:700} .cm{color:#6b7f6e;font-style:italic}
.str{color:#0b7a75} .op{color:#c25e00;font-weight:600} .num{color:#2a62b8}
.legend{display:flex;gap:16px;flex-wrap:wrap;font-size:12.5px;color:var(--muted);margin:6px 0 0}
.legend span{display:inline-flex;align-items:center;gap:6px}
.sw{font-family:"SF Mono",Menlo,monospace;font-size:12px;padding:1px 6px;border-radius:4px;background:#fff;border:1px solid var(--line)}
pre code .ln{display:block;padding:0 10px;margin:0 -10px}
.ln.b0,.ln.b2{background:#fff} .a .ln.b1,.a .ln.b3{background:#f4f8fd} .b .ln.b1,.b .ln.b3{background:#f0f9f5}
.ln.bandh{font:600 11.5px/2 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;color:var(--a);
letter-spacing:.04em;text-transform:uppercase;border-top:1px dashed #c9dbf3;margin-top:4px}
.ln.b1.bandh{border-top:none;margin-top:0} .b .ln.bandh{color:#13875e;border-top-color:#bfe6d6}
footer{margin-top:36px;padding-top:14px;border-top:1px solid var(--line);color:var(--muted);font-size:12.5px}
@media print{body{background:#fff;font-size:11.5pt} main{padding:0;max-width:none}
.step,.flow,table{break-inside:avoid} pre code{font-size:9pt} h2{break-after:avoid}}
"""

page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Supplementary Table 1 — iGEA algorithm</title><style>{CSS}</style></head>
<body><main>
<p class="kicker">Supplementary Table 1</p>
<h1>Iterative Gene Enrichment Analysis (iGEA): algorithm</h1>
<p class="lead">iGEA repeatedly identifies the most significant enriched term, removes its overlapping genes from the input list and repeats the enrichment analysis on the remaining genes. This is done separately for each gene-set collection, until termination criteria are met (Part 1). Random gene lists processed by the same procedure then give each selected term an empirical q-value (term level). Terms with q below the threshold form the cross-collection network, whose largest connected cluster is compared with those of the random lists (network level) (Part 2).</p>

<div class="flow">
 <div class="box a"><b>1 · Iterative selection</b><span>per collection; p &lt; 0.01; genes of each selected term removed</span></div>
 <div class="arrow">→</div>
 <div class="box b"><b>2 · Term level</b><span>empirical q-value from random lists; keep q &lt; 0.10</span></div>
 <div class="arrow">→</div>
 <div class="box c"><b>3 · Network level</b><span>largest cluster vs random lists; empirical p</span></div>
 <div class="arrow">→</div>
 <div class="box d"><b>4 · Communities</b><span>Leiden; hub and bridge genes (descriptive)</span></div>
</div>

<h2>Parameters</h2>
{ptable(PARAMS1 + PARAMS2)}

<h2>Part 1 · Iterative selection</h2>
<div class="legend"><span><span class="sw"><span class="kw">FOR</span></span> keyword</span><span><span class="sw"><span class="op">←</span></span> assignment / set operation</span><span><span class="sw"><span class="num">0.01</span></span> number</span><span><span class="sw"><span class="cm">// note</span></span> comment</span></div>
{block("Iterative selection, separately for each collection", "1", "a", P1, "Every gene set is restricted to the background before the size filter, so the input list, the gene sets and the random lists come from one gene universe. Worked example of the test: N = 350, M = 4,500, |T| = 50 and k = 15 give the table [[15, 35], [335, 4,115]].")}

<h2>Part 2 · Two-level empirical significance testing</h2>
{block("Random gene lists, analysed exactly like the input", "2a", "b", P2_RAND, "Two kinds of random lists are drawn, both independently from the background and both analysed exactly like the input. Lists of the input size are the reference for the network test and for the first q-value round. Reduced-size lists are needed because a q-value compares the input with chance, and chance selections come only from the input genes that are noise. Genes already explained by real terms cannot produce chance terms, so once terms are accepted the reference must be random lists the size of the unexplained part of the input (step 2b).")}
{block("Term level: empirical q-value within each collection", "2b", "b", P2_Q, "HIV example. Round 1 compares the HIV terms with the 265-gene lists, which select about 16.0 terms per list by chance across the ten collections. The rounds settle when the accepted terms explain 168 genes; the other 97 genes are matched by the 97-gene lists, which select about 7.2 terms by chance. For Reactome HIV transcription elongation (p = 0.0085; 9 HIV Reactome terms have p ≤ 0.0085), the 265-gene lists give 2.07 chance Reactome terms per list with p ≤ 0.0085, so q = 2.07/9 = 0.23 and the term would be removed; the 97-gene lists give 0.87, so q = 0.87/9 = 0.097 and the term is kept.")}
{block("Network from the kept terms", "2c", "c", P2_NET)}
{block("Network level: the same filter and network for every random list", "2d", "c", P2_TEST, "The metrics describe the same cluster and are not independent tests. With 10,000 random lists the smallest attainable p-value is 0.0001.")}
{block("Communities, hub genes and bridge genes", "2e", "d", P2_COMM)}

<footer>Implementation: <code>code/fast_igea.py</code> (Part 1), <code>code/per_run_null.py</code> (Part 2), <code>code/network_modules.py</code> (communities); iGEA source code (GPL-3.0).</footer>
</main></body></html>"""
open("out/Supp Table 01 - iGEA pseudocode.html", "w", encoding="utf8").write(page)
print("ok", len(page))
