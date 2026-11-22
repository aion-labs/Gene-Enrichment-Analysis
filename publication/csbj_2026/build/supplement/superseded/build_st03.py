"""Supplementary Table 3 (HTML): the standardized LLM prompt (v5) and the second-turn prompt."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import html, re
CSS = open(str(_PUB) + "/build/supplement/st16.css").read() + """
pre.prompt{white-space:pre-wrap;word-break:break-word;background:#fff;border:1px solid var(--line);border-radius:8px;padding:14px 16px;font:12.5px/1.5 ui-monospace,Menlo,Consolas,monospace}
pre.dot{max-height:none}.tag{display:inline-block;font-size:11px;padding:1px 7px;border-radius:9px;background:#eef2f7;color:#40506a;margin-left:6px}
"""
P = open(str(_PUB) + "/analyses/step09_end_to_end/hiv_go_only/HIV_dependency_factors.symbols_20261120_114431/llm_prompt.txt").read()
src = open(str(_REPO) + "/code/llm_export.py").read()
FOLLOW = re.search(r'FOLLOWUP = """(.*?)"""', src, re.S).group(1)
# split the example prompt: before the DOT, DOT, after
i = P.index("**NETWORK (DOT):**"); j = P.index("\n}\n", i) + 3
pre, dot, post = P[:i], P[i:j], P[j:]
CTX = [("HIV dependency factors (Brass, König and Zhou screens)", "an siRNA screen in human cells for host factors required for HIV infection (HIV dependency factors)"),
       ("Hypoxia signature", "a meta analysis of genes that are consistently regulated by mRNA expression when cells are exposed to hypoxia"),
       ("Influenza vaccination signature", "genes whose expression in blood 3 days after seasonal influenza vaccination correlated positively with the antibody response (hemagglutination-inhibition titer) at day 28, from a human systems-vaccinology study")]
e = lambda s: html.escape(s, quote=False)
doc = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Supplementary Table 3 — LLM prompt</title><style>{CSS}</style></head><body><main>
<p class="kicker">Supplementary Table 3</p><h1>Standardized prompt for LLM interpretation of an iGEA network</h1>
<p class="lead">The iGEA application generates the prompt from the analysis results (file <code>llm_prompt.txt</code> in every run folder; function
<code>llm_prompt</code> in <code>code/llm_export.py</code>). The same template was used for every network in this study (HIV, hypoxia,
vaccination, random control lists and the three HIV screens); only the statistics, the community table, the network and the list of collections
change with the data. For the ORA networks, only the sentences that describe the method were changed (ORA instead of iGEA; FDR instead of the
empirical q-value). Model: gpt-6.2-sol (OpenAI Responses API, default sampling settings, no stored conversation); each first prompt was sent three
times as independent requests.</p>
<div class="cards">
<div class="card a"><b>A. First prompt</b><span>Statistics, communities, network (DOT), interpretation guide, task. No information on the dataset or experiment.</span></div>
<div class="card b"><b>B. Second turn</b><span>Sent after each first answer; reveals how the gene list was generated and asks for a revised interpretation.</span></div>
<div class="card c"><b>C. Context sentences</b><span>The one sentence inserted into the second turn for each dataset.</span></div></div>
<h2>A. First prompt (example: HIV dependency factors, three Gene Ontology collections, the network of Figure 2)</h2>
<p class="note">Complete prompt as sent. Data-dependent parts: the statistics table, the community table, the network and the list of collections.
The network block is shown in full below the prompt.</p>
<pre class="prompt">{e(pre)}<span class="muted">[network in DOT format: see below]</span>
{e(post)}</pre>
<h3>Network block of this prompt</h3>
<pre class="prompt dot">{e(dot)}</pre>
<h2>B. Second turn</h2>
<p class="note">Sent in the same conversation as each first answer (not for the random control lists). <code>{{context}}</code> is replaced by the sentence in C.</p>
<pre class="prompt">{e(FOLLOW.strip())}</pre>
<h2>C. Context sentences</h2>
<div class="tbl"><table><tr><th>Dataset</th><th>Sentence inserted into the second turn</th></tr>
{''.join(f'<tr><td>{e(a)}</td><td>{e(b)}</td></tr>' for a, b in CTX)}</table></div>
</main></body></html>"""
open("out/Supp Table 03 - LLM prompt.html", "w").write(doc)
txt = ("Supplementary Table 3. Standardized prompt for LLM interpretation of an iGEA network.\n\n"
       "A. FIRST PROMPT (example: HIV dependency factors, three Gene Ontology collections)\n\n" + P +
       "\n\nB. SECOND TURN ({context} = sentence in C)\n\n" + FOLLOW.strip() +
       "\n\nC. CONTEXT SENTENCES\n" + "\n".join(f"- {a}: {b}" for a, b in CTX) + "\n")
open("out/Supp Table 03 - LLM prompt.txt", "w").write(txt)
