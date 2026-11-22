"""Rebuild llm_prompt.txt in each end-to-end run folder with the current prompt (code/llm_export.py),
keeping the statistics, communities and DOT already computed. The v2 (structured) prompt and its
responses are kept as llm_prompt_v2_structured.txt and llm_v2_structured/."""
import pathlib as _pl  # [publication] repository-relative paths (added when archived)
_REPO = next(p for p in _pl.Path(__file__).resolve().parents if (p / "pyproject.toml").exists())
_PUB = _REPO / "publication" / "csbj_2026"
import glob, os, re, shutil, sys
sys.path.insert(0, (str(_REPO) + "/code"))
import llm_export as L
HERE = os.path.dirname(os.path.abspath(__file__))
for f in sorted(glob.glob(f"{HERE}/**/llm_prompt.txt", recursive=True)):
    if "_to_delete" in f: continue
    d = os.path.dirname(f); old = open(f).read()
    v2 = os.path.join(d, "llm_prompt_v2_structured.txt")
    if not os.path.exists(v2): shutil.copy2(f, v2)
    if os.path.isdir(os.path.join(d, "llm")) and not os.path.exists(os.path.join(d, "llm_v2_structured")):
        os.rename(os.path.join(d, "llm"), os.path.join(d, "llm_v2_structured"))
    base = open(v2).read()
    dot = open(os.path.join(d, "network_llm.dot")).read()
    libs_present = list(dict.fromkeys(re.findall(r'library="([^"]+)"', dot)))
    order = list(L.LIBRARY_DESCRIPTIONS); libs_present.sort(key=lambda l: order.index(l) if l in order else len(order))
    libs = "\n".join(f"- **{l}**: {L.LIBRARY_DESCRIPTIONS.get(l, 'Gene set collection')}" for l in libs_present)
    comm_guide = ("\n- **Communities**: each `subgraph cluster_Cn` block is one community; term and gene nodes carry "
                  "`community`. Bridge genes are listed outside the blocks with `bridge=true` and a `participation` "
                  "value (0 = all terms in one community; higher = spread over communities)." if "subgraph cluster_" in dot else "")
    middle = base[base.index("**STATISTICAL CONTEXT:**"):base.index("**NETWORK INTERPRETATION GUIDE:**")]
    middle = re.sub(r"How to use the statistics:.*", L._stats_section.__code__.co_consts and
                    "How to use the statistics: the network is not distinguishable from what random gene lists produce if the empirical p for the number of nodes in the largest cluster is 0.05 or higher. In that case, say so plainly, keep the rest short and tentative, and do not propose an experimental context. Do not use p-values, q-values or iteration numbers to rank biological importance.", middle)
    new = f"{L.INTRO}\n\n{middle}{L._guide_and_task(libs, comm_guide)}"
    assert dot in new
    open(f, "w").write(new)
    print(os.path.relpath(d, HERE), len(base), "->", len(new))
