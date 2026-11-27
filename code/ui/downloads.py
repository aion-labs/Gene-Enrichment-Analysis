"""
Downloads of the iGEA results, with the statistics of this run.

Every iGEA table that can be downloaded carries the term q-value, whether the term is in the network
(q below the threshold) and its community, exactly as shown in the iteration tables. Links are labeled
with what they contain rather than with the file extension.
"""
import base64
import io
import tarfile
from typing import Dict, List, Optional, Set

import pandas as pd

from ui.utils import download_link


def bytes_link(data: bytes, filename: str, label: str) -> str:
    """Download link for binary content (e.g. a tar.gz archive) built in memory."""
    b64 = base64.b64encode(data).decode()
    return f'<a href="data:application/octet-stream;base64,{b64}" download="{filename}">{label}</a>'


def _safe(name: str) -> str:
    """File-name-safe version of a collection or list name (e.g. 'C5: Gene Ontology: X' -> 'C5_Gene_Ontology_X')."""
    import re
    return re.sub(r"[^A-Za-z0-9-]+", "_", name).strip("_")


def _full_stats(state):
    """Statistics computed for all collections (used for the q-values in the iteration tables)."""
    return state.get("per_run_full") or None


def kept_terms(state) -> Optional[Set[str]]:
    """Raw names of the terms that passed the q-value filter, or None if statistics are not available."""
    stored = _full_stats(state)
    if not stored or "terms" not in stored:
        return None
    t = stored["terms"]
    return set(t.loc[t["Kept"].astype(bool), "Term"])


def q_threshold(state) -> Optional[float]:
    stored = _full_stats(state)
    return stored["attrs"].get("q_threshold") if stored else None


def igea_table(state, library_name: str, it) -> pd.DataFrame:
    """One collection's iGEA table as shown in the app: iteration p-value, q-value, In network and Community."""
    from ui.per_run_benchmark import term_stats
    from ui.rendering import add_term_stats

    df = it.to_dataframe()
    if df.empty or "Iteration" not in df.columns:
        return df
    stats, _ = term_stats(state, library_name)
    df = add_term_stats(df.set_index("Iteration"), stats).reset_index()
    if "Kept" in df.columns:
        df = df.rename(columns={"Kept": "In network"})
        stored = _full_stats(state)
        mods = stored.get("modules") if stored else None
        comm = ({(r.Library, r.Term): r.Community for r in mods.terms.itertuples()}
                if mods is not None and mods.summary.get("n_terms") else {})
        pos = list(df.columns).index("In network") + 1
        df.insert(pos, "Community", [comm.get((library_name, t), "") if k else ""
                                     for t, k in zip(df["Term"], df["In network"])])
    if "Library" in df.columns:
        df = df.rename(columns={"Library": "Collection"})
        df = df[["Collection"] + [c for c in df.columns if c != "Collection"]]
    return df


def igea_combined(state, order: List[str]) -> pd.DataFrame:
    frames = [igea_table(state, lib, state.iter_enrich[lib]) for lib in order if lib in state.iter_enrich]
    frames = [f for f in frames if not f.empty]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _archive(files: Dict[str, str]) -> bytes:
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        for name, text in files.items():
            data = text.encode("utf-8")
            info = tarfile.TarInfo(name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
    return buf.getvalue()


def render_results_downloads(st, state, ora_results: Dict[str, List[dict]], ora_tsv: Optional[str], order: List[str]) -> None:
    """Download table for the ORA and iGEA results of all collections."""
    name = _safe(state.get("gene_set_name") or "gene_list")
    qt = q_threshold(state)
    igea = igea_combined(state, order)
    igea_link = download_link(igea.to_csv(sep="\t", index=False), f"{name}_iGEA_results", "tsv",
                              label="iGEA results (TSV)") if not igea.empty else "—"
    igea_files = {f"{_safe(lib)}_iGEA.tsv": igea_table(state, lib, state.iter_enrich[lib]).to_csv(sep="\t", index=False)
                  for lib in order if lib in state.iter_enrich and not state.iter_enrich[lib].to_dataframe().empty}
    igea_arch = bytes_link(_archive(igea_files), f"{name}_iGEA_per_collection.tar.gz",
                           "one file per collection (tar.gz)") if igea_files else "—"
    ora_link = download_link(ora_tsv, f"{name}_ORA_results", "tsv", label="ORA results (TSV)") if ora_tsv else "—"
    ora_files = {}
    for lib, recs in (ora_results or {}).items():
        if recs:
            ora_files[f"{_safe(lib)}_ORA.tsv"] = pd.DataFrame(recs).rename(columns={"Library": "Collection"}).to_csv(sep="\t", index=False)
    ora_arch = bytes_link(_archive(ora_files), f"{name}_ORA_per_collection.tar.gz",
                          "one file per collection (tar.gz)") if ora_files else "—"
    q_note = (f"q-value, whether the term is in the network (q < {qt:g}) and its community"
              if qt is not None and qt < 1 else "q-value and community (no q-value filter)") if kept_terms(state) is not None \
        else "statistics not computed: no q-values"
    genes = sorted(state.gene_set.genes) if state.get("gene_set") is not None and state.gene_set.genes else []
    genes_link = download_link("\n".join(genes), f"{name}_validated_genes", "txt", label="Validated genes (TXT)") if genes else "—"
    st.markdown("#### 📥 Download results")
    st.markdown(
        "| Results | All collections | Per collection | Contents |\n"
        "|---|---|---|---|\n"
        f"| **iGEA** (all iterations) | {igea_link} | {igea_arch} | Selected terms with iteration p-value, {q_note}, "
        "and the genes assigned to each term |\n"
        f"| **ORA** (iteration 1, every gene set) | {ora_link} | {ora_arch} | Every tested gene set with p-value and "
        "Benjamini–Hochberg FDR within its collection |\n"
        f"| **Input** | {genes_link} | | Validated gene symbols used in the analysis |",
        unsafe_allow_html=True,
    )
    st.caption("Network statistics, communities and hub genes are downloaded in the Network statistics section below; "
               "the network and the AI prompt in the network section.")
