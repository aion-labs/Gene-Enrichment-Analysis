"""
Gene Set Enrichment Utilities
-----------------------------
Load MSigDB gene set collections (GMT format) and run:
  - ORA  (Over-Representation Analysis, Fisher's exact test)
  - GSEA (Pre-ranked Gene Set Enrichment Analysis)

Collections included:
  Hallmark (50), KEGG Medicus (658), Reactome (1736), GO BP (7608)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

GENE_SETS_DIR = '/root/capsule/data/gene_sets'

COLLECTIONS = {
    'Hallmark': 'h.all.v2024.1.Hs.symbols.gmt',
    'KEGG': 'c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt',
    'Reactome': 'c2.cp.reactome.v2024.1.Hs.symbols.gmt',
    'GO_BP': 'c5.go.bp.v2024.1.Hs.symbols.gmt',
}


# ---------------------------------------------------------------------------
# GMT loading
# ---------------------------------------------------------------------------

def load_gmt(filepath):
    """Parse a GMT file into {set_name: [gene_symbols]}."""
    gene_sets = {}
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            name = parts[0]
            genes = parts[2:]          # parts[1] is description/URL
            gene_sets[name] = genes
    return gene_sets


def load_collections(gene_sets_dir=None):
    """
    Load all available collections.

    Returns {collection_name: {set_name: [genes]}}.
    """
    if gene_sets_dir is None:
        gene_sets_dir = GENE_SETS_DIR
    collections = {}
    for coll_name, filename in COLLECTIONS.items():
        path = os.path.join(gene_sets_dir, filename)
        if os.path.exists(path):
            collections[coll_name] = load_gmt(path)
    return collections


def filter_gene_sets(gene_sets, tested_genes, min_size=15, max_size=500):
    """
    Keep gene sets whose overlap with *tested_genes* falls within
    [min_size, max_size].  Returns {name: [overlapping genes]}.
    """
    tested = set(tested_genes)
    filtered = {}
    for name, genes in gene_sets.items():
        overlap = [g for g in genes if g in tested]
        if min_size <= len(overlap) <= max_size:
            filtered[name] = overlap
    return filtered


# ---------------------------------------------------------------------------
# ORA – Over-Representation Analysis (Fisher's exact test)
# ---------------------------------------------------------------------------

def run_ora(up_genes, down_genes, all_tested, gene_sets, collection_name=''):
    """
    Fisher's exact test for over-representation of DE gene lists in gene
    sets.  Tests both up-regulated and down-regulated lists.

    Returns a DataFrame with columns:
        collection, direction, term, overlap, term_size, query_size,
        p_value, fdr, genes
    """
    all_tested_set = set(all_tested)
    results = []

    for direction, query_genes, label in [
        ('up', up_genes, 'Upregulated in R'),
        ('down', down_genes, 'Downregulated in R'),
    ]:
        query_set = set(query_genes)
        for gs_name, gs_genes in gene_sets.items():
            gs_in_tested = set(gs_genes) & all_tested_set
            if len(gs_in_tested) < 3:
                continue
            overlap = query_set & gs_in_tested
            a = len(overlap)
            if a == 0:
                continue
            b = len(query_set) - a
            c = len(gs_in_tested) - a
            d = len(all_tested_set) - a - b - c

            _, pval = stats.fisher_exact([[a, b], [c, d]],
                                         alternative='greater')
            results.append({
                'collection': collection_name,
                'direction': label,
                'term': gs_name,
                'overlap': a,
                'term_size': len(gs_in_tested),
                'query_size': len(query_set),
                'p_value': pval,
                'genes': ', '.join(sorted(overlap)),
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, fdr, _, _ = multipletests(df['p_value'], method='fdr_bh')
        df['fdr'] = fdr
        df = df.sort_values('p_value')
    return df


# ---------------------------------------------------------------------------
# GSEA – Pre-ranked Gene Set Enrichment Analysis
# ---------------------------------------------------------------------------

def _compute_es(metrics_p, hit_indicator):
    """Enrichment score via weighted running-sum statistic."""
    hit_weights = metrics_p * hit_indicator
    hit_norm = hit_weights.sum()
    if hit_norm == 0:
        return 0.0

    miss = 1.0 - hit_indicator
    miss_norm = miss.sum()
    if miss_norm == 0:
        return 0.0

    running = np.cumsum(hit_weights / hit_norm) - np.cumsum(miss / miss_norm)
    max_pos = running.max()
    min_neg = running.min()
    return float(max_pos) if abs(max_pos) >= abs(min_neg) else float(min_neg)


def _perm_null_vectorized(metrics_p, N, set_size, nperm, rng):
    """
    Null distribution of ES for *nperm* random gene sets of *set_size*.
    Fully vectorised with numpy for speed.
    """
    # Random hit indicators: (nperm, N)
    perm_hits = np.zeros((nperm, N), dtype=np.float64)
    for i in range(nperm):
        idx = rng.choice(N, size=set_size, replace=False)
        perm_hits[i, idx] = 1.0

    # Hit running sum
    hw = metrics_p[np.newaxis, :] * perm_hits
    hw_norm = np.maximum(hw.sum(axis=1, keepdims=True), 1e-10)
    run_hit = np.cumsum(hw / hw_norm, axis=1)

    # Miss running sum
    ms = 1.0 - perm_hits
    ms_norm = np.maximum(ms.sum(axis=1, keepdims=True), 1e-10)
    run_miss = np.cumsum(ms / ms_norm, axis=1)

    running = run_hit - run_miss
    max_pos = running.max(axis=1)
    min_neg = running.min(axis=1)
    return np.where(np.abs(max_pos) >= np.abs(min_neg), max_pos, min_neg)


def run_gsea_prerank(de_df, gene_sets, collection_name='',
                     nperm=1000, min_size=15, max_size=500,
                     weighted_score_type=1, seed=42):
    """
    Pre-ranked GSEA.

    Ranking metric: sign(log2FC) * -log10(p_value).

    Returns DataFrame with columns:
        collection, term, es, nes, size, p_value, fdr, matched_genes
    """
    rng = np.random.default_rng(seed)

    # --- build ranking ---
    ranking = {}
    for _, row in de_df.iterrows():
        pval = max(row['p_value'], 1e-300)
        ranking[row['gene']] = np.sign(row['log2FC']) * (-np.log10(pval))

    sorted_genes = sorted(ranking, key=ranking.get, reverse=True)
    sorted_metrics = np.array([ranking[g] for g in sorted_genes])
    gene_to_idx = {g: i for i, g in enumerate(sorted_genes)}
    N = len(sorted_genes)
    metrics_p = np.abs(sorted_metrics) ** weighted_score_type

    # --- filter gene sets ---
    tested = set(ranking)
    filtered = {}
    for name, genes in gene_sets.items():
        overlap = [g for g in genes if g in tested]
        if min_size <= len(overlap) <= max_size:
            filtered[name] = overlap

    if not filtered:
        return pd.DataFrame()

    # --- group by set size (cache permutation nulls per size) ---
    size_groups = {}
    for name, genes in filtered.items():
        sz = len(genes)
        size_groups.setdefault(sz, []).append((name, genes))

    results = []
    null_cache = {}

    for sz in sorted(size_groups):
        if sz not in null_cache:
            null_cache[sz] = _perm_null_vectorized(metrics_p, N, sz,
                                                   nperm, rng)
        perm_null = null_cache[sz]

        # Mean of positive / negative null for NES
        pos_null = perm_null[perm_null >= 0]
        neg_null = perm_null[perm_null < 0]
        mean_pos = pos_null.mean() if len(pos_null) > 0 else 1.0
        mean_neg = np.abs(neg_null.mean()) if len(neg_null) > 0 else 1.0

        for name, genes in size_groups[sz]:
            hit = np.zeros(N)
            for g in genes:
                hit[gene_to_idx[g]] = 1.0

            es = _compute_es(metrics_p, hit)

            # NES
            if es >= 0:
                nes = es / mean_pos if mean_pos > 0 else 0.0
            else:
                nes = es / mean_neg if mean_neg > 0 else 0.0

            # Empirical p-value
            if es >= 0:
                pval = np.mean(perm_null >= es)
            else:
                pval = np.mean(perm_null <= es)
            pval = max(pval, 1.0 / (nperm + 1))

            results.append({
                'collection': collection_name,
                'term': name,
                'es': round(es, 4),
                'nes': round(nes, 4),
                'size': sz,
                'p_value': pval,
                'matched_genes': ', '.join(sorted(genes)),
            })

    df = pd.DataFrame(results)
    if len(df) > 0:
        _, fdr, _, _ = multipletests(df['p_value'], method='fdr_bh')
        df['fdr'] = fdr
        df = df.sort_values('p_value')
    return df
