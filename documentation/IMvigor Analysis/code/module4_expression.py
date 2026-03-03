"""
Module 4: Differential Expression & Pathway Analysis
IMvigor210 Biomarker Analysis
"""

import os
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests
from gene_set_utils import load_collections, filter_gene_sets, run_ora, run_gsea_prerank

RESULTS_DIR = '/results'



def filter_genes_for_de(log_expr, cohorts, min_mean_expr=1.0,
                        min_variance_quantile=0.25):
    """
    Pre-filter genes before differential expression testing to reduce the
    multiple-testing burden (FDR penalty).

    Both filters are independent of the test statistic so FDR control is
    preserved.

    Filters
    -------
    1. Minimum expression: drop genes whose mean log2 expression across the
       analysis cohort falls below *min_mean_expr*.
    2. Minimum variance: drop genes in the bottom *min_variance_quantile*
       of expression variance across cohort samples (genes that barely
       change across patients).

    Returns
    -------
    filtered_expr : DataFrame
        Expression matrix (genes x all patients) restricted to passing genes.
    filter_stats : dict
        Counts at each filtering step for logging.
    """
    # Restrict to patients in the binary analysis cohort
    cohort_ids = [p for p in cohorts['responders'] + cohorts['non_responders']
                  if p in log_expr.columns]
    expr_cohort = log_expr[cohort_ids]

    n_initial = expr_cohort.shape[0]

    # --- Filter 1: minimum mean expression ---
    gene_means = expr_cohort.mean(axis=1)
    pass_expr = gene_means >= min_mean_expr
    expr_cohort = expr_cohort[pass_expr]
    n_after_expr = expr_cohort.shape[0]

    # --- Filter 2: remove bottom variance quantile ---
    gene_vars = expr_cohort.var(axis=1)
    var_threshold = gene_vars.quantile(min_variance_quantile)
    pass_var = gene_vars >= var_threshold
    expr_cohort = expr_cohort[pass_var]
    n_final = expr_cohort.shape[0]

    # Return the full patient matrix (not just cohort columns) for the
    # genes that passed, so downstream code still has all patients.
    filtered_expr = log_expr.loc[expr_cohort.index]

    filter_stats = {
        'n_initial': n_initial,
        'n_after_expression_filter': n_after_expr,
        'n_after_variance_filter': n_final,
        'min_mean_expr_threshold': min_mean_expr,
        'variance_quantile_threshold': min_variance_quantile,
        'variance_cutoff_value': round(float(var_threshold), 4),
    }

    return filtered_expr, filter_stats


def differential_expression(log_expr, cohorts):
    """
    Mann-Whitney U test per gene: responders vs non-responders.
    Returns full DE table with log2FC, p-value, FDR.
    """
    responders = cohorts['responders']
    non_responders = cohorts['non_responders']

    # Ensure patients are in expression matrix
    r_ids = [p for p in responders if p in log_expr.columns]
    nr_ids = [p for p in non_responders if p in log_expr.columns]

    r_expr = log_expr[r_ids]
    nr_expr = log_expr[nr_ids]

    results = []
    for gene in log_expr.index:
        r_vals = r_expr.loc[gene].values
        nr_vals = nr_expr.loc[gene].values

        # Skip genes with no variance
        if np.std(r_vals) == 0 and np.std(nr_vals) == 0:
            continue

        log2fc = np.mean(r_vals) - np.mean(nr_vals)
        stat, p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        results.append({
            'gene': gene,
            'log2FC': round(log2fc, 4),
            'mean_responder': round(np.mean(r_vals), 4),
            'mean_nonresponder': round(np.mean(nr_vals), 4),
            'statistic': stat,
            'p_value': p,
        })

    de_df = pd.DataFrame(results)
    _, fdr, _, _ = multipletests(de_df['p_value'], method='fdr_bh')
    de_df['fdr'] = fdr
    de_df = de_df.sort_values('p_value')

    return de_df


def volcano_plot(de_df):
    """Volcano plot of differential expression results."""
    df = de_df.copy()
    df['neg_log10_fdr'] = -np.log10(df['fdr'].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(10, 8))

    # Classify: significant up, down, or NS
    sig_up = (df['fdr'] < 0.05) & (df['log2FC'] > 1)
    sig_down = (df['fdr'] < 0.05) & (df['log2FC'] < -1)
    ns = ~(sig_up | sig_down)

    ax.scatter(df.loc[ns, 'log2FC'], df.loc[ns, 'neg_log10_fdr'],
               c='grey', alpha=0.3, s=10, label=f'NS ({ns.sum()})')
    ax.scatter(df.loc[sig_up, 'log2FC'], df.loc[sig_up, 'neg_log10_fdr'],
               c='#e74c3c', alpha=0.6, s=15, label=f'Up in R ({sig_up.sum()})')
    ax.scatter(df.loc[sig_down, 'log2FC'], df.loc[sig_down, 'neg_log10_fdr'],
               c='#3498db', alpha=0.6, s=15, label=f'Down in R ({sig_down.sum()})')

    ax.axhline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5)
    ax.axvline(1, color='grey', linestyle='--', alpha=0.5)
    ax.axvline(-1, color='grey', linestyle='--', alpha=0.5)

    # Label top genes
    top_genes = pd.concat([
        df[sig_up].nlargest(5, 'neg_log10_fdr'),
        df[sig_down].nlargest(5, 'neg_log10_fdr'),
    ])
    for _, row in top_genes.iterrows():
        ax.annotate(row['gene'], (row['log2FC'], row['neg_log10_fdr']),
                     fontsize=7, alpha=0.8,
                     xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel('log2 Fold Change (R vs NR)')
    ax.set_ylabel('-log10(FDR)')
    ax.set_title('Differential Expression: Responders vs Non-Responders')
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_13_volcano.png'), dpi=150)
    plt.close()
    print("  Saved fig_13_volcano.png")


def de_heatmap(de_df, log_expr, cohorts):
    """Heatmap of top 50 DE genes x patients, clustered."""
    # Top 25 up + top 25 down by p-value among significant
    sig = de_df[(de_df['fdr'] < 0.05) & (de_df['log2FC'].abs() > 0.5)]
    up = sig[sig['log2FC'] > 0].head(25)
    down = sig[sig['log2FC'] < 0].head(25)
    top_genes = pd.concat([up, down])['gene'].tolist()

    if len(top_genes) < 4:
        # Relax threshold if too few significant
        top_genes = de_df.head(50)['gene'].tolist()

    # Patient order: responders then non-responders
    r_ids = [p for p in cohorts['responders'] if p in log_expr.columns]
    nr_ids = [p for p in cohorts['non_responders'] if p in log_expr.columns]
    patient_order = r_ids + nr_ids

    # Z-score per gene
    hm_data = log_expr.loc[log_expr.index.isin(top_genes), patient_order]
    hm_z = hm_data.subtract(hm_data.mean(axis=1), axis=0).div(hm_data.std(axis=1), axis=0)
    hm_z = hm_z.clip(-3, 3)

    # Annotation colors
    col_colors = ['#e74c3c'] * len(r_ids) + ['#3498db'] * len(nr_ids)

    fig_height = max(6, len(top_genes) * 0.25)
    g = sns.clustermap(hm_z, col_cluster=False, row_cluster=True,
                       cmap='RdBu_r', center=0, vmin=-3, vmax=3,
                       figsize=(14, fig_height),
                       col_colors=col_colors,
                       yticklabels=True, xticklabels=False,
                       cbar_kws={'label': 'Z-score'})
    g.ax_heatmap.set_xlabel('Patients (red=CR/PR, blue=SD/PD)')
    g.fig.suptitle('Top DE Genes Heatmap', y=1.01)
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_14_de_heatmap.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved fig_14_de_heatmap.png")



def plot_ora(ora_df, filepath):
    """Bar plot of top ORA hits across all collections."""
    if ora_df.empty:
        print("  No ORA results to plot")
        return
    plot_df = ora_df[ora_df['fdr'] < 0.25].head(25).copy()
    if plot_df.empty:
        print("  No significant ORA enrichments to plot (FDR < 0.25)")
        return

    plot_df['neg_log10_fdr'] = -np.log10(plot_df['fdr'].clip(lower=1e-20))
    plot_df['label'] = (plot_df['collection'] + ' | '
                        + plot_df['direction'].str[:4] + ': '
                        + plot_df['term'])

    fig, ax = plt.subplots(figsize=(11, max(5, len(plot_df) * 0.45)))
    colors = ['#e74c3c' if 'Up' in d else '#3498db'
              for d in plot_df['direction']]
    ax.barh(range(len(plot_df)), plot_df['neg_log10_fdr'], color=colors)
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df['label'], fontsize=7)
    ax.set_xlabel('-log10(FDR)')
    ax.axvline(-np.log10(0.05), color='grey', linestyle='--', alpha=0.5,
               label='FDR = 0.05')
    ax.legend(fontsize=8)
    ax.set_title('ORA – Over-Representation Analysis (Fisher\'s exact)')
    plt.tight_layout()
    plt.savefig(filepath, dpi=150)
    plt.close()
    print(f"  Saved {os.path.basename(filepath)}")


def plot_gsea(gsea_df, filepath):
    """Dot plot of top GSEA hits (by |NES|) across all collections."""
    if gsea_df.empty:
        print("  No GSEA results to plot")
        return
    plot_df = gsea_df[gsea_df['fdr'] < 0.25].copy()
    if plot_df.empty:
        print("  No significant GSEA enrichments to plot (FDR < 0.25)")
        return
    plot_df = plot_df.reindex(
        plot_df['nes'].abs().sort_values(ascending=False).index
    ).head(25)

    fig, ax = plt.subplots(figsize=(11, max(5, len(plot_df) * 0.45)))
    colors = ['#e74c3c' if n > 0 else '#3498db' for n in plot_df['nes']]
    sizes = (-np.log10(plot_df['fdr'].clip(lower=1e-20))) * 8
    ax.scatter(plot_df['nes'], range(len(plot_df)), c=colors, s=sizes,
               alpha=0.7, edgecolors='k', linewidths=0.3)
    ax.set_yticks(range(len(plot_df)))
    labels = plot_df['collection'] + ' | ' + plot_df['term']
    ax.set_yticklabels(labels, fontsize=7)
    ax.axvline(0, color='grey', linewidth=0.5)
    ax.set_xlabel('Normalized Enrichment Score (NES)')
    ax.set_title('GSEA – Pre-ranked Gene Set Enrichment Analysis')
    plt.tight_layout()
    plt.savefig(filepath, dpi=150)
    plt.close()
    print(f"  Saved {os.path.basename(filepath)}")


def run_module4(data):
    """Execute Module 4: Differential expression & pathway analysis."""
    print("\n" + "=" * 60)
    print("MODULE 4: Differential Expression & Pathway Analysis")
    print("=" * 60)

    log_expr = data['log_expr_filtered']
    cohorts = data['cohorts']

    # 4.1 Pre-filter genes to reduce multiple-testing burden
    print("\n[4.1] Pre-filtering genes for differential expression...")
    filtered_expr, fstats = filter_genes_for_de(log_expr, cohorts)
    print(f"  Initial genes: {fstats['n_initial']}")
    print(f"  After expression filter (mean log2 >= {fstats['min_mean_expr_threshold']}): "
          f"{fstats['n_after_expression_filter']}")
    print(f"  After variance filter (bottom {int(fstats['variance_quantile_threshold']*100)}% "
          f"removed, var >= {fstats['variance_cutoff_value']}): "
          f"{fstats['n_after_variance_filter']}")

    # 4.2 Differential expression on filtered gene set
    print("\n[4.2] Differential expression analysis (Mann-Whitney U per gene)...")
    de_df = differential_expression(filtered_expr, cohorts)
    de_df.to_csv(os.path.join(RESULTS_DIR, 'table_05_de_results.csv'), index=False)
    n_sig = ((de_df['fdr'] < 0.05) & (de_df['log2FC'].abs() > 1)).sum()
    print(f"  Saved table_05_de_results.csv ({len(de_df)} genes tested)")
    print(f"  Significant (|log2FC|>1, FDR<0.05): {n_sig}")

    # Top 50 up + 50 down
    up50 = de_df[de_df['log2FC'] > 0].head(50)
    down50 = de_df[de_df['log2FC'] < 0].head(50)
    top_de = pd.concat([up50, down50])
    top_de.to_csv(os.path.join(RESULTS_DIR, 'table_06_top_de_genes.csv'), index=False)
    print(f"  Saved table_06_top_de_genes.csv ({len(top_de)} genes)")

    # 4.3 Volcano plot
    print("\n[4.3] Volcano plot...")
    volcano_plot(de_df)

    # 4.4 DE heatmap
    print("\n[4.4] DE heatmap...")
    de_heatmap(de_df, log_expr, cohorts)

    # ------------------------------------------------------------------
    # 4.5  Gene-set enrichment (ORA + GSEA) across multiple collections
    # ------------------------------------------------------------------
    print("\n[4.5] Loading gene set collections (Hallmark, KEGG, Reactome, GO BP)...")
    collections = load_collections()
    for cname, gsets in collections.items():
        print(f"  {cname}: {len(gsets)} gene sets loaded")

    tested_genes = set(de_df['gene'])

    # --- 4.6  ORA ---
    print("\n[4.6] Over-Representation Analysis (ORA, Fisher's exact)...")
    de_sorted = de_df.sort_values('p_value')
    sig_genes = de_sorted[de_sorted['fdr'] < 0.1]
    up_genes = sorted(sig_genes[sig_genes['log2FC'] > 0].head(200)['gene'])
    down_genes = sorted(sig_genes[sig_genes['log2FC'] < 0].head(200)['gene'])

    # Save ORA input gene lists for traceability
    ora_input = pd.DataFrame({
        'direction': ['up'] * len(up_genes) + ['down'] * len(down_genes),
        'gene': up_genes + down_genes,
    })
    ora_input.to_csv(os.path.join(RESULTS_DIR, 'table_06b_ora_input_genes.csv'), index=False)
    print(f"  ORA input: {len(up_genes)} up-regulated, {len(down_genes)} down-regulated genes")
    print(f"  Saved table_06b_ora_input_genes.csv")

    ora_frames = []
    for cname, raw_gsets in collections.items():
        filt = filter_gene_sets(raw_gsets, tested_genes, min_size=15, max_size=500)
        print(f"  ORA {cname}: {len(filt)} sets after size filter")
        if filt:
            ora_res = run_ora(set(up_genes), set(down_genes), tested_genes,
                              filt, collection_name=cname)
            ora_frames.append(ora_res)

    ora_df = pd.concat(ora_frames, ignore_index=True) if ora_frames else pd.DataFrame()
    if not ora_df.empty:
        ora_df.to_csv(os.path.join(RESULTS_DIR, 'table_06c_ora_enrichment.csv'), index=False)
        n_sig_ora = (ora_df['fdr'] < 0.05).sum()
        print(f"  Saved table_06c_ora_enrichment.csv  ({len(ora_df)} tests, {n_sig_ora} FDR<0.05)")
        plot_ora(ora_df, os.path.join(RESULTS_DIR, 'fig_16a_ora_enrichment.png'))
    else:
        print("  No ORA results")

    # --- 4.7  GSEA (pre-ranked) ---
    print("\n[4.7] Pre-ranked Gene Set Enrichment Analysis (GSEA)...")

    # Save GSEA ranking for traceability
    gsea_ranking = de_df[['gene', 'log2FC', 'p_value']].copy()
    gsea_ranking['ranking_metric'] = (
        np.sign(gsea_ranking['log2FC'])
        * -np.log10(gsea_ranking['p_value'].clip(lower=1e-300))
    )
    gsea_ranking = gsea_ranking.sort_values('ranking_metric', ascending=False)
    gsea_ranking.to_csv(os.path.join(RESULTS_DIR, 'table_06b_gsea_input_ranking.csv'), index=False)
    print(f"  GSEA input: {len(gsea_ranking)} ranked genes")
    print(f"  Saved table_06b_gsea_input_ranking.csv")

    gsea_frames = []
    for cname, raw_gsets in collections.items():
        filt = filter_gene_sets(raw_gsets, tested_genes, min_size=15, max_size=500)
        print(f"  GSEA {cname}: {len(filt)} sets after size filter ...")
        if filt:
            gsea_res = run_gsea_prerank(de_df, filt, collection_name=cname,
                                        nperm=1000, min_size=15, max_size=500)
            gsea_frames.append(gsea_res)

    gsea_df = pd.concat(gsea_frames, ignore_index=True) if gsea_frames else pd.DataFrame()
    if not gsea_df.empty:
        gsea_df.to_csv(os.path.join(RESULTS_DIR, 'table_06d_gsea_enrichment.csv'), index=False)
        n_sig_gsea = (gsea_df['fdr'] < 0.05).sum()
        print(f"  Saved table_06d_gsea_enrichment.csv  ({len(gsea_df)} tests, {n_sig_gsea} FDR<0.05)")
        plot_gsea(gsea_df, os.path.join(RESULTS_DIR, 'fig_16b_gsea_enrichment.png'))
    else:
        print("  No GSEA results")

    print("\n[Module 4 Complete]")
    return {
        'de_results': de_df,
        'ora_enrichment': ora_df,
        'gsea_enrichment': gsea_df,
    }


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/code')
    from module1_data import run_module1
    os.makedirs(RESULTS_DIR, exist_ok=True)
    data = run_module1()
    run_module4(data)
