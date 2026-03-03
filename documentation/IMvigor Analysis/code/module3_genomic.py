"""
Module 3: Genomic Alteration Analysis
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
from sklearn.metrics import roc_curve, auc

RESULTS_DIR = '/results'


def mutation_frequency_analysis(mut_matrix, cohorts):
    """
    Per-gene mutation frequency in responders vs non-responders.
    Returns filtered DataFrame with Fisher's exact p-values.
    """
    binary_df = cohorts['binary']
    responders = cohorts['responders']
    non_responders = cohorts['non_responders']

    # Subset to patients with genomic data AND in binary cohort
    all_binary = responders + non_responders
    mut_binary = mut_matrix.reindex(all_binary).dropna(how='all')

    r_ids = [p for p in responders if p in mut_binary.index]
    nr_ids = [p for p in non_responders if p in mut_binary.index]

    results = []
    for gene in mut_binary.columns:
        r_mut = int(mut_binary.loc[r_ids, gene].sum())
        r_wt = len(r_ids) - r_mut
        nr_mut = int(mut_binary.loc[nr_ids, gene].sum())
        nr_wt = len(nr_ids) - nr_mut

        total_mut = r_mut + nr_mut
        if total_mut < 10:
            continue

        table = [[r_mut, r_wt], [nr_mut, nr_wt]]
        odds, p = stats.fisher_exact(table)

        results.append({
            'gene': gene,
            'resp_mutated': r_mut,
            'resp_total': len(r_ids),
            'resp_freq': round(r_mut / len(r_ids), 4),
            'nonresp_mutated': nr_mut,
            'nonresp_total': len(nr_ids),
            'nonresp_freq': round(nr_mut / len(nr_ids), 4),
            'total_mutated': total_mut,
            'odds_ratio': round(odds, 4),
            'p_value': p,
        })

    res_df = pd.DataFrame(results)
    if len(res_df) > 0:
        _, fdr, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['fdr'] = fdr
        res_df = res_df.sort_values('p_value')

    return res_df


def mutation_frequency_barplot(mut_assoc_df):
    """Horizontal bar chart comparing mutation frequency in R vs NR for top 20 genes."""
    top = mut_assoc_df.head(20).copy()
    top = top.sort_values('total_mutated', ascending=True)

    fig, ax = plt.subplots(figsize=(10, 8))
    y = np.arange(len(top))
    height = 0.35

    ax.barh(y - height / 2, top['resp_freq'] * 100, height, label='CR/PR', color='#e74c3c')
    ax.barh(y + height / 2, top['nonresp_freq'] * 100, height, label='SD/PD', color='#3498db')

    ax.set_yticks(y)
    ax.set_yticklabels(top['gene'], fontsize=9)
    ax.set_xlabel('Mutation Frequency (%)')
    ax.set_title('Mutation Frequency: Responders vs Non-Responders (Top 20 Genes)')
    ax.legend()

    # Mark significant genes
    for i, (_, row) in enumerate(top.iterrows()):
        if row.get('fdr', 1) < 0.05:
            ax.text(max(row['resp_freq'], row['nonresp_freq']) * 100 + 1,
                    i, '*', fontsize=12, color='red', va='center')

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_09_mutation_freq_comparison.png'), dpi=150)
    plt.close()
    print("  Saved fig_09_mutation_freq_comparison.png")


def oncoplot(mut_matrix, genomic_df, cohorts):
    """
    Oncoplot: top 20 most-mutated genes x patients, colored by alteration type.
    """
    binary_df = cohorts['binary']
    all_binary = cohorts['responders'] + cohorts['non_responders']
    mut_binary = mut_matrix.reindex(all_binary).dropna(how='all')

    # Top 20 genes by total mutation count
    gene_counts = mut_binary.sum().sort_values(ascending=False)
    top_genes = gene_counts.head(20).index.tolist()

    # Build alteration type matrix
    alt_records = genomic_df[genomic_df['patient_id'].isin(mut_binary.index) &
                             genomic_df['gene'].isin(top_genes)].copy()

    # Simplify alteration types
    alt_map = {}
    for _, row in alt_records.iterrows():
        key = (row['patient_id'], row['gene'])
        alt_map[key] = row['alteration_type']

    # Sort patients by response then mutation count
    patient_order = []
    for resp_status in ['CR/PR', 'SD/PD']:
        resp_patients = binary_df[binary_df['binaryResponse'] == resp_status].index
        sub = mut_binary.reindex([p for p in resp_patients if p in mut_binary.index])
        sub_sorted = sub[top_genes].sum(axis=1).sort_values(ascending=False)
        patient_order.extend(sub_sorted.index.tolist())

    # Build numeric matrix for plotting
    alt_type_colors = {
        'known_short': 0, 'likely_short': 1, 'known_rearrangement': 2,
        'likely_rearrangement': 3, 'known_copy_number': 4, 'ambiguous': 5,
    }
    n_types = len(alt_type_colors) + 1  # +1 for 'other'

    plot_matrix = np.full((len(top_genes), len(patient_order)), np.nan)
    for gi, gene in enumerate(top_genes):
        for pi, patient in enumerate(patient_order):
            if mut_binary.loc[patient, gene] == 1:
                key = (patient, gene)
                alt = alt_map.get(key, 'other')
                plot_matrix[gi, pi] = alt_type_colors.get(alt, n_types - 1)

    fig, (ax_top, ax_main) = plt.subplots(2, 1, figsize=(16, 10),
                                            gridspec_kw={'height_ratios': [1, 8]})

    # Top bar: response annotation
    resp_colors = []
    for p in patient_order:
        if p in cohorts['responders']:
            resp_colors.append('#e74c3c')
        else:
            resp_colors.append('#3498db')
    ax_top.bar(range(len(patient_order)), [1] * len(patient_order), color=resp_colors, width=1.0)
    ax_top.set_xlim(-0.5, len(patient_order) - 0.5)
    ax_top.set_xticks([])
    ax_top.set_yticks([])
    ax_top.set_title('Oncoplot: Top 20 Mutated Genes')

    # Legend for response
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#e74c3c', label='CR/PR'),
                       Patch(facecolor='#3498db', label='SD/PD')]
    ax_top.legend(handles=legend_elements, loc='upper right', fontsize=8)

    # Main heatmap
    cmap = plt.cm.get_cmap('Set3', n_types)
    masked = np.ma.masked_invalid(plot_matrix)
    ax_main.pcolormesh(masked, cmap=cmap, edgecolors='white', linewidth=0.5)
    ax_main.set_yticks(np.arange(len(top_genes)) + 0.5)
    ax_main.set_yticklabels(top_genes, fontsize=9)
    ax_main.set_xticks([])
    ax_main.set_xlabel(f'Patients (n={len(patient_order)})')
    ax_main.invert_yaxis()

    # Gene frequency annotation on the right
    for gi, gene in enumerate(top_genes):
        freq = gene_counts[gene]
        pct = 100 * freq / len(patient_order)
        ax_main.text(len(patient_order) + 1, gi + 0.5,
                     f'{int(freq)} ({pct:.0f}%)', va='center', fontsize=8)

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_10_oncoplot.png'), dpi=150)
    plt.close()
    print("  Saved fig_10_oncoplot.png")


def tmb_roc(binary_df):
    """ROC curve for TMB as continuous predictor of response."""
    sub = binary_df[['FMOne mutation burden per MB', 'response_binary']].dropna()
    y = sub['response_binary'].values
    tmb = sub['FMOne mutation burden per MB'].values

    fpr, tpr, thresholds = roc_curve(y, tmb)
    roc_auc = auc(fpr, tpr)

    # Youden's J
    j_scores = tpr - fpr
    best_idx = np.argmax(j_scores)
    best_threshold = thresholds[best_idx]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.plot(fpr, tpr, color='darkorange', lw=2, label=f'TMB (AUC = {roc_auc:.3f})')
    ax.plot([0, 1], [0, 1], 'k--', alpha=0.5)
    ax.scatter(fpr[best_idx], tpr[best_idx], color='red', s=100, zorder=5,
               label=f'Optimal threshold = {best_threshold:.1f}')
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Curve: TMB as Predictor of Response')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_11_tmb_roc.png'), dpi=150)
    plt.close()
    print(f"  Saved fig_11_tmb_roc.png (AUC={roc_auc:.3f}, threshold={best_threshold:.1f})")
    return roc_auc, best_threshold


def co_occurrence_matrix(mut_matrix, cohorts):
    """
    Co-occurrence / mutual exclusivity matrix of top 15 genes.
    Fisher's exact test pairwise.
    """
    all_binary = cohorts['responders'] + cohorts['non_responders']
    mut_binary = mut_matrix.reindex(all_binary).dropna(how='all')

    gene_counts = mut_binary.sum().sort_values(ascending=False)
    top_genes = gene_counts.head(15).index.tolist()
    sub = mut_binary[top_genes]

    n_genes = len(top_genes)
    log_or_matrix = np.zeros((n_genes, n_genes))
    p_matrix = np.ones((n_genes, n_genes))

    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            g1 = sub.iloc[:, i]
            g2 = sub.iloc[:, j]
            both = ((g1 == 1) & (g2 == 1)).sum()
            g1_only = ((g1 == 1) & (g2 == 0)).sum()
            g2_only = ((g1 == 0) & (g2 == 1)).sum()
            neither = ((g1 == 0) & (g2 == 0)).sum()

            table = [[both, g1_only], [g2_only, neither]]
            odds, p = stats.fisher_exact(table)
            log_or = np.log2(odds) if odds > 0 else 0
            log_or_matrix[i, j] = log_or
            log_or_matrix[j, i] = log_or
            p_matrix[i, j] = p
            p_matrix[j, i] = p

    fig, ax = plt.subplots(figsize=(10, 8))
    mask = np.eye(n_genes, dtype=bool)
    # Clip for visualization
    clipped = np.clip(log_or_matrix, -3, 3)
    sns.heatmap(clipped, mask=mask, xticklabels=top_genes, yticklabels=top_genes,
                cmap='RdBu_r', center=0, vmin=-3, vmax=3, ax=ax,
                annot=True, fmt='.1f', annot_kws={'fontsize': 7})

    # Mark significant pairs
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            if p_matrix[i, j] < 0.05:
                ax.text(j + 0.5, i + 0.85, '*', ha='center', va='center',
                        fontsize=10, color='black', fontweight='bold')

    ax.set_title('Co-occurrence (red) / Mutual Exclusivity (blue)\nlog2(OR), * = p < 0.05')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_12_co_occurrence.png'), dpi=150)
    plt.close()
    print("  Saved fig_12_co_occurrence.png")


def run_module3(data):
    """Execute Module 3: Genomic alteration analysis."""
    print("\n" + "=" * 60)
    print("MODULE 3: Genomic Alteration Analysis")
    print("=" * 60)

    mut_matrix = data['mut_matrix']
    cohorts = data['cohorts']
    binary_df = cohorts['binary']

    # 3.1 Mutation frequency and association
    print("\n[3.1] Mutation frequency per gene & Fisher's exact tests...")
    mut_assoc = mutation_frequency_analysis(mut_matrix, cohorts)
    mut_assoc.to_csv(os.path.join(RESULTS_DIR, 'table_04_mutation_association.csv'), index=False)
    print(f"  Saved table_04_mutation_association.csv ({len(mut_assoc)} genes tested)")
    if 'fdr' in mut_assoc.columns:
        sig = mut_assoc[mut_assoc['fdr'] < 0.05]
        print(f"  Significant after FDR: {len(sig)} genes")

    # 3.2 Mutation frequency bar plot
    print("\n[3.2] Mutation frequency comparison plot...")
    if len(mut_assoc) > 0:
        mutation_frequency_barplot(mut_assoc)

    # 3.3 Oncoplot
    print("\n[3.3] Oncoplot...")
    oncoplot(mut_matrix, data['genomic'], cohorts)

    # 3.4 TMB ROC curve
    print("\n[3.4] TMB ROC analysis...")
    tmb_auc, tmb_threshold = tmb_roc(binary_df)

    # 3.5 Co-occurrence matrix
    print("\n[3.5] Co-occurrence / mutual exclusivity analysis...")
    co_occurrence_matrix(mut_matrix, cohorts)

    print("\n[Module 3 Complete]")
    return {
        'mut_association': mut_assoc,
        'tmb_auc': tmb_auc,
        'tmb_threshold': tmb_threshold,
    }


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/code')
    from module1_data import run_module1
    os.makedirs(RESULTS_DIR, exist_ok=True)
    data = run_module1()
    run_module3(data)
