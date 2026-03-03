"""
Module 1: Data Loading, QC & Preprocessing
IMvigor210 Biomarker Analysis
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

DATA_DIR = '/root/capsule/data'
RESULTS_DIR = '/results'
SCRATCH_DIR = '/scratch'


def load_clinical(path=None):
    """Load clinical data."""
    if path is None:
        path = os.path.join(DATA_DIR, 'imvigor210_clinical_data.csv')
    df = pd.read_csv(path)
    df = df.set_index('patient_id')
    return df


def load_expression(path=None):
    """Load gene expression data. Returns genes x patients matrix."""
    if path is None:
        path = os.path.join(DATA_DIR, 'imvigor210_gene_expression.csv')
    df = pd.read_csv(path)
    # Store gene metadata
    gene_info = df[['gene_id', 'gene_symbol', 'entrez_id']].copy()
    # Set gene_symbol as index, drop metadata columns
    expr = df.drop(columns=['gene_id', 'gene_symbol', 'entrez_id'])
    expr.index = df['gene_symbol']
    expr.index.name = 'gene_symbol'
    return expr, gene_info


def load_genomic(path=None):
    """Load genomic alteration data in long format."""
    if path is None:
        path = os.path.join(DATA_DIR, 'imvigor210_genomic_alterations.csv')
    df = pd.read_csv(path)
    return df


def normalize_expression(expr_df, size_factors):
    """
    Normalize expression using size factors and log2 transform.
    Formula: log2(count / sizeFactor + 1)
    """
    # Align size factors to expression columns
    common = expr_df.columns.intersection(size_factors.index)
    expr_aligned = expr_df[common]
    sf_aligned = size_factors[common]

    # Divide each column by its size factor, then log2(x + 1)
    normalized = expr_aligned.div(sf_aligned, axis=1)
    log_expr = np.log2(normalized + 1)
    return log_expr


def filter_low_expression(expr_df, min_count=10, min_samples=15):
    """
    Remove lowly expressed genes.
    Keep genes with at least min_count in at least min_samples.
    Then deduplicate gene symbols by keeping the one with highest mean.
    """
    # Filter: gene must have >= min_count in >= min_samples
    mask = (expr_df >= min_count).sum(axis=1) >= min_samples
    filtered = expr_df[mask].copy()

    # Deduplicate gene symbols: keep row with highest mean expression
    filtered['_mean'] = filtered.mean(axis=1)
    filtered = filtered.sort_values('_mean', ascending=False)
    filtered = filtered[~filtered.index.duplicated(keep='first')]
    filtered = filtered.drop(columns=['_mean'])

    return filtered


def create_mutation_matrix(genomic_df, patient_ids):
    """
    Create binary patient x gene mutation matrix.
    1 = any alteration present, 0 = no alteration.
    Patients without genomic data get NaN.
    """
    # Pivot: for each patient-gene pair, mark as 1 if any alteration exists
    if genomic_df.empty:
        return pd.DataFrame(index=patient_ids)

    mut_pairs = genomic_df[['patient_id', 'gene']].drop_duplicates()
    mut_pairs['mutated'] = 1
    mut_matrix = mut_pairs.pivot_table(
        index='patient_id', columns='gene', values='mutated', fill_value=0
    )

    # Reindex to all patients - those without genomic data get NaN
    patients_with_data = set(genomic_df['patient_id'].unique())
    mut_matrix = mut_matrix.reindex(patient_ids)
    # Set patients without any genomic data to NaN
    for pid in patient_ids:
        if pid not in patients_with_data:
            mut_matrix.loc[pid] = np.nan

    return mut_matrix


def prepare_cohorts(clinical_df):
    """
    Define analysis cohorts.
    Returns dict with:
      - 'binary': DataFrame with responders (CR/PR) and non-responders (SD/PD), NE excluded
      - 'survival': Full DataFrame (all patients including NE) for survival analysis
      - 'responders': patient_ids of CR/PR
      - 'non_responders': patient_ids of SD/PD
    """
    binary_df = clinical_df[clinical_df['binaryResponse'].isin(['CR/PR', 'SD/PD'])].copy()
    binary_df['response_binary'] = (binary_df['binaryResponse'] == 'CR/PR').astype(int)

    responders = binary_df[binary_df['response_binary'] == 1].index.tolist()
    non_responders = binary_df[binary_df['response_binary'] == 0].index.tolist()

    return {
        'binary': binary_df,
        'survival': clinical_df.copy(),
        'responders': responders,
        'non_responders': non_responders,
    }


def run_module1():
    """Execute Module 1: load, preprocess, QC."""
    print("=" * 60)
    print("MODULE 1: Data Loading, QC & Preprocessing")
    print("=" * 60)

    # --- Load data ---
    print("\n[1.1] Loading clinical data...")
    clinical = load_clinical()
    print(f"  Clinical: {clinical.shape[0]} patients, {clinical.shape[1]} variables")

    print("[1.2] Loading gene expression data...")
    expr_raw, gene_info = load_expression()
    print(f"  Expression: {expr_raw.shape[0]} genes, {expr_raw.shape[1]} samples")

    print("[1.3] Loading genomic alteration data...")
    genomic = load_genomic()
    print(f"  Genomic: {genomic.shape[0]} records, {genomic['patient_id'].nunique()} patients")

    # --- Verify sample overlap ---
    print("\n[1.4] Verifying sample overlap...")
    clin_ids = set(clinical.index)
    expr_ids = set(expr_raw.columns)
    gen_ids = set(genomic['patient_id'].unique())
    print(f"  Clinical IDs: {len(clin_ids)}")
    print(f"  Expression IDs: {len(expr_ids)}")
    print(f"  Genomic IDs: {len(gen_ids)}")
    print(f"  Clinical ∩ Expression: {len(clin_ids & expr_ids)}")
    print(f"  Clinical ∩ Genomic: {len(clin_ids & gen_ids)}")

    # --- Normalize expression ---
    print("\n[1.5] Normalizing gene expression...")
    size_factors = clinical['sizeFactor']
    # Use only patients present in both clinical and expression
    common_patients = sorted(clin_ids & expr_ids)
    size_factors = size_factors.reindex(common_patients).dropna()

    log_expr = normalize_expression(expr_raw, size_factors)
    print(f"  Normalized expression: {log_expr.shape[0]} genes, {log_expr.shape[1]} samples")

    # --- Filter low expression ---
    print("[1.6] Filtering low-expression genes...")
    log_expr_filtered = filter_low_expression(log_expr)
    print(f"  After filtering: {log_expr_filtered.shape[0]} genes retained")

    # --- Mutation matrix ---
    print("[1.7] Creating mutation matrix...")
    mut_matrix = create_mutation_matrix(genomic, list(clinical.index))
    n_with_data = mut_matrix.dropna(how='all').shape[0]
    print(f"  Mutation matrix: {mut_matrix.shape[0]} patients, {mut_matrix.shape[1]} genes")
    print(f"  Patients with genomic data: {n_with_data}")

    # --- Prepare cohorts ---
    print("\n[1.8] Preparing analysis cohorts...")
    cohorts = prepare_cohorts(clinical)
    print(f"  Binary analysis: {len(cohorts['responders'])} responders, "
          f"{len(cohorts['non_responders'])} non-responders")
    print(f"  Survival analysis: {cohorts['survival'].shape[0]} patients")

    # --- Cohort summary table ---
    print("\n[1.9] Generating cohort summary...")
    summary_rows = []
    for var in ['Best Confirmed Overall Response', 'binaryResponse', 'Sex',
                'Enrollment IC', 'Immune phenotype', 'TCGA Subtype',
                'Baseline ECOG Score', 'Received platinum']:
        if var in clinical.columns:
            counts = clinical[var].value_counts(dropna=False)
            for val, cnt in counts.items():
                summary_rows.append({
                    'variable': var,
                    'value': str(val),
                    'count': cnt,
                    'percent': round(100.0 * cnt / len(clinical), 1)
                })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(RESULTS_DIR, 'table_01_cohort_summary.csv'), index=False)
    print(f"  Saved table_01_cohort_summary.csv")

    # --- Figure 1: Expression distribution pre/post normalization ---
    print("\n[1.10] Generating expression distribution plots...")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Pre-normalization (raw counts, log2 for visualization)
    raw_sample = expr_raw.iloc[:, :5]  # first 5 samples
    for col in raw_sample.columns:
        vals = np.log2(raw_sample[col] + 1)
        vals = vals[vals > 0]
        axes[0].hist(vals, bins=50, alpha=0.5, label=col[:10])
    axes[0].set_title('Raw Counts (log2 + 1)')
    axes[0].set_xlabel('log2(count + 1)')
    axes[0].set_ylabel('Frequency')
    axes[0].legend(fontsize=6)

    # Post-normalization
    norm_sample = log_expr_filtered.iloc[:, :5]
    for col in norm_sample.columns:
        vals = norm_sample[col]
        vals = vals[vals > 0]
        axes[1].hist(vals, bins=50, alpha=0.5, label=col[:10])
    axes[1].set_title('Normalized (log2(count/SF + 1))')
    axes[1].set_xlabel('log2(normalized count + 1)')
    axes[1].set_ylabel('Frequency')
    axes[1].legend(fontsize=6)

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_01_expression_distribution.png'), dpi=150)
    plt.close()
    print("  Saved fig_01_expression_distribution.png")

    # --- Figure 2: PCA of expression ---
    print("[1.11] Running PCA on normalized expression...")
    # Transpose: samples x genes
    expr_t = log_expr_filtered.T
    # Only keep patients in binary cohort for coloring
    binary_patients = cohorts['binary'].index
    expr_pca_input = expr_t.reindex(binary_patients).dropna()

    scaler = StandardScaler()
    expr_pca_input.columns = expr_pca_input.columns.astype(str)
    expr_scaled = scaler.fit_transform(expr_pca_input)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(expr_scaled)

    pca_df = pd.DataFrame({
        'PC1': pcs[:, 0],
        'PC2': pcs[:, 1],
        'Response': cohorts['binary'].loc[expr_pca_input.index, 'binaryResponse']
    })

    fig, ax = plt.subplots(figsize=(8, 6))
    colors = {'CR/PR': '#e74c3c', 'SD/PD': '#3498db'}
    for resp, color in colors.items():
        mask = pca_df['Response'] == resp
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                   c=color, label=resp, alpha=0.6, s=30)
    ax.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    ax.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    ax.set_title('PCA of Gene Expression by Response')
    ax.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_02_pca_expression.png'), dpi=150)
    plt.close()
    print("  Saved fig_02_pca_expression.png")

    # --- Return all processed data ---
    results = {
        'clinical': clinical,
        'expr_raw': expr_raw,
        'log_expr': log_expr,
        'log_expr_filtered': log_expr_filtered,
        'gene_info': gene_info,
        'genomic': genomic,
        'mut_matrix': mut_matrix,
        'cohorts': cohorts,
    }

    print("\n[Module 1 Complete]")
    return results


if __name__ == '__main__':
    os.makedirs(RESULTS_DIR, exist_ok=True)
    run_module1()
