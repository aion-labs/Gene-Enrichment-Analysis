"""
Module 5: Predictive Gene Expression Signature
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
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_curve, auc, roc_auc_score

RESULTS_DIR = '/results'


def build_feature_matrix(log_expr, de_df, cohorts, top_n_variance=5000, de_pval_thresh=0.05):
    """
    Feature preselection: top by variance, then filter to DE p < threshold.
    Returns X (patients x genes), y (binary response), gene list.
    """
    r_ids = [p for p in cohorts['responders'] if p in log_expr.columns]
    nr_ids = [p for p in cohorts['non_responders'] if p in log_expr.columns]
    all_ids = r_ids + nr_ids

    expr_sub = log_expr[all_ids].T  # patients x genes
    y = np.array([1] * len(r_ids) + [0] * len(nr_ids))

    return expr_sub, y, all_ids


def feature_preselection(X_train, y_train, gene_names, top_n_variance=5000, de_pval_thresh=0.05):
    """
    Recompute feature preselection within a CV fold to avoid leakage.
    1. Top by variance
    2. Then filter to genes with Mann-Whitney p < threshold (training set only)
    """
    # Variance filter
    var = X_train.var(axis=0)
    top_var_idx = var.nlargest(top_n_variance).index
    X_filt = X_train[top_var_idx]

    # DE filter on training set
    selected_genes = []
    for gene in X_filt.columns:
        r_vals = X_filt.loc[y_train == 1, gene].values
        nr_vals = X_filt.loc[y_train == 0, gene].values
        if np.std(r_vals) == 0 and np.std(nr_vals) == 0:
            continue
        _, p = stats.mannwhitneyu(r_vals, nr_vals, alternative='two-sided')
        if p < de_pval_thresh:
            selected_genes.append(gene)

    if len(selected_genes) < 10:
        # Fallback: use top variance genes
        selected_genes = top_var_idx[:500].tolist()

    return selected_genes


def nested_cv_elastic_net(X, y, gene_names, n_outer=5, n_inner=5):
    """
    Nested CV: outer loop for performance, inner loop for tuning C.
    Returns per-fold AUCs, mean AUC, fold predictions.
    """
    outer_cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=42)
    C_grid = [0.001, 0.01, 0.1, 1.0, 10.0]

    fold_aucs = []
    fold_fprs = []
    fold_tprs = []
    all_y_true = []
    all_y_prob = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Feature preselection on training set only
        selected = feature_preselection(X_train, y_train, gene_names)
        X_train_sel = X_train[selected]
        X_test_sel = X_test[selected]

        # Scale
        scaler = StandardScaler()
        X_train_sc = scaler.fit_transform(X_train_sel)
        X_test_sc = scaler.transform(X_test_sel)

        # Inner CV for hyperparameter tuning
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=fold_idx)
        best_c = 1.0
        best_auc = 0

        for c in C_grid:
            inner_aucs = []
            for itrain, ival in inner_cv.split(X_train_sc, y_train):
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    model = LogisticRegression(C=c, penalty='elasticnet', solver='saga',
                                               l1_ratio=0.5, max_iter=2000, random_state=42)
                    model.fit(X_train_sc[itrain], y_train[itrain])
                    probs = model.predict_proba(X_train_sc[ival])[:, 1]
                    try:
                        inner_aucs.append(roc_auc_score(y_train[ival], probs))
                    except ValueError:
                        inner_aucs.append(0.5)
            mean_inner = np.mean(inner_aucs)
            if mean_inner > best_auc:
                best_auc = mean_inner
                best_c = c

        # Refit on full training set with best C
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            model = LogisticRegression(C=best_c, penalty='elasticnet', solver='saga',
                                       l1_ratio=0.5, max_iter=2000, random_state=42)
            model.fit(X_train_sc, y_train)

        probs = model.predict_proba(X_test_sc)[:, 1]
        try:
            fold_auc = roc_auc_score(y_test, probs)
        except ValueError:
            fold_auc = 0.5
        fold_aucs.append(fold_auc)

        fpr, tpr, _ = roc_curve(y_test, probs)
        fold_fprs.append(fpr)
        fold_tprs.append(tpr)

        all_y_true.extend(y_test)
        all_y_prob.extend(probs)

        print(f"    Fold {fold_idx+1}: AUC={fold_auc:.3f} (C={best_c}, {len(selected)} features)")

    return {
        'fold_aucs': fold_aucs,
        'fold_fprs': fold_fprs,
        'fold_tprs': fold_tprs,
        'all_y_true': np.array(all_y_true),
        'all_y_prob': np.array(all_y_prob),
    }


def nested_cv_random_forest(X, y, gene_names, n_outer=5, n_inner=5):
    """
    Nested CV for random forest classifier.
    Returns per-fold AUCs, mean AUC, fold predictions.
    """
    outer_cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=42)

    param_grid = {
        'n_estimators': [100, 300],
        'max_depth': [5, 10, None],
        'min_samples_leaf': [5, 10],
    }

    fold_aucs = []
    fold_fprs = []
    fold_tprs = []
    all_y_true = []
    all_y_prob = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Feature preselection on training set only
        selected = feature_preselection(X_train, y_train, gene_names)
        X_train_sel = X_train[selected]
        X_test_sel = X_test[selected]

        # Scale (RF doesn't need it, but for consistency)
        scaler = StandardScaler()
        X_train_sc = pd.DataFrame(scaler.fit_transform(X_train_sel),
                                   columns=X_train_sel.columns, index=X_train_sel.index)
        X_test_sc = pd.DataFrame(scaler.transform(X_test_sel),
                                  columns=X_test_sel.columns, index=X_test_sel.index)

        # Inner CV
        inner_cv = StratifiedKFold(n_splits=n_inner, shuffle=True, random_state=fold_idx)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            gs = GridSearchCV(
                RandomForestClassifier(random_state=42, n_jobs=1),
                param_grid, cv=inner_cv, scoring='roc_auc', n_jobs=1
            )
            gs.fit(X_train_sc, y_train)

        probs = gs.predict_proba(X_test_sc)[:, 1]
        try:
            fold_auc = roc_auc_score(y_test, probs)
        except ValueError:
            fold_auc = 0.5
        fold_aucs.append(fold_auc)

        fpr, tpr, _ = roc_curve(y_test, probs)
        fold_fprs.append(fpr)
        fold_tprs.append(tpr)

        all_y_true.extend(y_test)
        all_y_prob.extend(probs)

        print(f"    Fold {fold_idx+1}: AUC={fold_auc:.3f} (best params: {gs.best_params_})")

    return {
        'fold_aucs': fold_aucs,
        'fold_fprs': fold_fprs,
        'fold_tprs': fold_tprs,
        'all_y_true': np.array(all_y_true),
        'all_y_prob': np.array(all_y_prob),
    }


def baseline_comparators(cohorts):
    """
    Compute AUC for TMB-only and PD-L1 IC-only as baseline comparators.
    Returns dict with ROC data.
    """
    binary_df = cohorts['binary']
    y = binary_df['response_binary'].values
    comparators = {}

    # TMB
    tmb_col = 'FMOne mutation burden per MB'
    if tmb_col in binary_df.columns:
        sub = binary_df[[tmb_col, 'response_binary']].dropna()
        fpr, tpr, _ = roc_curve(sub['response_binary'], sub[tmb_col])
        comparators['TMB'] = {
            'fpr': fpr, 'tpr': tpr, 'auc': auc(fpr, tpr),
            'y_true': sub['response_binary'].values, 'y_prob': sub[tmb_col].values,
        }

    # PD-L1 IC Level (ordinal)
    if 'IC Level' in binary_df.columns:
        ic_map = {'IC0': 0, 'IC1': 1, 'IC2+': 2}
        sub = binary_df[['IC Level', 'response_binary']].dropna()
        sub['ic_score'] = sub['IC Level'].map(ic_map)
        sub = sub.dropna(subset=['ic_score'])
        fpr, tpr, _ = roc_curve(sub['response_binary'], sub['ic_score'])
        comparators['PD-L1 IC'] = {
            'fpr': fpr, 'tpr': tpr, 'auc': auc(fpr, tpr),
            'y_true': sub['response_binary'].values, 'y_prob': sub['ic_score'].values,
        }

    return comparators


def extract_elastic_net_signature(X, y, gene_names, best_c=None, de_df=None):
    """
    Refit elastic net on full dataset with best C and extract non-zero coefficients.
    """
    # Feature preselection on full dataset
    selected = feature_preselection(X, y, gene_names)
    X_sel = X[selected]

    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X_sel)

    if best_c is None:
        best_c = 0.1  # default

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = LogisticRegression(C=best_c, penalty='elasticnet', solver='saga',
                                   l1_ratio=0.5, max_iter=2000, random_state=42)
        model.fit(X_sc, y)

    coefs = pd.DataFrame({
        'gene': selected,
        'coefficient': model.coef_[0],
    })
    coefs['abs_coef'] = coefs['coefficient'].abs()
    coefs = coefs[coefs['coefficient'] != 0].sort_values('abs_coef', ascending=False)
    coefs['rank'] = range(1, len(coefs) + 1)
    coefs = coefs.drop(columns='abs_coef')

    # Compute gene signature score for all patients
    sig_genes = coefs['gene'].tolist()
    sig_coefs = coefs.set_index('gene')['coefficient']

    return coefs, model, scaler, selected, sig_genes, sig_coefs


def extract_rf_importance(X, y, gene_names):
    """
    Refit RF on full dataset and extract top 30 feature importances.
    """
    selected = feature_preselection(X, y, gene_names)
    X_sel = X[selected]

    scaler = StandardScaler()
    X_sc = pd.DataFrame(scaler.fit_transform(X_sel), columns=X_sel.columns, index=X_sel.index)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = RandomForestClassifier(n_estimators=300, max_depth=10,
                                        min_samples_leaf=5, random_state=42, n_jobs=1)
        model.fit(X_sc, y)

    importances = pd.DataFrame({
        'gene': selected,
        'importance': model.feature_importances_,
    }).sort_values('importance', ascending=False)
    importances['rank'] = range(1, len(importances) + 1)
    top30 = importances.head(30)

    return top30, model


def plot_roc_comparison(en_results, rf_results, comparators):
    """Plot all ROC curves on same figure."""
    fig, ax = plt.subplots(figsize=(8, 7))

    # Mean ROC for elastic net (interpolated)
    mean_fpr = np.linspace(0, 1, 100)
    en_tprs = []
    for fpr, tpr in zip(en_results['fold_fprs'], en_results['fold_tprs']):
        en_tprs.append(np.interp(mean_fpr, fpr, tpr))
    mean_en_tpr = np.mean(en_tprs, axis=0)
    mean_en_auc = np.mean(en_results['fold_aucs'])
    std_en_auc = np.std(en_results['fold_aucs'])
    ax.plot(mean_fpr, mean_en_tpr, color='darkred', lw=2,
            label=f'Elastic Net (AUC={mean_en_auc:.3f}±{std_en_auc:.3f})')

    # Mean ROC for RF
    rf_tprs = []
    for fpr, tpr in zip(rf_results['fold_fprs'], rf_results['fold_tprs']):
        rf_tprs.append(np.interp(mean_fpr, fpr, tpr))
    mean_rf_tpr = np.mean(rf_tprs, axis=0)
    mean_rf_auc = np.mean(rf_results['fold_aucs'])
    std_rf_auc = np.std(rf_results['fold_aucs'])
    ax.plot(mean_fpr, mean_rf_tpr, color='darkgreen', lw=2,
            label=f'Random Forest (AUC={mean_rf_auc:.3f}±{std_rf_auc:.3f})')

    # Comparators
    colors = {'TMB': 'orange', 'PD-L1 IC': 'purple'}
    for name, data in comparators.items():
        ax.plot(data['fpr'], data['tpr'], color=colors.get(name, 'grey'), lw=1.5,
                linestyle='--', label=f'{name} (AUC={data["auc"]:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('ROC Comparison: Gene Signature vs Individual Biomarkers')
    ax.legend(loc='lower right')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_17_roc_comparison.png'), dpi=150)
    plt.close()
    print("  Saved fig_17_roc_comparison.png")


def plot_signature_coefficients(en_coefs):
    """Bar plot of elastic net signature coefficients."""
    if len(en_coefs) == 0:
        print("  No non-zero elastic net coefficients to plot")
        return
    top = en_coefs.head(30).copy()
    top = top.sort_values('coefficient')

    fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.3)))
    colors = ['#e74c3c' if c > 0 else '#3498db' for c in top['coefficient']]
    ax.barh(range(len(top)), top['coefficient'], color=colors)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['gene'], fontsize=8)
    ax.set_xlabel('Coefficient')
    ax.set_title('Elastic Net Gene Signature Coefficients')
    ax.axvline(0, color='black', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_18_signature_coefficients.png'), dpi=150)
    plt.close()
    print("  Saved fig_18_signature_coefficients.png")


def plot_rf_importance(rf_top):
    """Bar plot of RF feature importances."""
    top = rf_top.head(30).copy()
    top = top.sort_values('importance')

    fig, ax = plt.subplots(figsize=(8, max(4, len(top) * 0.3)))
    ax.barh(range(len(top)), top['importance'], color='#27ae60')
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top['gene'], fontsize=8)
    ax.set_xlabel('Importance (Mean Decrease in Impurity)')
    ax.set_title('Random Forest Top Feature Importances')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_18b_rf_feature_importance.png'), dpi=150)
    plt.close()
    print("  Saved fig_18b_rf_feature_importance.png")


def run_module5(data, m4_results):
    """Execute Module 5: Predictive gene expression signature."""
    print("\n" + "=" * 60)
    print("MODULE 5: Predictive Gene Expression Signature")
    print("=" * 60)

    log_expr = data['log_expr_filtered']
    cohorts = data['cohorts']
    de_df = m4_results['de_results']

    # Build feature matrix
    print("\n[5.1] Building feature matrix...")
    X, y, patient_ids = build_feature_matrix(log_expr, de_df, cohorts)
    gene_names = list(X.columns)
    print(f"  Samples: {X.shape[0]}, Starting genes: {X.shape[1]}")
    print(f"  Class balance: {y.sum()} responders, {(1-y).sum()} non-responders")

    # 5.2 Elastic net nested CV
    print("\n[5.2] Elastic net logistic regression (nested 5-fold CV)...")
    en_results = nested_cv_elastic_net(X, y, gene_names)
    mean_en_auc = np.mean(en_results['fold_aucs'])
    std_en_auc = np.std(en_results['fold_aucs'])
    print(f"  Mean AUC: {mean_en_auc:.3f} ± {std_en_auc:.3f}")

    # 5.3 Random forest nested CV
    print("\n[5.3] Random forest (nested 5-fold CV)...")
    rf_results = nested_cv_random_forest(X, y, gene_names)
    mean_rf_auc = np.mean(rf_results['fold_aucs'])
    std_rf_auc = np.std(rf_results['fold_aucs'])
    print(f"  Mean AUC: {mean_rf_auc:.3f} ± {std_rf_auc:.3f}")

    # 5.4 Baseline comparators
    print("\n[5.4] Baseline comparators (TMB, PD-L1)...")
    comparators = baseline_comparators(cohorts)
    for name, comp_data in comparators.items():
        print(f"  {name} AUC: {comp_data['auc']:.3f}")

    # 5.5 ROC comparison plot
    print("\n[5.5] ROC comparison plot...")
    plot_roc_comparison(en_results, rf_results, comparators)

    # 5.6 Extract elastic net signature
    print("\n[5.6] Extracting elastic net signature (full dataset refit)...")
    en_coefs, en_model, en_scaler, en_selected, en_sig_genes, en_sig_coefs = \
        extract_elastic_net_signature(X, y, gene_names, de_df=de_df)
    en_coefs.to_csv(os.path.join(RESULTS_DIR, 'table_07_signature_genes_elasticnet.csv'), index=False)
    print(f"  Saved table_07_signature_genes_elasticnet.csv ({len(en_coefs)} genes)")

    # 5.7 Extract RF importances
    print("\n[5.7] Extracting RF feature importances (full dataset refit)...")
    rf_top, rf_model = extract_rf_importance(X, y, gene_names)
    rf_top.to_csv(os.path.join(RESULTS_DIR, 'table_07b_signature_genes_rf.csv'), index=False)
    print(f"  Saved table_07b_signature_genes_rf.csv")

    # Overlap
    en_genes = set(en_coefs['gene'].tolist())
    rf_genes = set(rf_top['gene'].tolist())
    overlap = en_genes & rf_genes
    print(f"  Overlap between EN signature and RF top 30: {len(overlap)} genes")
    if overlap:
        print(f"    {sorted(overlap)}")

    # 5.8 Coefficient / importance plots
    print("\n[5.8] Signature coefficient and importance plots...")
    plot_signature_coefficients(en_coefs)
    plot_rf_importance(rf_top)

    # Model comparison table
    model_comp = {
        'model': ['Elastic Net', 'Random Forest'],
        'mean_auc': [mean_en_auc, mean_rf_auc],
        'std_auc': [std_en_auc, std_rf_auc],
    }
    for name, comp_data in comparators.items():
        model_comp['model'].append(name)
        model_comp['mean_auc'].append(comp_data['auc'])
        model_comp['std_auc'].append(0.0)

    model_comp_df = pd.DataFrame(model_comp)
    model_comp_df.to_csv(os.path.join(RESULTS_DIR, 'table_08_model_comparison.csv'), index=False)
    print(f"\n  Saved table_08_model_comparison.csv")

    # Compute gene signature score for all patients (for M6)
    # Use elastic net: score = sum(coef * z-scored expression)
    common_patients = [p for p in data['clinical'].index if p in log_expr.columns]
    sig_expr = log_expr.loc[log_expr.index.isin(en_sig_genes), common_patients]
    sig_z = sig_expr.subtract(sig_expr.mean(axis=1), axis=0).div(sig_expr.std(axis=1), axis=0)
    gene_sig_score = pd.Series(0.0, index=common_patients)
    for gene, coef in en_sig_coefs.items():
        if gene in sig_z.index:
            gene_sig_score += coef * sig_z.loc[gene]

    print("\n[Module 5 Complete]")
    return {
        'en_results': en_results,
        'rf_results': rf_results,
        'comparators': comparators,
        'en_coefs': en_coefs,
        'rf_importances': rf_top,
        'en_model': en_model,
        'en_scaler': en_scaler,
        'en_selected_genes': en_selected,
        'en_sig_genes': en_sig_genes,
        'en_sig_coefs': en_sig_coefs,
        'gene_sig_score': gene_sig_score,
        'model_comparison': model_comp_df,
    }


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/code')
    from module1_data import run_module1
    from module4_expression import run_module4
    os.makedirs(RESULTS_DIR, exist_ok=True)
    data = run_module1()
    m4 = run_module4(data)
    run_module5(data, m4)
