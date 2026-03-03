"""
Module 6: Integrated Classifier, Survival Model & Summary
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
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import StratifiedKFold, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (roc_curve, auc, roc_auc_score, accuracy_score,
                             precision_score, recall_score, f1_score,
                             confusion_matrix, ConfusionMatrixDisplay)
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test

RESULTS_DIR = '/results'


def build_integrated_features(data, m4_results, m5_results):
    """
    Build integrated feature set: gene signature score,
    TMB, IC Level, Immune phenotype, TCGA subtype, ECOG, Sex.
    Returns X, y, patient_ids for binary analysis.
    """
    binary_df = data['cohorts']['binary'].copy()
    gene_sig_score = m5_results['gene_sig_score']

    # Gene signature score
    binary_df['gene_sig_score'] = gene_sig_score.reindex(binary_df.index)

    # TMB
    tmb_col = 'FMOne mutation burden per MB'
    binary_df['tmb'] = binary_df[tmb_col]
    binary_df['tmb_missing'] = binary_df['tmb'].isna().astype(int)
    binary_df['tmb'] = binary_df['tmb'].fillna(binary_df['tmb'].median())

    # IC Level ordinal
    ic_map = {'IC0': 0, 'IC1': 1, 'IC2+': 2}
    binary_df['ic_ordinal'] = binary_df['IC Level'].map(ic_map)

    # Immune phenotype - dummies
    if 'Immune phenotype' in binary_df.columns:
        dummies = pd.get_dummies(binary_df['Immune phenotype'], prefix='immpheno', drop_first=True)
        binary_df = pd.concat([binary_df, dummies], axis=1)

    # TCGA Subtype - dummies
    if 'TCGA Subtype' in binary_df.columns:
        dummies = pd.get_dummies(binary_df['TCGA Subtype'], prefix='tcga', drop_first=True)
        binary_df = pd.concat([binary_df, dummies], axis=1)

    # ECOG
    binary_df['ecog'] = pd.to_numeric(binary_df['Baseline ECOG Score'], errors='coerce')

    # Sex
    binary_df['sex_male'] = (binary_df['Sex'] == 'M').astype(int)

    # Identify feature columns by category
    expr_features = ['gene_sig_score']
    clinical_features = ['tmb', 'tmb_missing', 'ic_ordinal', 'ecog', 'sex_male'] + \
                        [c for c in binary_df.columns if c.startswith('immpheno_') or c.startswith('tcga_')]

    all_features = expr_features + clinical_features
    y = binary_df['response_binary']

    # Drop rows with any NaN in features
    sub = binary_df[all_features + ['response_binary']].dropna()
    X = sub[all_features]
    y = sub['response_binary'].values

    return X, y, sub.index.tolist(), expr_features, clinical_features


def run_cv_model(X, y, model_class, param_grid=None, model_name='Model', n_outer=5):
    """
    Run 5-fold stratified CV for a classifier.
    Returns per-fold metrics and pooled predictions.
    """
    outer_cv = StratifiedKFold(n_splits=n_outer, shuffle=True, random_state=42)
    scaler = StandardScaler()

    fold_metrics = []
    all_y_true = []
    all_y_prob = []
    all_y_pred = []
    fold_fprs = []
    fold_tprs = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        X_train_sc = scaler.fit_transform(X_train)
        X_test_sc = scaler.transform(X_test)

        if param_grid is not None:
            inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=fold_idx)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                gs = GridSearchCV(model_class, param_grid, cv=inner_cv,
                                  scoring='roc_auc', n_jobs=1)
                gs.fit(X_train_sc, y_train)
                model = gs.best_estimator_
        else:
            model = model_class
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                model.fit(X_train_sc, y_train)

        probs = model.predict_proba(X_test_sc)[:, 1]
        preds = model.predict(X_test_sc)

        try:
            fold_auc = roc_auc_score(y_test, probs)
        except ValueError:
            fold_auc = 0.5

        fpr, tpr, _ = roc_curve(y_test, probs)
        fold_fprs.append(fpr)
        fold_tprs.append(tpr)

        fold_metrics.append({
            'fold': fold_idx + 1,
            'auc': fold_auc,
            'accuracy': accuracy_score(y_test, preds),
            'precision': precision_score(y_test, preds, zero_division=0),
            'recall': recall_score(y_test, preds, zero_division=0),
            'f1': f1_score(y_test, preds, zero_division=0),
        })

        all_y_true.extend(y_test)
        all_y_prob.extend(probs)
        all_y_pred.extend(preds)

    return {
        'fold_metrics': pd.DataFrame(fold_metrics),
        'all_y_true': np.array(all_y_true),
        'all_y_prob': np.array(all_y_prob),
        'all_y_pred': np.array(all_y_pred),
        'fold_fprs': fold_fprs,
        'fold_tprs': fold_tprs,
        'model_name': model_name,
    }


def run_integrated_classifiers(X, y, expr_features, clinical_features):
    """
    Run multiple classifiers: expression-only, clinical-only, integrated LR, integrated GB.
    """
    results = {}

    # Expression-only logistic regression
    print("  [a] Expression-only logistic regression...")
    X_expr = X[expr_features]
    results['Expression-only LR'] = run_cv_model(
        X_expr, y, LogisticRegression(C=1.0, max_iter=1000, random_state=42),
        model_name='Expression-only LR'
    )

    # Clinical-only logistic regression
    print("  [b] Clinical-only logistic regression...")
    X_clin = X[clinical_features]
    results['Clinical-only LR'] = run_cv_model(
        X_clin, y, LogisticRegression(C=1.0, max_iter=1000, random_state=42),
        model_name='Clinical-only LR'
    )

    # Integrated logistic regression
    print("  [c] Integrated logistic regression...")
    results['Integrated LR'] = run_cv_model(
        X, y, LogisticRegression(C=1.0, max_iter=1000, random_state=42),
        model_name='Integrated LR'
    )

    # Integrated gradient boosting
    print("  [d] Integrated gradient boosting (with tuning)...")
    gb_param_grid = {
        'n_estimators': [100, 200],
        'max_depth': [3, 5],
        'learning_rate': [0.05, 0.1],
    }
    results['Integrated GB'] = run_cv_model(
        X, y, GradientBoostingClassifier(random_state=42),
        param_grid=gb_param_grid, model_name='Integrated GB'
    )

    return results


def plot_classifier_rocs(all_results):
    """ROC curves for all integrated model comparisons."""
    fig, ax = plt.subplots(figsize=(8, 7))
    colors = {
        'Expression-only LR': '#e74c3c',
        'Clinical-only LR': '#3498db',
        'Integrated LR': '#2ecc71',
        'Integrated GB': '#9b59b6',
    }

    mean_fpr = np.linspace(0, 1, 100)
    for name, res in all_results.items():
        tprs = []
        for fpr, tpr in zip(res['fold_fprs'], res['fold_tprs']):
            tprs.append(np.interp(mean_fpr, fpr, tpr))
        mean_tpr = np.mean(tprs, axis=0)
        mean_auc = res['fold_metrics']['auc'].mean()
        std_auc = res['fold_metrics']['auc'].std()

        ax.plot(mean_fpr, mean_tpr, color=colors.get(name, 'grey'), lw=2,
                label=f'{name} (AUC={mean_auc:.3f}±{std_auc:.3f})')

    ax.plot([0, 1], [0, 1], 'k--', alpha=0.3)
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.set_title('Classifier Comparison: ROC Curves')
    ax.legend(loc='lower right', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_19_classifier_roc.png'), dpi=150)
    plt.close()
    print("  Saved fig_19_classifier_roc.png")


def plot_confusion_matrix(y_true, y_pred, model_name='Integrated GB'):
    """Confusion matrix for pooled out-of-fold predictions."""
    cm = confusion_matrix(y_true, y_pred)
    fig, ax = plt.subplots(figsize=(6, 5))
    disp = ConfusionMatrixDisplay(cm, display_labels=['SD/PD', 'CR/PR'])
    disp.plot(ax=ax, cmap='Blues', values_format='d')
    ax.set_title(f'Confusion Matrix ({model_name})')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_19b_confusion_matrix.png'), dpi=150)
    plt.close()
    print("  Saved fig_19b_confusion_matrix.png")


def integrated_survival(data, m4_results, m5_results):
    """
    Cox PH model with integrated features.
    Compare C-index: expression alone vs clinical alone vs integrated.
    KM by integrated risk score.
    """
    survival_df = data['cohorts']['survival'].copy()
    gene_sig_score = m5_results['gene_sig_score']

    # Add features
    survival_df['gene_sig_score'] = gene_sig_score.reindex(survival_df.index)

    tmb_col = 'FMOne mutation burden per MB'
    survival_df['tmb'] = survival_df[tmb_col].fillna(survival_df[tmb_col].median())
    survival_df['tmb_missing'] = survival_df[tmb_col].isna().astype(int)

    ic_map = {'IC0': 0, 'IC1': 1, 'IC2+': 2}
    survival_df['ic_ordinal'] = survival_df['IC Level'].map(ic_map)

    if 'Immune phenotype' in survival_df.columns:
        dummies = pd.get_dummies(survival_df['Immune phenotype'], prefix='immpheno', drop_first=True)
        survival_df = pd.concat([survival_df, dummies], axis=1)

    if 'TCGA Subtype' in survival_df.columns:
        dummies = pd.get_dummies(survival_df['TCGA Subtype'], prefix='tcga', drop_first=True)
        survival_df = pd.concat([survival_df, dummies], axis=1)

    survival_df['ecog'] = pd.to_numeric(survival_df['Baseline ECOG Score'], errors='coerce')
    survival_df['sex_male'] = (survival_df['Sex'] == 'M').astype(int)

    expr_features = ['gene_sig_score']
    clinical_features = ['tmb', 'tmb_missing', 'ic_ordinal', 'ecog', 'sex_male'] + \
                        [c for c in survival_df.columns if c.startswith('immpheno_') or c.startswith('tcga_')]
    all_features = expr_features + clinical_features

    c_indices = {}

    # Expression-only Cox
    expr_df = survival_df[['os', 'censOS'] + expr_features].dropna()
    try:
        cph_expr = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cph_expr.fit(expr_df, duration_col='os', event_col='censOS')
        c_indices['Expression-only'] = round(cph_expr.concordance_index_, 4)
    except Exception as e:
        print(f"  Expression-only Cox failed: {e}")
        c_indices['Expression-only'] = np.nan

    # Clinical-only Cox
    clin_df = survival_df[['os', 'censOS'] + clinical_features].dropna()
    try:
        cph_clin = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cph_clin.fit(clin_df, duration_col='os', event_col='censOS')
        c_indices['Clinical-only'] = round(cph_clin.concordance_index_, 4)
    except Exception as e:
        print(f"  Clinical-only Cox failed: {e}")
        c_indices['Clinical-only'] = np.nan

    # Integrated Cox
    full_df = survival_df[['os', 'censOS'] + all_features].dropna()
    try:
        cph_full = CoxPHFitter(penalizer=0.1)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            cph_full.fit(full_df, duration_col='os', event_col='censOS')
        c_indices['Integrated'] = round(cph_full.concordance_index_, 4)

        # Risk score for KM plot
        risk_scores = cph_full.predict_partial_hazard(full_df)
        median_risk = risk_scores.median()
        full_df['risk_group'] = np.where(risk_scores >= median_risk, 'High Risk', 'Low Risk')

        # KM by risk group
        fig, ax = plt.subplots(figsize=(8, 6))
        kmf = KaplanMeierFitter()
        for group, color in [('Low Risk', '#2ecc71'), ('High Risk', '#e74c3c')]:
            mask = full_df['risk_group'] == group
            kmf.fit(full_df.loc[mask, 'os'], event_observed=full_df.loc[mask, 'censOS'], label=group)
            kmf.plot_survival_function(ax=ax, color=color)

        m1 = full_df['risk_group'] == 'High Risk'
        m2 = full_df['risk_group'] == 'Low Risk'
        lr = logrank_test(full_df.loc[m1, 'os'], full_df.loc[m2, 'os'],
                          full_df.loc[m1, 'censOS'], full_df.loc[m2, 'censOS'])
        ax.text(0.6, 0.95, f'Log-rank p = {lr.p_value:.4f}', transform=ax.transAxes, fontsize=10)
        ax.set_title('Overall Survival by Integrated Risk Score')
        ax.set_xlabel('Time')
        ax.set_ylabel('Survival Probability')
        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, 'fig_19c_km_integrated_risk.png'), dpi=150)
        plt.close()
        print("  Saved fig_19c_km_integrated_risk.png")

    except Exception as e:
        print(f"  Integrated Cox failed: {e}")
        c_indices['Integrated'] = np.nan
        full_df = pd.DataFrame()

    surv_comp = pd.DataFrame([
        {'model': k, 'c_index': v} for k, v in c_indices.items()
    ])
    surv_comp.to_csv(os.path.join(RESULTS_DIR, 'table_09b_integrated_survival.csv'), index=False)
    print(f"  Saved table_09b_integrated_survival.csv")

    return c_indices, full_df


def summary_figure(data, all_classifier_results, c_indices):
    """Multi-panel summary figure with 4 key panels."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))

    # Panel A: Model AUC comparison (bar chart)
    ax = axes[0, 0]
    models = []
    aucs = []
    for name, res in all_classifier_results.items():
        models.append(name)
        aucs.append(res['fold_metrics']['auc'].mean())
    colors = ['#e74c3c', '#3498db', '#2ecc71', '#9b59b6'][:len(models)]
    ax.bar(range(len(models)), aucs, color=colors)
    ax.set_xticks(range(len(models)))
    ax.set_xticklabels(models, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Mean AUC')
    ax.set_title('A. Response Prediction: Model Comparison')
    ax.set_ylim(0.4, 0.9)
    for i, v in enumerate(aucs):
        ax.text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=9)

    # Panel B: Survival C-index comparison
    ax = axes[0, 1]
    names = list(c_indices.keys())
    vals = [c_indices[n] for n in names]
    valid_mask = [not (isinstance(v, float) and np.isnan(v)) for v in vals]
    valid_names = [n for n, m in zip(names, valid_mask) if m]
    valid_vals = [v for v, m in zip(vals, valid_mask) if m]
    ax.bar(range(len(valid_names)), valid_vals, color=['#e74c3c', '#3498db', '#2ecc71'][:len(valid_names)])
    ax.set_xticks(range(len(valid_names)))
    ax.set_xticklabels(valid_names, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('C-index')
    ax.set_title('B. Survival Prediction: C-index Comparison')
    ax.set_ylim(0.4, 0.8)
    for i, v in enumerate(valid_vals):
        ax.text(i, v + 0.01, f'{v:.3f}', ha='center', fontsize=9)

    # Panel C: Response rate by IC level
    ax = axes[1, 0]
    binary_df = data['cohorts']['binary']
    if 'IC Level' in binary_df.columns:
        grouped = binary_df.groupby('IC Level')['response_binary'].agg(['mean', 'count'])
        grouped = grouped.sort_index()
        ax.bar(range(len(grouped)), grouped['mean'] * 100,
               color=['#3498db', '#f39c12', '#e74c3c'][:len(grouped)])
        ax.set_xticks(range(len(grouped)))
        labels = [f"{idx}\n(n={int(row['count'])})" for idx, row in grouped.iterrows()]
        ax.set_xticklabels(labels, fontsize=9)
        ax.set_ylabel('Response Rate (%)')
        ax.set_title('C. Response Rate by PD-L1 IC Level')
        for i, (_, row) in enumerate(grouped.iterrows()):
            ax.text(i, row['mean'] * 100 + 1, f'{row["mean"]*100:.1f}%', ha='center', fontsize=9)

    # Panel D: Key biomarker effect sizes
    ax = axes[1, 1]
    # Show top significant associations as a summary
    ax.text(0.5, 0.9, 'Key Findings Summary', ha='center', va='top',
            transform=ax.transAxes, fontsize=12, fontweight='bold')

    best_auc = max(aucs)
    best_model = models[aucs.index(best_auc)]
    findings = [
        f'Best classifier: {best_model} (AUC={best_auc:.3f})',
    ]
    if valid_vals:
        best_ci = max(valid_vals)
        best_surv = valid_names[valid_vals.index(best_ci)]
        findings.append(f'Best survival model: {best_surv} (C-index={best_ci:.3f})')

    for i, finding in enumerate(findings):
        ax.text(0.1, 0.7 - i * 0.15, f'• {finding}', transform=ax.transAxes,
                fontsize=10, va='top')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis('off')
    ax.set_title('D. Summary')

    plt.suptitle('IMvigor210 Biomarker Analysis Summary', fontsize=14, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_20_summary_figure.png'), dpi=150)
    plt.close()
    print("  Saved fig_20_summary_figure.png")


def generate_tableau_outputs(data, m4_results, m5_results, classifier_results, survival_full_df):
    """Generate all Tableau-ready long-format CSV files."""
    print("\n[6.6] Generating Tableau-ready outputs...")
    clinical = data['clinical']
    log_expr = data['log_expr_filtered']
    genomic = data['genomic']
    cohorts = data['cohorts']
    binary_df = cohorts['binary']
    gene_sig_score = m5_results['gene_sig_score']
    de_df = m4_results['de_results']
    mut_matrix = data['mut_matrix']

    # --- 1. tableau_patient_features.csv ---
    rows = []
    for pid in clinical.index:
        base = {
            'patient_id': pid,
            'response': clinical.loc[pid, 'Best Confirmed Overall Response'] if 'Best Confirmed Overall Response' in clinical.columns else None,
            'binary_response': clinical.loc[pid, 'binaryResponse'] if 'binaryResponse' in clinical.columns else None,
            'os': clinical.loc[pid, 'os'] if 'os' in clinical.columns else None,
            'censOS': clinical.loc[pid, 'censOS'] if 'censOS' in clinical.columns else None,
        }
        # Clinical variables
        for var in ['Sex', 'IC Level', 'TC Level', 'Immune phenotype', 'TCGA Subtype',
                     'Baseline ECOG Score', 'Received platinum', 'Met Disease Status',
                     'Tobacco Use History', 'Lund']:
            if var in clinical.columns:
                rows.append({**base, 'variable_name': var,
                             'variable_value': str(clinical.loc[pid, var]),
                             'variable_type': 'clinical'})

        # TMB, neoantigen
        for var in ['FMOne mutation burden per MB', 'Neoantigen burden per MB']:
            if var in clinical.columns:
                val = clinical.loc[pid, var]
                rows.append({**base, 'variable_name': var,
                             'variable_value': str(val) if pd.notna(val) else '',
                             'variable_type': 'genomic'})

        # Gene signature score
        if pid in gene_sig_score.index:
            rows.append({**base, 'variable_name': 'gene_signature_score',
                         'variable_value': str(round(gene_sig_score[pid], 4)),
                         'variable_type': 'score'})

    tableau_patients = pd.DataFrame(rows)
    tableau_patients.to_csv(os.path.join(RESULTS_DIR, 'tableau_patient_features.csv'), index=False)
    print(f"  Saved tableau_patient_features.csv ({len(tableau_patients)} rows)")

    # --- 2. tableau_gene_expression.csv ---
    top_de_genes = de_df.head(100)['gene'].tolist()
    expr_rows = []
    for gene in top_de_genes:
        if gene not in log_expr.index:
            continue
        de_row = de_df[de_df['gene'] == gene].iloc[0]
        for pid in log_expr.columns:
            br = clinical.loc[pid, 'binaryResponse'] if pid in clinical.index and 'binaryResponse' in clinical.columns else None
            resp = clinical.loc[pid, 'Best Confirmed Overall Response'] if pid in clinical.index and 'Best Confirmed Overall Response' in clinical.columns else None
            expr_rows.append({
                'patient_id': pid,
                'gene_symbol': gene,
                'log2_expression': round(log_expr.loc[gene, pid], 4),
                'response': resp,
                'binary_response': br,
                'de_log2fc': round(de_row['log2FC'], 4),
                'de_fdr': round(de_row['fdr'], 6),
            })
    tableau_expr = pd.DataFrame(expr_rows)
    tableau_expr.to_csv(os.path.join(RESULTS_DIR, 'tableau_gene_expression.csv'), index=False)
    print(f"  Saved tableau_gene_expression.csv ({len(tableau_expr)} rows)")

    # --- 3. tableau_mutations.csv ---
    mut_rows = []
    for _, row in genomic.iterrows():
        pid = row['patient_id']
        br = clinical.loc[pid, 'binaryResponse'] if pid in clinical.index and 'binaryResponse' in clinical.columns else None
        resp = clinical.loc[pid, 'Best Confirmed Overall Response'] if pid in clinical.index and 'Best Confirmed Overall Response' in clinical.columns else None
        mut_rows.append({
            'patient_id': pid,
            'gene': row['gene'],
            'alteration_type': row['alteration_type'],
            'is_mutated': 1,
            'response': resp,
            'binary_response': br,
        })
    tableau_mut = pd.DataFrame(mut_rows)
    tableau_mut.to_csv(os.path.join(RESULTS_DIR, 'tableau_mutations.csv'), index=False)
    print(f"  Saved tableau_mutations.csv ({len(tableau_mut)} rows)")

    # --- 4. tableau_model_performance.csv ---
    perf_rows = []
    for model_name, res in classifier_results.items():
        for _, frow in res['fold_metrics'].iterrows():
            for metric in ['auc', 'accuracy', 'precision', 'recall', 'f1']:
                perf_rows.append({
                    'model_name': model_name,
                    'fold': int(frow['fold']),
                    'metric': metric,
                    'value': round(frow[metric], 4),
                })
    tableau_perf = pd.DataFrame(perf_rows)
    tableau_perf.to_csv(os.path.join(RESULTS_DIR, 'tableau_model_performance.csv'), index=False)
    print(f"  Saved tableau_model_performance.csv ({len(tableau_perf)} rows)")

    # --- 5. tableau_survival.csv ---
    surv_rows = []
    for pid in clinical.index:
        base = {
            'patient_id': pid,
            'os': clinical.loc[pid, 'os'] if 'os' in clinical.columns else None,
            'censOS': clinical.loc[pid, 'censOS'] if 'censOS' in clinical.columns else None,
            'response': clinical.loc[pid, 'Best Confirmed Overall Response'] if 'Best Confirmed Overall Response' in clinical.columns else None,
        }
        # IC Level stratification
        if 'IC Level' in clinical.columns:
            surv_rows.append({**base, 'stratification': 'IC_Level',
                              'stratum_value': str(clinical.loc[pid, 'IC Level'])})
        # Immune phenotype
        if 'Immune phenotype' in clinical.columns:
            surv_rows.append({**base, 'stratification': 'Immune_phenotype',
                              'stratum_value': str(clinical.loc[pid, 'Immune phenotype'])})
        # TMB group
        tmb_col = 'FMOne mutation burden per MB'
        if tmb_col in clinical.columns and pd.notna(clinical.loc[pid, tmb_col]):
            median_tmb = clinical[tmb_col].median()
            group = 'High' if clinical.loc[pid, tmb_col] >= median_tmb else 'Low'
            surv_rows.append({**base, 'stratification': 'TMB_group', 'stratum_value': group})

        # Integrated risk score group (if available)
        if len(survival_full_df) > 0 and 'risk_group' in survival_full_df.columns and pid in survival_full_df.index:
            surv_rows.append({**base, 'stratification': 'Risk_score_group',
                              'stratum_value': survival_full_df.loc[pid, 'risk_group']})

    tableau_surv = pd.DataFrame(surv_rows)
    tableau_surv.to_csv(os.path.join(RESULTS_DIR, 'tableau_survival.csv'), index=False)
    print(f"  Saved tableau_survival.csv ({len(tableau_surv)} rows)")


def analysis_summary_table(m2_results, m3_results, m4_results, m5_results,
                           classifier_results, c_indices):
    """Generate a key findings summary table."""
    findings = []

    # Module 2 findings
    if m2_results and 'univariate' in m2_results:
        uni = m2_results['univariate']
        sig_vars = uni[uni['fdr'] < 0.05]['variable'].tolist() if 'fdr' in uni.columns else []
        findings.append({
            'module': 'Clinical',
            'finding': f'{len(sig_vars)} clinical variables significantly associated with response (FDR<0.05)',
            'detail': ', '.join(sig_vars[:5]) if sig_vars else 'None',
        })
    if m2_results and 'cox_c_index' in m2_results:
        findings.append({
            'module': 'Clinical',
            'finding': f'Cox PH model C-index',
            'detail': f'{m2_results["cox_c_index"]:.4f}',
        })

    # Module 3 findings
    if m3_results:
        findings.append({
            'module': 'Genomic',
            'finding': 'TMB as response predictor',
            'detail': f'AUC={m3_results.get("tmb_auc", "N/A"):.3f}',
        })

    # Module 4 findings
    if m4_results and 'de_results' in m4_results:
        de = m4_results['de_results']
        n_sig = ((de['fdr'] < 0.05) & (de['log2FC'].abs() > 1)).sum()
        findings.append({
            'module': 'Expression',
            'finding': f'{n_sig} differentially expressed genes (|log2FC|>1, FDR<0.05)',
            'detail': '',
        })

    # Module 5 findings
    if m5_results and 'model_comparison' in m5_results:
        mc = m5_results['model_comparison']
        for _, row in mc.iterrows():
            findings.append({
                'module': 'Signature',
                'finding': f'{row["model"]} AUC',
                'detail': f'{row["mean_auc"]:.3f}±{row["std_auc"]:.3f}',
            })

    # Module 6 findings
    for name, res in classifier_results.items():
        mean_auc = res['fold_metrics']['auc'].mean()
        findings.append({
            'module': 'Integrated',
            'finding': f'{name} AUC',
            'detail': f'{mean_auc:.3f}',
        })

    for name, ci in c_indices.items():
        if not (isinstance(ci, float) and np.isnan(ci)):
            findings.append({
                'module': 'Integrated',
                'finding': f'{name} survival C-index',
                'detail': f'{ci:.4f}',
            })

    summary_df = pd.DataFrame(findings)
    summary_df.to_csv(os.path.join(RESULTS_DIR, 'table_10_analysis_summary.csv'), index=False)
    print(f"  Saved table_10_analysis_summary.csv")
    return summary_df


def run_module6(data, m2_results, m3_results, m4_results, m5_results):
    """Execute Module 6: Integrated classifier, survival model & summary."""
    print("\n" + "=" * 60)
    print("MODULE 6: Integrated Classifier, Survival Model & Summary")
    print("=" * 60)

    # 6.1 Build integrated features
    print("\n[6.1] Building integrated feature set...")
    X, y, patient_ids, expr_features, clinical_features = \
        build_integrated_features(data, m4_results, m5_results)
    print(f"  Features: {X.shape[1]} ({len(expr_features)} expression, {len(clinical_features)} clinical)")
    print(f"  Samples: {X.shape[0]} ({y.sum()} R, {(1-y).sum()} NR)")

    # 6.2 Run integrated classifiers
    print("\n[6.2] Running integrated classifiers (5-fold CV)...")
    classifier_results = run_integrated_classifiers(X, y, expr_features, clinical_features)

    # Print summary
    print("\n  Classifier performance summary:")
    perf_rows = []
    for name, res in classifier_results.items():
        fm = res['fold_metrics']
        perf_rows.append({
            'model': name,
            'mean_auc': fm['auc'].mean(),
            'std_auc': fm['auc'].std(),
            'mean_accuracy': fm['accuracy'].mean(),
            'mean_f1': fm['f1'].mean(),
        })
        print(f"    {name}: AUC={fm['auc'].mean():.3f}±{fm['auc'].std():.3f}, "
              f"Acc={fm['accuracy'].mean():.3f}, F1={fm['f1'].mean():.3f}")

    perf_df = pd.DataFrame(perf_rows)
    perf_df.to_csv(os.path.join(RESULTS_DIR, 'table_09_classifier_performance.csv'), index=False)
    print(f"  Saved table_09_classifier_performance.csv")

    # 6.3 ROC plot
    print("\n[6.3] Classifier ROC comparison plot...")
    plot_classifier_rocs(classifier_results)

    # 6.4 Confusion matrix (using best model: Integrated GB)
    print("\n[6.4] Confusion matrix...")
    best_model = 'Integrated GB'
    if best_model in classifier_results:
        res = classifier_results[best_model]
        plot_confusion_matrix(res['all_y_true'], res['all_y_pred'], best_model)

    # 6.5 Integrated survival model
    print("\n[6.5] Integrated survival model...")
    c_indices, survival_full_df = integrated_survival(data, m4_results, m5_results)
    for name, ci in c_indices.items():
        print(f"    {name} C-index: {ci}")

    # 6.6 Tableau outputs
    generate_tableau_outputs(data, m4_results, m5_results, classifier_results, survival_full_df)

    # 6.7 Summary figure
    print("\n[6.7] Summary figure...")
    summary_figure(data, classifier_results, c_indices)

    # 6.8 Analysis summary table
    print("\n[6.8] Analysis summary table...")
    analysis_summary_table(m2_results, m3_results, m4_results, m5_results,
                           classifier_results, c_indices)

    print("\n[Module 6 Complete]")
    return {
        'classifier_results': classifier_results,
        'c_indices': c_indices,
    }


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/code')
    from module1_data import run_module1
    from module2_clinical import run_module2
    from module3_genomic import run_module3
    from module4_expression import run_module4
    from module5_signature import run_module5
    os.makedirs(RESULTS_DIR, exist_ok=True)
    data = run_module1()
    m2 = run_module2(data)
    m3 = run_module3(data)
    m4 = run_module4(data)
    m5 = run_module5(data, m4)
    run_module6(data, m2, m3, m4, m5)
