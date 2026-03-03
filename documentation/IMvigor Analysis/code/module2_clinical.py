"""
Module 2: Clinical Feature Association Analysis
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
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

RESULTS_DIR = '/results'


def univariate_association(binary_df):
    """
    Test association of each clinical variable with binary response.
    Returns DataFrame with variable, test, statistic, p-value.
    """
    results = []

    # Categorical variables and their column names
    categorical_vars = [
        'IC Level', 'TC Level', 'Immune phenotype', 'Sex',
        'Baseline ECOG Score', 'TCGA Subtype', 'Lund',
        'Met Disease Status', 'Tobacco Use History',
        'Received platinum', 'Intravesical BCG administered',
    ]

    # Continuous variables
    continuous_vars = [
        'FMOne mutation burden per MB',
        'Neoantigen burden per MB',
        'sizeFactor',
    ]

    y = binary_df['response_binary']

    for var in categorical_vars:
        if var not in binary_df.columns:
            continue
        sub = binary_df[[var, 'response_binary']].dropna()
        if sub[var].nunique() < 2:
            continue
        ct = pd.crosstab(sub[var], sub['response_binary'])
        if ct.shape[0] < 2 or ct.shape[1] < 2:
            continue
        # Use chi-square for larger tables, Fisher's for 2x2
        if ct.shape == (2, 2):
            stat, p = stats.fisher_exact(ct)
            test_name = "Fisher's exact"
        else:
            stat, p, _, _ = stats.chi2_contingency(ct)
            test_name = "Chi-square"
        results.append({
            'variable': var, 'test': test_name,
            'statistic': round(stat, 4), 'p_value': p
        })

    for var in continuous_vars:
        if var not in binary_df.columns:
            continue
        sub = binary_df[[var, 'response_binary']].dropna()
        if len(sub) < 10:
            continue
        resp = sub.loc[sub['response_binary'] == 1, var]
        nonresp = sub.loc[sub['response_binary'] == 0, var]
        stat, p = stats.mannwhitneyu(resp, nonresp, alternative='two-sided')
        results.append({
            'variable': var, 'test': 'Mann-Whitney U',
            'statistic': round(stat, 4), 'p_value': p
        })

    res_df = pd.DataFrame(results)
    if len(res_df) > 0:
        _, fdr, _, _ = multipletests(res_df['p_value'], method='fdr_bh')
        res_df['fdr'] = fdr
    return res_df


def response_rate_barplots(binary_df):
    """Generate response rate bar plots for key stratifiers."""
    plot_vars = [
        ('IC Level', 'fig_03_response_by_ic_level.png'),
        ('Immune phenotype', 'fig_04_response_by_immune_phenotype.png'),
        ('TCGA Subtype', 'fig_05_response_by_tcga_subtype.png'),
    ]

    for var, fname in plot_vars:
        if var not in binary_df.columns:
            continue
        sub = binary_df[[var, 'response_binary']].dropna()
        if sub[var].nunique() < 2:
            continue

        # Compute response rate per group
        grouped = sub.groupby(var)['response_binary'].agg(['mean', 'count'])
        grouped.columns = ['response_rate', 'n']
        grouped = grouped.sort_values('response_rate', ascending=False)

        fig, ax = plt.subplots(figsize=(8, 5))
        bars = ax.bar(range(len(grouped)), grouped['response_rate'] * 100,
                       color=sns.color_palette('Set2', len(grouped)))
        ax.set_xticks(range(len(grouped)))
        labels = [f"{idx}\n(n={int(row['n'])})" for idx, row in grouped.iterrows()]
        ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
        ax.set_ylabel('Response Rate (%)')
        ax.set_title(f'Response Rate by {var}')
        ax.set_ylim(0, max(grouped['response_rate'] * 100) * 1.2 + 5)

        for bar, rate in zip(bars, grouped['response_rate']):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f'{rate*100:.1f}%', ha='center', va='bottom', fontsize=9)

        plt.tight_layout()
        plt.savefig(os.path.join(RESULTS_DIR, fname), dpi=150)
        plt.close()
        print(f"  Saved {fname}")


def tmb_boxplot(binary_df):
    """TMB distribution by response status."""
    sub = binary_df[['FMOne mutation burden per MB', 'binaryResponse']].dropna()
    fig, ax = plt.subplots(figsize=(6, 5))
    sns.boxplot(data=sub, x='binaryResponse', y='FMOne mutation burden per MB',
                order=['CR/PR', 'SD/PD'], palette={'CR/PR': '#e74c3c', 'SD/PD': '#3498db'},
                ax=ax, fliersize=3)
    sns.stripplot(data=sub, x='binaryResponse', y='FMOne mutation burden per MB',
                  order=['CR/PR', 'SD/PD'], color='black', alpha=0.3, size=3, ax=ax)
    ax.set_xlabel('Response')
    ax.set_ylabel('Tumor Mutation Burden (per MB)')
    ax.set_title('TMB by Response Status')

    # Add p-value
    resp_tmb = sub.loc[sub['binaryResponse'] == 'CR/PR', 'FMOne mutation burden per MB']
    nonresp_tmb = sub.loc[sub['binaryResponse'] == 'SD/PD', 'FMOne mutation burden per MB']
    _, p = stats.mannwhitneyu(resp_tmb, nonresp_tmb, alternative='two-sided')
    ax.text(0.5, 0.95, f'p = {p:.4f}', transform=ax.transAxes, ha='center', fontsize=10)

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_06_tmb_by_response.png'), dpi=150)
    plt.close()
    print("  Saved fig_06_tmb_by_response.png")


def multivariate_logistic(binary_df):
    """
    Multivariate logistic regression for response prediction.
    Returns DataFrame of coefficients, OR, 95% CI, p-value.
    """
    import statsmodels.api as sm

    df = binary_df.copy()

    # Prepare features
    features = {}

    # IC Level - ordinal encode
    ic_map = {'IC0': 0, 'IC1': 1, 'IC2+': 2}
    if 'IC Level' in df.columns:
        df['ic_ordinal'] = df['IC Level'].map(ic_map)
        features['IC Level'] = 'ic_ordinal'

    # Immune phenotype - dummy encode
    if 'Immune phenotype' in df.columns:
        dummies = pd.get_dummies(df['Immune phenotype'], prefix='immpheno', drop_first=True)
        df = pd.concat([df, dummies], axis=1)
        for col in dummies.columns:
            features[col] = col

    # TMB - median impute + missingness indicator
    if 'FMOne mutation burden per MB' in df.columns:
        tmb_col = 'FMOne mutation burden per MB'
        df['tmb_missing'] = df[tmb_col].isna().astype(int)
        df['tmb_imputed'] = df[tmb_col].fillna(df[tmb_col].median())
        features['TMB'] = 'tmb_imputed'
        features['TMB_missing'] = 'tmb_missing'

    # ECOG
    if 'Baseline ECOG Score' in df.columns:
        df['ecog'] = pd.to_numeric(df['Baseline ECOG Score'], errors='coerce')
        features['ECOG'] = 'ecog'

    # TCGA Subtype - dummy encode
    if 'TCGA Subtype' in df.columns:
        dummies = pd.get_dummies(df['TCGA Subtype'], prefix='tcga', drop_first=True)
        df = pd.concat([df, dummies], axis=1)
        for col in dummies.columns:
            features[col] = col

    # Sex
    if 'Sex' in df.columns:
        df['sex_male'] = (df['Sex'] == 'M').astype(int)
        features['Sex (Male)'] = 'sex_male'

    # Build design matrix
    feature_cols = list(features.values())
    X = df[feature_cols].dropna()
    X = X.astype(float)
    y = df.loc[X.index, 'response_binary'].astype(float)

    X = sm.add_constant(X)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        model = sm.Logit(y, X).fit(disp=0, maxiter=100)

    # Extract results
    results = []
    for i, name in enumerate(model.params.index):
        if name == 'const':
            display_name = 'Intercept'
        else:
            # Map back to readable name
            display_name = name
            for readable, col_name in features.items():
                if col_name == name:
                    display_name = readable
                    break
        coef = model.params.iloc[i]
        ci = model.conf_int().iloc[i]
        results.append({
            'variable': display_name,
            'coefficient': round(coef, 4),
            'odds_ratio': round(np.exp(coef), 4),
            'or_ci_lower': round(np.exp(ci[0]), 4),
            'or_ci_upper': round(np.exp(ci[1]), 4),
            'p_value': model.pvalues.iloc[i],
        })

    return pd.DataFrame(results), model


def forest_plot(logistic_df):
    """Forest plot of adjusted odds ratios from logistic regression."""
    # Exclude intercept
    plot_df = logistic_df[logistic_df['variable'] != 'Intercept'].copy()
    plot_df = plot_df.sort_values('odds_ratio', ascending=True).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(8, max(4, len(plot_df) * 0.5)))
    y_pos = range(len(plot_df))

    ax.errorbar(
        plot_df['odds_ratio'], y_pos,
        xerr=[plot_df['odds_ratio'] - plot_df['or_ci_lower'],
              plot_df['or_ci_upper'] - plot_df['odds_ratio']],
        fmt='o', color='darkblue', ecolor='steelblue', elinewidth=2,
        capsize=4, markersize=6
    )
    ax.axvline(x=1, color='red', linestyle='--', alpha=0.7)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(plot_df['variable'], fontsize=9)
    ax.set_xlabel('Odds Ratio (95% CI)')
    ax.set_title('Multivariate Logistic Regression - Forest Plot')
    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_07_forest_plot.png'), dpi=150)
    plt.close()
    print("  Saved fig_07_forest_plot.png")


def km_survival_analysis(survival_df):
    """
    Kaplan-Meier survival curves stratified by IC Level, Immune phenotype, TMB.
    Also runs Cox PH regression.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    stratifiers = [
        ('IC Level', axes[0]),
        ('Immune phenotype', axes[1]),
    ]

    for var, ax in stratifiers:
        if var not in survival_df.columns:
            continue
        sub = survival_df[['os', 'censOS', var]].dropna()
        groups = sub[var].unique()
        kmf = KaplanMeierFitter()

        for group in sorted(groups):
            mask = sub[var] == group
            kmf.fit(sub.loc[mask, 'os'], event_observed=sub.loc[mask, 'censOS'], label=str(group))
            kmf.plot_survival_function(ax=ax)

        # Log-rank test (pairwise between all groups)
        if len(groups) == 2:
            g1, g2 = groups
            m1, m2 = sub[var] == g1, sub[var] == g2
            lr = logrank_test(sub.loc[m1, 'os'], sub.loc[m2, 'os'],
                              sub.loc[m1, 'censOS'], sub.loc[m2, 'censOS'])
            ax.text(0.6, 0.95, f'p = {lr.p_value:.4f}', transform=ax.transAxes, fontsize=9)

        ax.set_title(f'OS by {var}')
        ax.set_xlabel('Time')
        ax.set_ylabel('Survival Probability')

    # TMB: median split
    tmb_col = 'FMOne mutation burden per MB'
    if tmb_col in survival_df.columns:
        sub = survival_df[['os', 'censOS', tmb_col]].dropna()
        median_tmb = sub[tmb_col].median()
        sub['TMB_group'] = np.where(sub[tmb_col] >= median_tmb, 'High TMB', 'Low TMB')
        kmf = KaplanMeierFitter()
        for group in ['High TMB', 'Low TMB']:
            mask = sub['TMB_group'] == group
            kmf.fit(sub.loc[mask, 'os'], event_observed=sub.loc[mask, 'censOS'], label=group)
            kmf.plot_survival_function(ax=axes[2])

        m1 = sub['TMB_group'] == 'High TMB'
        m2 = sub['TMB_group'] == 'Low TMB'
        lr = logrank_test(sub.loc[m1, 'os'], sub.loc[m2, 'os'],
                          sub.loc[m1, 'censOS'], sub.loc[m2, 'censOS'])
        axes[2].text(0.6, 0.95, f'p = {lr.p_value:.4f}', transform=axes[2].transAxes, fontsize=9)
        axes[2].set_title(f'OS by TMB (median split)')
        axes[2].set_xlabel('Time')
        axes[2].set_ylabel('Survival Probability')

    plt.tight_layout()
    plt.savefig(os.path.join(RESULTS_DIR, 'fig_08_km_survival.png'), dpi=150)
    plt.close()
    print("  Saved fig_08_km_survival.png")


def cox_regression(survival_df):
    """
    Cox PH regression with clinical covariates.
    Returns DataFrame with hazard ratios, CI, p-values.
    """
    df = survival_df.copy()

    # Prepare features
    ic_map = {'IC0': 0, 'IC1': 1, 'IC2+': 2}
    if 'IC Level' in df.columns:
        df['ic_ordinal'] = df['IC Level'].map(ic_map)

    if 'Immune phenotype' in df.columns:
        dummies = pd.get_dummies(df['Immune phenotype'], prefix='immpheno', drop_first=True)
        df = pd.concat([df, dummies], axis=1)

    tmb_col = 'FMOne mutation burden per MB'
    if tmb_col in df.columns:
        df['tmb_imputed'] = df[tmb_col].fillna(df[tmb_col].median())
        df['tmb_missing'] = df[tmb_col].isna().astype(int)

    if 'Baseline ECOG Score' in df.columns:
        df['ecog'] = pd.to_numeric(df['Baseline ECOG Score'], errors='coerce')

    if 'TCGA Subtype' in df.columns:
        dummies = pd.get_dummies(df['TCGA Subtype'], prefix='tcga', drop_first=True)
        df = pd.concat([df, dummies], axis=1)

    if 'Sex' in df.columns:
        df['sex_male'] = (df['Sex'] == 'M').astype(int)

    # Identify feature columns
    feature_cols = []
    for c in df.columns:
        if c in ['ic_ordinal', 'tmb_imputed', 'tmb_missing', 'ecog', 'sex_male'] or \
           c.startswith('immpheno_') or c.startswith('tcga_'):
            feature_cols.append(c)

    cox_df = df[['os', 'censOS'] + feature_cols].dropna()

    cph = CoxPHFitter()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cph.fit(cox_df, duration_col='os', event_col='censOS')

    summary = cph.summary.copy()
    summary = summary.rename(columns={
        'exp(coef)': 'hazard_ratio',
        'exp(coef) lower 95%': 'hr_ci_lower',
        'exp(coef) upper 95%': 'hr_ci_upper',
        'p': 'p_value',
    })
    summary = summary[['hazard_ratio', 'hr_ci_lower', 'hr_ci_upper', 'p_value']].round(4)
    c_index = cph.concordance_index_
    print(f"  Cox C-index: {c_index:.4f}")

    return summary, c_index, cph


def run_module2(data):
    """Execute Module 2: Clinical feature analysis."""
    print("\n" + "=" * 60)
    print("MODULE 2: Clinical Feature Association Analysis")
    print("=" * 60)

    cohorts = data['cohorts']
    binary_df = cohorts['binary']
    survival_df = cohorts['survival']

    # 2.1 Univariate associations
    print("\n[2.1] Univariate associations...")
    uni_df = univariate_association(binary_df)
    uni_df.to_csv(os.path.join(RESULTS_DIR, 'table_02_univariate_clinical.csv'), index=False)
    print(f"  Saved table_02_univariate_clinical.csv ({len(uni_df)} tests)")
    sig = uni_df[uni_df['fdr'] < 0.05] if 'fdr' in uni_df.columns else uni_df[uni_df['p_value'] < 0.05]
    print(f"  Significant after FDR correction: {len(sig)}")

    # 2.2 Response rate bar plots
    print("\n[2.2] Response rate bar plots...")
    response_rate_barplots(binary_df)

    # 2.3 TMB boxplot
    print("\n[2.3] TMB boxplot...")
    tmb_boxplot(binary_df)

    # 2.4 Multivariate logistic regression
    print("\n[2.4] Multivariate logistic regression...")
    logistic_df, logistic_model = multivariate_logistic(binary_df)
    logistic_df.to_csv(os.path.join(RESULTS_DIR, 'table_03_multivariate_logistic.csv'), index=False)
    print(f"  Saved table_03_multivariate_logistic.csv")

    # 2.5 Forest plot
    print("\n[2.5] Forest plot...")
    forest_plot(logistic_df)

    # 2.6 Kaplan-Meier survival curves
    print("\n[2.6] Kaplan-Meier survival analysis...")
    km_survival_analysis(survival_df)

    # 2.7 Cox PH regression
    print("\n[2.7] Cox proportional hazards regression...")
    cox_df, c_index, cph = cox_regression(survival_df)
    cox_df.to_csv(os.path.join(RESULTS_DIR, 'table_03b_cox_regression.csv'))
    print(f"  Saved table_03b_cox_regression.csv")

    print("\n[Module 2 Complete]")
    return {
        'univariate': uni_df,
        'logistic': logistic_df,
        'cox_summary': cox_df,
        'cox_c_index': c_index,
    }


if __name__ == '__main__':
    import sys
    sys.path.insert(0, '/code')
    from module1_data import run_module1
    os.makedirs(RESULTS_DIR, exist_ok=True)
    data = run_module1()
    run_module2(data)
