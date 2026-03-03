"""
Module 7: PDF Report Generator
IMvigor210 Biomarker Analysis
"""

import os
import datetime
import numpy as np
import pandas as pd
from fpdf import FPDF
from PIL import Image

RESULTS_DIR = '/results'

# --- Colours (R, G, B) ---
NAVY = (30, 50, 100)
DARK_GREY = (60, 60, 60)
LIGHT_FILL = (235, 240, 250)
WHITE = (255, 255, 255)
TABLE_HEADER_BG = (30, 50, 100)
TABLE_HEADER_FG = (255, 255, 255)


# =========================================================================
# Helpers
# =========================================================================

def load_csv(filename):
    """Load a CSV from RESULTS_DIR; return empty DataFrame if missing."""
    path = os.path.join(RESULTS_DIR, filename)
    if os.path.exists(path):
        return pd.read_csv(path)
    return pd.DataFrame()


def fmt_p(val):
    """Format a p-value for display."""
    try:
        v = float(val)
    except (TypeError, ValueError):
        return str(val)
    if v < 0.001:
        return f"{v:.2e}"
    return f"{v:.4f}"


def fmt_num(val, decimals=2):
    """Format a number for display."""
    try:
        v = float(val)
    except (TypeError, ValueError):
        return str(val)
    if abs(v) >= 1000:
        return f"{v:,.0f}"
    return f"{v:.{decimals}f}"


# =========================================================================
# IMvigorReport  (FPDF subclass)
# =========================================================================

class IMvigorReport(FPDF):
    """PDF report for the IMvigor210 biomarker analysis pipeline."""

    def __init__(self):
        super().__init__(orientation='P', unit='mm', format='Letter')
        self.set_auto_page_break(auto=True, margin=20)
        self.set_margins(15, 20, 15)
        self._on_title_page = False

    # --- header / footer ---
    def header(self):
        if self._on_title_page:
            return
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(*DARK_GREY)
        self.cell(0, 5, 'IMvigor210 Biomarker Analysis Report', align='L')
        self.ln(8)

    def footer(self):
        if self._on_title_page:
            return
        self.set_y(-15)
        self.set_font('Helvetica', 'I', 8)
        self.set_text_color(*DARK_GREY)
        self.cell(0, 10, f'Page {self.page_no()}', align='C')

    # --- building blocks ---

    def section_heading(self, number, title):
        """Add a numbered section heading."""
        self.set_font('Helvetica', 'B', 14)
        self.set_text_color(*NAVY)
        self.cell(0, 10, f'{number}. {title}', new_x="LMARGIN", new_y="NEXT")
        self.set_draw_color(*NAVY)
        self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
        self.ln(4)

    def narrative(self, text):
        """Add a block of body text."""
        self.set_font('Helvetica', '', 10)
        self.set_text_color(*DARK_GREY)
        self.multi_cell(0, 5, text)
        self.ln(3)

    def safe_image(self, filename, caption=None, max_w=170):
        """Embed an image from RESULTS_DIR with auto aspect-ratio scaling."""
        path = os.path.join(RESULTS_DIR, filename)
        if not os.path.exists(path):
            self.set_font('Helvetica', 'I', 9)
            self.set_text_color(180, 0, 0)
            self.cell(0, 6, f'[{filename} not found]',
                      new_x="LMARGIN", new_y="NEXT")
            self.set_text_color(*DARK_GREY)
            return

        img = Image.open(path)
        w_px, h_px = img.size
        aspect = h_px / w_px
        display_w = min(max_w, self.w - self.l_margin - self.r_margin)
        display_h = display_w * aspect

        # Page break if needed
        if self.get_y() + display_h + 12 > self.h - 20:
            self.add_page()

        x_offset = (self.w - display_w) / 2
        self.image(path, x=x_offset, w=display_w)
        if caption:
            self.set_font('Helvetica', 'I', 8)
            self.set_text_color(*DARK_GREY)
            self.cell(0, 5, caption, new_x="LMARGIN", new_y="NEXT", align='C')
        self.ln(4)

    def add_table(self, df, col_widths=None, max_rows=None, font_size=7):
        """Render a DataFrame as a formatted table."""
        if df.empty:
            self.narrative("No data available.")
            return
        if max_rows and len(df) > max_rows:
            df = df.head(max_rows)

        cols = list(df.columns)
        n_cols = len(cols)
        usable_w = self.w - self.l_margin - self.r_margin
        if col_widths is None:
            col_widths = [usable_w / n_cols] * n_cols

        # Estimate table height and page-break if needed
        est_h = (len(df) + 1) * 5 + 4
        if self.get_y() + est_h > self.h - 25:
            self.add_page()

        # Header
        self.set_font('Helvetica', 'B', font_size)
        self.set_fill_color(*TABLE_HEADER_BG)
        self.set_text_color(*TABLE_HEADER_FG)
        for i, col in enumerate(cols):
            self.cell(col_widths[i], 5, str(col), border=1, fill=True, align='C')
        self.ln()

        # Rows
        self.set_font('Helvetica', '', font_size)
        self.set_text_color(*DARK_GREY)
        for row_idx, (_, row) in enumerate(df.iterrows()):
            if row_idx % 2 == 0:
                self.set_fill_color(*LIGHT_FILL)
                fill = True
            else:
                self.set_fill_color(*WHITE)
                fill = True

            for i, col in enumerate(cols):
                val = row[col]
                text = str(val) if pd.notna(val) else ''
                # Truncate long strings
                if len(text) > 40:
                    text = text[:37] + '...'
                self.cell(col_widths[i], 5, text, border=1, fill=fill, align='C')
            self.ln()
        self.ln(4)


# =========================================================================
# Section builders
# =========================================================================

def build_title_page(pdf):
    """Section 1: Title page."""
    pdf._on_title_page = True
    pdf.add_page()
    pdf.ln(50)
    pdf.set_font('Helvetica', 'B', 28)
    pdf.set_text_color(*NAVY)
    pdf.cell(0, 15, 'IMvigor210', align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.set_font('Helvetica', '', 16)
    pdf.cell(0, 10, 'Biomarker Analysis Report', align='C',
             new_x="LMARGIN", new_y="NEXT")
    pdf.ln(10)
    pdf.set_font('Helvetica', '', 12)
    pdf.set_text_color(*DARK_GREY)
    pdf.cell(0, 8,
             'Atezolizumab (anti-PD-L1) in Metastatic Urothelial Carcinoma',
             align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.cell(0, 8,
             'Comprehensive Biomarker Discovery & Predictive Modelling',
             align='C', new_x="LMARGIN", new_y="NEXT")
    pdf.ln(20)
    pdf.set_font('Helvetica', 'I', 10)
    today = datetime.date.today().strftime('%B %d, %Y')
    pdf.cell(0, 6, f'Generated: {today}', align='C',
             new_x="LMARGIN", new_y="NEXT")
    pdf._on_title_page = False


def build_executive_summary(pdf):
    """Section 2: Executive summary."""
    pdf.add_page()
    pdf.section_heading(1, 'Executive Summary')

    summary = load_csv('table_10_analysis_summary.csv')
    cohort = load_csv('table_01_cohort_summary.csv')

    # Total patients
    n_patients = ''
    if not cohort.empty:
        total = cohort.loc[cohort['variable'] == 'binaryResponse', 'count']
        if len(total) > 0:
            n_patients = str(int(total.sum()))

    # Response rate
    resp_rate = ''
    if not cohort.empty:
        cr_pr = cohort.loc[(cohort['variable'] == 'binaryResponse') &
                           (cohort['value'] == 'CR/PR')]
        if len(cr_pr) > 0:
            resp_rate = f"{cr_pr['percent'].values[0]:.1f}%"

    # DE genes
    n_de = ''
    if not summary.empty:
        de_row = summary[summary['finding'].str.contains('differentially expressed',
                                                          case=False, na=False)]
        if len(de_row) > 0:
            n_de = de_row['finding'].values[0]

    # Best integrated AUC
    best_auc = ''
    if not summary.empty:
        int_rows = summary[summary['finding'].str.contains('Integrated.*AUC',
                                                            na=False)]
        if len(int_rows) > 0:
            details = int_rows['detail'].astype(str)
            try:
                aucs = [float(d) for d in details]
                best_auc = f"{max(aucs):.3f}"
            except ValueError:
                best_auc = details.values[0]

    # Build narrative
    parts = []
    parts.append(
        'This report summarises the IMvigor210 biomarker analysis pipeline, '
        'which integrates clinical, genomic, and transcriptomic data from a '
        'Phase-II trial of atezolizumab in metastatic urothelial carcinoma.'
    )
    if n_patients:
        line = f'The analysis cohort comprises {n_patients} evaluable patients'
        if resp_rate:
            line += f' with an objective response rate (CR/PR) of {resp_rate}'
        parts.append(line + '.')
    if n_de:
        parts.append(f'{n_de}.')
    if best_auc:
        parts.append(
            f'The best integrated classifier achieved an AUC of {best_auc} '
            'using a combination of expression and clinical features.'
        )

    pdf.narrative(' '.join(parts))

    # Summary figure
    pdf.safe_image('fig_20_summary_figure.png',
                   caption='Figure: Pipeline summary overview')


def build_cohort_overview(pdf):
    """Section 3: Cohort overview."""
    pdf.add_page()
    pdf.section_heading(2, 'Cohort Overview')

    cohort = load_csv('table_01_cohort_summary.csv')

    if not cohort.empty:
        n_total = int(cohort.loc[cohort['variable'] == 'binaryResponse',
                                  'count'].sum()) if len(
            cohort[cohort['variable'] == 'binaryResponse']) > 0 else '?'

        resp_rows = cohort[cohort['variable'] == 'Best Confirmed Overall Response']
        resp_dist = ', '.join(
            f"{r['value']}: {int(r['count'])} ({r['percent']}%)"
            for _, r in resp_rows.iterrows()
        )

        text = (
            f'The trial enrolled {n_total} evaluable patients '
            f'(excluding NE). Best confirmed overall response distribution: '
            f'{resp_dist}.'
        )
        pdf.narrative(text)

        # Table
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Cohort Summary', new_x="LMARGIN", new_y="NEXT")
        widths = [45, 50, 30, 30]
        pdf.add_table(cohort, col_widths=widths)
    else:
        pdf.narrative('Cohort summary data not available.')

    pdf.safe_image('fig_01_expression_distribution.png',
                   caption='Figure: Gene expression distributions before and after normalization')
    pdf.safe_image('fig_02_pca_expression.png',
                   caption='Figure: PCA of normalized gene expression coloured by response')


def build_clinical_associations(pdf):
    """Section 4: Clinical associations."""
    pdf.add_page()
    pdf.section_heading(3, 'Clinical Feature Associations')

    univar = load_csv('table_02_univariate_clinical.csv')
    logistic = load_csv('table_03_multivariate_logistic.csv')
    cox = load_csv('table_03b_cox_regression.csv')

    # Narrative
    parts = []
    if not univar.empty:
        n_tested = len(univar)
        n_sig = (univar['fdr'] < 0.05).sum() if 'fdr' in univar.columns else 0
        parts.append(
            f'{n_tested} clinical variables were tested for association with '
            f'response; {n_sig} were significant at FDR < 0.05.'
        )
        if n_sig > 0:
            top = univar.sort_values('p_value').head(1)
            parts.append(
                f'The strongest association was {top["variable"].values[0]} '
                f'(p = {fmt_p(top["p_value"].values[0])}).'
            )

    if not logistic.empty:
        log_excl = logistic[logistic['variable'] != 'Intercept']
        sig_log = log_excl[log_excl['p_value'] < 0.05] if 'p_value' in log_excl.columns else pd.DataFrame()
        parts.append(
            f'Multivariate logistic regression identified {len(sig_log)} '
            f'significant predictors (p < 0.05) from {len(log_excl)} variables.'
        )

    if not cox.empty:
        sig_cox = cox[cox['p_value'] < 0.05] if 'p_value' in cox.columns else pd.DataFrame()
        parts.append(
            f'Cox proportional-hazards regression found {len(sig_cox)} '
            f'significant covariates for overall survival.'
        )

    pdf.narrative(' '.join(parts))

    # Univariate table
    if not univar.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Univariate Associations',
                 new_x="LMARGIN", new_y="NEXT")
        uv_display = univar.copy()
        uv_display['p_value'] = uv_display['p_value'].apply(fmt_p)
        uv_display['fdr'] = uv_display['fdr'].apply(fmt_p)
        uv_display['statistic'] = uv_display['statistic'].apply(
            lambda x: fmt_num(x, 2))
        widths = [40, 30, 30, 30, 30]
        pdf.add_table(uv_display, col_widths=widths)

    # Logistic table
    if not logistic.empty:
        log_display = logistic[logistic['variable'] != 'Intercept'].copy()
        for c in ['coefficient', 'odds_ratio', 'or_ci_lower', 'or_ci_upper']:
            if c in log_display.columns:
                log_display[c] = log_display[c].apply(lambda x: fmt_num(x, 3))
        if 'p_value' in log_display.columns:
            log_display['p_value'] = log_display['p_value'].apply(fmt_p)
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Multivariate Logistic Regression',
                 new_x="LMARGIN", new_y="NEXT")
        widths = [35, 22, 22, 22, 22, 25]
        pdf.add_table(log_display, col_widths=widths)

    # Figures
    pdf.safe_image('fig_07_forest_plot.png',
                   caption='Figure: Forest plot of univariate clinical associations')
    pdf.safe_image('fig_08_km_survival.png',
                   caption='Figure: Kaplan-Meier overall survival by response')


def build_genomic_landscape(pdf):
    """Section 5: Genomic landscape."""
    pdf.add_page()
    pdf.section_heading(4, 'Genomic Alteration Landscape')

    mut = load_csv('table_04_mutation_association.csv')
    summary = load_csv('table_10_analysis_summary.csv')

    parts = []
    if not mut.empty:
        n_genes = len(mut)
        n_sig = (mut['fdr'] < 0.05).sum() if 'fdr' in mut.columns else 0
        parts.append(
            f'{n_genes} genes were tested for mutation-response association; '
            f'{n_sig} reached FDR < 0.05.'
        )
        top = mut.sort_values('p_value').head(1)
        parts.append(
            f'The top gene was {top["gene"].values[0]} '
            f'(OR = {fmt_num(top["odds_ratio"].values[0], 2)}, '
            f'p = {fmt_p(top["p_value"].values[0])}).'
        )

    if not summary.empty:
        tmb_row = summary[summary['finding'].str.contains('TMB.*predictor',
                                                           na=False, case=False)]
        if len(tmb_row) > 0:
            parts.append(f'Tumour mutational burden as a response predictor: '
                         f'{tmb_row["detail"].values[0]}.')

    pdf.narrative(' '.join(parts))

    # Table: top 15 mutations
    if not mut.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Top 15 Mutated Genes',
                 new_x="LMARGIN", new_y="NEXT")
        mut_display = mut.head(15)[['gene', 'total_mutated', 'resp_freq',
                                     'nonresp_freq', 'odds_ratio',
                                     'p_value', 'fdr']].copy()
        mut_display['resp_freq'] = mut_display['resp_freq'].apply(
            lambda x: fmt_num(x, 3))
        mut_display['nonresp_freq'] = mut_display['nonresp_freq'].apply(
            lambda x: fmt_num(x, 3))
        mut_display['odds_ratio'] = mut_display['odds_ratio'].apply(
            lambda x: fmt_num(x, 2))
        mut_display['p_value'] = mut_display['p_value'].apply(fmt_p)
        mut_display['fdr'] = mut_display['fdr'].apply(fmt_p)
        widths = [22, 22, 22, 26, 22, 26, 26]
        pdf.add_table(mut_display, col_widths=widths)

    pdf.safe_image('fig_09_mutation_freq_comparison.png',
                   caption='Figure: Mutation frequency in responders vs non-responders')
    pdf.safe_image('fig_10_oncoplot.png',
                   caption='Figure: Oncoplot of top mutated genes')
    pdf.safe_image('fig_11_tmb_roc.png',
                   caption='Figure: TMB as a response predictor (ROC)')


def build_differential_expression(pdf):
    """Section 6: Differential expression."""
    pdf.add_page()
    pdf.section_heading(5, 'Differential Gene Expression')

    de = load_csv('table_05_de_results.csv')
    top_de = load_csv('table_06_top_de_genes.csv')

    parts = []
    if not de.empty:
        n_tested = len(de)
        sig_mask = (de['fdr'] < 0.05) & (de['log2FC'].abs() > 1)
        n_sig = sig_mask.sum()
        n_up = ((de['fdr'] < 0.05) & (de['log2FC'] > 1)).sum()
        n_down = ((de['fdr'] < 0.05) & (de['log2FC'] < -1)).sum()

        parts.append(
            f'Mann-Whitney U tests were performed on {n_tested:,} genes '
            f'(after independent pre-filtering) comparing responders (CR/PR) '
            f'vs non-responders (SD/PD).'
        )
        parts.append(
            f'{n_sig} genes were significant at |log2FC| > 1 and FDR < 0.05 '
            f'({n_up} up-regulated in responders, {n_down} down-regulated).'
        )
        if n_sig > 0:
            top_up = de[de['log2FC'] > 0].head(3)['gene'].tolist()
            top_down = de[de['log2FC'] < 0].head(3)['gene'].tolist()
            if top_up:
                parts.append(f'Top up-regulated: {", ".join(top_up)}.')
            if top_down:
                parts.append(f'Top down-regulated: {", ".join(top_down)}.')

    pdf.narrative(' '.join(parts))

    # Table: top 10 up + 10 down
    if not de.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Top Differentially Expressed Genes',
                 new_x="LMARGIN", new_y="NEXT")
        up = de[de['log2FC'] > 0].head(10)
        down = de[de['log2FC'] < 0].head(10)
        de_display = pd.concat([up, down])[
            ['gene', 'log2FC', 'mean_responder', 'mean_nonresponder',
             'p_value', 'fdr']].copy()
        de_display['log2FC'] = de_display['log2FC'].apply(
            lambda x: fmt_num(x, 3))
        de_display['mean_responder'] = de_display['mean_responder'].apply(
            lambda x: fmt_num(x, 2))
        de_display['mean_nonresponder'] = de_display['mean_nonresponder'].apply(
            lambda x: fmt_num(x, 2))
        de_display['p_value'] = de_display['p_value'].apply(fmt_p)
        de_display['fdr'] = de_display['fdr'].apply(fmt_p)
        widths = [25, 22, 30, 35, 28, 28]
        pdf.add_table(de_display, col_widths=widths, max_rows=20)

    pdf.safe_image('fig_13_volcano.png',
                   caption='Figure: Volcano plot of differential expression')
    pdf.safe_image('fig_14_de_heatmap.png',
                   caption='Figure: Heatmap of top differentially expressed genes')


def build_pathway_enrichment(pdf):
    """Section 7: Pathway enrichment (ORA + GSEA)."""
    pdf.add_page()
    pdf.section_heading(6, 'Pathway Enrichment Analysis')

    ora = load_csv('table_06c_ora_enrichment.csv')
    gsea = load_csv('table_06d_gsea_enrichment.csv')

    parts = []

    # ORA narrative
    if not ora.empty:
        n_sig_ora = (ora['fdr'] < 0.05).sum()
        collections_ora = ora['collection'].nunique()
        parts.append(
            f"Over-Representation Analysis (ORA, Fisher's exact test) was run "
            f"across {collections_ora} gene set collection(s). "
            f"{n_sig_ora} gene sets were enriched at FDR < 0.05."
        )
        if n_sig_ora > 0:
            top_ora = ora[ora['fdr'] < 0.05].head(3)
            top_terms = ', '.join(
                f'{r["term"]} ({r["direction"]})'
                for _, r in top_ora.iterrows()
            )
            parts.append(f'Top ORA hits: {top_terms}.')

    # GSEA narrative
    if not gsea.empty:
        n_sig_gsea = (gsea['fdr'] < 0.05).sum()
        collections_gsea = gsea['collection'].nunique()
        parts.append(
            f'Pre-ranked GSEA (1,000 permutations) was run across '
            f'{collections_gsea} collection(s). '
            f'{n_sig_gsea} gene sets reached FDR < 0.05.'
        )
        if n_sig_gsea > 0:
            top_gsea = gsea[gsea['fdr'] < 0.05].reindex(
                gsea['nes'].abs().sort_values(ascending=False).index
            ).head(3)
            top_g = ', '.join(
                f'{r["term"]} (NES={fmt_num(r["nes"], 2)})'
                for _, r in top_gsea.iterrows()
            )
            parts.append(f'Top GSEA hits: {top_g}.')

    if not parts:
        parts.append('No pathway enrichment results available.')
    pdf.narrative(' '.join(parts))

    # ORA table
    if not ora.empty and (ora['fdr'] < 0.25).any():
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Top ORA Enrichments',
                 new_x="LMARGIN", new_y="NEXT")
        ora_display = ora[ora['fdr'] < 0.25].head(10)[
            ['collection', 'direction', 'term', 'overlap',
             'term_size', 'p_value', 'fdr']].copy()
        ora_display['p_value'] = ora_display['p_value'].apply(fmt_p)
        ora_display['fdr'] = ora_display['fdr'].apply(fmt_p)
        widths = [22, 25, 42, 18, 20, 22, 22]
        pdf.add_table(ora_display, col_widths=widths)

    # GSEA table
    if not gsea.empty and (gsea['fdr'] < 0.25).any():
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Top GSEA Enrichments',
                 new_x="LMARGIN", new_y="NEXT")
        gsea_display = gsea[gsea['fdr'] < 0.25].copy()
        gsea_display = gsea_display.reindex(
            gsea_display['nes'].abs().sort_values(ascending=False).index
        ).head(10)[['collection', 'term', 'es', 'nes', 'size',
                     'p_value', 'fdr']].copy()
        gsea_display['es'] = gsea_display['es'].apply(
            lambda x: fmt_num(x, 3))
        gsea_display['nes'] = gsea_display['nes'].apply(
            lambda x: fmt_num(x, 3))
        gsea_display['p_value'] = gsea_display['p_value'].apply(fmt_p)
        gsea_display['fdr'] = gsea_display['fdr'].apply(fmt_p)
        widths = [22, 42, 18, 18, 16, 22, 22]
        pdf.add_table(gsea_display, col_widths=widths)

    pdf.safe_image('fig_16a_ora_enrichment.png',
                   caption='Figure: ORA enrichment bar plot')
    pdf.safe_image('fig_16b_gsea_enrichment.png',
                   caption='Figure: GSEA enrichment dot plot')


def build_predictive_signature(pdf):
    """Section 8: Predictive gene expression signature."""
    pdf.add_page()
    pdf.section_heading(7, 'Predictive Gene Expression Signature')

    en_genes = load_csv('table_07_signature_genes_elasticnet.csv')
    rf_genes = load_csv('table_07b_signature_genes_rf.csv')
    model_cmp = load_csv('table_08_model_comparison.csv')

    parts = []

    if not en_genes.empty:
        n_en = len(en_genes)
        top_pos = en_genes[en_genes['coefficient'] > 0].head(3)['gene'].tolist()
        top_neg = en_genes[en_genes['coefficient'] < 0].head(3)['gene'].tolist()
        parts.append(
            f'Elastic Net regularised logistic regression selected '
            f'{n_en} signature genes via nested cross-validation.'
        )
        if top_pos:
            parts.append(
                f'Top positive predictors (higher in R): {", ".join(top_pos)}.'
            )
        if top_neg:
            parts.append(
                f'Top negative predictors (higher in NR): {", ".join(top_neg)}.'
            )

    if not rf_genes.empty:
        n_rf = len(rf_genes)
        parts.append(f'Random Forest identified {n_rf} important features.')

    if not en_genes.empty and not rf_genes.empty:
        overlap = set(en_genes['gene']) & set(rf_genes['gene'])
        if overlap:
            parts.append(
                f'{len(overlap)} genes overlap between Elastic Net and '
                f'Random Forest selections.'
            )

    if not model_cmp.empty:
        parts.append('Model comparison (5-fold nested CV):')
        for _, row in model_cmp.iterrows():
            parts.append(
                f'  {row["model"]}: AUC = {fmt_num(row["mean_auc"], 3)} '
                f'\u00b1 {fmt_num(row["std_auc"], 3)}'
            )

    pdf.narrative(' '.join(parts))

    # Model comparison table
    if not model_cmp.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Model Comparison',
                 new_x="LMARGIN", new_y="NEXT")
        mc_display = model_cmp.copy()
        mc_display['mean_auc'] = mc_display['mean_auc'].apply(
            lambda x: fmt_num(x, 3))
        mc_display['std_auc'] = mc_display['std_auc'].apply(
            lambda x: fmt_num(x, 3))
        widths = [55, 45, 45]
        pdf.add_table(mc_display, col_widths=widths, font_size=9)

    # EN signature genes table (top 15)
    if not en_genes.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Top 15 Elastic Net Signature Genes',
                 new_x="LMARGIN", new_y="NEXT")
        en_display = en_genes.head(15).copy()
        en_display['coefficient'] = en_display['coefficient'].apply(
            lambda x: fmt_num(x, 4))
        widths = [50, 50, 30]
        pdf.add_table(en_display, col_widths=widths, font_size=9)

    pdf.safe_image('fig_17_roc_comparison.png',
                   caption='Figure: ROC comparison across predictive models')
    pdf.safe_image('fig_18_signature_coefficients.png',
                   caption='Figure: Elastic Net signature gene coefficients')


def build_integrated_model(pdf):
    """Section 9: Integrated model & survival."""
    pdf.add_page()
    pdf.section_heading(8, 'Integrated Model & Survival Analysis')

    clf_perf = load_csv('table_09_classifier_performance.csv')
    surv = load_csv('table_09b_integrated_survival.csv')

    parts = []
    if not clf_perf.empty:
        best = clf_perf.sort_values('mean_auc', ascending=False).iloc[0]
        parts.append(
            f'The best-performing classifier was {best["model"]} with a '
            f'mean AUC of {fmt_num(best["mean_auc"], 3)}.'
        )
        # Compare expression-only vs integrated
        expr_only = clf_perf[clf_perf['model'].str.contains('Expression',
                                                             case=False, na=False)]
        clin_only = clf_perf[clf_perf['model'].str.contains('Clinical',
                                                             case=False, na=False)]
        if len(expr_only) > 0 and len(clin_only) > 0:
            parts.append(
                f'Expression-only: AUC = {fmt_num(expr_only["mean_auc"].values[0], 3)}, '
                f'Clinical-only: AUC = {fmt_num(clin_only["mean_auc"].values[0], 3)}.'
            )

    if not surv.empty:
        parts.append('Survival model C-indices:')
        for _, row in surv.iterrows():
            parts.append(f'  {row["model"]}: C-index = {fmt_num(row["c_index"], 4)}')

    pdf.narrative(' '.join(parts))

    # Classifier performance table
    if not clf_perf.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Classifier Performance',
                 new_x="LMARGIN", new_y="NEXT")
        clf_display = clf_perf.copy()
        for c in ['mean_auc', 'std_auc', 'mean_accuracy', 'mean_f1']:
            if c in clf_display.columns:
                clf_display[c] = clf_display[c].apply(
                    lambda x: fmt_num(x, 3))
        widths = [40, 28, 28, 32, 28]
        pdf.add_table(clf_display, col_widths=widths, font_size=8)

    # Survival C-index table
    if not surv.empty:
        pdf.set_font('Helvetica', 'B', 10)
        pdf.set_text_color(*NAVY)
        pdf.cell(0, 7, 'Table: Survival Model C-Indices',
                 new_x="LMARGIN", new_y="NEXT")
        surv_display = surv.copy()
        surv_display['c_index'] = surv_display['c_index'].apply(
            lambda x: fmt_num(x, 4))
        widths = [70, 50]
        pdf.add_table(surv_display, col_widths=widths, font_size=9)

    pdf.safe_image('fig_19_classifier_roc.png',
                   caption='Figure: Integrated classifier ROC curves')
    pdf.safe_image('fig_19c_km_integrated_risk.png',
                   caption='Figure: Kaplan-Meier survival by integrated risk score')


def build_methods(pdf):
    """Section 10: Methods summary."""
    pdf.add_page()
    pdf.section_heading(9, 'Methods Summary')

    text = (
        "Data source: IMvigor210 Phase-II clinical trial of atezolizumab "
        "(anti-PD-L1) in metastatic urothelial carcinoma. "
        "Gene expression was normalised using DESeq2-style size factors "
        "(log2(count/sizeFactor + 1)) and low-expression genes were filtered "
        "(minimum 10 counts in at least 15 samples). "
        "Prior to differential expression testing, an independent pre-filter "
        "removed genes with mean log2 expression < 1.0 and genes in the "
        "bottom 25% of expression variance, reducing the multiple-testing "
        "burden while preserving FDR control."
    )
    pdf.narrative(text)

    text = (
        "Differential expression was assessed using Mann-Whitney U tests "
        "(two-sided) comparing responders (CR/PR) vs non-responders (SD/PD). "
        "P-values were corrected for multiple testing using the "
        "Benjamini-Hochberg FDR procedure."
    )
    pdf.narrative(text)

    text = (
        "Over-Representation Analysis (ORA) used Fisher's exact test "
        "(one-sided, greater) to test whether up- or down-regulated gene "
        "lists were enriched in curated gene sets from MSigDB v2024.1 "
        "(Hallmark, KEGG Medicus, Reactome, GO Biological Process). "
        "Gene sets were filtered to 15-500 genes overlapping the tested "
        "gene universe. FDR correction was applied within each analysis."
    )
    pdf.narrative(text)

    text = (
        "Pre-ranked Gene Set Enrichment Analysis (GSEA) used the ranking "
        "metric sign(log2FC) x -log10(p-value). Enrichment scores were "
        "computed via a weighted running-sum statistic, normalised against "
        "1,000 permutation-based null distributions (cached per gene-set "
        "size for efficiency). Empirical p-values and FDR were calculated "
        "from the permutation distributions."
    )
    pdf.narrative(text)

    text = (
        "Predictive gene expression signatures were built using Elastic Net "
        "regularised logistic regression and Random Forest classifiers within "
        "a nested 5-fold cross-validation framework. Model performance was "
        "compared by mean AUC across outer folds."
    )
    pdf.narrative(text)

    text = (
        "Integrated classifiers combined expression-derived signature scores, "
        "clinical features, and genomic features (TMB, key mutations) using "
        "logistic regression and gradient boosting. Survival analysis used "
        "Cox proportional-hazards models, with concordance index (C-index) "
        "as the performance metric."
    )
    pdf.narrative(text)


# =========================================================================
# Orchestrator
# =========================================================================

def run_module7():
    """Execute Module 7: Generate PDF report."""
    print("\n" + "=" * 60)
    print("MODULE 7: PDF Report Generation")
    print("=" * 60)

    pdf = IMvigorReport()

    print("\n[7.1] Building title page...")
    build_title_page(pdf)

    print("[7.2] Building executive summary...")
    build_executive_summary(pdf)

    print("[7.3] Building cohort overview...")
    build_cohort_overview(pdf)

    print("[7.4] Building clinical associations...")
    build_clinical_associations(pdf)

    print("[7.5] Building genomic landscape...")
    build_genomic_landscape(pdf)

    print("[7.6] Building differential expression...")
    build_differential_expression(pdf)

    print("[7.7] Building pathway enrichment...")
    build_pathway_enrichment(pdf)

    print("[7.8] Building predictive signature...")
    build_predictive_signature(pdf)

    print("[7.9] Building integrated model...")
    build_integrated_model(pdf)

    print("[7.10] Building methods summary...")
    build_methods(pdf)

    out_path = os.path.join(RESULTS_DIR, 'report_imvigor210.pdf')
    pdf.output(out_path)
    size_mb = os.path.getsize(out_path) / (1024 * 1024)
    print(f"\n  Saved {out_path} ({size_mb:.1f} MB, {pdf.page_no()} pages)")
    print("\n[Module 7 Complete]")


if __name__ == '__main__':
    run_module7()
