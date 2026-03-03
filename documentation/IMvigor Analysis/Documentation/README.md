# IMvigor210 Biomarker Analysis Pipeline

Comprehensive biomarker discovery and predictive modelling pipeline for the IMvigor210 clinical trial dataset. The pipeline integrates clinical, genomic, and transcriptomic data from a Phase II trial of atezolizumab (anti-PD-L1) in 348 patients with metastatic urothelial carcinoma.

## Dataset

The pipeline operates on three linked datasets:

| Dataset | Description | Dimensions |
|---------|-------------|------------|
| Clinical outcomes | Demographics, response, survival, PD-L1 status | 348 patients x 27 variables |
| RNA-seq gene expression | Raw counts | 31,287 genes x 348 patients |
| Genomic alterations | Long-format mutation records | 1,647 records across 264 patients, 222 genes |

Response is classified as CR (complete response), PR (partial response), SD (stable disease), PD (progressive disease), or NE (not evaluable). For binary analysis, CR/PR are responders and SD/PD are non-responders; NE patients are excluded.

## Pipeline Structure

The pipeline is organised into seven sequential modules orchestrated by `run.py`. Each module reads shared data objects from upstream modules and writes results (CSV tables, PNG figures, Tableau exports) to `/results`.

```
run (bash entrypoint)
 -> run.py (orchestrator)
     -> module1_data.py           Data loading, QC, preprocessing
     -> module2_clinical.py       Clinical feature association analysis
     -> module3_genomic.py        Genomic alteration analysis
     -> module4_expression.py     Differential expression & pathway analysis
     -> module5_signature.py      Predictive gene expression signature
     -> module6_integrated.py     Integrated classifier, survival model & summary
     -> module7_report.py         PDF report generation
```

Supporting utility:

- `gene_set_utils.py` -- GMT file loading, ORA (Fisher's exact), and pre-ranked GSEA used by Module 4.

### Execution order and dependencies

```
Module 1 (data) -> shared data dict passed to all downstream modules
   |-- Module 2 (clinical)     -- independent
   |-- Module 3 (genomic)      -- independent
   |-- Module 4 (expression)   -- independent
   |      \-- Module 5 (signature) -- needs DE results from Module 4
   |             \-- Module 6 (integrated) -- needs results from Modules 2-5
   \-- Module 7 (report)       -- reads CSV/PNG outputs from all modules
```

`run.py` runs them sequentially: M1 -> M2 -> M3 -> M4 -> M5 -> M6 -> M7.

---

## Module Details

### Module 1 -- Data Loading, QC & Preprocessing (`module1_data.py`)

Loads the three source datasets and prepares them for analysis:

- **Expression normalisation**: DESeq2-style size-factor normalisation followed by log2 transform: `log2(count / sizeFactor + 1)`.
- **Low-expression gene filtering**: Removes genes with fewer than 10 counts in fewer than 15 samples. Deduplicates gene symbols by keeping the highest-mean row.
- **Mutation matrix**: Converts long-format genomic alteration records into a binary patient-by-gene matrix (1 = any alteration present, NaN for patients without genomic data).
- **Cohort definition**: Defines binary analysis cohort (CR/PR vs SD/PD, NE excluded) and survival cohort (all patients).

**Outputs**:
- `table_01_cohort_summary.csv` -- Patient counts and demographics
- `fig_01_expression_distribution.png` -- Pre/post normalisation expression distributions
- `fig_02_pca_expression.png` -- PCA of normalised expression coloured by response

### Module 2 -- Clinical Feature Association Analysis (`module2_clinical.py`)

Tests associations between clinical variables and treatment response:

- **Univariate analysis**: Chi-square or Fisher's exact test for categorical variables (IC Level, TC Level, Immune phenotype, Sex, ECOG, TCGA Subtype, Lund, etc.); Mann-Whitney U for continuous variables (TMB, neoantigen burden, sizeFactor). All p-values are FDR-corrected (Benjamini-Hochberg).
- **Multivariate logistic regression**: Predicts binary response from IC Level (ordinal), Immune phenotype (dummy-coded), TMB (median-imputed with missingness indicator), ECOG, TCGA Subtype, and Sex. Reports odds ratios with 95% CIs.
- **Survival analysis**: Kaplan-Meier curves stratified by IC Level, Immune phenotype, and TMB (median split) with log-rank tests. Cox proportional hazards regression with the same covariates; reports hazard ratios and concordance index.

**Outputs**:
- `table_02_univariate_clinical.csv` -- Univariate association results
- `table_03_multivariate_logistic.csv` -- Logistic regression coefficients, ORs, CIs
- `table_03b_cox_regression.csv` -- Cox regression hazard ratios
- `fig_03` through `fig_08` -- Response rate bar plots, TMB boxplot, forest plot, KM curves

### Module 3 -- Genomic Alteration Analysis (`module3_genomic.py`)

Analyses genomic alterations in the 264-patient subset with mutation data:

- **Mutation frequency**: Per-gene mutation frequency in responders vs non-responders, filtered to genes altered in at least 10 patients. Fisher's exact test with FDR correction.
- **Oncoplot**: Top 20 most-mutated genes by patient, coloured by alteration type (short variant, rearrangement, copy number, ambiguous), with response status annotation.
- **TMB ROC analysis**: ROC curve for TMB as a continuous predictor of response, with optimal threshold via Youden's J statistic.
- **Co-occurrence analysis**: Pairwise Fisher's exact tests on the top 15 genes to identify co-occurring or mutually exclusive mutation patterns. Displayed as a log2 odds ratio heatmap.

**Outputs**:
- `table_04_mutation_association.csv` -- Per-gene mutation-response association
- `fig_09_mutation_freq_comparison.png` -- Mutation frequency bar chart (R vs NR)
- `fig_10_oncoplot.png` -- Oncoplot
- `fig_11_tmb_roc.png` -- TMB ROC curve
- `fig_12_co_occurrence.png` -- Co-occurrence heatmap

### Module 4 -- Differential Expression & Pathway Analysis (`module4_expression.py`)

Identifies genes and pathways differentially active between responders and non-responders:

- **Gene pre-filtering**: Removes genes with mean log2 expression < 1.0 and genes in the bottom 25% of expression variance. These independent filters reduce the multiple-testing burden while preserving FDR control.
- **Differential expression**: Mann-Whitney U test per gene on log2-normalised expression. Computes log2 fold-change (mean R - mean NR). Significance threshold: |log2FC| > 1 and FDR < 0.05.
- **Volcano plot and heatmap**: Volcano plot with labelled top hits; clustered heatmap of top DE genes by patient.
- **ORA (Over-Representation Analysis)**: Fisher's exact test (one-sided) on top up- and down-regulated gene lists against MSigDB collections (Hallmark, KEGG Medicus, Reactome, GO Biological Process). Gene sets filtered to 15-500 genes overlapping the tested universe.
- **Pre-ranked GSEA**: Ranking metric is `sign(log2FC) * -log10(p-value)`. Weighted running-sum enrichment scores normalised against 1,000 permutation-based null distributions (cached per gene-set size for efficiency). Empirical p-values and FDR correction.

**Outputs**:
- `table_05_de_results.csv` -- Full DE table (gene, FC, p, FDR)
- `table_06_top_de_genes.csv` -- Top 50 up + top 50 down genes
- `table_06b_ora_input_genes.csv`, `table_06b_gsea_input_ranking.csv` -- Traceability inputs
- `table_06c_ora_enrichment.csv`, `table_06d_gsea_enrichment.csv` -- Enrichment results
- `fig_13_volcano.png` -- Volcano plot
- `fig_14_de_heatmap.png` -- DE heatmap
- `fig_16a_ora_enrichment.png` -- ORA bar plot
- `fig_16b_gsea_enrichment.png` -- GSEA dot plot

### Module 5 -- Predictive Gene Expression Signature (`module5_signature.py`)

Builds and evaluates predictive models of treatment response using gene expression:

- **Feature preselection**: Top 5,000 genes by variance, further filtered to those with Mann-Whitney p < 0.05. Preselection is recomputed within each CV fold to prevent data leakage.
- **Elastic Net logistic regression**: Nested 5-fold stratified CV (outer loop for performance estimation, inner loop for tuning regularisation parameter C). L1 ratio = 0.5.
- **Random Forest classifier**: Same nested CV framework. Inner CV tunes n_estimators, max_depth, and min_samples_leaf via grid search.
- **Baseline comparators**: TMB-only and PD-L1 IC Level-only predictors for benchmarking.
- **Signature extraction**: Elastic Net refit on full dataset yields non-zero gene coefficients; Random Forest refit yields top 30 features by mean decrease in impurity. Overlap between the two is reported.
- **Gene signature score**: For each patient, the score is the weighted sum of z-scored expression of signature genes (using Elastic Net coefficients). This score is passed to Module 6.

**Outputs**:
- `table_07_signature_genes_elasticnet.csv` -- Gene, coefficient, rank
- `table_07b_signature_genes_rf.csv` -- Gene, importance, rank
- `table_08_model_comparison.csv` -- AUC for each predictor
- `fig_17_roc_comparison.png` -- ROC curves for all models
- `fig_18_signature_coefficients.png` -- Elastic Net coefficient bar plot
- `fig_18b_rf_feature_importance.png` -- Random Forest importance bar plot

### Module 6 -- Integrated Classifier, Survival Model & Summary (`module6_integrated.py`)

Combines expression, clinical, and genomic features into integrated models:

- **Integrated feature set**: Gene signature score (from M5), TMB (median-imputed with missingness indicator), IC Level (ordinal), Immune phenotype (dummy-coded), TCGA Subtype (dummy-coded), ECOG, and Sex.
- **Classifiers**: Expression-only LR, Clinical-only LR, Integrated LR, and Integrated Gradient Boosting (with hyperparameter tuning). All evaluated via 5-fold stratified CV reporting AUC, accuracy, precision, recall, and F1.
- **Confusion matrix**: Pooled out-of-fold predictions from the best model.
- **Integrated survival model**: Cox PH models comparing expression-only, clinical-only, and integrated feature sets by concordance index. Kaplan-Meier curve stratified by integrated risk score (median split).
- **Tableau-ready exports**: Five long-format CSV files for interactive visualisation:
  - `tableau_patient_features.csv` -- One row per patient per variable
  - `tableau_gene_expression.csv` -- Top 100 DE genes per patient
  - `tableau_mutations.csv` -- One row per patient per gene mutation
  - `tableau_model_performance.csv` -- One row per model per fold per metric
  - `tableau_survival.csv` -- One row per patient per biomarker stratification
- **Summary figure**: Four-panel overview (model AUC comparison, survival C-index comparison, response rate by IC Level, key findings).

**Outputs**:
- `table_09_classifier_performance.csv` -- Per-model AUC, accuracy, F1
- `table_09b_integrated_survival.csv` -- Survival model C-index comparison
- `table_10_analysis_summary.csv` -- Key findings across all modules
- `fig_19_classifier_roc.png`, `fig_19b_confusion_matrix.png`, `fig_19c_km_integrated_risk.png`
- `fig_20_summary_figure.png` -- Multi-panel summary
- Five `tableau_*.csv` files

### Module 7 -- PDF Report (`module7_report.py`)

Generates a formatted multi-page PDF report (`report_imvigor210.pdf`) that assembles narrative text, tables, and figures from all modules into sections:

1. Title page
2. Executive summary (with summary figure)
3. Cohort overview
4. Clinical feature associations
5. Genomic alteration landscape
6. Differential gene expression
7. Pathway enrichment (ORA + GSEA)
8. Predictive gene expression signature
9. Integrated model and survival analysis
10. Methods summary

### Gene Set Utilities (`gene_set_utils.py`)

Shared utility module for pathway enrichment analyses:

- **GMT loading**: Parses MSigDB GMT files (Hallmark, KEGG Medicus, Reactome, GO Biological Process, v2024.1).
- **ORA**: Fisher's exact test (one-sided, greater) for gene list overlap against gene sets.
- **Pre-ranked GSEA**: Weighted running-sum enrichment statistic with vectorised permutation-based null distribution generation. Null distributions are cached by gene-set size for efficiency.

---

## How to Run

```bash
# From the capsule root
python code/run.py
```

Or via the Code Ocean "Reproducible Run" entrypoint:

```bash
./code/run
```

All outputs are written to `/results`.

## Dependencies

- Python 3.12+
- pandas, numpy, scipy, scikit-learn, statsmodels
- matplotlib, seaborn
- lifelines (survival analysis)
- fpdf2 (PDF generation)
- Pillow (image handling for PDF)
- adjustText

## Input Data

Expected in `/root/capsule/data/`:

| File | Description |
|------|-------------|
| `imvigor210_clinical_data.csv` | Clinical outcomes (patient_id indexed, 27 columns including response, survival, PD-L1 status, TMB) |
| `imvigor210_gene_expression.csv` | RNA-seq counts (gene_id, gene_symbol, entrez_id + patient columns) |
| `imvigor210_genomic_alterations.csv` | Long-format genomic records (patient_id, gene, alteration_type) |
| `gene_sets/*.gmt` | MSigDB gene set collections (Hallmark, KEGG, Reactome, GO BP) |
