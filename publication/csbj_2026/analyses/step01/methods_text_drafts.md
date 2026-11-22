# Draft text for v7 (from Step 1) — not yet in the manuscript

Author decisions, 14 Nov 2026: keep Fisher's exact test and mention the hypergeometric option; report the reproducible z-scores and add empirical p-values.

---

## Methods — enrichment test (replaces the `hypergeom.sf` paragraph)

For each gene-set collection, gene sets were first filtered by size (10–600 genes by default). Let G be the set of genes annotated to at least one retained gene set in the collection. The input list L and the background B were each intersected with G, giving the collection-specific background size M = |B ∩ G| and input size n = |L ∩ G|. For a gene set t with K genes, the overlap k = |L ∩ G ∩ t| was tested with Fisher's exact test (two-sided; SciPy `fisher_exact`) on the 2 × 2 table [[k, K − k], [n − k, M − K − n + k]]. Within each collection, p-values were adjusted with the Benjamini–Hochberg procedure for the ORA comparison. iGEA also offers a one-sided hypergeometric test, P(X ≥ k) with X ~ Hypergeometric(M, K, n) (SciPy `hypergeom.sf(k − 1, M, K, n)`). The two tests gave very similar results for the HIV dataset: 252 versus 254 ORA terms at p < 0.01 (the hypergeometric test added one GO Cellular Component and one GO Molecular Function term), the same 87 terms at FDR < 0.05, and the same 62 iGEA terms. All results reported here use Fisher's exact test, which was also used to generate the random-list null distributions.

*Note for the author:* K is the full size of the gene set as given in the GMT file (it is not restricted to the background). With the default all-genes background this has little effect; say so if Reviewer 5 (minor 1) asks.

## Methods — statistical benchmarking (replaces the z-score description)

To assess the iGEA network against chance, iGEA was applied to random gene lists of 50 to 1,000 genes (in steps of 50; 1,000 lists per size; 20,000 lists in total). Each list was sampled uniformly without replacement from the default background and analysed with the same settings as the real data (Fisher's exact test, p < 0.01, minimum overlap 3, gene-set size 10–600, at most 30 iterations), separately for each collection. For an input list of size s, the null distribution was taken from the 1,000 random lists of size s rounded up to the next multiple of 50 (300 for the 272-gene HIV list). For each random list, the iGEA network was assembled from the same collections as the real analysis, and the size, number of genes, number of terms, number of edges and edges per gene of its largest connected cluster were recorded.

Each observed metric is reported with two statistics. The z-score is (x_obs − mean_null) / SD_null. Because the null distributions are counts with right-skewed distributions (skewness 0.6–1.3 for the HIV comparison), a z-score cannot be converted to a probability under a normal assumption, so we also report an empirical p-value, p = (r + 1) / (n + 1), where n is the number of random lists (1,000) and r is the number of random lists whose largest cluster had a value at least as large as the observed one [Phipson & Smyth 2010]. Adding 1 to both counts treats the observed list as one member of the reference set, so the p-value is never zero; with 1,000 random lists the smallest attainable value is 1/1,001 ≈ 0.001. The five metrics are correlated (all describe the same cluster) and are reported as complementary descriptions of network structure, not as independent tests.

Reference to add: Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. *Stat Appl Genet Mol Biol.* 2010;9:Article 39. doi:10.2202/1544-6115.1585

## Results — HIV network (replaces the z-score sentence)

The largest connected cluster contained 179 genes and 62 terms connected by 412 edges (2.30 edges per gene). None of the 1,000 random lists of matched size produced a cluster as large in genes, terms or edges (largest random cluster: 98 genes, 33 terms, 159 edges; empirical p ≤ 0.001 for each; z = 7.38, 8.6 and 11.45). Connectivity per gene was also higher than expected by chance (2.30 versus a random mean of 1.46; z = 3.6; empirical p = 0.005, 4 of 1,000 random lists at or above the observed value).

## Table for the main text or Supp Table 5

| Metric | Observed | Random mean ± SD | Random max | z | Random ≥ observed | Empirical p |
|---|---|---|---|---|---|---|
| Cluster size (genes + terms) | 241 | 46.4 ± 25.0 | 124 | 7.79 | 0 / 1,000 | ≤ 0.001 |
| Genes | 179 | 36.0 ± 19.4 | 98 | 7.38 | 0 / 1,000 | ≤ 0.001 |
| Terms | 62 | 10.4 ± 6.0 | 33 | 8.6 | 0 / 1,000 | ≤ 0.001 |
| Edges | 412 | 53.5 ± 31.3 | 159 | 11.45 | 0 / 1,000 | ≤ 0.001 |
| Edges per gene | 2.30 | 1.46 ± 0.23 | 3.31 | 3.6 | 4 / 1,000 | 0.005 |

*Rounding:* terms and edges-per-gene z are given to one decimal because the null mean varies in the third decimal between runs (8.58–8.59; 3.60–3.63) — when two random clusters tie for largest, the code picks one in an order that is not fixed. The empirical p-values do not change. A deterministic tie-break is a candidate fix for Step 12.
