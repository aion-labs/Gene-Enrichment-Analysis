# Draft text for v7 — per-run random-list statistics (supersedes the statistics part of `step01/methods_text_drafts.md`)

14 Nov 2026. Numbers from `step02_per_run_null/` (10,000 random lists, seed 2026). The enrichment-test paragraph in `step01/methods_text_drafts.md` (Fisher's exact test, hypergeometric option) is unchanged.

---

## Methods — Statistical evaluation of iGEA networks

Because terms are selected adaptively, with each iteration conditioned on the genes removed before it, nominal p-values are used only to select terms. Statistical support is assessed for the resulting network by comparison with random gene lists processed in exactly the same way.

For each analysis, iGEA draws n random gene lists (n = 10,000 by default; user-adjustable). Each list contains as many genes as the validated input list and is sampled uniformly without replacement from the background used in the analysis. Every random list is analysed with the same gene-set collections and the same parameters as the input (statistical test, p-value threshold, minimum overlap, gene-set size limits and maximum number of iterations), and its cross-collection network is assembled. For the largest connected cluster of each network we record the number of nodes (genes + terms), genes, terms and edges, and the mean number of edges per gene. Because the random lists match the input in size, background and every parameter, no approximation by precomputed or size-binned reference distributions is needed.

For each metric we report the z-score, (x_obs − mean_random) / SD_random, and the empirical p-value p = (r + 1) / (n + 1), where r is the number of random lists whose largest cluster has a value at least as large as the observed one [Phipson & Smyth 2010]. Counting the observed list as one member of the reference set keeps p above zero; with 10,000 random lists the smallest attainable value is 1/10,001 ≈ 0.0001. The random-list distributions are right-skewed counts, so z-scores are given for reference and inference is based on the empirical p-values. The five metrics describe the same cluster, are correlated, and are reported as complementary descriptions rather than independent tests. Network-level significance indicates that the input genes form a more connected cross-collection network than random lists; it does not establish the significance of individual terms. Random lists are generated from a fixed seed (default 2026), so results are exactly reproducible.

**Implementation and run time.** To make per-analysis random lists practical, iGEA includes a vectorised implementation of the iterative procedure: each collection is held as a sparse gene-set × gene matrix, so the overlaps of all gene sets with the remaining genes are obtained in one matrix–vector product; p-values are computed only for gene sets that meet the minimum overlap; and the selected term is recomputed with the same SciPy function as the reference implementation. Its output is identical to the reference implementation (selected terms, removed genes and p-values; verified by automated tests). On a 4-core Linux virtual machine on a MacBook Pro, 10,000 random lists took about 10 s for 75 genes, 1 min for 272 genes and 4 min for 426 genes; 1,000 lists of 1,000 genes took about 3 min. The original implementation required about 6 s per random list of 300 genes. The hosted web application uses 1,000 random lists by default because of its limited computing resources.

## Results — HIV (replace the z-score sentence)

The largest connected cluster contained 179 genes and 62 terms from all 10 collections, connected by 412 edges (2.30 edges per gene). Compared with 10,000 random lists of 272 genes analysed in the same way, none produced a cluster with as many genes, terms or edges (largest random cluster: 112 genes, 41 terms, 195 edges; empirical p ≤ 0.0001 for each; z = 8.0, 9.0 and 12.2). Edges per gene were also higher than in random lists (2.30 versus 1.45; z = 3.4; p = 0.008).

## Results — hypoxia and vaccination (replace "z-scores > 3 across metrics of genes and edges")

For both datasets the largest cluster exceeded every random list in size, genes, terms and edges (10,000 random lists each; empirical p ≤ 0.0001 for each metric). Hypoxia (75 genes): 64 genes, 32 terms, 171 edges from 9 collections (z = 10.0–16.2); edges per gene 2.67 versus 1.38 (z = 3.3; p = 0.011). Vaccination (426 genes): 272 genes, 75 terms, 702 edges from 9 collections (z = 8.5–15.1); edges per gene 2.58 versus 1.44 (z = 5.6; p = 0.0013).

## Table (main text or Supplementary Table 5)

| Dataset (genes) | Metric | Observed | Random mean ± SD (max) | z | r / 10,000 | Empirical p |
|---|---|---|---|---|---|---|
| HIV (272) | Cluster size | 241 | 42.9 ± 23.6 (148) | 8.39 | 0 | ≤ 0.0001 |
| | Genes | 179 | 33.1 ± 18.2 (112) | 8.04 | 0 | ≤ 0.0001 |
| | Terms | 62 | 9.7 ± 5.8 (41) | 8.97 | 0 | ≤ 0.0001 |
| | Edges | 412 | 49.1 ± 29.7 (195) | 12.24 | 0 | ≤ 0.0001 |
| | Edges per gene | 2.30 | 1.45 ± 0.25 (5.0) | 3.41 | 80 | 0.008 |
| Hypoxia (75) | Cluster size | 96 | 12.9 ± 7.9 (57) | 10.59 | 0 | ≤ 0.0001 |
| | Genes | 64 | 9.3 ± 5.4 (39) | 10.04 | 0 | ≤ 0.0001 |
| | Terms | 32 | 3.6 ± 2.6 (19) | 11.10 | 0 | ≤ 0.0001 |
| | Edges | 171 | 13.6 ± 9.7 (74) | 16.22 | 0 | ≤ 0.0001 |
| | Edges per gene | 2.67 | 1.38 ± 0.39 (5.75) | 3.29 | 106 | 0.011 |
| Vaccination (426) | Cluster size | 347 | 63.1 ± 32.8 (214) | 8.65 | 0 | ≤ 0.0001 |
| | Genes | 272 | 50.1 ± 26.1 (175) | 8.49 | 0 | ≤ 0.0001 |
| | Terms | 75 | 13.0 ± 7.2 (45) | 8.66 | 0 | ≤ 0.0001 |
| | Edges | 702 | 73.4 ± 41.5 (270) | 15.14 | 0 | ≤ 0.0001 |
| | Edges per gene | 2.58 | 1.44 ± 0.20 (3.36) | 5.64 | 12 | 0.0013 |

Settings for all three: 10 collections (MSigDB v2025.1: Hallmark; Reactome, KEGG Legacy, KEGG MEDICUS, BioCarta, PID, WikiPathways; GO BP, CC, MF), all_genes background (23,391 validated symbols), Fisher's exact test, p < 0.01, minimum overlap 3, gene-set size 10–600, no iteration cap; seed 2026.

Reference: Phipson B, Smyth GK. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn. Stat Appl Genet Mol Biol. 2010;9:Article 39.

## Text that becomes obsolete
- "1000 random lists … 20,000 random iGEA runs", "Null distributions were generated for 11 collections …", "random gene lists of similar size, in this case 300", "the user-selected raw p-value threshold was matched provided it was ≤ 0.01 … capped at 30", "some parameters were not matched during benchmarking, including minimum overlap and maximum gene set size", and the hosted-app sentence about not "recomputing empirical null distributions".
- Supplementary Figure 6 (cluster sizes of random networks vs list size) can be regenerated with the new engine in ~30 min if kept; it still supports the list-size recommendation (R5-m5, Step 4).
