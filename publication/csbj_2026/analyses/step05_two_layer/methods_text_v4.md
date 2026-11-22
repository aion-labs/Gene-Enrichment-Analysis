# Draft text — iGEA: redundancy reduction, empirical multiple-testing correction and network significance (supersedes v3)

15 Nov 2026. Author decisions:
- The iterative step keeps the manuscript's view: it **reduces within-collection redundancy in preparation for the next steps**. It is not presented as a signal or error-control step.
- Statistics at **two levels**, named for what they are: (1) **term level — empirical multiple-testing correction** (empirical FDR / q-value within each collection, default q < 0.10), the counterpart of ORA's Benjamini–Hochberg FDR; (2) **network level — empirical significance of the cross-collection network**. Umbrella term: **two-level empirical significance testing**. The network test is a significance test of the combined network, not a multiple-testing correction, and is not called one.
- ORA is shown at the standard BH FDR < 0.05.

Numbers: 10,000 random lists, seed 2026 (`step05_two_layer/`, `02 Figures/networks/`).

---

## Methods — Iterative gene-set selection (existing paragraph, unchanged in substance)

iGEA is performed separately for each gene-set collection. The process begins with ORA on the full input gene list. The most significant term is selected, the input genes overlapping it are assigned to that term and removed, and ORA is repeated on the remaining genes until no term reaches the selection threshold (nominal p < 0.01 by default). Each gene is therefore attributed to at most one term per collection. The purpose of this step is to reduce within-collection redundancy before the collections are integrated, so that the resulting network is not dominated by repeated, overlapping terms and its structure reflects convergence across collections. The nominal p-value is used only as a selection criterion; it is not interpreted as evidence for individual terms, which is assessed in the next two steps.

## Methods — Two-level empirical significance testing

Because each iterative test is conditional on earlier selections, analytical multiple-testing corrections such as Benjamini–Hochberg do not apply to iGEA terms. iGEA therefore corrects for multiple testing empirically at the term level, within each collection, and then tests the significance of the combined network across collections. Both levels use random gene lists that are processed in exactly the same way as the input.

*Random gene lists.* For each analysis, n random gene lists (default 10,000) with as many genes as the validated input list are drawn without replacement from the analysis background. Each random list is analysed with the same collections and parameters (statistical test, selection threshold, minimum overlap, gene-set size limits and iteration cap) and passes through both levels below exactly as the input does. List i is generated from a fixed seed (default 2026), so all results are reproducible.

*Term level — empirical multiple-testing correction (q-value).* The q-value estimates the false discovery rate among terms at least as strong as a given term, the counterpart of the Benjamini–Hochberg FDR used for ORA, but estimated from random lists rather than from analytical assumptions. For a term selected in collection c with nominal p-value p, the empirical q-value is the expected number of terms selected in c from random lists with p-value ≤ p (mean over random lists) divided by the number of input terms in c with p-value ≤ p, made monotone in p as in the Benjamini–Hochberg procedure. Random lists consist entirely of noise, whereas part of an input list is explained by genuine terms; to avoid over-estimating chance selections, the expected number of chance terms is re-estimated from random lists with as many genes as remain unexplained after removing the genes of accepted terms, repeating until the set of accepted terms is stable. The number of unexplained genes is rounded up to the next point of a size grid, which can only increase the estimate. Terms with q ≥ 0.10 (default) are removed before the network is assembled. Each random list is scored and filtered in the same way, leaving itself out of the reference.

*Network level — significance of the cross-collection network.* The retained terms of all collections are joined through their genes into a gene–term network, and its largest connected cluster is compared with the largest clusters of the random lists' networks, built after the same filtering. For each metric (nodes, genes, terms, edges, edges per gene) we report an empirical p-value, p = (r + 1)/(n + 1), where r is the number of random lists whose largest cluster is at least as large [Phipson & Smyth 2010]; z-scores are given descriptively, because the random distributions are skewed counts. The five metrics describe the same cluster and are not independent tests. A significant result indicates that genes attributed to terms in one collection link to terms in other collections more than chance allows, i.e. convergence across independently curated annotations. The statistics always refer to the collections shown: when a figure displays a subset of collections (for example the GO collections in Figure 2), the random lists are analysed with the same subset.

*Choice of the default q threshold.* We use q < 0.10 rather than the 0.05 conventional for ORA because the term-level correction is followed by the network-level test, which independently requires cross-collection support; the random lists show how little survives both levels (below). Thresholds between 0.05 and 0.25 are common in exploratory enrichment analysis (the GSEA documentation recommends FDR < 0.25 [Subramanian et al. 2005]). The threshold is a user parameter; setting it to 1 keeps all selected terms while still reporting q-values.

References to add: Phipson B, Smyth GK. Stat Appl Genet Mol Biol. 2010;9:Article 39. Subramanian A, et al. PNAS 2005;102:15545–50.

## Results — table (HIV, 10 collections)

| | Selected terms (iterative, p < 0.01) | Term level: retained (q < 0.10) | Network level: largest cluster |
|---|---|---|---|
| HIV dependency factors (272 genes) | 62 | 49 | 49 terms, 162 genes, 342 edges, 9 collections |
| Random lists (mean per list) | 16.2 | 1.3 | 0.9 terms, 4.2 genes |
| Network test | | | p ≤ 0.0001 for nodes, genes, terms and edges; edges per gene p = 0.009 |

Text: "Of the 62 HIV terms selected by iGEA, 49 had q < 0.10 and all 49 formed a single connected cluster of 162 genes spanning 9 collections. Random lists of the same size retained on average 1.3 terms after the term-level correction, and their largest clusters contained 0.9 terms and 4.2 genes; none approached the HIV cluster (p ≤ 0.0001)."

Other datasets: hypoxia 32 of 32 retained, cluster 64 genes / 32 terms, p ≤ 0.0001, edges per gene p = 0.005 (random 5.7 → 1.0 terms per list); vaccination 66 of 75, cluster 267 genes / 66 terms, p ≤ 0.0001, edges per gene p = 0.003 (random 21.7 → 1.3). In all three datasets every retained term belonged to the largest cluster.

Figure 2 (GO only, GO-only random lists): 25 of 33 retained; cluster 130 genes, 25 terms, 196 edges; p ≤ 0.0001 (nodes, genes, terms, edges), edges per gene p = 0.005.

Kept for the response letter / supplement only (not part of the framing): random 272-gene lists give 31.4 ORA terms vs 16.2 iGEA terms at p < 0.01 across the 10 collections.

## Figure 1 legend (R5-m2)

**Figure 1. ORA and iGEA results for GO Biological Process, HIV dependency factors.** (A) ORA: 14 terms significant at Benjamini–Hochberg FDR < 0.05 (57 genes). (B) iGEA: 5 terms retained at empirical q < 0.10 (55 genes). Terms found by both methods are shown in yellow; terms found only by ORA in blue (A) and only by iGEA in red (B); gene symbols label the input genes linked to each term. The two methods use different error measures (analytical FDR over all tested gene sets for ORA; empirical q-value from 10,000 random gene lists processed by the full iterative procedure for iGEA; see Methods), so the panels compare each method at its default rather than at an identical error level. iGEA covers almost the same genes as ORA with about a third of the terms, because each gene is attributed to one term.

(Matched 10 % thresholds, response letter only, R4-1: ORA FDR < 0.10 gives 20 terms / 73 genes, including all 5 iGEA terms.)

## Response-letter paragraph (R2-1, R5-M1, R3-M2)

We agree that network-level significance does not by itself establish that individual terms are reliable, and we have revised both the method and the manuscript. The iterative step is now described explicitly for what it is: redundancy reduction within each collection, with the nominal p-value used only to select terms. Statistical significance is then assessed at two levels, both using random gene lists of the same size, drawn from the same background and processed by the complete iterative procedure with the same collections and parameters. At the term level, multiple testing is corrected empirically: every selected term receives an empirical q-value (an estimate of the false discovery rate); because the reference is generated by the same adaptive procedure, it does not depend on the independence assumptions that make analytical corrections inapplicable to iterative selection. Terms with q ≥ 0.10 are removed before the network is built. At the network level, the largest connected cluster of the network is compared with the clusters of the random lists, which have been filtered identically, using empirical p-values in place of z-scores. For the HIV list, random lists retain on average 1.3 terms after the term-level correction, whereas 49 of 62 HIV terms are retained and all form a single cluster that no random list approaches (p ≤ 0.0001). The 10 % default is justified in the Methods and is a user parameter.
