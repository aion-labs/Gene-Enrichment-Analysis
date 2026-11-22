> **Superseded by `methods_text_v4.md` (15 Nov 2026): the iterative step is redundancy reduction, not a signal-refinement stage.**

# Draft text — iGEA's three-stage signal refinement (supersedes the statistics text in `step02_per_run_null/methods_text_v2.md`)

15 Nov 2026. Author decisions: ORA shown at the standard BH FDR < 0.05; iGEA default q < 0.10, justified explicitly; describe the three stages as signal-to-noise refinement rather than "statistical correction". Numbers: 10,000 random lists, seed 2026 (`step05_two_layer/`, `02 Figures/networks/`).

Proposed name: **three-stage signal refinement** (alternatives: "layered noise reduction", "three levels of evidence"). Stages: (1) iterative gene attribution, (2) collection-level empirical filtering, (3) cross-collection network test.

---

## Methods — Signal refinement in iGEA

iGEA does not apply a single multiple-testing correction to a fixed list of tests. Its terms are selected adaptively, and the evidence that a selected term reflects the input rather than chance accumulates over three stages. Stage 1 reduces redundancy and limits chance selections within each collection; stage 2 estimates, for every selected term, how often terms of that strength arise from random gene lists; stage 3 asks whether the terms that remain, taken across collections, form a network that random lists do not produce. All three stages are applied in exactly the same way to the input list and to every random list used as reference.

**Stage 1 — Iterative gene attribution (within each collection).** In each iteration, the gene set with the lowest p-value is selected and its overlapping genes are removed before the next iteration, so each gene is attributed to at most one term per collection and each new term must be supported by genes not already explained. This is not an error-control procedure and the nominal p-value (p < 0.01) serves only as a selection criterion. It does, however, reduce chance selections as well as redundancy: random gene lists of the HIV input size (272 genes) yielded on average 31.4 ORA terms at p < 0.01 across the ten collections, but 16.2 iGEA terms.

**Stage 2 — Collection-level empirical q-value.** Because each test is conditional on earlier selections, analytical corrections such as Benjamini–Hochberg do not apply to iGEA terms. Instead, n random gene lists (default 10,000) with as many genes as the validated input are drawn without replacement from the analysis background and processed by the complete iterative procedure with the same collections and parameters. For a term selected in collection c with nominal p-value p, the empirical q-value is the expected number of terms selected by chance in c with p-value ≤ p (mean over random lists) divided by the number of input terms in c with p-value ≤ p, made monotone in p as in the Benjamini–Hochberg procedure. Because the random lists consist entirely of noise, whereas part of the input list is explained by genuine terms, this estimate is conservative. We therefore re-estimate the expected number of chance terms from random lists with as many genes as remain unexplained after removing the genes of accepted terms (q < threshold), repeating until the set of accepted terms is stable; the number of unexplained genes is rounded up to the next point of a size grid, which can only increase the estimated number of chance terms. Terms with q ≥ 0.10 (default) are removed before the network is assembled. Each random list is scored and filtered in the same way, leaving itself out of the reference.

**Stage 3 — Cross-collection network test.** The retained terms of all collections are joined through their genes into a gene–term network, and the largest connected cluster is compared with the largest clusters of the random lists, which have passed through stages 1 and 2 identically. For each cluster metric (nodes, genes, terms, edges, edges per gene) we report an empirical p-value, p = (r + 1)/(n + 1), where r is the number of random lists whose largest cluster is at least as large [Phipson & Smyth 2010]; z-scores are given descriptively. Network-level significance reflects convergence across independently curated collections: genes attributed to terms in one collection link to terms in others more often than chance allows. Because collections differ in content, the statistics always refer to the collections actually shown; when a figure displays a subset (for example the GO collections in Figure 2), the random lists are analysed with the same subset.

**Choice of the default q threshold.** We use q < 0.10 rather than the 0.05 conventional for ORA because stage 2 is not the final filter: a 10% expected false-term rate at the collection level is followed by a network-level test that independently requires cross-collection support, and the random lists show how little chance material survives both stages. In random lists of 272 genes the q filter reduced selected terms from 16.2 to 1.3 per list, and their largest clusters contained 0.9 terms on average, whereas 49 of the 62 HIV terms passed and all 49 formed a single connected cluster. Thresholds of 0.05–0.25 are common in exploratory enrichment analysis (for example, the GSEA documentation recommends FDR < 0.25), and the threshold is a user parameter; setting it to 1 retains all selected terms while still reporting q-values.

Reference to add: Phipson B, Smyth GK. Stat Appl Genet Mol Biol. 2010;9:Article 39. (GSEA FDR 0.25: Subramanian A, et al. PNAS 2005;102:15545–50.)

## Results — one table that shows the three stages (HIV, 10 collections)

| | ORA (p < 0.01) | Stage 1: iGEA terms | Stage 2: q < 0.10 | Stage 3: largest cluster |
|---|---|---|---|---|
| HIV dependency factors (272 genes) | 252 terms | 62 terms | 49 terms | 49 terms, 162 genes, 342 edges, 9 collections |
| Random lists, mean per list | 31.4 terms | 16.2 terms | 1.3 terms | 0.9 terms, 4.2 genes |
| Network test | | | | p ≤ 0.0001 for nodes, genes, terms and edges; edges per gene p = 0.009 |

Note for the text: term counts at stage 1 should not be read as a signal-to-noise ratio — ORA's 252 terms are highly redundant (mean Jaccard 0.33 among Reactome FDR terms); iGEA's are not.

Other datasets (stages 2 and 3): hypoxia 32 of 32 terms retained, cluster 64 genes / 32 terms, p ≤ 0.0001, edges per gene p = 0.005 (random: 5.7 → 1.0 terms per list); vaccination 66 of 75, cluster 267 genes / 66 terms, p ≤ 0.0001, edges per gene p = 0.003 (random: 21.7 → 1.3). In all three datasets every retained term belonged to the largest cluster.

Figure 2 (GO collections only, statistics with GO-only random lists): 25 of 33 terms retained; cluster 130 genes, 25 terms, 196 edges; p ≤ 0.0001 (nodes, genes, terms, edges), edges per gene p = 0.005.

## Figure 1 legend (R5-m2: state the different criteria)

**Figure 1. ORA and iGEA results for GO Biological Process, HIV dependency factors.** (A) ORA: 14 terms significant at Benjamini–Hochberg FDR < 0.05 (57 genes). (B) iGEA: 5 terms retained at empirical q < 0.10 (55 genes). Terms found by both methods are shown in yellow; terms found only by ORA in blue (A) and only by iGEA in red (B); gene symbols label the input genes linked to each term. The two methods use different error measures (analytical FDR over all tested gene sets for ORA; empirical q-value from 10,000 random gene lists processed by the full iterative procedure for iGEA, see Methods), so the panels compare outputs at each method's default, not at an identical error level. iGEA covers almost the same genes as ORA with about a third of the terms, because each gene is attributed to one term.

(Counts at matched 10 % thresholds for the response letter, R4-1: ORA FDR < 0.10 gives 20 terms / 73 genes, which include all 5 iGEA terms.)

## Response-letter paragraph (R2-1, R5-M1, R3-M2)

We agree that network-level significance does not by itself establish that individual terms are reliable, and we have revised both the method and the manuscript. iGEA now assesses every selected term with an empirical q-value obtained by running the complete iterative procedure — including gene removal — on random gene lists of the same size, drawn from the same background and analysed with the same collections and parameters. Because the reference is generated by the same adaptive procedure, the q-value does not rely on the independence assumptions that make analytical corrections inapplicable to iterative selection. Terms with q ≥ 0.10 are removed before the network is built, and the same filter is applied to every random list, so the network-level test that follows compares like with like. We describe the result as three stages of signal refinement (iterative attribution, collection-level q-values and a cross-collection network test; new Methods section and Table X). For the HIV list, random lists of the same size yield 16.2 chance terms on average after stage 1 but 1.3 after stage 2, whereas 49 of 62 HIV terms are retained and all form a single cluster that no random list approaches (p ≤ 0.0001). Empirical p-values replace z-scores for inference. We justify the 10 % default explicitly and make it a user parameter.
