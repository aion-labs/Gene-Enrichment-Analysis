# iGEA pipeline as implemented (code commit 4eefc2f)

Read-only audit, 14 Nov 2026. File/function for each stage. Notes marked **≠ text** differ from the v6 manuscript.

| # | Stage | Implementation | What it does |
|---|---|---|---|
| 1 | Input parsing | `cli.py` / `streamlit_app.py` → `GeneConverter.validate_and_map_symbol`, `GeneSet` | Upper-cases, maps aliases to current NCBI symbols (`Homo_sapiens.gene_info`, 9 Aug 2025), drops duplicates. Hard cap 800 genes. HIV: 272 lines → 272 genes. |
| 2 | Background | `BackgroundGeneSet` | Default `data/backgrounds/all_genes.txt` (23,590 lines); symbol validation keeps 23,391. |
| 3 | Collection preparation | `GeneSetLibrary`, `Enrichment._compute_enrichment` | Keep sets with 10 ≤ size ≤ 600. Collection genes G_c = union of kept sets. Collection background M = |background ∩ G_c|. |
| 4 | Per-term test | `enrichment.compute_pvalue` | Input restricted to collection genes: n = |input ∩ G_c|. k = |input ∩ G_c ∩ term|. K = full term size (**not** restricted to background). 2×2 table [[k, K−k],[n−k, M−K−n+k]]. **Default test: Fisher's exact, two-sided** (`scipy.stats.fisher_exact`) — **≠ text**, which says `hypergeom.sf`. BH-FDR over all tested sets in the collection. |
| 5 | Iteration | `IterativeEnrichment._compute_enrichment` | Each iteration re-runs stage 4 on the remaining genes (M unchanged). Keep terms with k ≥ min overlap (3); take the lowest p; ties keep GMT order (stable sort). Stop if p ≥ threshold (0.01), no genes left, no eligible term, or max iterations reached. Remove the selected term's overlapping genes. Collections are processed independently. |
| 6 | Network | `NetworkConnectivityAnalyzer.add_igea_results` | Bipartite, unweighted: term — gene edges for the genes removed at that term's iteration. Within a collection each gene has ≤ 1 edge; a gene's degree = number of collections that assigned it. |
| 7 | Largest cluster | `_compute_connected_components` (BFS) | Connected components; largest by node count. |
| 8 | Null | `parallel_null_distribution.compute_null_distribution_from_raw_permutations` | Stored random runs (`permutations/merged_permutation_results.tsv`), filtered to: list size = input size **rounded up** to next 50 (**≠ text**, "similar size"/"nearest"), selected collections, p ≤ threshold, iteration ≤ 30. Rebuilds each random network, takes its largest cluster. |
| 9 | Statistics | `network_connectivity_benchmark.benchmark_cluster` | z = (obs − mean)/SD. "Percentile" = Φ(z), i.e. **normal approximation, not empirical**. |
| 10 | Export | `cli.py`, `ui/rendering.py`, `ui/dot_utils.py` | TSV/JSON/DOT, statistical report (Supp Table 5 format). |

## Random lists (stored null) — `tools/generate_permutation_distribution.py`
- 20 sizes (50–1000, step 50) × 1,000 lists; uniform sampling **without replacement** from `all_genes.txt` (23,590 symbols, unvalidated; `skip_validation=True`).
- Seed = permutation index (1–1000), the same seed reused for every size. `random.sample(list(background.genes), size)` samples from a Python **set**, whose order depends on hash randomisation — exact lists are probably not regenerable from the seed alone (to confirm).
- Run parameters: Fisher's exact, min overlap 3, sizes 10–600, max 30 iterations. Script comment says p threshold 0.05, but stored rows all have p < 0.01 and no run has a gap in its iteration sequence, so the stored runs are equivalent to running at p < 0.01.
- Stored per run: Library, Iteration, Term, iteration p, genes removed, list size, file name. Collections stored: the 10 network collections + the aggregate C2:CP (not used for the 10-collection null). No HPO, no STRING.

## Hypergeometric formula for the Methods (R5-m1)
For term t in collection c: P(X ≥ k), X ~ Hypergeometric(M, K, n), with
M = |B ∩ G_c|, K = |t| (full term size), n = |L ∩ G_c|, k = |L ∩ G_c ∩ t|, where B = background, L = remaining input list, G_c = genes in sets of c within the size limits.
p = Σ_{i=k}^{min(K,n)} C(K,i) C(M−K, n−i) / C(M, n).
Current code default is the two-sided Fisher's exact test on the same 2×2 table; switching the HIV analysis to one-sided hypergeometric gives 254 / 87 / 62 (vs 252 / 87 / 62).
