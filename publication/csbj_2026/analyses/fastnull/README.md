# Fast random-list null — prototype (14 Nov 2026)

Prototype only; not in the iGEA repo. Files: `fast_igea.py` (engine), `validate_speed.py` (equivalence + timing vs iGEA code), `fast_null.py` (null generator), CSV outputs.

## Idea
Same algorithm and statistics as iGEA, computed differently:
1. Each collection is loaded once as a sparse gene-set × gene matrix.
2. Each iteration computes the overlap with every gene set in one sparse matrix–vector product.
3. Exact two-sided Fisher p-values are computed only for sets with overlap ≥ the minimum overlap (vectorised, log-space); the best candidate and any near-ties (within 1e-6) are recomputed with `scipy.stats.fisher_exact`, so the selected term and its p-value are exactly what the iGEA code produces.
4. Ties resolved by GMT order (as the stable sort in iGEA). Gene-set size K = unique genes; size filter on raw GMT length (as iGEA; 46 GO sets have duplicate entries in MSigDB v2025.1).
5. Random lists drawn with a seeded NumPy generator from a sorted background, so every list is exactly reproducible (the stored null's lists are not).

## Validation
- HIV + 6 random lists (300 and 1,000 genes) × 10 collections = 70 iGEA runs: **identical terms, order, removed genes, and bit-identical p-values** (max relative difference 0).
- New 1,000-list null at size 300 vs the stored null: distributions indistinguishable (KS p = 0.94, 0.83, 0.91, 0.89, 0.47 for cluster size, genes, terms, edges, edges/gene); HIV z 7.74 / 7.32 / 8.62 / 11.37 / 3.55 vs stored 7.79 / 7.38 / 8.59 / 11.45 / 3.62.

## Speed (Linux VM on the author's MacBook Pro, 4 cores)
| Task | Current code | Fast engine |
|---|---|---|
| One random list, 300 genes, 10 collections (1 core) | 6.3 s | ≈ 0.03 s |
| One random list, 1,000 genes (1 core) | 24 s | ≈ 0.8 s |
| 1,000 lists of 300 (4 cores) | ≈ 26 min (est.) | 6.9 s |
| 1,000 lists of 1,000 (4 cores) | ≈ 100 min (est.) | ≈ 3 min |
| 10,000 lists of 272 (exact HIV size, 4 cores) | ≈ 4–5 h (est.) | 55 s |
| Full 20,000-list reference null | ≈ 80 core-hours (est.) | ≈ 30 min on 4 cores (est.) |

## HIV against an exact-size null (272 genes, 10,000 lists)
| Metric | Random mean ± SD (max) | z | Random ≥ HIV | Empirical p |
|---|---|---|---|---|
| Cluster size | 42.8 ± 23.7 (133) | 8.36 | 0 | ≤ 0.0001 |
| Genes | 33.0 ± 18.2 (104) | 8.01 | 0 | ≤ 0.0001 |
| Terms | 9.8 ± 5.9 (35) | 8.93 | 0 | ≤ 0.0001 |
| Edges | 49.3 ± 29.9 (180) | 12.13 | 0 | ≤ 0.0001 |
| Edges per gene | 1.46 ± 0.25 (3.75) | 3.34 | 84 | 0.0085 |
