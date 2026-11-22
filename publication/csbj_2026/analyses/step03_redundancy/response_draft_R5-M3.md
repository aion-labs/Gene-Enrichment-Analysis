# Draft — response to Reviewer 5, major comment 3 (collection redundancy), and Methods/Results text

## Response
We thank the reviewer for this point. The aggregate MSigDB Canonical Pathways collection (C2:CP) was not part of the network analysis; it appeared only in the description of the previous precomputed reference set, which we have now replaced (see response to comment 2). The 10 collections in the network are Hallmark, six C2:CP sub-collections (Reactome, KEGG Legacy, KEGG MEDICUS, BioCarta, PID, WikiPathways) and the three GO collections; no gene set is counted twice through an aggregate.

As requested, we quantified redundancy between the collections that remain. Of 10,445 gene sets within the size limits, 74 (0.7%) are identical to, or have a Jaccard index ≥ 0.9 with, a gene set in another collection (new Supplementary Table X). In the HIV network, four cross-collection term pairs were highly similar (Jaccard 0.81–0.92; GPI-anchor biosynthesis in KEGG MEDICUS, KEGG Legacy and WikiPathways; N-glycan biosynthesis in KEGG Legacy and WikiPathways), sharing 3–4 network genes.

We then repeated the analysis (i) after removing every gene set with Jaccard ≥ 0.9 to a set in a higher-priority collection (37 sets removed), and (ii) leaving out each collection in turn, each time recomputing the random-list reference with the same collections. With deduplicated collections the largest cluster was essentially unchanged (181 genes, 63 terms, 418 edges vs 179, 62, 412) and remained significant (10,000 random lists: p ≤ 0.0001 for cluster size, genes, terms and edges; p = 0.007 for edges per gene). Leaving out any single collection retained a cluster of 145–179 genes, with p ≤ 0.001 (1,000 random lists) for cluster size, genes, terms and edges in all ten cases and p = 0.007–0.033 for edges per gene. Network connectivity therefore does not depend on duplicated annotations or on any single collection.

## Methods sentence
Redundancy between collections was assessed as the number of gene sets (within the size limits) that are identical to, or have a Jaccard index ≥ 0.9 with, a gene set in another collection. For sensitivity analyses, deduplicated collections were built by removing such sets from the lower-priority collection (order: Hallmark, Reactome, KEGG Legacy, KEGG MEDICUS, BioCarta, PID, WikiPathways, GO BP, GO CC, GO MF), and the network and its random-list reference were recomputed.

## Numbers (from `s03_redundancy_matrix_J90.csv`, `s03_sensitivity_dedup.csv`, `s03_sensitivity_loo.csv`)
| Configuration | Cluster genes / terms / edges | p (size, genes, terms, edges) | p (edges/gene) | Random lists |
|---|---|---|---|---|
| All 10 (as published) | 179 / 62 / 412 | ≤ 0.0001 | 0.008 | 10,000 |
| All 10, deduplicated | 181 / 63 / 418 | ≤ 0.0001 | 0.007 | 10,000 |
| Without Hallmark | 178 / 61 / 403 | ≤ 0.001 | 0.013 | 1,000 |
| Without Reactome | 170 / 53 / 348 | ≤ 0.001 | 0.011 | 1,000 |
| Without KEGG Legacy | 178 / 57 / 392 | ≤ 0.001 | 0.007 | 1,000 |
| Without KEGG MEDICUS | 179 / 59 / 398 | ≤ 0.001 | 0.009 | 1,000 |
| Without BioCarta | 179 / 59 / 402 | ≤ 0.001 | 0.011 | 1,000 |
| Without PID | 179 / 60 / 404 | ≤ 0.001 | 0.013 | 1,000 |
| Without WikiPathways | 174 / 56 / 376 | ≤ 0.001 | 0.015 | 1,000 |
| Without GO BP | 149 / 49 / 305 | ≤ 0.001 | 0.033 | 1,000 |
| Without GO CC | 169 / 51 / 346 | ≤ 0.001 | 0.025 | 1,000 |
| Without GO MF | 145 / 50 / 323 | ≤ 0.001 | 0.013 | 1,000 |
