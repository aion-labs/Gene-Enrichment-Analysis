# Reviewer analyses — results log

Code: `iterative_enrichment` @ 4eefc2f (read-only). Annotations: MSigDB v2025.1.Hs. Symbol mapping: NCBI gene_info 9 Aug 2025.
Environment: Linux VM on author's Mac, Python 3.10.12, pandas 2.2.3, SciPy 1.15.3. Scripts and outputs in `stepNN/`.

---

## Step 1 — Code audit, reproduction, collection inventory (14 Nov 2026)
Answers: R4-3, R5-M9, R5-M10, R5-m1, R2-1 (sampling), R3-m3 (max possible terms). Files: `step01/PIPELINE.md`, `step01/results/`.

### 1. Reproduction of the HIV analysis
| Quantity | Paper (v5/v6) | Reproduced | Match |
|---|---|---|---|
| Input genes | 272 | 272 | ✓ |
| ORA terms p < 0.01 (overlap ≥ 3), 10 collections | 252 | 252 | ✓ |
| … of which FDR < 0.05 | 87 | 87 | ✓ |
| iGEA terms | 62 | 62 | ✓ |
| Largest cluster | 241 nodes, 179 genes, 62 terms, 412 edges, 2.3017 edges/gene | identical | ✓ |
| Null list size | 300 | 300 | ✓ |
| z: cluster size / genes / terms / edges / edges per gene | 7.74 / 7.39 / 8.44 / 11.30 / **3.00** | **7.79 / 7.38 / 8.58 / 11.45 / 3.60** | ✗ |

**Discrepancy — z-scores.** The observed network reproduces exactly. The null does not: the shipped code with the stored random runs, restricted to the 10 network collections, gives null means 46.43 / 36.00 / 10.42 / 53.54 / 1.459 — identical to the **v3** report (Feb 2026; z 7.79 / 7.38 / 8.58 / 11.45 / 3.63), which already used exactly these 10 collections. The v5/v6 Supp Table 5 (regenerated 30 Aug 2026) has null means 37.95 / 29.38 / 8.57 / 43.51 / 1.424, which no tested configuration reproduces (closest: dropping GO CC or GO MF from the null, 37.99 / 29.49). The 30 Aug change report's explanation ("nulls rebuilt from 10-collection random networks") is not borne out — the null was already 10-collection. **Recommendation:** use the reproducible values in v7, report empirical p-values (Step 3), and note the correction in the response letter. Conclusions are unaffected; edges/gene becomes stronger (3.60 vs 3.00).

### 2. Corrections the Methods need (code ≠ text)
1. **Test:** code (app, CLI, and random runs) uses two-sided Fisher's exact test; text says `hypergeom.sf`. One-sided hypergeometric gives 254 / 87 / 62 — negligible change. Either describe Fisher, or switch and re-run (would also change the null). Recommend: describe what was run.
2. **"Percentile"** in the statistical report is Φ(z), not an empirical percentile. Replace with empirical p (Step 3).
3. **Null bin:** list size is rounded **up** to the next multiple of 50 (272 → 300), not "nearest".
4. **Term size K** is the full term size, not restricted to the background; the input list is restricted to the collection's genes, not to the background (R5-m1, R5-M10).
5. **Random lists:** uniform without replacement from `all_genes.txt` (23,590 symbols, unvalidated), whereas the real run uses the validated background (23,391). Seed = list index, reused across sizes; exact lists likely not regenerable (set ordering). Stored runs are nonetheless complete and fixed — reproducibility rests on the stored file, which should be deposited.
6. **Max iterations:** app/CLI default is 10; the HIV run reached 12 (GO BP), so the paper run used ≥ 12. State the value used (null cap 30).

### 3. Stored-null integrity
611,110 rows; 11 collections (10 network + aggregate C2:CP); no HPO, no STRING; all p < 0.01; max 30 iterations; **no run has a gap** in its iteration sequence, so filtering by p is equivalent to the real stopping rule. The aggregate C2:CP is **excluded** from the 10-collection null (filtered by name) — supports the R5-M3 rebuttal.

### 4. Collection inventory (R4-3, R5-M9) → `s01_collection_inventory.csv`
All 10 collections: MSigDB v2025.1.Hs; used in the ORA comparison, the iGEA network and the random-list null. Max possible terms (sets 10–600 with ≥ 3 HIV genes) per collection: 22, 228, 47, 11, 10, 33, 113, 793, 136, 115 — **1,508 in total** (R3-m3), vs 252 ORA (p < 0.01), 87 ORA (FDR < 0.05), 62 iGEA. Figure/LLM-usage columns to be completed with the author.

### 5. ORA vs iGEA recount (replaces "5 identical / 4 additional / 2 limited") → `s01_recount_ora_vs_igea.csv`
| Category | n | Collections |
|---|---|---|
| Identical | 1 | Hallmark |
| All FDR terms captured, iGEA adds terms | 3 | KEGG MEDICUS, GO CC, GO MF |
| No FDR-significant ORA term; iGEA selects terms | 4 | KEGG Legacy, BioCarta, PID, WikiPathways |
| Limited overlap | 2 | Reactome (65 FDR vs 9 iGEA, 7 shared), GO BP (14 vs 12, 4 shared) |
The old "5 identical" cannot be reproduced on 10 collections under any reading; the text must change.

### 6. Parameter derivability (R5-M10, R3-M1) → `s01_parameter_derivability.csv`
Derivable from stored nulls: collection subsets; **stricter** p thresholds (truncate each run — not yet implemented; corrects the earlier assumption that this was impossible); max iterations ≤ 30. Not derivable: looser p, min overlap, size limits, background, test, exact list size.

### 7. Decisions and empirical p-values (14 Nov 2026, after review)
- **Test:** keep Fisher's exact (two-sided), mention the hypergeometric option. Verified: hypergeometric gives the same 62 iGEA terms (identical in 9 collections; WikiPathways same 6 terms in a different order), 254 vs 252 ORA terms (+GOCC CILIARY BASAL BODY, +GOMF TRANSCRIPTION COREGULATOR ACTIVITY), the same 87 at FDR < 0.05.
- **Statistics:** report the reproducible z-scores plus empirical p = (r+1)/(n+1) (`s01c_empirical_p.py` → `s01_hiv_z_and_empirical_p.csv`; per-run null values in `s01_null_runs_size300.csv`, captured from the code's own null builder).

| Metric | Obs | Null mean ± SD (max) | z | r / 1000 | Emp. p | Skew |
|---|---|---|---|---|---|---|
| Cluster size | 241 | 46.4 ± 25.0 (124) | 7.79 | 0 | 0.001 | 0.59 |
| Genes | 179 | 36.0 ± 19.4 (98) | 7.38 | 0 | 0.001 | 0.62 |
| Terms | 62 | 10.4 ± 6.0 (33) | 8.59 | 0 | 0.001 | 0.67 |
| Edges | 412 | 53.5 ± 31.3 (159) | 11.45 | 0 | 0.001 | 0.72 |
| Edges/gene | 2.30 | 1.46 ± 0.23 (3.31) | 3.62 | 4 | 0.005 | 1.29 |

- **Minor non-determinism:** null means for terms and edges/gene vary in the 3rd decimal between runs (z terms 8.58–8.59, edges/gene 3.60–3.63) because ties for the largest random cluster are broken by set order. Report those z to one decimal; fix in Step 12.
- Draft Methods/Results text: `step01/methods_text_drafts.md`.

---

## Fast null — prototype (14 Nov 2026) → `fastnull/README.md`
Vectorised re-implementation of the iGEA selection step. Bit-identical to the iGEA code on 70 test runs; a fresh 1,000-list null at size 300 matches the stored null (KS p 0.47–0.94). 1,000 random lists of 300 genes take 7 s on 4 cores instead of ~26 min; 10,000 exact-size lists for HIV take 55 s. Makes per-analysis, parameter-matched, exact-size, any-background nulls practical (R3-M1, R5-M2, R5-M10, R2-1). Exact-size (272) null, 10,000 lists: HIV empirical p ≤ 0.0001 for size/genes/terms/edges; edges per gene p = 0.0085.

---

## Per-run null — implemented (14 Nov 2026)
Branch `csbj-revision` (from tag `csbj-v6-submitted` = 4eefc2f):
- `f41ccf9` fix: deterministic choice of the largest cluster (ties no longer depend on hash order)
- `b9bd84b` feat: `code/fast_igea.py` (engine), `code/per_run_null.py` (null + statistics), `test/test_fast_igea.py` (4 tests pass: HIV Fisher + hypergeometric, random 300/600 lists, reproducibility n_jobs 1 vs 2)
- `8fb0b6b` feat(cli): `--benchmark` now uses the per-run null by default; `--null-lists`, `--null-seed`, `--null-jobs`, `--null-source per-run|stored`
Not yet: Streamlit app, README/CLI docs, `click<8.2` pin in requirements.

HIV via the CLI (272 validated genes, background all_genes validated = 23,391, 10 collections, Fisher, p < 0.01, min overlap 3, 10–600, no iteration cap) → `step02_per_run_null/`:

| Metric | Obs | 1,000 lists: mean ± SD (max) | z | r | p | 10,000 lists: mean ± SD (max) | z | r | p |
|---|---|---|---|---|---|---|---|---|---|
| Cluster size | 241 | 42.5 ± 23.4 (125) | 8.47 | 0 | 0.001 | 42.9 ± 23.6 (148) | 8.39 | 0 | 0.0001 |
| Genes | 179 | 32.7 ± 18.0 (99) | 8.12 | 0 | 0.001 | 33.1 ± 18.2 (112) | 8.04 | 0 | 0.0001 |
| Terms | 62 | 9.7 ± 5.8 (29) | 8.96 | 0 | 0.001 | 9.7 ± 5.8 (41) | 8.97 | 0 | 0.0001 |
| Edges | 412 | 48.7 ± 29.6 (160) | 12.29 | 0 | 0.001 | 49.1 ± 29.7 (195) | 12.24 | 0 | 0.0001 |
| Edges per gene | 2.30 | 1.45 ± 0.26 (3.38) | 3.31 | 11 | 0.012 | 1.45 ± 0.25 (5.0) | 3.41 | 80 | 0.0081 |
Runtime: 7.7 s (1,000) and 61 s (10,000) on 4 cores.

### App and defaults (14 Nov 2026)
- `4f6e129` (csbj-revision): Streamlit section replaced by per-run statistics (button, progress bar, n lists default **10,000**, seed, empirical p explanation, downloads); per_run_null/CLI default 10,000; `click<8.2` pinned; README section. AppTest: app renders and runs end-to-end without exceptions; n_jobs 1 and 2 give working results.
- `ba15721` (new branch `csbj-revision-cloud`, from `origin/streamlit-cloud`): same feature for the hosted app, configured for Community Cloud limits (0.078–2 CPUs, 690 MB–2.7 GB; docs.streamlit.io): in-process (n_jobs = 1), default 1,000 lists, max 10,000. AppTest end-to-end OK (example list: 1,000 lists in 2 s). **Not pushed** — pushing `csbj-revision-cloud` onto `streamlit-cloud` redeploys the hosted app.
- Note: running the app rewrites `data/libraries/alias.json` (formatting only); restored, not committed.

## Hypoxia and vaccination with per-run null (14 Nov 2026) → `step02_per_run_null/`
- `ca47acc` fix(cli): libraries with no selected terms were left out of the null (hypoxia: BioCarta). Found because chunked and CLI nulls disagreed on 3/200 lists; after the fix they agree 200/200. HIV unaffected (all 10 collections had terms). App was already correct.
- Hypoxia (75 genes; CLI, 10,000 lists, 10 s): cluster 96 nodes / 64 genes / 32 terms / 171 edges / 9 collections; z 10.59/10.04/11.10/16.22/3.29; p ≤ 0.0001 ×4; edges/gene p = 0.011.
- Vaccination (426 genes, MSigDB C7 NAKAYA … HAI 28DY RESPONSE AT 3DY POSITIVE; 10,000 lists in 3 chunks via `null_chunk.py`, ~4 min): cluster 347 / 272 / 75 / 702 / 9 collections; z 8.65/8.49/8.66/15.14/5.64; p ≤ 0.0001 ×4; edges/gene p = 0.0013.
- The old ">3 z-scores" claim is now documented and stronger. Draft text: `step02_per_run_null/methods_text_v2.md`.

## Step 2 — Collection redundancy (R5-M3) (14 Nov 2026) → `step03_redundancy/`
- 10,445 gene sets (10–600) across the 10 collections; **74 (0.7%) have Jaccard ≥ 0.9 with a set in another collection, 12 identical**. Matrix: `s03_redundancy_matrix_J90.csv`. Largest overlaps: GO MF↔GO BP (13), WikiPathways↔KEGG Legacy (5).
- HIV network: 4 cross-collection term pairs with J ≥ 0.8 (GPI anchor ×3, N-glycan), 3–4 shared network genes each (`s03_hiv_similar_term_pairs.csv`).
- Deduplicated collections (37 sets removed by priority order): cluster 181 / 63 / 418, p ≤ 0.0001 ×4, edges/gene p = 0.0072 (10,000 lists).
- Leave-one-collection-out (1,000 lists each): cluster 145–179 genes; p ≤ 0.001 for size/genes/terms/edges in all 10; edges/gene p 0.007–0.033.
- Conclusion: connectivity does not depend on duplicated annotations or any single collection. Draft response + Methods: `step03_redundancy/response_draft_R5-M3.md`.

## Step 4 (in progress) — term-level error, HIV × Reactome worked example (14 Nov 2026) → `step04_term_fdr/`
Two calibrations, both running the iterative procedure itself on random lists (no analytical correction):
- **A. Pooled empirical q** (whole-run null, same 10,000 random lists as the network null): q = mean random terms with p' ≤ p ÷ real terms with p' ≤ p, made monotone. Reactome q: 0, 0, 0, 0.034, 0.047, 0.082, 0.082, 0.099, 0.225. Random lists select 2.27 Reactome terms on average.
- **B. Step-wise conditional p**: at iteration k, 10,000 random lists of n_k genes (= remaining input genes) from background minus genes already removed; statistic = best p over all sets; p_step = (r+1)/(N+1). Reactome: 0.0001, 0.0002, 0.0003, 0.122, 0.200, 0.384, 0.388, 0.488, 0.806. Median random best p ≈ 0.002–0.003 at every step.
- ForwardStop (G'Sell et al. 2016) on B: α = 0.05 → first 4 terms; α = 0.10 → first 5. A at q < 0.05 → 5 terms.
- **ORA comparison (same random lists, HIV × Reactome)** → `step04_term_fdr/hiv_reactome_ora_calibration.tsv`: 85 ORA terms at p < 0.01 (overlap ≥ 3). BH FDR < 0.05: 65; empirical q < 0.05: 83; minP/Westfall–Young FWER < 0.05: 46. Random lists give 4.41 ORA terms per list (vs 2.27 iGEA terms). The step-wise iGEA test is the iterative analogue of minP (family-wise), which explains its strictness; ORA loses 39 of 85 terms under the same family-wise logic. Per-term FDR is inflated by redundancy in ORA (85 overlapping terms enlarge the denominator). HIV transcription elongation is not significant in ORA (p > 0.01); it reaches p = 0.0085 in iGEA only after removal.
- **All 10 collections, HIV: ORA FDR < 0.05 vs iGEA q < 0.05** (10,000 random lists, 62 s) → `hiv_all_ora_vs_igea_summary.tsv`, `hiv_all_ora_vs_igea_terms.tsv`. Terms: ORA 87; iGEA 62 total, 25 at q < 0.05 (14 at 0.05–0.10, 18 at 0.10–0.25, 5 > 0.25); both 16; iGEA-only 9; ORA-only 71. Distinct genes: ORA 123, iGEA q<0.05 101, shared 98; 25 ORA-only genes (18 sit in iGEA terms with q ≥ 0.05; 7 in none, incl. CXCR4, NMT1). ORA FDR terms with < 50% of genes in kept iGEA terms: 8 within the same collection, **2 across the network** (Reactome IL-12 signaling ×2); at iGEA q < 0.10: 39 terms, 137 genes, 110/123 ORA genes, 0 ORA terms uncovered. Gained (iGEA q<0.05, not ORA FDR): KEGG pancreatic cancer (EGF/EGFR/AKT1/RELA…), BioCarta cardiac EGF, BioCarta Rab, PID ARF3 (AP2M1, ARF1, CD4, GOSR2), GOBP DNA-templated transcription elongation (Mediator), GOCC Golgi transport complex (COG2–4), GOCC ER protein complex (OST/GPI), GOCC ERGIC membrane, GOCC nuclear exosome. Mean pairwise Jaccard among ORA FDR terms: Reactome 0.33, GO BP 0.15, GO CC 0.15 (iGEA assigned genes are disjoint by design — recompute on full overlaps for R5-m4).
- **Chance terms per random list (272 genes, 10 collections, p < 0.01):** ORA 31.4 (median 28, 95th pct 60) vs iGEA 16.2 (median 16, 95th 25) → `chance_terms_ora_vs_igea.tsv`. GO BP 16.9 vs 6.4; Reactome 4.4 vs 2.2. Iteration roughly halves chance terms but does not remove them.
- **Signal-adjusted q (π0-style), HIV:** genes assigned to accepted terms treated as signal; chance terms re-estimated from random lists of the remaining genes (background minus signal), iterated to convergence (3,000 lists/round). Converged at 142 noise genes, 9.8 chance terms per list (vs 16.2). Terms: q < 0.05: 36 (vs 25 standard), q < 0.10: 50 (vs 39), q < 0.25: 62 (vs 57). Heuristic — needs validation by planted-signal simulation before adoption. → `hiv_q_signal_adjusted_alpha0.05.tsv`.

## Two-layer statistics implemented (15 Nov 2026) → `step05_two_layer/`
Author decision: no planted-signal simulation; signal-adjusted q < 0.10 applied as a term filter to the input AND to every random list before networks are built (collection level, then network level).
- Code: `f0851ca` (csbj-revision), `7c9fb55` (csbj-revision-cloud; hosted defaults 1,000 lists, 200 per size step, in-process). Tests: 5 pass (incl. reproducibility, monotone q, filter applied to random lists).
- Cost: HIV 10,000 lists + adjustment 97 s vs 61 s without (≈1.6×); hypoxia 18 s; vaccination ≈5 min (chunked, `run_two_layer_chunked.py`).

| Dataset | Terms kept (q < 0.10) | Noise genes (final) | Chance terms/random list (before → after filter) | Cluster genes / terms / edges / coll. | p size, genes, terms, edges | Edges/gene p |
|---|---|---|---|---|---|---|
| HIV (272) | 49 / 62 | 110 | 16.2 → 1.26 | 162 / 49 / 342 / 9 | ≤ 0.0001 | 0.009 |
| Hypoxia (75) | 32 / 32 | 11 | 5.7 → 0.97 | 64 / 32 / 171 / 9 | ≤ 0.0001 | 0.005 |
| Vaccination (426) | 66 / 75 | 159 | 21.7 → ~0.9 | 267 / 66 / 662 / 9 | ≤ 0.0001 | 0.003 |

HIV dropped terms: all 6 WikiPathways (q = 0.101, just above threshold), Reactome HIV transcription elongation (0.118), GO BP positive regulation by host of viral process (0.117), protein polyubiquitination, RNP complex biogenesis, positive regulation of IFN-β production, N-linked glycosylation via asparagine, Arp2/3 actin nucleation (0.219). Z-scores become very large because filtered random networks are tiny (mean < 1 term in the largest cluster) — report empirical p as the inference; z only as description.

## Network from filtered terms by default (15 Nov 2026)
- `8cac536` (csbj-revision), `dc687cf` (cloud): CLI runs statistics by default before assembling the network; DOT + AI prompt contain only q < 0.10 terms; combined table lists all terms with q and "In network"; q = 1 keeps all (q still reported); `--no-benchmark` = unfiltered, no statistics. App: statistics auto-run after each analysis and rebuild `iter_dot` from kept terms.
- `5b030f5` fix: full app "Generate Network" crashed (render_network signature mismatch, present since csbj-v6-submitted).
- Note: with 1,000 random lists (hosted default) HIV keeps 55/62 terms vs 49/62 with 10,000 — q near 0.10 (WikiPathways at 0.101) is sensitive to the number of random lists. Paper uses 10,000.

## Framing decisions (15 Nov 2026)
- Figure 1A: ORA at the standard BH FDR < 0.05 (`02 Figures/networks/Fig1A_ORA_GOBP_FDR0.05.dot`); the FDR 0.10 version is for the response letter only.
- iGEA described as **three-stage signal refinement** (not "statistical correction"): (1) iterative gene attribution, (2) collection-level empirical q (default 0.10, justified), (3) cross-collection network test. Funnel (HIV / random per list): ORA 252 / 31.4 → iGEA 62 / 16.2 → q<0.10 49 / 1.3 → largest cluster 49 terms / 0.9. All retained terms are in the largest cluster in HIV, hypoxia and vaccination. Draft: `step05_two_layer/methods_text_v3_three_stages.md`.
- **Framing revised (15 Nov 2026, author):** the iterative step is **redundancy reduction in preparation for the next steps** (manuscript's existing view), not a signal stage. Signal is assessed at **two levels**: within each collection (empirical q, default 0.10) and across collections (network test). Chance-term reduction by iteration (31.4 ORA vs 16.2 iGEA) kept as data only. Draft: `step05_two_layer/methods_text_v4.md` (v3 marked superseded).
- **Naming (15 Nov 2026, author):** term level = empirical multiple-testing correction (empirical FDR / q-value); network level = empirical significance of the cross-collection network; umbrella = two-level empirical significance testing. Updated in `methods_text_v4.md`.

## Triage (15 Nov 2026, author)
Keep: LLM experiments (re-run), stability (90 % subsamples), light communities/hub analysis, small code fixes + Supp Tables 3–5. Text only: 800-gene advice (per-run statistics are calibrated at any size; limit becomes a run-time note). Decline with reasons: topGO / clusterProfiler / REVIGO / EnrichmentMap benchmarks (GO-hierarchy or GO-semantic tools; would test within-GO redundancy reduction, not iGEA's contribution per R2-3), Bushman benchmark, annotation-matched null, alternative assignment schemes, planted-signal simulation. **Cite every tool reviewers named** (topGO elim/weight, parent-child, iGA, clusterProfiler simplify, REVIGO, EnrichmentMap, DAVID, + relevant CSBJ papers from R3) in Introduction/Discussion, with a short method-category table. Annotation-loss (R5-M5) answered from the ORA comparison (25 genes; IL-12; CXCR4).
- LLM context check: filtered GO network (25 terms) has no context-suggestive labels; full 10-collection network keeps REACTOME IL-12 signaling and KEGG MEDICUS EGF-EGFR-PI3K-NFKB (immune, not viral).

## LLM demonstration package (15 Nov 2026) → `llm_package/`
Author decision: middle path — LLM shown as a demonstration of turning the network into natural language, not as an accuracy benchmark; provenance (model, date, full prompt/response), repeats, traceability check, random-list safeguard. Provider: **OpenAI API** (author's key, run in Terminal; `run_llm.py`, stdlib only, Responses API, saves full JSON + exact returned model).
17 calls: HIV GO iGEA ×3 (25 terms), HIV GO ORA FDR<0.05 ×3 (19 terms, same prompt/format), HIV 10-coll iGEA ×3 (49), random 272 genes unfiltered ×3 (16 terms) and filtered ×3 (2 terms; network p = 0.27, not significant), hypoxia ×1, vaccination ×1. Prompt = Supp Table 2 updated to 10 collections + traceability instruction + "say so if too weak". Graph named `enrichment_network` (no dataset names; hypoxia contains genuine HALLMARK: HYPOXIA term). Random control list: seed 999 (`work/random_list_stats.py`).
Note for code tidy-up: the app's AI prompt still lists HPO/Protein Interaction and derives `library` from node colour (root of the "Unknown Library" bug) — align with this prompt and write the library attribute directly.

## Stability under input perturbation (R1-1, R5-M5) (15 Nov 2026) → `step06_stability/`
100 random 90 % subsamples of HIV (245 genes, seed 2026), full two-level pipeline (q < 0.10), shared reference of 5,000 random 245-gene lists (+500 per size step).
- **Network level: stable.** 100/100 subsamples p ≤ 0.0002 (cluster genes). Cluster 142 genes median (77–173); gene Jaccard with the full-list cluster median 0.74 (0.43–0.87).
- **Term level: moderately stable.** Kept terms per subsample 47 (33–62). Of the 49 full-list terms, 34 recovered on average (70 %); kept-set Jaccard vs full median 0.55. Term recurrence: median 0.73; 43/49 in ≥ 50 % of subsamples, 16 in ≥ 80 %, 6 in < 50 %.
- **Near ties:** 34 of 62 full-run iterations had a runner-up within 2× of the selected p; near-tie terms recur less (median 0.54 vs 0.74).
- **Terms only in subsamples** (13 per subsample, 117 distinct): mostly the six WikiPathways terms that sit at q = 0.101 in the full run (threshold effect), plus alternative ER–Golgi/vesicle terms (intra-Golgi transport, ERGIC, vesicle tethering complex). Gene-set Jaccard with the nearest full-list term median 0.14.
- Reading: path dependence affects which specific terms are reported (esp. near ties and terms near q = 0.10), not the significance or gene composition of the network.

## LLM interpretation run (15 Nov 2026) — llm_package/

**Run.** gpt-6-astra via OpenAI Responses API (returned model = requested), 17 calls, full request/response JSON + manifest in `llm_package/responses/`. Input networks are neutral (graph name `enrichment_network`, no dataset names).

**Traceability** (`work/traceability.py` → `work/traceability.tsv`; gene symbols = NCBI human symbols, abbreviations/collection codes excluded):
| Job | Reps | Genes mentioned in network | Input genes covered | Input terms cited | [outside] markers |
|---|---|---|---|---|---|
| hiv_go_igea | 3 | 100% (88–89) | 68% | 96–100% | 3 |
| hiv_go_ora | 3 | 100% (48–75) | 61–95% | 100% | 4–6 |
| hiv_all_igea | 3 | 100% (73–89) | 45–55% | 69–88% | 4 |
| random filtered | 3 | 100% (6) | 100% | 100% | 6 |
| random unfiltered | 3 | 100%* (46–48) | 96–100% | 100% | 2–4 |
| hypoxia | 1 | 100% (35) | 55% | 53% | 3 |
| vaccination | 1 | 100% (73) | 27% | 52% | 5 |
*Only out-of-network symbols across all 17 responses: ESR1/ESR2 (2 responses), each mentioned explicitly to say they are NOT in the network. No invented term names.

**Context inference (neutral HIV networks).** 9/9 responses inferred a loss-of-function (CRISPR/RNAi) screen for a trafficking/cell-surface phenotype — consistent with the true source (host-factor screens); none named HIV/virus; all flagged this as [outside the network]. Themes were stable across repeats (Golgi transport/COG, nuclear pore/NUP, glycosylation, clathrin endocytosis).

**ORA vs iGEA (GO).** Same headline themes; ORA text bridged via broad terms (protein localization to organelle, complex disassembly); iGEA text covered more distinct modules (glycolipid biosynthesis, clathrin vesicle, MED transcription, autophagy).

**Positive controls.** Hypoxia → low-oxygen transcriptomics (input contains HALLMARK: HYPOXIA, partly given). Vaccination → innate phagocyte (monocyte/macrophage) programme; did not name vaccination.

**Negative control (key caveat).** Random 272-gene list: unfiltered (16 terms) → "epithelial development + nuclear-receptor" story; filtered (2 terms/6 genes, network p≈0.27) → "steroid-hormone transcription", described as small but coherent. The model never declared the random network uninterpretable → the LLM does not detect chance structure; the q-filter and network significance must precede LLM interpretation (filter reduced the random story from 16 terms/multiple modules to 2 terms).

**Blind scoring.** `llm_package/blind_scoring/` — R01–R17 shuffled (seed 20261115), scoring sheet (5 criteria + guessed context), key in separate file.

## Step 07 — Communities, hub genes, bridges (R5-M8) — step07_modules/

**Method (code commit 7702f2b, `code/network_modules.py`).** Leiden (modularity, igraph) on the term–term projection of the q<0.10 network (terms linked by shared genes, weight = number of shared genes); best of 100 seeds (seed 2026), stability = mean ARI vs reported partition; genes assigned to the community of their terms; gene degree (= collections placing the gene; ≤ 10), participation coefficient (Guimerà & Amaral 2005), bridge = terms in >1 community. <6 terms → connected components. Descriptive, no test. CLI writes `network_modules_*.tsv` and a Community column; app shows an expander.

| Network | Terms/genes | Communities | Modularity | Stability (ARI; identical) | Bridge genes | Max degree |
|---|---|---|---|---|---|---|
| HIV all (q<0.10) | 49/162 | 7 | 0.66 | 1.00; 100% | 23 | 6 |
| HIV GO-only (Fig 2) | 25/130 | 7 | 0.70 | 0.85; 58% | 6 | 3 |
| Hypoxia | 32/64 | 4 | 0.30 | 0.92; 83% | 27 | 9 |
| Vaccination | 66/267 | 5 | 0.37 | 0.98; 92% | 91 | 8 |

HIV communities: C1 Golgi vesicle transport (GOSR2, RAB2A, COG2, ARF1, MGAT1); C2 EGFR/PI3K/NF-kB signalling + endocytosis (AKT1, AP2M1, EGF, EGFR = degree 6, RELA); C3 Mediator transcription (MED14/4/6/7); C4 nuclear pore/export (NUP107/133/155/160/85); C5 RNA exosome/processing (EXOSC3/5); C6 ER glycosylation/GPI anchor (PIGY 6, PIGH, PIGK, DPM1, STT3A); C7 cullin-RING ligase. Hypoxia: HIF/hypoxia core, oxygen response/cell cycle, glycolysis (PGK1 degree 9, GAPDH 8), 2-OG dioxygenases/collagen (P4HA1/2, PLOD1, LOX). Vaccination: neutrophil degranulation/lysosome, phagocytosis/NADPH oxidase, cytokine/inflammasome, amide metabolism, sulfotransferases.

**Hubs vs annotation (HIV; `hiv_hub_annotation_bias.tsv`).** Spearman degree vs number of size-filtered sets containing the gene: iGEA ρ=0.35 (162 genes), ORA FDR<0.05 ρ=0.35 (123 genes) — same rank bias (inevitable: more annotated genes have more chances). The difference is magnitude and meaning: ORA max degree 54 (NUP genes repeated across 54 nested terms; top-7 hubs all nuclear pore), iGEA max 6 (= collections agreeing); iGEA hubs span modules (AKT1 864 sets, but also PIGY with only 24 sets at degree 6). Do NOT claim iGEA removes annotation bias; claim it bounds degree and turns it into cross-collection support.

**LLM modules vs Leiden (`llm_vs_communities.tsv`).** Term pairs cited on the same line of an LLM response: within-community co-mention 0.48–0.60 vs between 0.07–0.08 for HIV GO (ratio 6.6–7.7, 3/3 repeats); HIV all 2.0–6.6; hypoxia 2.0; vaccination 2.4. The LLM's thematic grouping follows the network's community structure (conservative: responses are partly organised by collection, which forces between-community co-mentions). LLM hub claims checked: hiv_go_igea r1 listed 15 degree-3 genes — exactly correct.

## Step 08 — LLM export v2: communities + statistics in DOT and prompt (code commit 136e13e)

**Code.** `code/llm_export.py`: two DOT files per network — display DOT (drawing; term nodes now carry `library`) and semantic DOT `network_llm.dot` (no styling): one `subgraph cluster_Cn` per community; term attrs library/iteration/p/q/community/genes; gene attrs degree/community, bridge genes outside clusters with bridge=true + participation; graph attrs input_genes/random_lists/q_threshold/network_p; edges one line per term (`"term" -- {"G1" "G2"}`). Valid for Graphviz (fdp draws clusters) and pydot. Prompt `llm_prompt.txt`: statistical context (term filter + 5 network metrics with empirical p; rule: if p(cluster nodes) ≥ 0.05 say so first, stay tentative, no experimental context), community table, DOT, guide, request (per-community naming with merge/split allowed, hubs & bridges, hypothesis, context), traceability rule. CLI writes both files (old combined_ai_analysis_prompt_enhanced.txt retired); app offers both downloads when statistics exist. Tests: 6 pass (test_llm_export, test_network_modules).

**Paper inputs** (`llm_package/v2/inputs`, `work/build_llm_inputs_v2.py`; statistics from the saved 10,000-list reports):
| Job | Terms/genes | Communities | Network p (cluster nodes) | ~tokens (v1) |
|---|---|---|---|---|
| hiv_go_igea | 25/130 | 7 | <0.0001 | 5.6k (6.4k) |
| hiv_all_igea | 49/162 | 7 | <0.0001 | 8.1k (10.3k) |
| hypoxia_all_igea | 32/64 | 4 | <0.0001 | 5.2k |
| vaccination_all_igea | 66/267 | 5 | <0.0001 | 11.6k (16.8k) |
| random filtered (q<0.10) | 2/6 | 1 (components) | 0.27 | 1.7k |
| random unfiltered (q=1) | 16/48 | 7 | 0.60 | 3.7k |
Random q=1 statistics: `work/random_list_stats_q1.py` (q threshold 1.0, 10,000 lists, reads the saved control list) → largest cluster 33 nodes/24 genes/9 terms vs random mean 42.9, p = 0.60.
Note (resolved): the gene background depends on the working directory. BackgroundGeneSet validates symbols with NCBI gene_info, which GeneConverter looks for (uncompressed) in data/recent_release of the repo or of the current directory; the repo ships only the .gz. Run from ~/s01wd (decompressed copy) → 23,391 validated genes; run elsewhere → silently 23,590 unvalidated genes. All paper analyses ran from ~/s01wd (23,391; stated in every report). Two runs today had started elsewhere and were redone from ~/s01wd: the q=1 random statistics (p 0.60 → 0.60; report now says 23,391) and step07 annotation_bias.py (unchanged, ρ = 0.35 both). The control list regenerates exactly from its seed with the 23,391-gene background. Repo fix pending: read gene_info.gz directly / fail loudly.

**Run:** `python3 run_llm.py --set v2 --model gpt-6-astra` (14 calls → `v2/responses`). v1 responses kept as the unguided comparison.

## Step 09 — Code review + end-to-end runs (code commits a3dbdbb, 49c921f) — step09_end_to_end/

**Review findings fixed.** (1) NCBI gene data was not in the repo (gitignored) and a missing file silently left the background unvalidated (23,590 vs 23,391 genes, depending on working directory) → gene_info.gz + human gene_history now versioned, read directly, error if missing. (2) Input genes outside the background were kept (HIV: 7 of 272 — FAM226A, IGHV1OR21-1, KRABD5, LINC00208, LINC01994, LOC400627, RPL32P3) while random lists drew all genes from the background; gene sets counted 95 non-background genes in K (1,575 of 659,814 memberships) → one universe: input and gene sets restricted to the background (size filter after restriction). (3) Defaults: 10 collections; max 30 iterations (no dataset exceeds 12); >800 genes warns in CLI. (4) `--null-cache` makes long runs resumable (identical results, verified). (5) `code/llm_client.py` + `--llm-model`.

**End-to-end CLI runs (defaults, 10,000 lists, run from /tmp to prove cwd independence):**
| Run | Input (in bg) | Terms kept / selected | Largest cluster (genes/terms/edges, collections) | Network p (nodes; edges/gene) | Communities (modularity) |
|---|---|---|---|---|---|
| HIV, 10 collections | 265 (272 validated) | 56/62 (was 49/62) | 168/56/381, 10 | <0.0001; 0.006 | 7 (0.65) |
| HIV, GO only (Fig 2) | 265 | 26/33 (was 25/33) | 138/26/212, 3 | <0.0001; 0.005 | 6 (0.68) |
| Hypoxia | 75 | 32/32 (unchanged) | 64/32/171, 9 | <0.0001; 0.005 | 4 (0.30) |
| Vaccination | 426 | 66/76 (was 66/75) | 267/66/662, 9 | <0.0001; 0.003 | 5 (0.37) |
| Random control (seed 999, 265 genes from bg; 264 validated), q<0.10 | 264 | 1/15 | 4/1/4, 1 | 0.49 | — |
| Random control, q=1 | 264 | 15/15 | 23/7/32, 5 | 0.63 (genes) | 9 |

HIV: same 62 selected terms; p-values shift ≤ 0.06 log10 units; 7 borderline terms (6 WikiPathways at q 0.101, Reactome HIV transcription elongation 0.118) now q 0.088–0.097 → kept. Smaller n (265) → random lists select fewer chance terms (16.0 vs 16.2) → slightly lower q. Vaccination: one extra non-kept GO BP term selected. Hypoxia identical.

**Consequences for the paper:** HIV numbers (Results table, Fig 1/2 counts, Supp Fig 2/3 DOTs, stability, community table, LLM inputs v2) must be regenerated from step09 outputs; hypoxia/vaccination unchanged. The earlier `llm_package/v2/inputs` are superseded by the step09 `llm_prompt.txt` files.

**LLM:** `bash step09_end_to_end/run_llm_all.sh gpt-6-astra` (Terminal, key exported) → 14 calls into each run's `llm/` folder.

## Step 10 — LLM prompt iterations (HIV, 15 Nov) → v5 chosen

| Version | Kept in HIV folder | What changed | Verdict |
|---|---|---|---|
| v2 structured | llm_prompt_v2_structured.txt, llm_v2_structured/ (all 6 runs, 14 calls) | stats + communities + semantic DOT; per-section analysis; exact term/gene citation per statement | mechanical, reads the network back |
| v3 narrative | llm_prompt_v3_narrative.txt, llm_v3_narrative/ | natural prose, genes in parentheses, no term quotes | lost communities |
| v4 communities | llm_prompt_v4_communities.txt, llm_v4_communities/ | exec summary, per-community biology, combined story, questions | lost experimental-context guess; questions not useful |
| **v5 (code 9201258)** | llm_prompt.txt, llm/ | exec summary, communities, putting it together (hub/bridge roles), **estimated experimental context**, closing invitation; **second turn** (llm_client --context) revealing the experiment: revised summary, fits, does not fit, may be new, priorities | **chosen** |

HIV v5 r1: guessed "cell-based loss-of-function screen (CRISPR or RNAi)", alternative "host-factor screen for infection or toxin entry"; second turn with "an siRNA screen in human cells for host factors required for HIV infection (HIV dependency factors)".
Requests are stateless and store=false. Final run: `bash step09_end_to_end/run_llm_all.sh gpt-6-astra` — 3 repeats per run folder (HIV, HIV GO, random q<0.10, random q=1, hypoxia, vaccination) + second turn for HIV, HIV GO, hypoxia ("a meta analysis of genes that are consistently regulated by mRNA expression when cells are exposed to hypoxia") and vaccination (day-3 blood genes correlated with day-28 HAI response after seasonal influenza vaccination). No second turn for random controls.

## Step 11 — Final LLM run (prompt v5, gpt-6-astra, 16 Nov) — checks in step09_end_to_end/llm_v5_checks.tsv
30 answers (18 first turns: 6 runs × 3; 12 second turns: HIV, HIV GO, hypoxia, vaccination × 3). One provider refusal on the first attempt (OpenAI bio_policy, HIV second turn r3); the identical request succeeded when resumed — no answers missing.
- **Invented genes:** none. Only out-of-network symbols: HIF1A (hypoxia, 4 answers; explicitly as the upstream regulator not in the list) and CDK1 (hypoxia r1; as CCNB1's partner) — both legitimate outside knowledge.
- **Communities:** every answer covers every community (HIV 7/7, HIV GO 6/6, hypoxia 4/4, vaccination 5/5, random q=1 9/9).
- **Random controls:** all 6 answers open by stating the network is not distinguishable from random lists (p ≈ 0.49 and 0.63), are short (350–690 words vs 1,130–1,250), and omit the experimental-context section, as instructed.
- **Experimental context (first turn, blind):** HIV 6/6 (incl. GO-only) → cell-based loss-of-function screen (CRISPR/RNAi), with "host-factor screen for a pathogen/infection" as alternative or within the hypothesis (r1, GO r1); hypoxia 3/3 → transcriptomics of cells under low oxygen; vaccination 3/3 → myeloid/monocyte-macrophage expression signature (vaccination not inferred; context-free list of blood genes).
- **Second turn (experiment revealed):** HIV: guess "partly held"; nuclear pore (NUP153, RANBP2), CD4 and transcription become central; new candidates: glycosylation–trafficking coupling, GPI machinery, RNA surveillance. Hypoxia: "matches extremely well"; new: 2-OG dioxygenase chromatin–collagen coupling, secretory stress, nuclear-receptor modulation. Vaccination: myeloid estimate held; key reinterpretation — genes mark future high responders, not vaccine-induced mechanisms; candidates: phagolysosomal competence, lipid turnover, regulated responsiveness.

## Step 12 — HIV results regenerated with one universe (16 Nov) — old outputs in `_superseded_2026-11-15_before_one_universe/` of each folder
All scripts now: input restricted to the background (265 genes), gene sets restricted before the 10–600 filter, max 30 iterations, q < 0.10 two-level statistics, gene data from the repo.

| Result | Before (272 genes) | Now (265 genes) | Source |
|---|---|---|---|
| HIV terms kept / selected | 49/62 | **56/62** | step09 CLI |
| HIV largest cluster (genes/terms/edges/collections) | 162/49/342/9 | **168/56/381/10**; p < 0.0001 (edges/gene 0.006) | step09 |
| Fig 2 GO network (kept, genes/terms/edges) | 25/33; 130/25/196 | **26/33; 138/26/212**; p < 0.0001 (edges/gene 0.005) | 02 Figures/networks/hiv_go_* |
| Fig 1 GO BP: ORA FDR<0.05 vs iGEA q<0.10 | 14 terms/57 genes vs 5/55, 4 shared | **unchanged** (same terms) | Fig1_counts.tsv |
| GO BP-only network p | — | cluster 20 genes/1 term, p = 0.0003 (terms p = 0.10) | hiv_gobp_statistical_report.txt |
| Supp Fig 2 / 3 | 49 terms/162 genes; 25/130 | **56/168; 26/138** | DOTs + previews regenerated |
| Communities HIV (10 coll.) | 7, Q 0.66, 23 bridges, max degree 6 | **7, Q 0.65, 31 bridges, max degree 7**; same 7 themes; 100% seed agreement | step07_modules/hiv_* |
| Communities HIV GO | 7, Q 0.70, 58% identical | **6, Q 0.68, 100% identical** (clathrin vesicle merged into C1) | step07 |
| Hub vs annotation (Spearman) | iGEA 0.35, ORA 0.35 | **iGEA 0.41 (168 genes), ORA 0.35**; max degree 7 vs 54 | step07 annotation_bias |
| ORA FDR<0.05 vs final iGEA (10 coll.) | (earlier vs iGEA q<0.05: 87 vs 25) | **87 ORA vs 56 iGEA terms, 19 shared; genes: ORA 123, iGEA 168, both 113, only ORA 10 (incl. CXCR4), only iGEA 55; 2 ORA terms mostly uncovered (Notch signalling; protein localization to organelle)** | step04 s04_final_ora_vs_igea.py |
| Chance terms per random list (ORA vs iGEA) | 31.4 vs 16.2 | **31.1 vs 16.0** (GO BP 16.5 vs 6.3) | step04 chance_terms |
| Stability (100 × 90% subsamples, 238 genes, 5,000 ref lists) | p ≤ 0.0002 100/100; cluster gene J 0.74; 34/49 recovered; median recurrence 0.73; near ties 34/62 | **p ≤ 0.001 100/100; cluster gene J 0.73 (median 0.75); 38/56 recovered (median 39); kept-term recurrence median 0.71 (46/56 ≥ 0.5); near ties 35/62, recurrence 0.58 vs 0.74** | step06 |
| Redundancy between collections (J ≥ 0.9) | 74 of 10,445 sets, 12 identical | **82 of 10,439, 12 identical** (sets restricted to background) | step03 |
| Deduplicated collections | cluster 181/63/418 (no q filter) | **57 kept, cluster 170/57/387, p = 0.001 (1,000 lists, q < 0.10)** | step03 |
| Leave-one-collection-out (q < 0.10, 1,000 lists) | p ≤ 0.001 all | **p = 0.001 (minimum) for cluster size in all 10; largest drop without GO MF (78 genes/24 terms)** | step03 |
Hypoxia and vaccination: unchanged networks (step09). Not regenerated (historical, superseded by the above): step01 reproduction, step02, step04 Reactome examples (q < 0.05 exploratory), step05 HIV files, llm_package v1/v2.

## Step 13 — Repository branches (16 Nov, local only, nothing pushed)
- `main` fast-forwarded to `csbj-revision` (edfdf1c, 17 commits over v6); tags `csbj-v7-revision` and `v2.0.0` on it; `csbj-v6-submitted` kept.
- `csbj-revision-cloud` (c17d549) = origin/streamlit-cloud + 14 commits: the 4 cloud statistics commits, the author's 2 unpushed March UI commits (statistics now above the network), and all 16 Nov changes ported (communities, LLM export/prompt v5, 10-collection defaults, one universe, max 30 iterations, gene data .gz support + error if missing, llm_client). Cloud keeps its own layout (no CLI, no tests, plain gene files, AI-prompt text box with platform buttons, now showing the v5 prompt and the semantic DOT download).
- Checked in the cloud worktree: background 23,391, HIV 265; IterativeEnrichment = FastCollection; communities + LLM export from app state; app starts (health ok); AppTest: example → iterative → Submit → Generate Network, no exceptions, prompt contains statistics, communities and experimental-context section.
- To publish: `git push origin main csbj-v7-revision v2.0.0` (GitHub code) and `git push origin csbj-revision-cloud:streamlit-cloud` (redeploys the public app). Both fast-forward.

## Step 14 — Supplementary Tables 3–5 (16 Nov) → `v7/03 Supplementary Material/` (builder: `build/build_supp_tables_3_5.py`)
- **Bug found and fixed (code 9531a05, csbj-revision; after the pushed v7 tag):** CLI regular (ORA) mode used the first/last GMT set sizes as the size filter (e.g. [200, 200] for Hallmark), so most sets were never tested. App unaffected; paper ORA numbers came from scripts (FastCollection) and are unaffected. Fixed CLI now matches the scripts exactly: 254 terms with p < 0.01, 87 with FDR < 0.05.
- Table 3 = CLI regular run (`step09_end_to_end/ora_regular/`), HIV 265 genes, 10 collections, ranked by p per collection. Table 4 = CLI iterative run: 62 terms with iteration p, q-value, in-network (56), community, plus full-list overlap/p/BH FDR (19 of 62 have full-list FDR < 0.05). Table 5 = statistical reports: Part A 10 collections, Part B GO network (Figure 2).

## Step 15 — Hosted app vs Streamlit Community Cloud (16 Nov)
Limits (docs.streamlit.io, Feb 2024 figures): 0.078–2 CPUs, 690 MB–2.7 GB RAM, ≤50 GB storage; sleep after 12 h without traffic; requirements.txt in root (uv, pip fallback), packages.txt = apt (Debian); Python chosen at deployment and fixed until redeploy; forum (Aug 2026): new deployments default to Python 3.14.
Found and fixed:
1. Enrichment pool size `cpu_count()-2` (≤0 on small machines → Pool error → **no results returned silently**; huge on big hosts) → bounded workers + serial fallback (main c09b8da); cloud serial by default (9670397).
2. Exact pins (pandas 2.2.3, scipy 1.15.2, Pillow 11.1.0) have no Python 3.14 wheels → cloud requirements as ranges (a1710c2); resolves with wheels only on x86_64 for 3.12/3.13/3.14; app flow + HIV terms/p identical (max rel diff 1.7e-13) on 3.12 and 3.14.
3. Reference data cached with st.cache_data (copy per call) → st.cache_resource (ab07b0a): peak RSS over 3 HIV sessions 489–521 MB (was 536–705).
4. Speed: fast_igea skips sets whose table probability ≥ threshold (exact lower bound of two-sided Fisher) — identical results (tests; HIV/hypoxia random lists and q identical) — ~20× faster for 800-gene lists (main 663d2ff/cc5adce). Hosted defaults (1,000 lists, 200 per grid size, 1 process, 10 collections, VM core): HIV 12.5 s iGEA + 17 s statistics, 536 MB; vaccination 19 + 26 s, 567 MB; random 800 genes 16 + 51 s, 433 MB.
Open: app default collections are 3 (Hallmark, Reactome, GO BP) in both app versions vs 10 in paper/CLI.

## Step 16 — Hosted app = public repo ybenita/iterative-gene-enrichment-app (17 Nov)
- The live app deploys from the public repo (folder "Iterative enrichment cloud version", branch main c9bdda7, entry streamlit_app.py → code/streamlit_app.py); it is aion origin/streamlit-cloud + the author's Feb–Aug UI commits; v6 statistics (stored permutations), 11 collections incl. HPO.
- New branch **v7** (6588894) in that repo, one commit on main: v7 analysis code (statistics, one universe, communities, AI export, 10 collections, 30 iterations, serial enrichment, shared caches, dependency ranges incl. streamlit>=1.48,<2) merged into the author's UI (docs button, default collections, no Select All, inline downloads, colours). llm_client not included. README updated.
- Tests (Python 3.12 and 3.14, Streamlit 1.67, pandas 2.3, SciPy 1.18): HIV iterative → 265 genes, 56/62 kept, cluster 168/56/381, 7 communities, v5 prompt, no exceptions; regular mode OK; peak memory 607–774 MB over 3 sessions (live version ≈800 MB for one).
- Publish: `git push origin v7` then `git push origin v7:main` (Community Cloud redeploys from main).
- (17 Nov) Author pushed v7 → GitHub main = 6588894 (live). Slimming commit 1c5fa96 on v7 (local, not pushed): no packages.txt (Graphviz unused), statsmodels → scipy BH (FDR identical to 1e-16; Table 3 reproduced), graphviz/streamlit-option-menu removed, old permutation benchmark + unused files + HPO/STRING GMTs removed, gene data as .gz (39 → 6 MB), results/ ignored; working tree 16 MB; HIV/regular-mode checks unchanged on 3.12 and 3.14.

## Step 17 — Terminology "collection"; Supp Tables 3–5 headers (17 Nov)
- Code (csbj-revision 5cd9a13, 612d060; public v7 8cecb59, c369b6c): app, downloads, statistical report, READMEs and LLM prompts say "collection"; CLI `--collections/-c` (old `--libraries/-l` still accepted); semantic DOT attribute `collection=`. LLM package of the paper (v5, 15–16 Nov) not rerun: its prompts said "library source".
- Supp Tables 3–5 rebuilt from the same 16 Nov runs with the column header "Collection" (Table 3, Table 4, Table 5 TERMS headers); diff against the previous files = header lines only. Previous files in `03 Supplementary Material/_superseded_2026-11-17_library_header/`.
- Licence: GPL-3.0 (LICENSE file; pyproject updated in both repos, 6bc613a / 7cadc0f).

## Step 18 — Supplementary Table 1 (pseudocode) for v7 (20 Nov)
- v6 pseudocode checked against the code (fast_igea.py, per_run_null.py, network_modules.py). Errors in v6: remaining genes and iteration counter reset inside the WHILE loop; size filter inverted (`max_gene_set_size < term_gene_size`); no minimum gene-set size; max iterations 10 (default now 30); no handling of "no term tested"; "Libraries".
- v7 file (`03 Supplementary Material/Supp Table 01 - iGEA pseudocode.txt`): Part 1 = iterative selection with one universe (input ∩ background; sets restricted to background before the 10–600 filter; M = collection universe, N = input in it), 2×2 table, Fisher two-sided / hypergeometric option, ties by GMT order, cap 30; Part 2 = random lists (seed (2026, n, i)), size grid (n/20 steps to n/4, max(200, lists/10) per size), signal-adjusted q (≤ 10 rounds, rounded-up size, leave-one-out, monotone per collection, accept level 0.10 when q = 1), q-filtered network, empirical p, z descriptive, Leiden communities. v6 untouched.

## Step 19 — Pseudocode check; two-sided Fisher selects depleted sets in large lists (20 Nov)
- Supp Table 1 (HTML + txt) checked line by line against fast_igea.py, per_run_null.py, network_connectivity_benchmark.py, network_modules.py. Fixed: exact reduced sizes (n − s, n − 2s, … ≥ n/4, s = ceil(n/20)); seeding (seed, size, j); q_threshold = 1 keeps all terms.
- **Finding:** selection uses the two-sided Fisher p, so an under-represented (depleted) gene set can be selected (k < N·K/M). Counted on the selected terms (k vs expected): HIV 0/62, hypoxia 0/32, vaccination 0/76; ORA p < 0.01: HIV 0/254, hypoxia 0/120, vaccination 0/497. Random lists (200 each, 10 collections): 265 genes 0/3,327; 426 genes 1/4,457; **800 genes 1,485/7,172 (21%)**. Paper results unaffected; large input lists are. Script /tmp/depl.py (VM). Awaiting author decision (enrichment-only selection).

## Step 20 — Switch to one-sided Fisher's exact test; everything rerun (20 Nov)
Author decision after Step 19. Code: `fisher_exact(alternative="greater")` = hypergeometric upper tail (csbj-revision 46af50e; public v7 42ea4fd); random-list cache key versioned (67911f7 / 55b3cc6) — the old cache had the same key and was silently reused by subset_stats.py and s03b at first; both rerun with fresh lists; report states "(one-sided)" (64ea03d / b2c306b). Test added (17 tests pass); hosted app HIV 56/62, same cluster.
Old outputs kept in `_superseded_2026-11-20_two_sided_fisher/` in step09_end_to_end, 02 Figures/networks, step03, step04 (3 files), step06, step07, 03 Supplementary Material (Tables 3–5).

| Result | Two-sided | One-sided |
|---|---|---|
| HIV 10 coll.: selected / kept; cluster genes/terms/edges/coll.; p | 62 / 56; 168/56/381/10; < 0.0001 (e/g 0.0063) | **identical** (p up to 0.19 log10 lower; q change ≤ 0.007; noise genes 97) |
| HIV GO only (Fig 2); GO BP only | 26/33; 138/26/212; 6 comm. | **identical** (fresh random lists) |
| Hypoxia | 32/32; 64/32/171; 4 comm. | **identical** |
| Vaccination | 66/76; 267/66/662/9; 5 comm. (Q 0.37) | **70/80; 273/70/691/9; 6 comm. (Q 0.37); p < 0.0001, e/g 0.0029**. +GOMF cadherin binding, GTPase binding, inward rectifier K channel, oxidoreductase (sulfur donors), sulfotransferase; −GOMF SNARE binding |
| Random control q < 0.10 / q = 1 | 1/15, p 0.49 / p 0.63 | identical (q = 1: p 0.63) |
| ORA HIV p < 0.01 / FDR < 0.05 | 254 / 87 | **255** (+GOMF transcription coregulator activity) / 87 (same set) |
| ORA FDR vs iGEA | 87 vs 56, 19 shared; genes 123 vs 168 | identical |
| Chance terms per random list ORA / iGEA | 31.1 / 16.0 | **31.3 / 16.1** |
| Figures 1, 2, Supp Fig 3 DOTs | — | byte-identical; Supp Fig 2 same network (one WikiPathways term reordered) |
| Dedup collections | 57 kept, 170/57/387, p 0.001 | identical |
| Leave-one-out (1,000 lists) | p 0.001 all | p 0.001 all; larger networks without Reactome (47 vs 38 kept), GO BP (48 vs 36), GO CC (44 vs 35); without GO MF 78/24 unchanged |
| Stability | 100/100 p ≤ 0.001; J 0.73; recurrence 0.71; near ties 35/62 | same (J 0.731, kept 49.1 vs 49.0) |
| Hubs vs annotation | 0.41 vs 0.35 | identical |
| Supp Tables 3–5 | — | rebuilt (p-values; Table 3 255 rows) |
| Supp Table 1 | — | test described as one-sided |
Pending: LLM rerun (author, Terminal: `bash step09_end_to_end/run_llm_all.sh gpt-6-astra`; prompts carry the new p-values and vaccination has a new network), then llm_v5 checks + LLM supp PDFs; manuscript Implementation text (two-sided → one-sided) when the docx is closed.

## Step 21 — New Supplementary Table 2 (gene-set collections); supplementary tables renumbered (20 Nov)
- New Supp Table 2 (`03 Supplementary Material/Supp Table 02 - Gene-set collections.html/.txt`; data `build/collection_table_data.py/.tsv`, page `build/build_supp_table_2.py`): per collection gene sets in file / used (10–600 after background restriction), genes covered (M), median size; HIV: input genes covered, possible terms (≥ 3 genes), ORA p < 0.01, ORA FDR < 0.05, iGEA selected, kept. Totals: 14,534 / 10,439 sets; 18,506 genes; HIV 246 of 265 covered; 1,508 possible; 255; 87; 62; 56. Answers R5-M9, R4-3, R3-m3.
- Renumbering: old 2–12 → 3–13 (Table 1 unchanged). Files: Supp Table 03/04/05 → 04 ORA results / 05 iGEA results / 06 Network statistical report; builder renamed build_supp_tables_4_6.py; build_llm_supp.py cites the prompt as Supplementary Table 3 (HTML/PDF regenerate with the LLM rerun).
- Manuscript: 15 citations renumbered as tracked changes + placeholder replaced by "Supplementary Table 2"; comment added. Reject-all view = submitted v6 text; 25 Zotero fields intact.
- Open: "Supplementary Tables X-Y" (LLM responses) — set when the LLM supplement is final.

## Step 22 — LLM rerun with one-sided-Fisher networks and "collection" prompt wording (20 Nov)
Author ran `run_llm_all.sh gpt-6-astra` (prompt v5, wording "collection", DOT attribute collection=): 30 answers (18 first, 12 second turns), no refusals. Old answers in `step09_end_to_end/_superseded_2026-11-20_two_sided_fisher/`. Checks: `analyze_llm_v5.py` (glob fixed to folders with network_llm.dot) → `llm_v5_checks.tsv`, `llm_v5_experimental_context.tsv`.
- Invented genes: none. Only out-of-network symbol: HIF1A (hypoxia, 5 of 6 answers; named as upstream regulator).
- Communities: every answer covers every community (HIV 7/7, HIV GO 6/6, hypoxia 4/4, vaccination 6/6, random q=1 9/9, random q<0.10 1/1).
- Random controls (6 answers): all open with "not distinguishable from random gene lists"; 336–681 words (real datasets 1,036–1,238); no experimental-context guess (one q=1 answer has a context heading that states there is insufficient support to nominate one).
- Experimental context, blind (first turn): HIV 6/6 genetic/loss-of-function screen with a cell-surface/uptake readout (no answer names HIV or infection; the three 10-collection answers mention virus-named terms only to say they are not evidence of infection); hypoxia 3/3 low-oxygen transcriptomics; vaccination 3/3 myeloid (monocyte/macrophage) signature, vaccination not inferred.
- Second turn: all 12 have the five sections. **HIV: 4 of 6 second turns add a provider safety limit** ("I can't nominate new HIV dependency targets or provide prioritized experiments…"; offers non-infectious follow-up instead) — previous run 1 of 6 + 1 refusal. Hypoxia: "matches the experiment very well". Vaccination: myeloid estimate held; genes mark future high responders rather than vaccine-induced mechanisms.
- LLM supplementary documents rebuilt (`03 Supplementary Material/LLM responses/html|pdf`, prompt cited as Supplementary Table 3; PDFs rendered in the cloud workspace — no Playwright on the VM). Old documents in `LLM responses/_superseded_2026-11-20_two_sided_fisher/`.

## Step 23 — Final LLM run with gpt-6.2-sol (20 Nov)
Model changed from gpt-6-astra (Step 22) to **gpt-6.2-sol**, because the astra second turns for HIV carried explicit provider safety limits that reviewers would question (author decision after the Sol HIV test). Author ran `bash run_llm_all.sh gpt-6.2-sol` (prompt v5, same networks): answers in `<run folder>/llm_gpt-6.2-sol/`; the astra answers stay in `llm/`. All 30 calls: requested = returned model gpt-6.2-sol, status completed, no refusals, no truncation. Checks: `LLM_DIR=llm_gpt-6.2-sol analyze_llm_v5.py` → `llm_v5_checks_gpt-6.2-sol.tsv`, `llm_v5_experimental_context_gpt-6.2-sol.tsv`.
- Invented genes: none. Out-of-network symbols, each explicitly marked as outside knowledge: HIF1A (hypoxia, 2 answers; "This is not in the network"), CDK1 (hypoxia, 1), P2RY12 and TMEM119 (vaccination, 2; named as microglia markers to rule a tissue in or out).
- Communities: every answer covers every community (HIV 7/7, HIV GO 6/6, hypoxia 4/4, vaccination 6/6, random q=1 9/9, random q<0.10 1/1).
- Length: first turns of real datasets 1,027–1,185 words; random controls 304–698; second turns 935–1,147.
- Random controls (6 answers): all open with "not distinguishable from random gene lists" and give the empirical p (0.49 at q < 0.10; 0.63 at q = 1, largest cluster 30 vs random mean 42); interpretation kept tentative; no experimental context nominated (one answer states the evidence is insufficient to favour an assay, perturbation or cell type).
- Experimental context, blind (first turn): HIV 6/6 (10 collections and GO) a loss-of-function / perturbation screen in cultured mammalian cells with a cell-surface, secretory-cargo or receptor-response readout (HIV and infection not named); hypoxia 3/3 expression comparison under sustained low oxygen (HIF-like response; matrix-producing cells); vaccination 3/3 bulk expression from blood or inflamed tissue with a myeloid (neutrophil and/or monocyte) contribution; vaccination not inferred.
- Second turn: all 12 have the five sections. HIV (6): no explicit refusal or safety statement; the models say the blind estimate was right about the perturbation screen but missed infection as the readout, and discuss CD4 (entry receptor), NUP153/RANBP2 (nuclear import), Tat/P-TEFb and elongation, reporter-readout artefacts. **All proposed follow-up experiments are in uninfected cells or with non-viral readouts** (6/6); one answer says so ("These can be tested without virus or viral material; they would clarify the screen without assigning or optimizing new HIV dependencies"). To be disclosed in the Results (author agreed). Vaccination: myeloid estimate held; follow-up with non-infectious fluorescent vaccine antigen.
- LLM supplementary documents rebuilt from the Sol answers: `build_llm_supp.py` now reads the answer folder from env `LLM_DIR`; `LLM_DIR=llm_gpt-6.2-sol python build_llm_supp.py "../../07 Reviewer Analyses/step09_end_to_end" html`; PDFs rendered in the cloud (HIV 10 collections, HIV GO, hypoxia, vaccination 17 pages each; random control 12). Astra documents in `LLM responses/_superseded_2026-11-20_gpt-6-astra/`.
- Manuscript Methods: model "gpt-6.2-sol … on 20 November 2026" (inside the existing tracked insertion); comment 81 notes the model change.

## Step 24 — ORA LLM baseline; Results rewritten (20 Nov)
- **ORA baseline (R2-4):** `step09_end_to_end/build_ora_llm_input.py` → `ora_go_fdr0.05/`: HIV, 3 GO collections, ORA BH FDR < 0.05 = 19 terms, each linked to all its genes → 79 genes, 5 Leiden communities (modularity 0.36; iGEA GO 0.68), NUP genes degree 8–9. Prompt = llm_export prompt with only the method sentences replaced (ORA introduction, FDR instead of network statistics, gene degree, q attribute; each replacement asserted). Author ran `run_llm_ora.sh` (gpt-6.2-sol × 3, first turn only): 1,067–1,116 words; same five themes as iGEA GO (pore/RNA export, Golgi, glycolipid/GPI, Mediator, RNA/complex disassembly) and same context (cell-based perturbation screen, protein abundance/localization readout); genes named 33–42 of 79 (iGEA GO answers 51–55 of 138, incl. CD4 and EGF/AKT1, absent from the ORA network); no out-of-network genes.
- LLM supplement: `build_llm_supp.py` handles the ORA folder (no statistics) → `LLM_HIV_GO_ORA_baseline.html/.pdf` (11 pages).
- Supplementary numbering for the LLM documents (unchanged from the renumbered v6 slots): 7 HIV GO, 8 HIV 10 collections, 9 ORA baseline, 10 DAVID (v6 report), 11 hypoxia, 12 vaccination, 13 random control. Methods placeholder → "Supplementary Tables 3, 7–9 and 11–13".
- **Manuscript Results rewritten** (tracked, author Claude, against v6; reject view = v6): subsections ORA and iGEA; two-level significance (+ random control); robustness (stability, collection redundancy); communities, hubs, bridges; LLM-assisted interpretation (Sol answers, ORA baseline, random control, DAVID shortened); additional datasets. 15 new margin comments with reviewer IDs. Still v6 and flagged: Supplementary Figures 1, 4, 5, 6 (to regenerate), HIV life-cycle sentence; figure legends.
- Note: the 15:0x commit of the model-name edit was briefly reverted on the device (file reappeared with the earlier content, same mtime pattern as an iCloud sync); re-committed, then the Results version committed on top. Check that the open copy in Word is the latest before editing.

## Step 25 — ORA vs iGEA as collections are added (20 Nov)
- **Structure** (`step09_end_to_end/scaling_ora_vs_igea.py` → `scaling_ora_vs_igea.tsv`; collections added GO → +Reactome → +KEGG ×2 → +WikiPathways/BioCarta/PID → +Hallmark): ORA (FDR < 0.05) 19→84→87 terms, 79→122→123 genes, 200→763→783 edges, max gene degree 9→53→54, Leiden communities 5→3→2 (modularity 0.13 at 10); iGEA (q < 0.10) 27→36→56 terms, 144→159→168 genes, 220→284→381 edges, max degree 3→4→7, communities 6→6→7. Small collections add no ORA terms at FDR < 0.05. Figure: `02 Figures/networks/ORA_vs_iGEA_10_collections_communities.png` (ORA C1 33 terms Golgi/glycosylation/signalling; C2 54 terms nuclear pore + Reactome viral/SUMO/mitosis/glucokinase terms carried by 6 NUP genes of degree 48–54).
- **LLM, 10-collection ORA** (`ora_all_fdr0.05/`, 87 terms, prompt 42.5k chars, 16,230 input tokens vs 13,070 for iGEA; author ran `run_llm_ora.sh gpt-6.2-sol ora_all_fdr0.05`, 3 repeats). Predefined scoring (`score_themes.py`, themes = 7 iGEA communities, marker = hub genes): themes discussed ORA 5–6 of 6 present (cullin-RING absent from ORA network) vs iGEA 7/7; genes named 30–35 vs 43–46; community sections 2 vs 7; sentences per theme (mean of 3) ORA vs iGEA: EGFR/AKT 3.0 vs 11.0, Mediator 1.0 vs 2.0, exosome 0.7 vs 2.0, GPI 1.7 vs 3.7, cullin 0 vs 1.7, Golgi 3.7 vs 3.7, pore 2.0 vs 2.0. Blind context: same class (cell-based perturbation screen, protein localization/abundance readout; r1 "low confidence"). All three ORA answers spend text explaining away redundant labels ("the same pore proteins appear under seemingly unrelated headings, including glucose metabolism, stress responses, and infection"); none inferred infection.
- Reading: the model recognises ORA redundancy and keeps the main themes; the difference is granularity (two sections; minor themes reduced to one sentence, Mediator dismissed as an "adipocyte" label) and coverage (fewer genes). GO-only ORA vs iGEA: 5 vs 6–7 sections, similar themes.

## Step 26 — ORA count check: BH denominator (20 Nov)
Author queried "only 87". Verified by two implementations (CLI regular mode and s04_final_ora_vs_igea.py): 87 = BH within each collection over **all** gene sets in the 10–600 range after background restriction (10,439, incl. sets with no input gene; code: enrichment.py multipletests over every tested set), reported if overlap ≥ 3. Variants (`step04_term_fdr/s04_ora_fdr_variants.py` → `hiv_ora_fdr_variants.tsv`), FDR < 0.05 / 0.10 / 0.25: all sets per collection 87 / 112 / 169 (pooled 83 / 99 / 172); sets with ≥ 1 input gene (4,895) 116 / 146 / 367 (pooled 101 / 156 / 366); sets with ≥ 3 input genes (1,508, "testable") 219 / 383 / 780 (pooled 228 / 403 / 798); nominal p < 0.01: 255. Decision pending (author): keep 87 (v6 convention) and state the denominator, or restrict BH to testable sets.

## Step 27 — ORA vs iGEA as collections are added, three datasets; journal limits (20 Nov)
- `step09_end_to_end/scaling_all.py` → `scaling_all.tsv` (= Supplementary Table 14, `03 Supplementary Material/Supp Table 14 - ORA vs iGEA as collections are added.txt`); networks per dataset in `scaling_networks/`. ORA computed with the vectorised one-sided Fisher (same as CLI; HIV reproduces 19/87 terms). iGEA = 10-collection kept terms restricted to each step.
- At 10 collections, ORA vs iGEA: HIV 87 vs 56 terms, 123 vs 168 genes, edges/gene 6.4 vs 2.3, max degree 54 vs 7, communities 2 vs 7, modularity 0.13 vs 0.65; hypoxia 52 vs 32, 50 vs 64 genes, 6.7 vs 2.7, 34 vs 9, 2 vs 4, 0.19 vs 0.30; vaccination 237 vs 70, 318 vs 273 genes, 10.8 vs 2.5, 76 vs 8, 4 vs 6, 0.25 vs 0.37. ORA with all set sizes (no 10–600): HIV 124 terms/194 genes, hypoxia 62/56, vaccination 372/375 (edges/gene 9.1, 8.4, 21.0).
- Trend with added collections: HIV ORA modularity 0.36 → 0.13 (communities 5 → 2) vs iGEA 0.68 → 0.65 (6 → 7); hypoxia ORA 0.29 → 0.19 (3 → 2) vs iGEA 0.36 → 0.30 (4 → 4); vaccination ORA already dense with GO alone (214 terms, 10 edges/gene, max degree 73; 0.24 → 0.25) vs iGEA 0.42 → 0.37 (6 → 6). iGEA modularity also falls slightly in hypoxia/vaccination (more cross-collection bridges), but stays above ORA at every step in all three datasets.
- Figures: `02 Figures/networks/ORA_vs_iGEA_GO_and_10_collections.png/.pdf` (HIV, app layout) and `ORA_vs_iGEA_hypoxia_vaccination_10_collections.png/.pdf`.
- CSBJ limits (spj.science.org/page/csbj/for-authors, read 20 Nov): Method articles ≤ 8 pages / 5,000 words incl. references and captions, ≤ 5 displayed items; Software/Web server articles < 15 pages / 10,000 words, ≤ 10 displayed items, < 50 references; Research articles ~15,000 words, 10 figures, 5 tables. Current v7 manuscript ≈ 8,500 words (v6 ≈ 6,000).

## Step 28 — Results: iGEA vs ORA section rewritten (20 Nov)
Subsection renamed "Comparison of iGEA and ORA": like-for-like ORA settings sentence (+ conventional ORA 124 terms), counts per collection, gene coverage ("with the same gene sets"), Figure 1, new paragraph on networks as collections are added (HIV: structure, density 6.4 vs 2.3 edges per gene, communities and the context added term by term), new paragraph on hypoxia and vaccination (incl. iGEA modularity caveat); LLM ORA paragraph rewritten (GO and 10-collection ORA answers). Supplementary Figure 7 (HIV, 4 panels) and 8 (hypoxia/vaccination) = PDFs in `03 Supplementary Material` + legends file; Supplementary Table 14 note on the GO rows (27/144 vs Figure 2 26/138). Built by rerunning /home/claude/ms4/results_edit.py (results_text.py) on the model-name version; reject view = v6; 32 comments, 27 fields.

## Step 29 — Additional reviewer analyses (20 Nov) → `step10_extras/`
Scripts share `common.py` (CLI data, background 23,391, 10 collections, defaults). All run on the VM in ≤ 3 min calls (background processes do not survive a call).
- **(1) Parameter sensitivity** (`s10a_sensitivity.py` → `s10a_parameter_sensitivity.tsv`; one parameter changed at a time; random lists generated with the changed parameters, 2,000 lists): HIV network significant in all 8 settings (cluster size p = 0.0005, the minimum; edges/gene p 0.005–0.012); kept terms 38–96 (p < 0.005: 46; p < 0.05: 96; overlap 2: 40; overlap 5: 38; sizes 15–500: 53; 5–1000: 52; ≤ 10 iterations: 54); Jaccard with the default kept terms 0.54–0.96, with the default genes 0.81–0.99. Random control: never significant (cluster p 0.46–1.0; 0–1 kept terms in every setting). → calibration holds by design (random lists matched to the parameters); settings change which terms represent the genes, not the genes or the network result.
- **(2) Chance selections by iteration and list size** (`s10b_iteration_fdr.py` → `s10b_chance_by_iteration.tsv`, `s10b_chance_by_list_size.tsv`): chance terms concentrate in early iterations (HIV, expected per random list: 5.9 at iteration 1, 3.6, 2.3, 1.5 …; observed 10, 9, 8, 6 …) and grow sub-linearly with list size (1.6 at 25 genes, 7.3 at 100, 17.5 at 300, 29.5 at 800 genes). The q-value already accounts for both (it is computed from the p-values of the same procedure on lists of the same size), so per-iteration rates are descriptive.
- **(3) Annotation-matched random lists** (`s10c_annotation_matched.py` + `matched_lib.py` → `s10c_matched_<dataset>.tsv`, `_terms.tsv`): input genes are better annotated than the background (mean gene sets per gene: HIV 45.9, hypoxia 58.3, vaccination 50.1; background 28.1). Lists matched gene-by-gene on annotation count (bin 0 + 9 quantile bins) select more chance terms (HIV 28.2 vs 16.1 per list; hypoxia 15.4 vs 5.7; vaccination 43.2 vs 22.0), so fewer terms pass q < 0.10 (HIV 35 vs 56; hypoxia 16 vs 32; vaccination 37 vs 70). The network remains significant in all three (cluster size/genes/edges r = 0 of 10,000, p < 0.0001; hypoxia terms r = 2; edges per gene r = 72, 98, 14). HIV terms lost under matching: all 6 WikiPathways, Reactome IL-12, retrograde TGN, GPI-anchored proteins, HIV transcription elongation, GOCC nuclear pore, clathrin-coated vesicle, spindle, autophagosome, cullin-RING, GOMF G protein activity, GTPase binding, hexosyltransferase, and the 6 GOBP terms. Caveat: matching on annotation also removes real biology (well-studied processes are well annotated), so this is a conservative sensitivity analysis, not a replacement for the uniform null.
- **(4) Uncertainty** (`s10d_uncertainty.py` → `s10d_q_bootstrap.tsv`, `s10d_network_p_intervals.tsv`): 200 bootstrap resamples of the random lists reproduce the kept sets (56/32/70); terms whose 95% q interval crosses 0.10: HIV 2 (Reactome HIV transcription elongation, q 0.097 [0.091–0.103]; GOBP positive regulation by host of viral process, 0.104), hypoxia 0, vaccination 6. Network p, exact 95% interval of the tail probability: HIV cluster size r = 0/10,000 → ≤ 0.00037; edges/gene r = 62 → 0.0048–0.0079; GO network edges/gene 0.0034–0.0061.
- (20 Nov, author) Parameter sensitivity → one sentence in the Results (robustness) + Supplementary Table 15 (`03 Supplementary Material/Supp Table 15 - Sensitivity to analysis settings.txt`). Sentence to insert at the end of the 90%-subsample paragraph when the docx is closed (the author was editing it in Word at 19:42).
- (20 Nov, 20:0x) Sentence inserted (tracked, Claude) at the end of the 90%-subsample paragraph, on top of the author's edited version (19:42); comment 224; reject view unchanged.

## Step 30 — Three HIV siRNA screens with iGEA (Brass, König, Zhou) (20 Nov) → `step11_three_screens/`
- Lists (`s11a_gene_lists.py`, from the papers' supplements in `v7/Konig and Zhou/`): König 2008 Table S3, 295 Gene IDs → 292 current symbols (285 validated in the background); Zhou 2008 Table S1, 390 primary hits minus those not expressed in activated T cells/macrophages = 312 → 311 symbols (307 validated). Inputs saved in `inputs/`. Brass = paper list (265).
- iGEA, defaults (10 collections, one-sided Fisher, 10,000 random lists, q < 0.10), CLI runs in step11: König 52/69 kept, cluster 148 genes/48 terms/473 edges/9 collections, p < 0.0001 (edges/gene 0.001), 8 communities (splicing; proteasome; transcription coregulators/Mediator; nuclear pore/import; RNA pol II/DNA repair; nucleases; tyrosine phosphatases; perinuclear ER). Zhou 32/62 kept, cluster 128/32/269/9, p < 0.0001 (edges/gene 0.0097), 5 communities (NF-κB/AKT signalling + autophagy; viral life cycle/RNA pol II + CD4/CXCR4; Mediator; mitochondrial complex I; K48 ubiquitination).
- Overlap (`s11b_compare.py`): genes pairwise Jaccard 0.017–0.029 (10–16 shared; 3 in all: MED6, MED7, RELA); kept terms 0.035–0.080 (3–8 shared; 2 in all: GOCC core Mediator, GOMF transcription coregulator activity).
- Communities (`s11b_community_similarity.tsv`, `s11b_community_best_matches.tsv`; footprint = union of the full gene sets of a community's terms; Jaccard vs 2,000 random term groups of the same sizes, `s11c_footprint_null.py`): 131 cross-screen pairs, 9 with p ≤ 0.01. Shared by all three: Mediator/transcription coregulators (J 0.56–0.78, p ≤ 0.0005); ubiquitin–proteasome (Brass cullin-RING–König proteasome 0.16, p 0.0025; –Zhou K48 ubiquitination 0.08, p 0.0065); RNA processing / RNA polymerase (Brass RNA exosome–König splicing 0.12, –Zhou RNA pol II 0.12; König RNA pol II–Zhou 0.10). Two screens: nuclear pore (Brass–König 0.35, 4 shared terms, p 0.001); NF-κB/AKT signalling (Brass–Zhou 0.095, p 0.027, hubs AKT1, RELA). Screen-specific: Brass Golgi transport and ER glycosylation/GPI; König nucleases, tyrosine phosphatases, perinuclear ER; Zhou mitochondrial complex I.
- **LLM, two turns, König and Zhou networks** (author ran `step11_three_screens/run_llm_screens.sh gpt-6.2-sol`; same prompt v5 and identical second-turn sentence as Brass; checks `s11d_llm_checks.py` → `s11d_llm_checks.tsv`, sections `s11d_llm_sections.json`): 12 answers, all completed; every answer covers every community (König 8/8, Zhou 5/5); no invented genes (one outside-knowledge mention: XPO1, König second turn, as the Rev export receptor). Blind turn: Brass 3/3 loss-of-function screen with a cell-surface/receptor readout (HIV not inferred); König 3/3 pooled loss-of-function screen with a growth/fitness or reporter readout (essential complexes; HIV not inferred); Zhou 3/3 functional screen with an infection-related readout, r1 "possibly HIV-related" — the Zhou network contains CD4, CXCR4, VPS4A and the term GOBP viral life cycle. Second turn, common to all three: the network fits an HIV dependency-factor screen; the model separates established HIV factors (Brass: CD4, NUP153/RANBP2, Tat/P-TEFb elongation; König: nuclear pore/export NUP214, RANBP2, NXF1; Zhou: CD4/CXCR4, VPS4A) from core complexes that may reflect cell fitness or the assay readout (Mediator, spliceosome, proteasome, complex I); all priorities use existing screen data or uninfected cells. Screen-specific "may be new": Brass Golgi/glycosylation surface composition and GPI anchors; König RNA processing–export handoff (NXF1, NUP98) and proteasome cell state; Zhou respiratory complex I / energy sensing (PRKAA1, AKT1) and Mediator selectivity.
- **Same comparison with ORA** (`s11e_ora_vs_igea_screens.py` → `s11e_ora_community_similarity.tsv`; FDR < 0.05, same gene sets/background/overlap): ORA terms Brass 87, König 246, Zhou 38; edges per gene 6.4 / 17.3 / 5.0; communities 2 / 2 / 3 (König: two blobs of 124 and 122 terms). Term overlap is higher than iGEA (Jaccard 0.05–0.14; 6 terms in all three) because redundant terms repeat across screens, but only 1 of 16 cross-screen community pairs matches at p ≤ 0.01 (Brass nuclear pore/protein localization blob ↔ König RNA processing blob, 54 vs 124 terms), so shared themes cannot be read from the ORA networks. iGEA: 9 of 131 pairs at p ≤ 0.01, resolving five specific shared themes (three in all screens). Bushman 2009 (pooled DAVID categories, ≥2 of 3 screens) already reported shared categories (Mediator, proteasome, RNA polymerase/splicing, nuclear pore/export, ER/Golgi transport, DNA repair …): iGEA recovers the same themes from each screen separately, as communities, without pooling or manual binning; ER/Golgi appears only in the Brass iGEA network.
- (20 Nov, 21:0x) Supplementary Table 16 (`03 Supplementary Material/Supp Table 16 - Three independent HIV screens.html/.pdf`; builder `build/build_supp_table_16.py`): networks, theme × screen grid, gene/term overlap, community matches (10 of 131 at p ≤ 0.01 — 9 strictly < 0.01), ORA comparison, LLM turn 1/turn 2 comparison. LLM answer documents 16B–C: `LLM responses/pdf/LLM_Konig_screen.pdf`, `LLM_Zhou_screen.pdf`. Manuscript: new Results subsection "Agreement between independent HIV screens" at the end of the Results (tracked, Claude; comment 225; citations as yellow [CITE] placeholders for the existing refs König, Zhou, Bushman).

## Step 31 — Missing supplementary items (21 Nov)
- Supp Table 3 (prompt v5): `03 Supplementary Material/Supp Table 03 - LLM prompt.html/.pdf/.txt` — full HIV GO prompt, network block, second-turn template, context sentences (builder `build/build_st03.py`).
- Supp Table 13 (DAVID): copied unchanged from v6 (`Supp Table 09 - DAVID cluster report.pdf`).
- Supp Figs 1–6 regenerated with v7 data at print size (builders `03 Supplementary Material/build/sf*.py`; legends `Supp Figures - legends.txt`).
  - SF1 per-collection both/iGEA-only/ORA-only: p<0.01 50/12/205; FDR<0.05 vs q<0.10 19/37/68 (matches Results text).
  - SF5 A–C: three random lists of 265 genes, no term filter (`step12_supp_figures/`: seeds 999 [= random control], 1000, 1001; CLI, 1,000 random lists, --q-threshold 1): largest clusters 30, 27, 67 nodes. D: HIV vs 10,000-list null (q<0.10): nodes 224 vs max 81; edges/gene 2.27, p = 0.0063.
  - SF6 (`step12_supp_figures/s12a_size_trend.py` → `s12a_size_trend.tsv`; 1,000 lists per size 50–1,000, engine `per_run_null.generate_null`): no filter, median largest cluster 7→147 nodes, edges/gene 1.2–1.47 (peak at 200, slight decline to 1.36 at 1,000); after q<0.10 median 0–7 nodes, edges/gene median ≤1. The v6 sentence "decrease in connectivity" does not hold as stated → text to update.
