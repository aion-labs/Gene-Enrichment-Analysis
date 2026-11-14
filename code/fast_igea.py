"""
Fast iterative enrichment engine.

Implements exactly the same selection procedure as ``IterativeEnrichment`` (iter_enrichment.py):
for each collection, gene sets are filtered by size, the input list is restricted to the genes of
the retained sets, every set is tested against the remaining genes, sets with overlap below
``min_overlap`` are skipped, the lowest p-value is selected (ties broken by GMT order), its
overlapping genes are removed, and the procedure repeats until no set reaches ``p_threshold``.

It differs only in how the work is done:
  * each collection is stored once as a sparse gene-set x gene matrix, so the overlaps of all
    sets with the remaining genes are a single matrix-vector product;
  * p-values are computed (vectorised) only for sets that meet the minimum overlap;
  * the best candidate and any near-ties are recomputed with the same SciPy function as
    ``enrichment.compute_pvalue``, so the selected term and its p-value are identical.

It returns only what is needed to build the network (term, iteration, p-value, genes removed),
which makes it suitable for generating random-list null distributions inside each run.
"""
from typing import Iterable, List, Optional, Tuple

import numpy as np
from scipy import sparse
from scipy.special import gammaln
from scipy.stats import chi2_contingency, fisher_exact, hypergeom

FISHER = "Fisher's Exact Test"
HYPERGEOMETRIC = "Hypergeometric Test"
CHI_SQUARED = "Chi-squared Test"


def _lchoose(a, b):
    return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)


def _fisher_two_sided_fast(k, K, n, M):
    """Vectorised two-sided Fisher p-values (log-space). Used only to rank candidates."""
    lo = np.maximum(0, n + K - M)
    hi = np.minimum(K, n)
    width = int((hi - lo).max()) + 1
    x = lo[:, None] + np.arange(width)[None, :]
    valid = x <= hi[:, None]
    logden = _lchoose(M, n)
    lp = _lchoose(K[:, None], x) + _lchoose(M - K[:, None], n - x) - logden
    lp = np.where(valid, lp, -np.inf)
    lobs = _lchoose(K, k) + _lchoose(M - K, n - k) - logden
    keep = lp <= lobs[:, None] + np.log1p(1e-10)
    return np.minimum(np.where(keep, np.exp(lp), 0.0).sum(1), 1.0)


def exact_pvalue(method: str, k: int, K: int, n: int, M: int) -> float:
    """p-value computed exactly as in enrichment.compute_pvalue."""
    table = [[k, K - k], [n - k, M - K - n + k]]
    if method == FISHER:
        return float(fisher_exact(table)[1])
    if method == HYPERGEOMETRIC:
        return float(hypergeom.sf(k - 1, M, K, n))
    if method == CHI_SQUARED:
        return float(chi2_contingency(table)[1])
    raise ValueError(f"Unsupported p_value_method: {method}")


class FastCollection:
    """One gene-set collection prepared for repeated iterative enrichment."""

    def __init__(self, terms: List[dict], background_genes: Iterable[str], name: str = "",
                 min_term_size: int = 10, max_term_size: int = 600):
        # Same filter as Enrichment: size = raw GMT length
        kept = [t for t in terms if min_term_size <= t["size"] <= max_term_size]
        self.name = name
        self.term_names = [t["name"] for t in kept]
        usets = [sorted(set(t["genes"])) for t in kept]
        self.genes = sorted(set().union(*map(set, usets))) if usets else []
        self.index = {g: i for i, g in enumerate(self.genes)}
        # Same K as compute_pvalue: number of unique genes in the set
        self.K = np.array([len(s) for s in usets], dtype=float)
        rows = np.repeat(np.arange(len(usets)), [len(s) for s in usets])
        cols = np.array([self.index[g] for s in usets for g in s], dtype=int)
        self.A = sparse.csr_matrix((np.ones(len(cols)), (rows, cols)),
                                   shape=(len(usets), len(self.genes)))
        self.members = [np.array([self.index[g] for g in s], dtype=int) for s in usets]
        # Same M as Enrichment: background genes present in the retained sets
        self.M = len(set(background_genes) & set(self.genes))

    @classmethod
    def from_library(cls, library, background_genes, min_term_size=10, max_term_size=600):
        """Build from a ``GeneSetLibrary``."""
        return cls(library.library, background_genes, library.name, min_term_size, max_term_size)

    def _pvalues(self, method, k, K, n):
        if method == FISHER:
            return _fisher_two_sided_fast(k, K, n, self.M)
        if method == HYPERGEOMETRIC:
            return hypergeom.sf(k - 1, self.M, K, n)
        if method == CHI_SQUARED:
            return np.array([exact_pvalue(method, int(a), int(b), int(n), self.M) for a, b in zip(k, K)])
        raise ValueError(f"Unsupported p_value_method: {method}")

    def iterate(self, gene_list: Iterable[str], p_threshold: float = 0.01, min_overlap: int = 3,
                max_iterations: Optional[int] = None,
                method: str = FISHER) -> List[Tuple[int, str, float, List[str]]]:
        """Run iterative enrichment on ``gene_list``.

        Returns a list of (iteration, term name, iteration p-value, genes removed).
        """
        if not self.term_names:
            return []
        x = np.zeros(len(self.genes))
        for g in gene_list:
            i = self.index.get(g)
            if i is not None:
                x[i] = 1.0
        out = []
        iteration = 1
        while max_iterations is None or iteration <= max_iterations:
            n = x.sum()
            if n == 0:
                break
            k = self.A @ x
            cand = np.flatnonzero(k >= min_overlap)
            if cand.size == 0:
                break
            p = self._pvalues(method, k[cand], self.K[cand], n)
            # Recompute the best candidate and near-ties exactly; first minimum = GMT order
            near = np.flatnonzero(p <= p.min() * (1 + 1e-6))
            exact = [exact_pvalue(method, int(k[cand[i]]), int(self.K[cand[i]]), int(n), self.M) for i in near]
            j = int(np.argmin(exact))
            p_best = exact[j]
            if p_best >= p_threshold:
                break
            t = cand[near[j]]
            removed_idx = self.members[t][x[self.members[t]] == 1]
            out.append((iteration, self.term_names[t], p_best, [self.genes[i] for i in removed_idx]))
            x[removed_idx] = 0.0
            iteration += 1
        return out
