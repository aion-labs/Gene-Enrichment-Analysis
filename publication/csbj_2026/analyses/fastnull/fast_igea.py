"""Prototype: vectorised iGEA engine (same statistics as the iGEA code, much faster).
Per collection: sparse set x gene matrix built once. Per iteration: overlaps for all sets in one
sparse mat-vec; exact two-sided Fisher p-values computed (vectorised) only for sets with
overlap >= min_overlap; pick lowest p (ties -> GMT order); remove its genes; repeat."""
import numpy as np
from scipy import sparse
from scipy.special import gammaln
from scipy.stats import fisher_exact

def _lchoose(a, b):
    return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)

def fisher_two_sided(k, K, n, M):
    """Fast exact two-sided Fisher p, vectorised over k, K (log-space hypergeometric pmf).
    Tables with probability <= P(observed) * (1 + 1e-10) are summed (SciPy uses 1 + 1e-14;
    1e-10 sits just above the log-space rounding error). Used to rank candidates; the selected
    term and any near-ties are then recomputed with scipy.stats.fisher_exact (see Collection.igea)."""
    k = np.asarray(k, float); K = np.asarray(K, float)
    lo = np.maximum(0, n + K - M); hi = np.minimum(K, n)
    width = int((hi - lo).max()) + 1
    x = lo[:, None] + np.arange(width)[None, :]
    valid = x <= hi[:, None]
    logden = _lchoose(M, n)
    lp = _lchoose(K[:, None], x) + _lchoose(M - K[:, None], n - x) - logden
    lp = np.where(valid, lp, -np.inf)
    lobs = _lchoose(K, k) + _lchoose(M - K, n - k) - logden
    keep = lp <= lobs[:, None] + np.log1p(1e-10)
    p = np.where(keep, np.exp(lp), 0.0).sum(1)
    return np.minimum(p, 1.0)

def scipy_fisher(k, K, n, M):
    return fisher_exact([[int(k), int(K - k)], [int(n - k), int(M - K - n + k)]])[1]

class Collection:
    def __init__(self, gmt_path, background, min_size=10, max_size=600):
        names, sets = [], []
        for line in open(gmt_path):
            parts = line.rstrip("\n").split("\t")
            g = parts[2:]
            if min_size <= len(g) <= max_size:
                names.append(parts[0]); sets.append(g)
        self.names = names
        self.genes = sorted(set().union(*map(set, sets)))          # G_c
        self.idx = {g: i for i, g in enumerate(self.genes)}
        self.K = np.array([len(set(s)) for s in sets], float)      # unique genes in the set, as in iGEA code (size filter uses raw GMT length)
        usets = [sorted(set(s)) for s in sets]
        rows = np.repeat(np.arange(len(usets)), [len(s) for s in usets])
        cols = np.array([self.idx[g] for s in usets for g in s])
        self.A = sparse.csr_matrix((np.ones(len(cols)), (rows, cols)), shape=(len(sets), len(self.genes)))
        self.M = len(set(background) & set(self.genes))
        self.members = [np.array(sorted(self.idx[g] for g in set(s))) for s in sets]

    def igea(self, gene_list, p_threshold=0.01, min_overlap=3, max_iterations=30):
        x = np.zeros(len(self.genes)); 
        for g in gene_list:
            i = self.idx.get(g)
            if i is not None: x[i] = 1
        out = []
        for it in range(1, max_iterations + 1):
            n = x.sum()
            if n == 0: break
            k = self.A @ x
            cand = np.flatnonzero(k >= min_overlap)
            if cand.size == 0: break
            p = fisher_two_sided(k[cand], self.K[cand], n, self.M)
            # exact SciPy p for the best candidate and any near-ties (within 1e-6 relative)
            near = np.flatnonzero(p <= p.min() * (1 + 1e-6))
            exact = np.array([scipy_fisher(k[cand[i]], self.K[cand[i]], n, self.M) for i in near])
            j = int(near[np.argmin(exact)]); pj = float(exact.min())   # first minimum = GMT order
            if pj >= p_threshold: break
            t = cand[j]
            removed = [self.genes[i] for i in self.members[t] if x[i] == 1]
            out.append((it, self.names[t], pj, removed))
            x[self.members[t]] = 0
        return out
