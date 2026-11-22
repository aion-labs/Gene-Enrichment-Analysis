"""Annotation-matched random-list generator (importable, so it works with spawned worker processes)."""
import numpy as np
_S = {}

def init(collections, background, need, bins, n_input, params):
    _S.update(collections=collections, background=background, need=need, bins=bins, n=n_input, params=params)

def gen(task):
    size, i = task
    rng = np.random.default_rng([2026, 7, _S["n"], i])        # matched list i (independent of size)
    lst = []
    for b, k in sorted(_S["need"].items()):
        pool = _S["bins"][b]
        lst += [pool[j] for j in rng.choice(len(pool), k, replace=False)]
    if size < len(lst):                                        # reduced-size lists: random subsets of matched lists
        lst = [lst[j] for j in np.random.default_rng([2026, 8, size, i]).choice(len(lst), size, replace=False)]
    out = []
    for c, col in enumerate(_S["collections"]):
        for it, term, p, removed in col.iterate(lst, **_S["params"]):
            out.append((c, term, it, p, removed))
    return out
