import numpy as np, networkx as nx
def pull_in(pos, terms):
    G = nx.Graph(); G.add_nodes_from(pos)
    G.add_edges_from((("T", i), ("G", g)) for i, x in enumerate(terms) for g in x["genes"])
    comps = sorted(nx.connected_components(G), key=len, reverse=True)
    P = {n: np.array(v, float) for n, v in pos.items()}
    main = np.array([P[n] for n in comps[0]]); c0 = main.mean(0); R = np.max(np.linalg.norm(main - c0, axis=1))
    for comp in comps[1:]:
        pts = np.array([P[n] for n in comp]); cc = pts.mean(0); r = np.max(np.linalg.norm(pts - cc, axis=1)) if len(pts) > 1 else 0
        d = cc - c0; d = d / (np.linalg.norm(d) or 1); new = c0 + d * (R + r + 0.12 * R)
        for n in comp: P[n] = P[n] - cc + new
    return {n: tuple(v) for n, v in P.items()}
