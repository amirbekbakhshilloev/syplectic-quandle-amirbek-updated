from collections import defaultdict


def padd(p, q, shift=0, c=1):
    for e, v in q.items():
        p[e + shift] = p.get(e + shift, 0) + c * v
    return {e: v for e, v in p.items() if v}


def pmul(p, q):
    r = defaultdict(int)
    for e1, v1 in p.items():
        for e2, v2 in q.items():
            r[e1 + e2] += v1 * v2
    return {e: v for e, v in r.items() if v}


D_LOOP = {2: -1, -2: -1}


def bracket(pd, order=None):
    crossings = list(pd)
    if order is None:
        remaining = list(range(len(crossings))); order = []; open_edges = set()
        while remaining:
            best = max(remaining, key=lambda c: (sum(1 for e in crossings[c] if e in open_edges), -c))
            order.append(best); remaining.remove(best)
            for e in crossings[best]:
                if e in open_edges: open_edges.remove(e)
                else: open_edges.add(e)
    states = {frozenset(): {0: 1}}
    for ci in order:
        i, j, k, l = crossings[ci]
        new = defaultdict(dict)
        for pairing, poly in states.items():
            partner = {}
            for a, b in pairing:
                partner[a] = b; partner[b] = a
            for (pairs, shift) in (([(i, j), (k, l)], 1), ([(i, l), (j, k)], -1)):
                part = dict(partner)
                loops = 0
                for x, y in pairs:
                    ex = part.pop(x, None); ey = part.pop(y, None)
                    if ex is None and ey is None:
                        if x == y: loops += 1
                        else: part[x] = y; part[y] = x
                    elif ex is not None and ey is None:
                        if ex == y: loops += 1
                        else:
                            part.pop(ex, None); part[ex] = y; part[y] = ex
                    elif ex is None and ey is not None:
                        if ey == x: loops += 1
                        else:
                            part.pop(ey, None); part[ey] = x; part[x] = ey
                    else:
                        if ex == y:
                            loops += 1
                        else:
                            part.pop(ex, None); part.pop(ey, None)
                            part[ex] = ey; part[ey] = ex
                key = frozenset(frozenset((a, b)) for a, b in part.items() if a < b)
                key = frozenset(tuple(sorted(p)) for p in key)
                term = {e + shift: v for e, v in poly.items()}
                for _ in range(loops): term = pmul(term, D_LOOP)
                new[key] = padd(dict(new[key]), term)
        states = {k: v for k, v in new.items() if v}
    assert list(states.keys()) == [frozenset()] or not states
    return states.get(frozenset(), {})


def span(p):
    return max(p) - min(p)
