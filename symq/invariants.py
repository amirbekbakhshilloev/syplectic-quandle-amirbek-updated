from collections import Counter


def closure(Q, S):
    op = Q.op
    S = set(S); frontier = list(S)
    while frontier:
        new = []
        cur = list(S)
        for x in frontier:
            for y in cur:
                for z in (int(op[x, y]), int(op[y, x])):
                    if z not in S:
                        S.add(z); new.append(z)
        frontier = new
    return S


def orbit(Q, x, gens):
    op = Q.op; inv = Q.inv
    O = {x}; frontier = [x]
    while frontier:
        new = []
        for y in frontier:
            for g in gens:
                for z in (int(op[y, g]), int(inv[y, g])):
                    if z not in O:
                        O.add(z); new.append(z)
        frontier = new
    return O


def poly_str(counter, var="q"):
    if not counter: return "0"
    terms = []
    for m in sorted(counter):
        a = counter[m]
        if m == 0: terms.append(f"{a}")
        elif m == 1: terms.append(f"{a}{var}")
        else: terms.append(f"{a}{var}^{{{m}}}")
    return " + ".join(terms)


def naive_polynomial(sols):
    return Counter(len(set(s)) for s in sols)


def enhanced_polynomial(Q, sols, cache=None):
    cache = {} if cache is None else cache
    C = Counter()
    for s in sols:
        key = frozenset(s)
        if key not in cache:
            cache[key] = len(closure(Q, key))
        C[cache[key]] += 1
    return C


def enhanced_polynomial_stream(Q, sol_iter, naive=False):
    cache = {}
    C = Counter(); Cn = Counter(); total = 0
    for s in sol_iter:
        key = frozenset(s)
        g = cache.get(key)
        if g is None:
            g = len(closure(Q, key)); cache[key] = g
        C[g] += 1; total += 1
        if naive: Cn[len(key)] += 1
    return (C, Cn, total) if naive else (C, total)
