import collections
import sympy as sp
from sympy.polys.matrices import DomainMatrix


def components(rels):
    rel = {k: (a, b, c) for k, a, b, c in rels}
    ends = collections.defaultdict(list)
    for k, a, b, c in rels:
        ends[a].append(k); ends[c].append(k)
    used = set(); comps = []
    for k0, a0, b0, c0 in rels:
        if k0 in used: continue
        steps = []; cur, k = a0, k0
        while True:
            a, b, c = rel[k]
            nxt = c if cur == a else a
            sign = +1 if cur == a else -1
            steps.append((cur, k, nxt, sign)); used.add(k)
            others = [kk for kk in ends[nxt] if kk != k]
            if not others:
                break
            k, cur = others[0], nxt
            if k == k0: break
        comps.append(steps)
    return comps


def component_of_arc(rels):
    m = {}
    for i, st in enumerate(components(rels)):
        for (x, k, y, s) in st:
            m[x] = i; m[y] = i
    return m


def linking_numbers(rels):
    comps = components(rels)
    comp = component_of_arc(rels)
    sign = {}
    for st in comps:
        for (x, k, y, s) in st: sign[k] = s
    under_over = collections.Counter()
    for k, a, b, c in rels:
        i, j = comp[a], comp[b]
        if i != j: under_over[(i, j)] += sign[k]
    ncomp = len(comps)
    lk = {}
    for i in range(ncomp):
        for j in range(i + 1, ncomp):
            lk[(i, j)] = (under_over[(i, j)], under_over[(j, i)])
    return ncomp, lk


def alexander(rels, multivariable=False):
    comp = component_of_arc(rels)
    arcs = sorted({x for r in rels for x in r[1:]})
    idx = {a: i for i, a in enumerate(arcs)}
    ncomp = max(comp.values()) + 1
    T = sp.symbols(f"t1:{ncomp + 1}") if multivariable else (sp.Symbol('t'),) * ncomp
    gens = tuple(dict.fromkeys(T))
    dom = sp.ZZ[gens]
    rows = []
    for k, a, b, c in rels:
        row = [dom.zero] * len(arcs)
        ti, tj = T[comp[a]], T[comp[b]]
        row[idx[a]] += dom.convert(1)
        row[idx[b]] += dom.convert(ti - 1)
        row[idx[c]] += dom.convert(-tj)
        rows.append(row)
    A = [r[1:] for r in rows[1:]]
    dm = DomainMatrix(A, (len(A), len(A)), dom)
    return sp.factor(dom.to_sympy(dm.det())), comp[arcs[0]], T
