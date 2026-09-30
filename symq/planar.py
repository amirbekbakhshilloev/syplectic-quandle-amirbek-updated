import itertools, math
import numpy as np


def traverse_components(rels):
    ends = {}
    for k, a, b, c in rels:
        ends.setdefault(a, []).append(k); ends.setdefault(c, []).append(k)
    rel = {k: (a, b, c) for k, a, b, c in rels}
    seen = set(); comps = []
    for k0, a0, b0, c0 in rels:
        if k0 in seen: continue
        seq = []; cur = a0; k = k0
        while True:
            seq.append((cur, k)); seen.add(k)
            a, b, c = rel[k]
            nxt = c if cur == a else a
            other = [kk for kk in ends[nxt] if kk != k]
            k2 = other[0] if other else k
            cur, k = nxt, k2
            if k == k0 and cur == a0: break
        comps.append(seq)
    return comps


def build_instance(rels, orient_bits):
    comps = traverse_components(rels)
    rel = {k: (a, b, c) for k, a, b, c in rels}
    oriented = []
    for i, seq in enumerate(comps):
        if (orient_bits >> i) & 1:
            arcs = [x for x, _ in seq]; ks = [k for _, k in seq]
            n = len(seq)
            rseq = [(arcs[(j) % n], ks[(j - 1) % n]) for j in range(n, 0, -1)]
            oriented.append(rseq)
        else:
            oriented.append(seq)
    start = {}; end = {}
    for seq in oriented:
        n = len(seq)
        for j, (x, k) in enumerate(seq):
            end[x] = k
            start[seq[(j + 1) % n][0]] = k
    sign = {}
    for k, (a, b, c) in rel.items():
        if a == c: sign[k] = 0
        else: sign[k] = +1 if (end[a] == k and start[c] == k) else -1
    overs = {}
    for k, (a, b, c) in rel.items():
        overs.setdefault(b, []).append(k)
    return oriented, start, end, sign, overs


def faces_count(rels, oriented, sign, order):
    ks = sorted({k for k, *_ in rels}); kid = {k: i for i, k in enumerate(ks)}
    V = len(ks)
    def port(k, role):
        s = sign[k]
        if role == 'uin': return 0
        if role == 'uout': return 2
        if role == 'oin': return 3 if s > 0 else 1
        if role == 'oout': return 1 if s > 0 else 3
    other = {}
    for seq in oriented:
        visits = []
        for j, (x, k_end) in enumerate(seq):
            k_start = seq[j - 1][1]
            visits.append((k_start, 'uout'))
            for ko in order.get(x, []):
                visits.append((ko, 'oin')); visits.append((ko, 'oout'))
            visits.append((k_end, 'uin'))
        dep = None
        for (k, role) in visits:
            if role in ('uout', 'oout'):
                dep = (kid[k], port(k, role))
            else:
                arr = (kid[k], port(k, role))
                other[dep] = arr; other[arr] = dep
    seen = set(); F = 0
    for v in range(V):
        for p in range(4):
            h = (v, p)
            if h in seen: continue
            F += 1
            while h not in seen:
                seen.add(h)
                w, q = other[h]
                h = (w, (q + 1) % 4)
    E = 2 * V
    return V - E + F


def search(rels, verbose=True):
    sols = []
    for ob in range(8):
        oriented, start, end, sign, overs = build_instance(rels, ob)
        arcs_multi = [x for x, ks in overs.items() if len(ks) > 1]
        choices = [list(itertools.permutations(overs[x])) for x in arcs_multi]
        base = {x: list(ks) for x, ks in overs.items() if len(ks) == 1}
        total = math.prod(len(c) for c in choices)
        cnt = 0
        for combo in itertools.product(*choices):
            order = dict(base)
            for x, perm in zip(arcs_multi, combo): order[x] = list(perm)
            chi = faces_count(rels, oriented, sign, order)
            if chi == 2:
                sols.append((ob, order, sign, oriented))
                cnt += 1
        if verbose: print(f"orientation {ob:03b}: {total} orderings tried, {cnt} planar", flush=True)
    return sols


def search_fast(rels, orientations=range(8), verbose=True, max_keep=50):
    from numba import njit

    @njit(cache=False)
    def run(nV, arc_start, arc_end, arc_nover, arc_over_in, arc_over_out, arc_perm_off,
            perms, perm_len, multi_arcs, radix, total, keep):
        nA = arc_start.shape[0]
        other = np.empty(4 * nV, dtype=np.int64)
        seen = np.empty(4 * nV, dtype=np.uint8)
        digits = np.zeros(multi_arcs.shape[0], dtype=np.int64)
        found = np.full((keep, multi_arcs.shape[0]), -1, dtype=np.int64)
        nfound = 0
        for it in range(total):
            r = it
            for j in range(multi_arcs.shape[0]):
                digits[j] = r % radix[j]; r //= radix[j]
            for x in range(nA):
                m = arc_nover[x]
                prev = arc_start[x]
                if m == 0:
                    pass
                else:
                    pidx = 0
                    if m > 1:
                        j = arc_perm_off[x]
                        pidx = digits[j]
                    for t in range(m):
                        if m > 1:
                            s = perms[x, pidx, t]
                        else:
                            s = 0
                        hin = arc_over_in[x, s]; hout = arc_over_out[x, s]
                        other[prev] = hin; other[hin] = prev
                        prev = hout
                other[prev] = arc_end[x]; other[arc_end[x]] = prev
            for h in range(4 * nV): seen[h] = 0
            F = 0
            for h0 in range(4 * nV):
                if seen[h0]: continue
                F += 1
                h = h0
                while not seen[h]:
                    seen[h] = 1
                    w = other[h]
                    h = (w // 4) * 4 + ((w % 4) + 1) % 4
            if nV - 2 * nV + F == 2:
                if nfound < keep:
                    for j in range(multi_arcs.shape[0]): found[nfound, j] = digits[j]
                nfound += 1
        return nfound, found

    results = []
    for ob in orientations:
        oriented, start, end, sign, overs = build_instance(rels, ob)
        ks = sorted({k for k, *_ in rels}); kid = {k: i for i, k in enumerate(ks)}
        nV = len(ks)
        def port(k, role):
            s = sign[k]
            return {'uin': 0, 'uout': 2, 'oin': 3 if s > 0 else 1, 'oout': 1 if s > 0 else 3}[role]
        arcs = sorted(start)
        aid = {x: i for i, x in enumerate(arcs)}
        nA = len(arcs)
        maxo = max(len(v) for v in overs.values())
        arc_start = np.array([4 * kid[start[x]] + 2 for x in arcs], dtype=np.int64)
        arc_end = np.array([4 * kid[end[x]] + 0 for x in arcs], dtype=np.int64)
        arc_nover = np.array([len(overs.get(x, [])) for x in arcs], dtype=np.int64)
        arc_over_in = np.zeros((nA, maxo), dtype=np.int64); arc_over_out = np.zeros((nA, maxo), dtype=np.int64)
        for x in arcs:
            for t, k in enumerate(overs.get(x, [])):
                arc_over_in[aid[x], t] = 4 * kid[k] + port(k, 'oin')
                arc_over_out[aid[x], t] = 4 * kid[k] + port(k, 'oout')
        multi = [x for x in arcs if len(overs.get(x, [])) > 1]
        maxp = max([math.factorial(len(overs[x])) for x in multi] + [1])
        perms = np.zeros((nA, maxp, maxo), dtype=np.int64)
        radix = np.array([math.factorial(len(overs[x])) for x in multi], dtype=np.int64)
        arc_perm_off = np.full(nA, -1, dtype=np.int64)
        for j, x in enumerate(multi):
            arc_perm_off[aid[x]] = j
            for pi, perm in enumerate(itertools.permutations(range(len(overs[x])))):
                perms[aid[x], pi, :len(perm)] = perm
        total = int(np.prod(radix)) if len(radix) else 1
        nfound, found = run(nV, arc_start, arc_end, arc_nover, arc_over_in, arc_over_out, arc_perm_off,
                            perms, np.zeros(1, dtype=np.int64), np.array([aid[x] for x in multi], dtype=np.int64),
                            radix, total, max_keep)
        if verbose: print(f"orientation {ob:03b}: {total} orderings, {nfound} planar", flush=True)
        for row in found[:min(nfound, max_keep)]:
            order = {x: list(v) for x, v in overs.items()}
            for j, x in enumerate(multi):
                perm = list(itertools.permutations(range(len(overs[x]))))[row[j]]
                order[x] = [overs[x][s] for s in perm]
            results.append((ob, order, sign, oriented, nfound))
    return results


def pd_code(rels, sol):
    ob, order, sign, oriented = sol[:4]
    def port(k, role):
        s = sign[k]
        return {'uin': 0, 'uout': 2, 'oin': 3 if s > 0 else 1, 'oout': 1 if s > 0 else 3}[role]
    slots = {}
    label = 0
    comp_edges = []
    for seq in oriented:
        visits = []
        for j, (x, k_end) in enumerate(seq):
            k_start = seq[j - 1][1]
            visits.append((k_start, 'uout'))
            for ko in order.get(x, []):
                visits.append((ko, 'oin')); visits.append((ko, 'oout'))
            visits.append((k_end, 'uin'))
        first = label
        dep = None
        for (k, role) in visits:
            if role in ('uout', 'oout'):
                dep = (k, port(k, role))
            else:
                slots[dep] = label; slots[(k, port(k, role))] = label; label += 1
        comp_edges.append((first, label - 1))
    ks = sorted({k for k, *_ in rels})
    return [tuple(slots[(k, p)] for p in range(4)) for k in ks], comp_edges


def search_local(rels, orientations=range(8), restarts=200, seed=1, verbose=True):
    import random
    rng = random.Random(seed)
    for ob in orientations:
        oriented, start, end, sign, overs = build_instance(rels, ob)
        multi = [x for x, ks in overs.items() if len(ks) > 1]
        perms = {x: list(itertools.permutations(overs[x])) for x in multi}
        best_chi = -10**9
        for r in range(restarts):
            order = {x: list(ks) for x, ks in overs.items()}
            for x in multi: order[x] = list(rng.choice(perms[x]))
            chi = faces_count(rels, oriented, sign, order)
            improved = True
            while improved:
                improved = False
                xs = multi[:]; rng.shuffle(xs)
                for x in xs:
                    cur = order[x]; bestp = cur; bc = chi
                    for p in perms[x]:
                        order[x] = list(p)
                        c = faces_count(rels, oriented, sign, order)
                        if c > bc: bc, bestp = c, list(p)
                    order[x] = bestp
                    if bc > chi: chi = bc; improved = True
            best_chi = max(best_chi, chi)
            if chi == 2:
                if verbose: print(f"orientation {ob:03b}: planar realisation found at restart {r}", flush=True)
                return (ob, order, sign, oriented)
        if verbose: print(f"orientation {ob:03b}: best chi {best_chi}", flush=True)
    return None
