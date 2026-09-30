H = [(1, 1, 3, 1), (2, 3, 1, 2), (3, 2, 4, 3), (4, 4, 3, 4)]

L1 = [
    (1, 1, 26, 2), (2, 27, 1, 26), (3, 2, 3, 1), (4, 4, 1, 3), (5, 5, 3, 4),
    (6, 6, 7, 5), (7, 7, 6, 8), (8, 28, 6, 29), (9, 6, 28, 3), (10, 8, 7, 9),
    (11, 10, 9, 7), (12, 9, 10, 11), (13, 12, 11, 10), (14, 11, 12, 13), (15, 14, 13, 12),
    (16, 15, 14, 19), (17, 15, 13, 22), (18, 16, 14, 20), (19, 16, 13, 21), (20, 17, 21, 14),
    (21, 18, 21, 13), (22, 17, 22, 20), (23, 18, 22, 19), (24, 23, 20, 21), (25, 24, 20, 22),
    (26, 23, 19, 26), (27, 24, 19, 25), (28, 28, 29, 30), (29, 31, 30, 29), (30, 30, 31, 32),
    (31, 33, 32, 31), (32, 32, 33, 34), (33, 35, 34, 33), (34, 36, 35, 40), (35, 36, 34, 38),
    (36, 37, 35, 41), (37, 37, 34, 39), (38, 42, 39, 35), (39, 43, 39, 34), (40, 42, 38, 41),
    (41, 43, 38, 40), (42, 44, 40, 25), (43, 44, 41, 39), (44, 45, 40, 27), (45, 45, 41, 38),
]

L2 = [
    (1, 1, 26, 2), (2, 27, 1, 26), (3, 2, 3, 1), (4, 4, 1, 3), (5, 5, 3, 4),
    (6, 6, 7, 5), (7, 7, 6, 8), (8, 28, 6, 29), (9, 6, 28, 47), (10, 8, 7, 9),
    (11, 10, 9, 7), (12, 9, 10, 11), (13, 12, 11, 10), (14, 11, 12, 13), (15, 14, 13, 12),
    (16, 15, 14, 19), (17, 15, 13, 22), (18, 16, 14, 20), (19, 16, 13, 21), (20, 17, 21, 14),
    (21, 18, 21, 13), (22, 17, 22, 20), (23, 18, 22, 19), (24, 23, 20, 21), (25, 24, 20, 22),
    (26, 23, 19, 26), (27, 24, 19, 25), (28, 28, 29, 30), (29, 31, 30, 29), (30, 30, 31, 32),
    (31, 33, 32, 31), (32, 32, 33, 34), (33, 35, 34, 33), (34, 36, 35, 40), (35, 36, 34, 38),
    (36, 37, 35, 41), (37, 37, 34, 39), (38, 42, 39, 35), (39, 43, 39, 34), (40, 42, 38, 41),
    (41, 43, 38, 40), (42, 44, 40, 25), (43, 44, 41, 39), (44, 45, 40, 46), (45, 45, 41, 38),
    (46, 67, 48, 47), (47, 48, 67, 49), (48, 68, 67, 69), (49, 67, 68, 3), (50, 49, 48, 50),
    (51, 51, 50, 48), (52, 50, 51, 52), (53, 53, 52, 51), (54, 52, 53, 54), (55, 55, 54, 53),
    (56, 56, 55, 58), (57, 56, 54, 61), (58, 57, 55, 59), (59, 57, 54, 60), (60, 62, 60, 55),
    (61, 63, 60, 54), (62, 62, 61, 59), (63, 63, 61, 58), (64, 64, 58, 46), (65, 64, 59, 60),
    (66, 65, 58, 66), (67, 65, 59, 61), (68, 68, 69, 70), (69, 71, 70, 69), (70, 70, 71, 72),
    (71, 73, 72, 71), (72, 72, 73, 74), (73, 75, 74, 73), (74, 76, 75, 78), (75, 76, 74, 81),
    (76, 77, 75, 79), (77, 77, 74, 80), (78, 82, 80, 75), (79, 83, 80, 74), (80, 82, 81, 79),
    (81, 83, 81, 78), (82, 84, 78, 66), (83, 84, 79, 80), (84, 85, 78, 27), (85, 85, 79, 81),
]

LINKS = {"H": H, "L1": L1, "L2": L2}


def allen_swenberg(n):
    if n == 1:
        return list(L1)
    B = [r for r in L2 if r[0] >= 46]
    rels = [r for r in L2 if r[0] < 46]
    top_in = {47: 47, 46: 46}
    next_arc = 86; next_cross = 86
    for copy in range(n - 2):
        ren = {}
        for x in range(48, 86): ren[x] = next_arc; next_arc += 1
        new_in47, new_in46 = next_arc, next_arc + 1; next_arc += 2
        ren[47] = new_in47; ren[46] = new_in46
        ren[3] = top_in[47]; ren[27] = top_in[46]
        blk = []
        for (k, a, b, c) in B:
            blk.append((next_cross, ren.get(a, a), ren.get(b, b), ren.get(c, c))); next_cross += 1
        rels = [(k, a, b, new_in47) if k == 9 else ((k, a, b, new_in46) if k == 44 else (k, a, b, c))
                for (k, a, b, c) in rels]
        rels += blk
        top_in = {47: new_in47, 46: new_in46}
    rels += B
    rels.sort()
    arcs = sorted({x for r in rels for x in r[1:]})
    amap = {a: i + 1 for i, a in enumerate(arcs)}
    return [(i + 1, amap[a], amap[b], amap[c]) for i, (k, a, b, c) in enumerate(rels)]


def presentation_isomorphism(rels_a, rels_b, fixed):
    A = [tuple(r[1:]) for r in rels_a]; B = set(tuple(r[1:]) for r in rels_b)
    if len(A) != len(B): return None
    arcsA = sorted({x for r in A for x in r}); arcsB = {x for r in B for x in r}
    phi = dict(fixed); used = set(phi.values())
    byarc = {}
    for r in A:
        for x in r: byarc.setdefault(x, []).append(r)

    def consistent():
        for r in A:
            if all(x in phi for x in r) and tuple(phi[x] for x in r) not in B:
                return False
        return True

    def extend():
        changed = True
        while changed:
            changed = False
            for r in A:
                known = [x in phi for x in r]
                if sum(known) == 2:
                    cands = [s for s in B if all((not known[i]) or s[i] == phi[r[i]] for i in range(3))]
                    i = known.index(False)
                    vals = {s[i] for s in cands}
                    if len(vals) == 0: return False
                    if len(vals) == 1:
                        v = vals.pop()
                        if r[i] in phi:
                            if phi[r[i]] != v: return False
                        else:
                            if v in used: return False
                            phi[r[i]] = v; used.add(v); changed = True
        return consistent()

    def search():
        snapshot = dict(phi)
        if not extend():
            phi.clear(); phi.update(snapshot); used.clear(); used.update(phi.values()); return False
        free = [x for x in arcsA if x not in phi]
        if not free: return True
        x = free[0]
        for v in sorted(arcsB - used):
            phi[x] = v; used.add(v)
            if search(): return True
            phi.clear(); phi.update(snapshot); used.clear(); used.update(phi.values())
            if not extend(): break
        phi.clear(); phi.update(snapshot); used.clear(); used.update(phi.values())
        return False

    return dict(phi) if search() else None
