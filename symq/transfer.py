import time, collections
import numpy as np
from symq.links import L2
from symq.coloring import Presentation, colorings_c

R_RELS = [r for r in L2 if r[0] <= 45]
B_RELS = [r for r in L2 if r[0] >= 46]


def boundary_tallies(Q):
    PR = Presentation(R_RELS, 'R'); PB = Presentation(B_RELS, 'B')
    t = time.time()
    solR = colorings_c(PR, Q, select=[47, 46, 3, 27])
    solB = colorings_c(PB, Q, select=[47, 46, 3, 27])
    Rt = collections.Counter(solR); Bt = collections.Counter(solB)
    return Rt, Bt, time.time() - t


def family_counts(Q, nmax):
    Rt, Bt, sec = boundary_tallies(Q)
    N = Q.N
    import scipy.sparse as sp
    rows, cols, vals = [], [], []
    for (a, b, c, d), m in Bt.items():
        rows.append(a * N + b); cols.append(c * N + d); vals.append(m)
    M = sp.csr_matrix((vals, (rows, cols)), shape=(N * N, N * N), dtype=np.int64)
    Rrows, Rcols, Rvals = [], [], []
    for (a, b, c, d), m in Rt.items():
        Rrows.append(a * N + b); Rcols.append(c * N + d); Rvals.append(m)
    Rm = sp.csr_matrix((Rvals, (Rrows, Rcols)), shape=(N * N, N * N), dtype=np.int64)
    counts = []
    P = sp.identity(N * N, dtype=np.int64, format='csr')
    Pf = sp.identity(N * N, dtype=np.float64, format='csr'); Mf = M.astype(np.float64)
    for n in range(1, nmax + 1):
        if Pf.max() * max(1, Rm.max()) > 4e18: break
        counts.append(int(Rm.multiply(P).sum()))
        P = P @ M; Pf = Pf @ Mf
    diag = sum(1 for w in range(N * N) if M[w, w] > 0)
    return counts, dict(R_colourings=sum(Rt.values()), B_colourings=sum(Bt.values()),
                        boundary_pairs_with_identity_colouring=diag, sec=round(sec, 1)), Rt, Bt


def sl2_perms(Q, p):
    V = Q.elements; idx = {v: i for i, v in enumerate(V)}
    perms = []
    for a in range(p):
        for b in range(p):
            for c in range(p):
                for d in range(p):
                    if (a * d - b * c) % p != 1: continue
                    perms.append([idx[((a * x + b * y) % p, (c * x + d * y) % p)] for (x, y) in V])
    op = Q.op
    for g in perms[:5] + perms[-5:]:
        G = np.array(g)
        assert np.array_equal(G[op], op[G[:, None], G[None, :]]), "not an automorphism"
    return perms


def pair_orbits(N, perms):
    seen = {}; reps = []
    for x in range(N):
        for y in range(N):
            if (x, y) in seen: continue
            rep = (x, y); reps.append(rep)
            for g in perms:
                w = (g[x], g[y])
                if w not in seen: seen[w] = (rep, g)
    return reps, seen


def family_counts_sym(Q, p, nmax, verbose=False):
    import scipy.sparse as sp
    perms = sl2_perms(Q, p)
    N = Q.N
    reps, where = pair_orbits(N, perms)
    PR = Presentation(R_RELS, 'R'); PB = Presentation(B_RELS, 'B')
    rowsB = {}; rowsR = {}
    t0 = time.time()
    for (x, y) in reps:
        sb = colorings_c(PB, Q, fixed={47: {x}, 46: {y}}, select=[3, 27])
        sr = colorings_c(PR, Q, fixed={47: {x}, 46: {y}}, select=[3, 27])
        rowsB[(x, y)] = collections.Counter(sb); rowsR[(x, y)] = collections.Counter(sr)
        if verbose: print('  rep', (x, y), 'B', len(sb), 'R', len(sr), round(time.time() - t0, 1), flush=True)
    def expand(rows):
        r_, c_, v_ = [], [], []
        for w, (rep, g) in where.items():
            for (u1, u2), m in rows[rep].items():
                r_.append(w[0] * N + w[1]); c_.append(g[u1] * N + g[u2]); v_.append(m)
        return sp.csr_matrix((v_, (r_, c_)), shape=(N * N, N * N), dtype=np.int64)
    M = expand(rowsB); Rm = expand(rowsR)
    Mi = M; Ri = Rm
    counts = []
    P = sp.identity(N * N, dtype=np.int64, format='csr')
    Pf = sp.identity(N * N, dtype=np.float64, format='csr'); Mf = Mi.astype(np.float64)
    for n in range(1, nmax + 1):
        if Pf.max() > 1e17: break
        counts.append(int(Ri.multiply(P).sum()))
        P = P @ Mi; Pf = Pf @ Mf
    diag_ok = all(Mi[w, w] > 0 for w in range(N * N))
    info = dict(orbit_reps=len(reps), R_colourings=int(Ri.sum()), B_colourings=int(Mi.sum()),
                identity_block_colouring_for_all_boundary_data=diag_ok, sec=round(time.time() - t0, 1))
    return counts, info, Mi, Ri


def subquandles(Q):
    op = Q.op.tolist(); N = Q.N
    assert N <= 20
    subs = []
    for m in range(1, 1 << N):
        ok = True
        for x in range(N):
            if not (m >> x) & 1: continue
            row = op[x]
            for y in range(N):
                if (m >> y) & 1 and not (m >> row[y]) & 1:
                    ok = False; break
            if not ok: break
        if ok: subs.append(m)
    return subs


def restricted_tallies(Q, mask):
    allowed = {i for i in range(Q.N) if (mask >> i) & 1}
    PR = Presentation(R_RELS, 'R'); PB = Presentation(B_RELS, 'B')
    fR = {lab: allowed for lab in PR.labels}; fB = {lab: allowed for lab in PB.labels}
    solR = colorings_c(PR, Q, fixed=fR, select=[47, 46, 3, 27])
    solB = colorings_c(PB, Q, fixed=fB, select=[47, 46, 3, 27])
    return collections.Counter(solR), collections.Counter(solB)


def hom_counts_restricted(Q, mask, nmax):
    from symq.closedform import minimal_polynomial, c_values
    N = Q.N
    Rt, Bt = restricted_tallies(Q, mask)
    def mat(T):
        A = np.zeros((N * N, N * N), dtype=object)
        for (a, b, c, d), m in T.items(): A[a * N + b, c * N + d] += m
        return A
    import scipy.sparse as sp
    M = mat(Bt); Rm = sp.csr_matrix(mat(Rt).astype(np.int64))
    states = sorted({a * N + b for (a, b, c, d) in Bt} | {c * N + d for (a, b, c, d) in Bt} |
                    {a * N + b for (a, b, c, d) in Rt} | {c * N + d for (a, b, c, d) in Rt})
    Ms = M[np.ix_(states, states)]
    Rs = sp.csr_matrix(Rm.toarray()[np.ix_(states, states)])
    mu, _ = minimal_polynomial(Ms, maxdeg=30)
    return c_values(Rs, Ms, nmax), mu


CLASP = [r for r in L2 if r[0] <= 5]
T0_COPY = [r for r in L2 if 6 <= r[0] <= 45]


def clasp_check(Q):
    import numpy as np
    N = Q.N
    PC = Presentation(CLASP, 'clasp')
    solC = colorings_c(PC, Q, select=[3, 27, 5, 26])
    AC = np.zeros((N * N, N * N), dtype=object)
    for (a, b, c, d) in solC: AC[a * N + b, c * N + d] += 1
    Rt, Bt, _ = boundary_tallies(Q)
    M = np.zeros((N * N, N * N), dtype=object); Rm = np.zeros((N * N, N * N), dtype=object)
    for (a, b, c, d), m in Bt.items(): M[a * N + b, c * N + d] += m
    for (a, b, c, d), m in Rt.items(): Rm[a * N + b, c * N + d] += m
    PT = Presentation(T0_COPY, 'T0')
    solT = colorings_c(PT, Q, select=[5, 26, 47, 46])
    AT = np.zeros((N * N, N * N), dtype=object)
    for (a, b, c, d) in solT: AT[a * N + b, c * N + d] += 1
    from symq.links import H
    Hc = colorings_c(Presentation(H), Q, count_only=True)
    return dict(T0_copy_transfer_equals_block=bool((AT == M).all()),
                R_transpose_equals_AC_times_M=bool((Rm.T == AC.dot(M)).all()),
                trace_AC=int(np.trace(AC)), H=Hc, trace_AC_equals_H=int(np.trace(AC)) == Hc)
