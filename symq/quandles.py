import itertools
import numpy as np


class Zn:
    def __init__(self, n):
        self.n = n; self.size = n; self.name = f"Z{n}"
        self.add = [[(a + b) % n for b in range(n)] for a in range(n)]
        self.mul = [[(a * b) % n for b in range(n)] for a in range(n)]
        self.neg = [(-a) % n for a in range(n)]
        self.char = None

    def elt(self, v):
        return v % self.n


class GF:
    IRRED = {(2, 2): [1, 1, 1], (2, 3): [1, 1, 0, 1], (3, 2): [2, 2, 1],
             (2, 4): [1, 1, 0, 0, 1], (5, 2): [2, 0, 1], (7, 2): [1, 0, 1]}

    def __init__(self, p, k):
        self.p, self.k = p, k; q = p ** k
        self.size = q; self.name = f"GF{q}"
        f = self.IRRED[(p, k)]
        def dec(x): return [(x // p ** i) % p for i in range(k)]
        def enc(c): return sum((ci % p) * p ** i for i, ci in enumerate(c))
        def pmul(a, b):
            prod = [0] * (2 * k - 1)
            for i, ai in enumerate(a):
                for j, bj in enumerate(b):
                    prod[i + j] += ai * bj
            for d in range(2 * k - 2, k - 1, -1):
                c = prod[d] % p
                if c:
                    for i in range(k + 1):
                        prod[d - k + i] -= c * f[i]
            return [x % p for x in prod[:k]]
        self.add = [[enc([x + y for x, y in zip(dec(a), dec(b))]) for b in range(q)] for a in range(q)]
        self.mul = [[enc(pmul(dec(a), dec(b))) for b in range(q)] for a in range(q)]
        self.neg = [enc([-x for x in dec(a)]) for a in range(q)]
        for a in range(1, q):
            assert any(self.mul[a][b] == 1 for b in range(q)), "not a field"

    def elt(self, v):
        return v


class Quandle:
    def __init__(self, op, name, elements=None, meta=None):
        self.op = np.array(op, dtype=np.int32)
        N = self.op.shape[0]
        self.N = N; self.name = name
        self.elements = elements; self.meta = meta or {}
        inv = np.zeros_like(self.op)
        for y in range(N):
            col = self.op[:, y]
            assert len(set(col.tolist())) == N, "right multiplication not bijective"
            inv[col, y] = np.arange(N)
        self.inv = inv

    def check_axioms(self):
        op = self.op; N = self.N
        assert all(op[x, x] == x for x in range(N)), "idempotency fails"
        for z in range(N):
            lhs = op[op, z]
            rhs = op[op[:, z][:, None], op[:, z][None, :]]
            assert np.array_equal(lhs, rhs), "self-distributivity fails"
        return True


def vectors(R, r):
    return list(itertools.product(range(R.size), repeat=r))


def symplectic_quandle(R, M, name=None):
    r = len(M)
    for i in range(r):
        assert M[i][i] == 0
        for j in range(r):
            assert R.add[M[i][j]][M[j][i]] == 0, "M not antisymmetric"
    V = vectors(R, r)
    index = {v: i for i, v in enumerate(V)}
    add, mul = R.add, R.mul

    def form(x, y):
        s = 0
        for i in range(r):
            if x[i] == 0: continue
            for j in range(r):
                if M[i][j] and y[j]:
                    s = add[s][mul[mul[x[i]][M[i][j]]][y[j]]]
        return s

    N = len(V)
    op = [[0] * N for _ in range(N)]
    for i, x in enumerate(V):
        for j, y in enumerate(V):
            c = form(x, y)
            z = tuple(add[x[k]][mul[c][y[k]]] for k in range(r))
            op[i][j] = index[z]
    Q = Quandle(op, name or f"Symp({R.name}^{r})", elements=V,
                meta={"ring": R.name, "rank": r, "M": [list(m) for m in M]})
    Q.form = form
    Q.zero = index[tuple([0] * r)]
    return Q


def standard_J(R, r, lam=1):
    M = [[0] * r for _ in range(r)]
    for b in range(0, r, 2):
        M[b][b + 1] = R.elt(lam) if isinstance(R, Zn) else lam
        M[b + 1][b] = R.neg[M[b][b + 1]]
    return M


def conjugation_quandle_from_perms(gens, name):
    def comp(p, q):
        return tuple(q[i] for i in p)
    idp = tuple(range(len(gens[0])))
    G = {idp}; frontier = [idp]
    while frontier:
        new = []
        for g in frontier:
            for s in gens:
                h = comp(g, s)
                if h not in G:
                    G.add(h); new.append(h)
        frontier = new
    G = sorted(G); idx = {g: i for i, g in enumerate(G)}
    def invp(p):
        q = [0] * len(p)
        for i, pi in enumerate(p): q[pi] = i
        return tuple(q)
    op = [[idx[comp(comp(invp(y), x), y)] for y in G] for x in G]
    return Quandle(op, name, elements=G)
