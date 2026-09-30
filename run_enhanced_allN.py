import sys, json, time, collections
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import sympy as sp
from symq.quandles import Zn, GF, symplectic_quandle, standard_J
from symq.transfer import subquandles, hom_counts_restricted

x = sp.symbols('x'); n = sp.symbols('n', integer=True, nonnegative=True)


def run(desc, Q):
    t0 = time.time()
    op = Q.op.tolist(); N = Q.N
    subs = subquandles(Q)
    def trivial(m):
        els = [i for i in range(N) if (m >> i) & 1]
        return all(op[a][b] == a for a in els for b in els)
    K = 12
    f = {}; mus = []
    for m in subs:
        if trivial(m):
            f[m] = [bin(m).count('1') ** 3] * K; mus.append([-1, 1])
        else:
            vals, mu = hom_counts_restricted(Q, m, K); f[m] = vals; mus.append(mu)
    order = sorted(subs, key=lambda m: bin(m).count('1'))
    e = {}
    for m in order:
        e[m] = [f[m][k] - sum(e[s][k] for s in order if s != m and (s & m) == s and s in e) for k in range(K)]
    coeff = collections.defaultdict(lambda: [0] * K)
    for m in subs:
        sz = bin(m).count('1')
        for k in range(K): coeff[sz][k] += e[m][k]
    L = sp.Integer(1)
    for mu in mus: L = sp.lcm(L, sum(c * x**i for i, c in enumerate(mu)))
    L = sp.Poly(L, x)
    d = L.degree()
    assert d + 2 <= K, "need more terms"
    roots = sp.roots(L.as_expr(), x)
    out = {}
    for sz in sorted(coeff):
        seq = coeff[sz]
        if all(v == 0 for v in seq): continue
        unknowns = []; expr = 0
        for r, mult in roots.items():
            for j in range(mult):
                a = sp.Symbol(f"a{len(unknowns)}"); unknowns.append(a)
                expr += a * n**j * r**n
        sol = sp.solve([sp.Eq(expr.subs(n, k + 1), seq[k]) for k in range(len(unknowns))], unknowns, dict=True)[0]
        closed = sp.simplify(expr.subs(sol))
        assert all(closed.subs(n, k + 1) == seq[k] for k in range(K)), "closed form fails"
        out[sz] = str(closed)
    rec = dict(quandle=desc, subquandles=len(subs), recurrence=str(sp.factor(L.as_expr())),
               coefficients={str(k): v for k, v in out.items()},
               values_n1_to_4={str(sz): coeff[sz][:4] for sz in sorted(coeff) if any(coeff[sz])},
               sec=round(time.time() - t0, 1))
    return rec


if __name__ == "__main__":
    todo = {"Z4^2, J": lambda: symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1)),
            "GF(4)^2, J": lambda: symplectic_quandle(GF(2, 2), standard_J(GF(2, 2), 2, 1)),
            "Z2^4, J+J": lambda: symplectic_quandle(Zn(2), standard_J(Zn(2), 4, 1))}
    for desc in (sys.argv[1:] or todo):
        rec = run(desc, todo[desc]())
        print(json.dumps(rec), flush=True)
        open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'results', 'enhanced_allN.jsonl'), 'a').write(json.dumps(rec) + "\n")
