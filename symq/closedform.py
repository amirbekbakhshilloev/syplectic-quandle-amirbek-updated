import numpy as np
import sympy as sp


def minimal_polynomial(A, maxdeg=40):
    n = A.shape[0]
    powers = [np.identity(n, dtype=object)]
    for d in range(1, maxdeg + 1):
        powers.append(powers[-1].dot(A))
        V = sp.Matrix([[int(P.flat[i]) for P in powers[:d]] for i in range(n * n)])
        rhs = sp.Matrix([int(powers[d].flat[i]) for i in range(n * n)])
        rows = [i for i in range(n * n) if any(V[i, j] != 0 for j in range(d)) or rhs[i] != 0]
        Vr = V.extract(rows, list(range(d))); rr = rhs.extract(rows, [0])
        try:
            sol, params = Vr.gauss_jordan_solve(rr)
        except ValueError:
            continue
        if params.shape[0] > 0:
            sol = sol.subs({p: 0 for p in params})
        coeffs = [-sol[j] for j in range(d)] + [1]
        acc = sum((coeffs[i] * powers[i] for i in range(d + 1)), np.zeros((n, n), dtype=object))
        assert all(v == 0 for v in acc.flat), "verification failed"
        assert all(sp.Rational(c).q == 1 for c in coeffs)
        return [int(c) for c in coeffs], powers
    raise RuntimeError("degree bound exceeded")


def c_values(Rm, A, nmax):
    R = Rm.toarray().astype(object)
    P = np.identity(A.shape[0], dtype=object); out = []
    for n in range(1, nmax + 1):
        out.append(int((R * P).sum()))
        P = P.dot(A)
    return out
