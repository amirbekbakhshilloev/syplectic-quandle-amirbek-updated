import sys, json, time
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from symq.links import LINKS
from symq.quandles import Zn, GF, symplectic_quandle, standard_J
from symq.coloring import Presentation, colorings_c, iter_colorings_c
from symq.invariants import enhanced_polynomial_stream, poly_str


def catalogue():
    Q = []
    for n in (4, 6, 8, 9):
        R = Zn(n)
        Q.append((f"Z{n}^2, J", R, standard_J(R, 2, 1)))
        for lam in range(2, n):
            if n % lam == 0:
                Q.append((f"Z{n}^2, {lam}J (degenerate)", R, standard_J(R, 2, lam)))
    for (p, k) in ((2, 2), (2, 3), (3, 2)):
        R = GF(p, k); Q.append((f"GF({p**k})^2, J", R, standard_J(R, 2, 1)))
    for p in (2, 3):
        R = Zn(p); Q.append((f"Z{p}^4, J+J", R, standard_J(R, 4, 1)))
    for p in (2, 3, 5):
        R = Zn(p); M = [[0, 1, 0], [R.neg[1], 0, 0], [0, 0, 0]]
        Q.append((f"Z{p}^3, J+0 (degenerate)", R, M))
        M2 = [[0, 1, 1], [R.neg[1], 0, 1], [R.neg[1], R.neg[1], 0]]
        Q.append((f"Z{p}^3, generic alternating (degenerate)", R, M2))
    return Q


if __name__ == "__main__":
    out = open(sys.argv[1], 'a')
    only = sys.argv[2].split(';') if len(sys.argv) > 2 else None
    names = ['H', 'L1', 'L2']
    for desc, R, M in catalogue():
        if only and desc not in only: continue
        Qd = symplectic_quandle(R, M, name=desc)
        Qd.check_axioms()
        rec = dict(quandle=desc, size=Qd.N)
        for nm in names:
            P = Presentation(LINKS[nm], nm)
            t = time.time()
            total = colorings_c(P, Qd, count_only=True)
            enh = None
            if total <= 2_000_000:
                C, total2 = enhanced_polynomial_stream(Qd, iter_colorings_c(P, Qd))
                assert total2 == total
                enh = poly_str(C)
            rec[nm] = dict(count=total, enhanced=enh, sec=round(time.time() - t, 1))
            print(desc, nm, rec[nm], flush=True)
        key = lambda nm: (rec[nm]['count'], rec[nm]['enhanced'])
        rec['separates_H_L1'] = key('H') != key('L1')
        rec['separates_H_L2'] = key('H') != key('L2')
        rec['separates_L1_L2'] = key('L1') != key('L2')
        out.write(json.dumps(rec) + "\n"); out.flush()
