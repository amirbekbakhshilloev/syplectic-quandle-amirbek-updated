from symq.links import H
from symq.coloring import Presentation, colorings
from symq.invariants import naive_polynomial, enhanced_polynomial, poly_str

H_R2 = [(1, 5, 3, 1), (2, 3, 1, 2), (3, 2, 4, 3), (4, 4, 3, 4), (5, 1, 4, 6), (6, 5, 4, 6)]


def demo(Q):
    out = {}
    for name, rels in (("H", H), ("H after R2", H_R2)):
        S = colorings(Presentation(rels, name), Q)
        out[name] = dict(count=len(S), naive=poly_str(naive_polynomial(S)),
                         corrected=poly_str(enhanced_polynomial(Q, S)))
    return out


if __name__ == "__main__":
    from symq.quandles import Zn, symplectic_quandle, standard_J
    for p in (3, 5):
        R = Zn(p); Q = symplectic_quandle(R, standard_J(R, 2, 1))
        for k, v in demo(Q).items():
            print(f"p={p}  {k:11s} |Hom|={v['count']:5d}  arc-count 'polynomial': {v['naive']:40s} corrected: {v['corrected']}")
