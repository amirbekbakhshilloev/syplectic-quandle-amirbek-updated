#!/usr/bin/env python3
import sys, os, json, time, collections
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.makedirs(os.path.join(HERE, "results"), exist_ok=True)

from symq.links import LINKS, H, L1, L2, allen_swenberg
from symq.quandles import Zn, GF, symplectic_quandle, standard_J
from symq.coloring import Presentation, colorings, colorings_c, iter_colorings_c, verify
from symq.invariants import (closure, naive_polynomial, enhanced_polynomial,
                             enhanced_polynomial_stream, poly_str)

FULL = "--full" in sys.argv
AXIS = {"H": [2, 3], "L1": [1, 2], "L2": [1, 2]}


def save(name, obj):
    with open(os.path.join(HERE, "results", name + ".json"), "w") as f:
        json.dump(obj, f, indent=1, default=str)
    print(f"  -> results/{name}.json", flush=True)


def T8():
    R = Zn(2); return symplectic_quandle(R, [[0, 1, 0], [1, 0, 0], [0, 0, 0]], name="(Z/2)^3, J+0")


def Zp2(p):
    R = Zn(p); return symplectic_quandle(R, standard_J(R, 2, 1), name=f"(Z_{p})^2")


def stage_links():
    import warnings; warnings.filterwarnings("ignore")
    from symq.topology import components, linking_numbers, alexander
    from symq.planar import search_fast, search_local, pd_code
    from symq.jones import bracket, span
    out = {}
    links = {"H": H, "L1": L1, "L2": L2, "L3": allen_swenberg(3)}
    out["L2_generated_equals_notebook"] = sorted(allen_swenberg(2)) == sorted(L2)
    from symq.links import presentation_isomorphism
    T0 = [r for r in L2 if 6 <= r[0] <= 45]; B = [r for r in L2 if r[0] >= 46]
    phi = presentation_isomorphism(T0, B, {5: 47, 26: 46, 47: 3, 46: 27})
    out["block_B_isomorphic_to_T0_copy_c6_c45"] = (phi is not None and
        set(tuple(phi[x] for x in r[1:]) for r in T0) == set(tuple(r[1:]) for r in B))
    print("B isomorphic to the copy of T0 (c6..c45):", out["block_B_isomorphic_to_T0_copy_c6_c45"], flush=True)
    for name, R in links.items():
        rec = {}
        comps = components(R)
        rec["crossings"] = len(R)
        rec["arcs_per_component"] = [len(c) for c in comps]
        nc, lk = linking_numbers(R)
        rec["linking_numbers(i_under_j, j_under_i)"] = {f"{i}-{j}": v for (i, j), v in lk.items()}
        rec["alexander_1var_minor"] = str(alexander(R)[0])
        d, c0, T = alexander(R, True)
        rec["alexander_multivar_minor"] = str(d)
        rec["multivar_minor_dropped_component"] = c0
        if name == "L1":
            sols = search_fast(R, verbose=False)
            rec["planar_realisations_exhaustive"] = {"orientation_choices": 8,
                                                     "orders_per_orientation": 589824,
                                                     "planar": [s[4] for s in sols], "orientations": [s[0] for s in sols]}
            sol = sols[0]
        elif name != "H":
            sol = search_local(R, restarts=300, verbose=False)
        else:
            sol = search_local(R, restarts=50, verbose=False)
        if sol is not None:
            pd, _ = pd_code(R, sol)
            rec["planar_found"] = True
            rec["kauffman_bracket_span"] = span(bracket(pd)) - 4
            try:
                import spherogram
                L = spherogram.Link(pd)
                rec["snappy_linking_matrix"] = L.linking_matrix()
                subl = []
                for keep in ([0, 1], [0, 2], [1, 2]):
                    S = L.sublink(keep); S.simplify("global")
                    subl.append((keep, len(S.crossings)))
                rec["sublink_crossings_after_simplify"] = subl
            except Exception as e:
                rec["snappy"] = f"not available ({e.__class__.__name__})"
            if name in ("L1", "L2"):
                rec["pd_code"] = pd
        else:
            rec["planar_found"] = False
        out[name] = rec
        print(name, {k: v for k, v in rec.items() if k != "pd_code"}, flush=True)
    save("links", out)


def stage_v1():
    from symq.reidemeister import H_R2
    out = {"table": {}, "R2": {}, "H_formula": {}}
    for p in (3, 5):
        Q = Zp2(p); Q.check_axioms()
        for name in ("H", "L1", "L2"):
            P = Presentation(LINKS[name], name)
            S = colorings_c(P, Q); verify(P, Q, S)
            if p == 3:
                assert sorted(S) == sorted(colorings(P, Q)), "C and Python solvers disagree"
            rec = dict(count=len(S), arc_colour_count=poly_str(naive_polynomial(S)),
                       corrected=poly_str(enhanced_polynomial(Q, S)))
            out["table"][f"p={p} {name}"] = rec; print(p, name, rec, flush=True)
        for nm, rels in (("H", H), ("H after R2", H_R2)):
            S = colorings(Presentation(rels, nm), Q)
            out["R2"][f"p={p} {nm}"] = dict(count=len(S), arc_colour_count=poly_str(naive_polynomial(S)),
                                            corrected=poly_str(enhanced_polynomial(Q, S)))
        print(out["R2"], flush=True)
    for p in (2, 3, 5, 7, 11):
        Q = Zp2(p); S = colorings_c(Presentation(H), Q)
        C = enhanced_polynomial(Q, S)
        formula = collections.Counter({1: p * p, 2: 3 * p * (p * p - 1), 3: p * (p - 2) * (p * p - 1),
                                       p * p: p * (p - 1) * (p * p - 1)})
        out["H_formula"][p] = dict(count=len(S), formula_count=2 * p**4 - p**2,
                                   polynomial=poly_str(C), matches=(C == formula and len(S) == 2 * p**4 - p**2))
    print(out["H_formula"], flush=True)
    out["python_equals_C"] = {}
    for desc, Q in (("(Z/4)^2", symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))), ("T8", T8())):
        for name in ("H", "L1", "L2"):
            P = Presentation(LINKS[name], name)
            out["python_equals_C"][f"{desc} {name}"] = sorted(colorings(P, Q)) == sorted(colorings_c(P, Q))
    print(out["python_equals_C"], flush=True)
    save("v1_vs_corrected", out)


def axis0_sym(rels, p):
    Q = Zp2(p); z = Q.zero; e1 = Q.elements.index((1, 0))
    P = Presentation(rels)
    c0 = colorings_c(P, Q, fixed={1: {z}, 2: {z}, 7: {z}}, count_only=True)
    c1 = colorings_c(P, Q, fixed={1: {z}, 2: {z}, 7: {e1}}, count_only=True)
    return c0 + (p * p - 1) * c1


def axis_nonzero_sym(rels, p):
    Q = Zp2(p); e1 = Q.elements.index((1, 0)); P = Presentation(rels)
    S = colorings_c(P, Q, fixed={1: {e1}})
    collinear = all(len(closure(Q, set(s))) <= p for s in S)
    return (p * p - 1) * len(S), collinear


def stage_primes():
    out = {}
    plist_L1 = [2, 3, 5, 7, 11, 13, 17, 19, 23] + ([29, 31] if FULL else [])
    plist_L2 = [2, 3, 5, 7, 11] + ([13] if FULL else [])
    nz_L1 = [p for p in plist_L1 if p <= 13]
    nz_L2 = [p for p in plist_L2 if p <= 7]
    for name, plist, nzl in (("L1", plist_L1, nz_L1), ("L2", plist_L2, nz_L2)):
        for p in plist:
            t = time.time()
            a0 = axis0_sym(LINKS[name], p)
            rec = dict(H=2 * p**4 - p**2, H_axis0=p**4, axis0=a0, extra_sky_colourings=a0 - p**4)
            if p in nzl:
                nz, col = axis_nonzero_sym(LINKS[name], p)
                rec.update(axis_nonzero=nz, axis_nonzero_equals_H=(nz == (p * p - 1) * p * p),
                           axis_nonzero_collinear=col, total=a0 + nz)
            rec["separates_from_H"] = (a0 != p**4) if p in nzl else (a0 > 2 * p**4 - p**2 or None)
            rec["sec"] = round(time.time() - t, 1)
            out[f"{name} p={p}"] = rec; print(name, p, rec, flush=True)
    p = 7; Q = Zp2(p)
    SH = colorings_c(Presentation(H), Q); CH = enhanced_polynomial(Q, SH)
    collinear_H = sum(v for k, v in CH.items() if k <= p)
    for name in ("L1", "L2"):
        tot = out[f"{name} p=7"]["total"]
        C = collections.Counter({k: v for k, v in CH.items() if k <= p}); C[p * p] = tot - collinear_H
        out[f"corrected (Z7)^2 {name}"] = poly_str(C)
    out["corrected (Z7)^2 H"] = poly_str(CH)
    S = colorings_c(Presentation(L1), Q)
    out["corrected (Z7)^2 L1 direct"] = poly_str(enhanced_polynomial(Q, S))
    save("primes", out)


def stage_family():
    from symq.transfer import family_counts, family_counts_sym
    out = {}
    for p in (2, 3, 5, 7):
        Q = Zp2(p)
        counts, info, M, Rm = family_counts_sym(Q, p, 8)
        out[f"(Z_{p})^2"] = dict(H=2 * p**4 - p**2, counts=counts, **info)
        print(p, out[f"(Z_{p})^2"], flush=True)
    for desc, Q in (("(Z/4)^2", symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))),
                    ("(F_4)^2", symplectic_quandle(GF(2, 2), standard_J(GF(2, 2), 2, 1))),
                    ("(Z/2)^4", symplectic_quandle(Zn(2), standard_J(Zn(2), 4, 1))),
                    ("(Z/2)^3, J+0", T8())):
        counts, info, Rt, Bt = family_counts(Q, 8)
        Hc = colorings_c(Presentation(H), Q, count_only=True)
        out[desc] = dict(H=Hc, counts=counts, R_colourings=info["R_colourings"], B_colourings=info["B_colourings"],
                         identity_block_colouring_for_all_boundary_data=(info["boundary_pairs_with_identity_colouring"] == Q.N ** 2))
        print(desc, out[desc], flush=True)
    Q = symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))
    out["direct check (Z/4)^2 n=1,2,3"] = [colorings_c(Presentation(allen_swenberg(n)), Q, count_only=True) for n in (1, 2, 3)]
    save("family", out)


def stage_closed():
    import sympy as sp, numpy as np, scipy.sparse as sps
    from symq.transfer import boundary_tallies
    from symq.closedform import minimal_polynomial, c_values
    x = sp.symbols("x")
    out = {}
    for desc, Q, formula in (("(Z/4)^2", symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1)), lambda n: 640 + 96 * 16**n),
                             ("(Z/2)^4", symplectic_quandle(Zn(2), standard_J(Zn(2), 4, 1)), lambda n: 736 + 480 * 100**n),
                             ("(Z/2)^3, J+0", T8(), lambda n: 176 + 48 * 16**n),
                             ("(F_4)^2", symplectic_quandle(GF(2, 2), standard_J(GF(2, 2), 2, 1)), None)):
        N = Q.N
        Rt, Bt, sec = boundary_tallies(Q)
        def mat(T):
            A = np.zeros((N * N, N * N), dtype=object)
            for (a, b, c, d), m in T.items(): A[a * N + b, c * N + d] += m
            return A
        M = mat(Bt); Rm = sps.csr_matrix(mat(Rt).astype(np.int64))
        mu, _ = minimal_polynomial(M, maxdeg=16)
        cv = c_values(Rm, M, len(mu) + 3)
        rec = dict(minimal_polynomial=str(sp.factor(sum(c * x**i for i, c in enumerate(mu)))), c=cv)
        if formula:
            rec["formula_matches_first_terms"] = all(cv[k] == formula(k + 1) for k in range(len(cv)))
            rec["H_equals_formula_at_n0"] = colorings_c(Presentation(H), Q, count_only=True) == formula(0)
        out[desc] = rec; print(desc, rec, flush=True)
    import run_enhanced_allN as E
    for desc, Q in (("Z4^2, J", symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))),
                    ("Z2^4, J+J", symplectic_quandle(Zn(2), standard_J(Zn(2), 4, 1))),
                    ("Z2^3, J+0", T8())):
        out[f"enhanced all n {desc}"] = E.run(desc, Q); print(out[f"enhanced all n {desc}"], flush=True)
    from symq.transfer import clasp_check
    for desc, Q in (("(Z/4)^2", symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))), ("T8", T8()),
                    ("(Z/2)^4", symplectic_quandle(Zn(2), standard_J(Zn(2), 4, 1))),
                    ("(F_4)^2", symplectic_quandle(GF(2, 2), standard_J(GF(2, 2), 2, 1))),
                    ("(Z_3)^2", Zp2(3)), ("(Z_5)^2", Zp2(5))):
        out[f"clasp check {desc}"] = clasp_check(Q); print(desc, out[f"clasp check {desc}"], flush=True)
    Z4 = symplectic_quandle(Zn(4), standard_J(Zn(4), 2, 1))
    for desc, Q, formula in (("(Z/4)^2", Z4, lambda n: {1: 16, 2: 216, 3: 192, 4: 120, 5: 48, 6: 48, 13: 96 * 16**n}),
                             ("T8", T8(), lambda n: {1: 8, 2: 96, 3: 72, 4: 48, 7: 48 * 16**n - 48})):
        res = {}
        for n in (0, 1, 2, 3):
            rels = H if n == 0 else allen_swenberg(n)
            C, tot = enhanced_polynomial_stream(Q, iter_colorings_c(Presentation(rels), Q))
            f = {k: v for k, v in formula(n).items() if v}
            res[n] = dict(direct=poly_str(C), matches_formula=(dict(C) == f))
        out[f"enhanced direct n<=3 {desc}"] = res; print(desc, res, flush=True)
    save("closed_forms", out)


def stage_other():
    import scan_other
    wanted = {"Z4^2, J", "Z4^2, 2J (degenerate)", "Z6^2, J", "Z6^2, 2J (degenerate)", "Z6^2, 3J (degenerate)",
              "GF(4)^2, J", "Z2^4, J+J", "Z2^3, J+0 (degenerate)", "Z3^3, J+0 (degenerate)"}
    if FULL:
        wanted |= {"GF(8)^2, J", "GF(9)^2, J", "Z8^2, J", "Z9^2, J", "Z3^4, J+J"}
    out = {}
    for desc, R, M in scan_other.catalogue():
        if desc not in wanted: continue
        Q = symplectic_quandle(R, M, name=desc); Q.check_axioms()
        rec = dict(size=Q.N)
        for nm in ("H", "L1", "L2"):
            rec[nm] = colorings_c(Presentation(LINKS[nm]), Q, count_only=True)
        rec["separates_H_L1"] = rec["H"] != rec["L1"]; rec["separates_H_L2"] = rec["H"] != rec["L2"]
        rec["separates_L1_L2"] = rec["L1"] != rec["L2"]
        out[desc] = rec; print(desc, rec, flush=True)
    save("other", out)


STAGES = dict(links=stage_links, v1=stage_v1, primes=stage_primes, family=stage_family,
              closed=stage_closed, other=stage_other)

if __name__ == "__main__":
    todo = [a for a in sys.argv[1:] if not a.startswith("--")] or list(STAGES)
    for s in todo:
        print(f"=== stage {s} ===", flush=True); t = time.time()
        STAGES[s]()
        print(f"=== stage {s} done in {time.time() - t:.0f}s ===", flush=True)
