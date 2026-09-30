from collections import defaultdict


class Presentation:
    def __init__(self, rels, name=""):
        self.name = name
        arcs = sorted({x for r in rels for x in r[1:]})
        self.labels = arcs
        self.idx = {a: i for i, a in enumerate(arcs)}
        self.n = len(arcs)
        self.raw = list(rels)
        self.rels = [(self.idx[a], self.idx[b], self.idx[c]) for (_, a, b, c) in rels]
        self.by_arc = defaultdict(list)
        for ri, (a, b, c) in enumerate(self.rels):
            for x in {a, b, c}:
                self.by_arc[x].append(ri)


def _bits(m):
    while m:
        low = m & -m
        yield low.bit_length() - 1
        m ^= low


def colorings(P, Q, fixed=None):
    op = Q.op.tolist(); N = Q.N
    rels = P.rels; by_arc = P.by_arc; n = P.n
    FULL = (1 << N) - 1
    dom = [FULL] * n
    if fixed:
        for lab, allowed in fixed.items():
            m = 0
            for v in allowed: m |= 1 << v
            dom[P.idx[lab]] = m
    out = []

    def revise(ri, dom):
        a, b, c = rels[ri]
        Da, Db, Dc = dom[a], dom[b], dom[c]
        nA = nB = nC = 0
        la = list(_bits(Da))
        for y in _bits(Db):
            hit = False
            for x in la:
                z = op[x][y]
                if (Dc >> z) & 1:
                    nC |= 1 << z; nA |= 1 << x; hit = True
            if hit: nB |= 1 << y
        changed = []
        if nA != Da:
            dom[a] = nA; changed.append(a)
        if nB != Db:
            dom[b] = nB; changed.append(b)
        if nC != Dc:
            dom[c] = nC; changed.append(c)
        return changed

    def propagate(dom, queue):
        inq = set(queue); queue = list(queue)
        while queue:
            ri = queue.pop(); inq.discard(ri)
            a, b, c = rels[ri]
            fulls = (dom[a] == FULL) + (dom[b] == FULL) + (dom[c] == FULL)
            if fulls >= 2: continue
            for x in revise(ri, dom):
                if dom[x] == 0: return False
                for rj in by_arc[x]:
                    if rj not in inq:
                        inq.add(rj); queue.append(rj)
        return True

    def dfs(dom):
        best = -1; bsize = N + 1
        for i in range(n):
            s = bin(dom[i]).count("1")
            if 1 < s < bsize:
                best, bsize = i, s
                if s == 2: break
        if best < 0:
            sol = tuple(d.bit_length() - 1 for d in dom)
            for a, b, c in rels:
                if op[sol[a]][sol[b]] != sol[c]:
                    return
            out.append(sol); return
        for v in _bits(dom[best]):
            d2 = dom[:]
            d2[best] = 1 << v
            if propagate(d2, by_arc[best]):
                dfs(d2)

    if propagate(dom, range(len(rels))):
        dfs(dom)
    return out


def verify(P, Q, sols):
    op = Q.op
    for s in sols:
        for a, b, c in P.rels:
            assert op[s[a], s[b]] == s[c]
    return True


import os, subprocess, shutil
_HERE = os.path.dirname(os.path.abspath(__file__))
_BIN = os.path.join(_HERE, "colorsolve")


def _ensure_binary():
    src = os.path.join(_HERE, "colorsolve.c")
    if not os.path.exists(_BIN) or os.path.getmtime(_BIN) < os.path.getmtime(src):
        cc = shutil.which("gcc") or shutil.which("cc") or shutil.which("clang")
        subprocess.check_call([cc, "-O3", "-o", _BIN, src])


def colorings_c(P, Q, fixed=None, count_only=False, select=None, branch_first=None):
    _ensure_binary()
    N = Q.N
    parts = [f"{N} {P.n} {len(P.rels)}", " ".join(map(str, Q.op.flatten().tolist()))]
    parts += [f"{a} {b} {c}" for a, b, c in P.rels]
    fixed = fixed or {}
    parts.append(str(len(fixed)))
    for lab, allowed in fixed.items():
        allowed = sorted(allowed)
        parts.append(f"{P.idx[lab]} {len(allowed)} " + " ".join(map(str, allowed)))
    args = [_BIN] + (["-c"] if count_only else [])
    if branch_first is not None:
        args += ["-b", str(len(branch_first))] + [str(P.idx[lab]) for lab in branch_first]
    if select is not None:
        args += ["-s", str(len(select))] + [str(P.idx[lab]) for lab in select]
    res = subprocess.run(args, input="\n".join(parts) + "\n", capture_output=True, text=True, check=True)
    lines = res.stdout.strip().split("\n")
    count = int(lines[-1].split()[1])
    if count_only:
        return count
    sols = [tuple(map(int, l.split())) for l in lines[:-1]]
    assert len(sols) == count
    return sols


def iter_colorings_c(P, Q, fixed=None, select=None):
    _ensure_binary()
    N = Q.N
    parts = [f"{N} {P.n} {len(P.rels)}", " ".join(map(str, Q.op.flatten().tolist()))]
    parts += [f"{a} {b} {c}" for a, b, c in P.rels]
    fixed = fixed or {}
    parts.append(str(len(fixed)))
    for lab, allowed in fixed.items():
        allowed = sorted(allowed)
        parts.append(f"{P.idx[lab]} {len(allowed)} " + " ".join(map(str, allowed)))
    args = [_BIN]
    if select is not None:
        args += ["-s", str(len(select))] + [str(P.idx[lab]) for lab in select]
    import tempfile
    with tempfile.TemporaryFile(mode="w+") as fin:
        fin.write("\n".join(parts) + "\n"); fin.flush(); fin.seek(0)
        proc = subprocess.Popen(args, stdin=fin, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        n = 0
        for line in proc.stdout:
            if line.startswith("#"):
                assert int(line.split()[1]) == n, "solver count mismatch"
                continue
            n += 1
            yield tuple(map(int, line.split()))
        proc.wait()
        assert proc.returncode == 0
